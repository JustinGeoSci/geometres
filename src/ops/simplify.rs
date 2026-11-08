// src/ops/simplify.rs

//! Ramer-Douglas-Peucker vertex reduction, vouched.
//!
//! # The contract
//!
//! RDP is **not topology-preserving**. Dropping vertices from a valid simple
//! polygon can make it cross itself, can push a hole outside its exterior,
//! can turn a spike into a sliver wound the wrong way, and can make the
//! parts of a multipolygon overlap. None of this is exotic: it is routine at
//! tolerances near the polygon's local feature size, which is where simplify
//! is used.
//!
//! So `simplify` takes the shape every constructive operation takes:
//! **construct, vouch exactly, and offer an unvouched result only through an
//! entry point that says so.** Every vouch decision is one the crate already
//! makes exactly:
//!
//! | Check | Decided by |
//! |---|---|
//! | ring construction (≥ 3 vertices, no repeated consecutive vertex) | [`LinearRing::new`] |
//! | simplicity: no crossing, pinch, fold-back or vertex on an edge | the exact self-contact test (`exact::first_self_contact`), `orient2d` on stored coordinates |
//! | winding and degeneracy | [`Polygon::new`], via the certified sign of [`LinearRing::signed_area`] |
//! | every hole inside the exterior; holes interior-disjoint | [`Polygon::new`], exact ring placement |
//! | multipolygon parts meet at most at points | [`MultiPolygon::new`] |
//!
//! This is assembly, not new numerics. The only inexact arithmetic is the
//! RDP distance test itself, *which* decides *which* vertices to drop; it
//! cannot make a result wrong, only different, because every result is then
//! judged exactly.
//!
//! # What the unvouched entry point can and cannot hand back
//!
//! An [`Unvouched<Polygon>`] holds a `Polygon`, and `Polygon::new` enforces
//! winding, degeneracy and hole placement. The one defect it does **not**
//! enforce is simplicity. So `simplify_unvouched` returns
//! `Ok(Unvouched(..))` exactly when the result is a constructible polygon
//! whose only defect is a ring that touches or crosses itself, and `Err`
//! otherwise. A hole pushed outside its exterior, holes that now overlap,
//! an inverted or collapsed ring, or (for a multipolygon) parts that now
//! overlap are not polygons at all, and are reported as the constructor's
//! error from both entry points. Nothing is repaired, dropped or reversed to
//! make a value fit.
//!
//! # Tolerance
//!
//! `tolerance` is a distance in the coordinates' units. A vertex is dropped
//! when its distance from the chord replacing it is **not greater than**
//! `tolerance`. Zero is the identity: the input comes back unchanged, which
//! is what "move nothing by more than 0" must mean for a vertex exactly on
//! the chord, where the two readings differ. A negative or non-finite
//! tolerance is an error rather than a request to keep everything.
//!
//! **The distance test is a selection rule, evaluated in `f64`.** A vertex
//! whose true distance lies within rounding of `tolerance` may go either
//! way. That choice is never a correctness claim: whichever vertices are
//! chosen, the result is vouched exactly or refused, so the `f64` test
//! decides *which* valid answer is returned, not whether it is valid.
//!
//! # Pinched holes
//!
//! A hole that passes through one of its vertices twice (a pinch) is
//! refused here as a self-contact. A caller that splits such a ring first
//! with [`ops::split_pinches`](crate::ops::split_pinches) hands over two or
//! more simple holes touching at the pinch, and those simplify normally.
//!
//! # No retry
//!
//! A refused result is refused. There is no retry at a reduced tolerance:
//! that would answer a question the caller did not ask, with a bound that
//! would be a tuned constant. A caller who wants that writes the loop.
//!
//! # Rings
//!
//! A ring has no start, so RDP on a ring cannot anchor a chord at "the first
//! and last vertex". The closing step is the degenerate-chord step of RDP
//! itself: the chord from vertex 0 to vertex 0 has zero length, so the
//! vertex farthest from vertex 0 is found and kept, and the two arcs between
//! them are simplified as open chains. Vertex 0 and its antipode are always
//! kept, so a ring keeps at least two vertices; whether it keeps a third is
//! for the vouch. Running RDP on the open chain from vertex 0 to vertex
//! `n − 1` instead would never consider the closing edge for removal, and
//! the result would depend on where the ring happens to start.
//!
//! # Lines
//!
//! A `LineString` is simplified as an open chain, endpoints kept; a closed
//! line (first vertex equal to last) is treated as a ring, so its closing
//! point is not privileged. Lines are vouched by their constructor only: a
//! line string may cross itself, so there is no simplicity to check.

use crate::GeometryError;
use crate::exact::kernel::first_self_contact;
use crate::types::{
    Coord2d, CoordSeq, LineString, LinearRing, MaybeVouched, MultiLineString, MultiPolygon,
    Polygon, Unvouched, UnvouchedReason,
};

// ===========================================================================
// The distance test
// ===========================================================================

/// Distance from `p` to the closed segment `ab`, or `None` where the
/// arithmetic did not produce a finite value.
///
/// `hypot` for the lengths (no spurious overflow or underflow in the
/// squares); the perpendicular distance is `|cross| / |ab|`,
/// which is one rounding of each factor. A non-finite result can only come
/// from an overflowing product at coordinates above about `1e154`; the
/// caller keeps such a vertex, which is the conservative reading.
pub(super) fn segment_distance(p: Coord2d, a: Coord2d, b: Coord2d) -> Option<f64> {
    let (abx, aby) = (b.x - a.x, b.y - a.y);
    let (apx, apy) = (p.x - a.x, p.y - a.y);
    let len = abx.hypot(aby);
    let d = if len == 0.0 {
        apx.hypot(apy)
    } else {
        // Projection parameter in units of |ab|², compared without dividing.
        // An overflow in either side would make the comparison meaningless,
        // so it is a "cannot measure" rather than a branch.
        let dot = apx * abx + apy * aby;
        let len2 = len * len;
        if !dot.is_finite() || !len2.is_finite() {
            return None;
        }
        if dot <= 0.0 {
            apx.hypot(apy)
        } else if dot >= len2 {
            (p.x - b.x).hypot(p.y - b.y)
        } else {
            (abx * apy - aby * apx).abs() / len
        }
    };
    d.is_finite().then_some(d)
}

// ===========================================================================
// RDP on an open chain, iteratively
// ===========================================================================

/// Marks in `keep` the vertices of `pts[lo..=hi]` that RDP retains, with
/// `pts[lo]` and `pts[hi]` always kept.
///
/// Iterative with an explicit stack: the recursion depth of textbook RDP is
/// the vertex count in the worst case (a monotone spiral), and a parcel ring
/// can have thousands of vertices.
///
/// The farthest vertex is the first one attaining the maximum, so the result
/// is deterministic and independent of evaluation order.
fn mark_chain(pts: &[Coord2d], lo: usize, hi: usize, tolerance: f64, keep: &mut [bool]) {
    keep[lo] = true;
    keep[hi] = true;
    let mut stack = vec![(lo, hi)];
    while let Some((lo, hi)) = stack.pop() {
        if hi <= lo + 1 {
            continue;
        }
        let (a, b) = (pts[lo], pts[hi]);
        let mut best: Option<(usize, f64)> = None;
        let mut unknown: Option<usize> = None;
        for (i, &p) in pts.iter().enumerate().take(hi).skip(lo + 1) {
            match segment_distance(p, a, b) {
                None => {
                    unknown.get_or_insert(i);
                }
                Some(d) => {
                    if best.is_none_or(|(_, bd)| d > bd) {
                        best = Some((i, d));
                    }
                }
            }
        }
        // A vertex whose distance could not be computed is kept and splits
        // the chain, **before** any finite maximum: it may be the farthest.
        let split = match (unknown, best) {
            (Some(i), _) => Some(i),
            (None, Some((i, d))) if d > tolerance => Some(i),
            _ => None,
        };
        if let Some(i) = split {
            keep[i] = true;
            stack.push((lo, i));
            stack.push((i, hi));
        }
    }
}

/// RDP over an open chain: the retained vertices, endpoints included.
fn simplify_chain(pts: &[Coord2d], tolerance: f64) -> Vec<Coord2d> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let mut keep = vec![false; pts.len()];
    mark_chain(pts, 0, pts.len() - 1, tolerance, &mut keep);
    pts.iter()
        .zip(&keep)
        .filter(|(_, k)| **k)
        .map(|(p, _)| *p)
        .collect()
}

/// RDP over a closed ring given without a closing vertex: the retained
/// vertices, in ring order, starting from vertex 0.
///
/// Vertex 0 is the anchor. Its antipode is the vertex farthest from it, by
/// the degenerate-chord step of RDP: the chord `0 → 0` has zero length, and
/// the distance to it is the distance to vertex 0. Both arcs are then open
/// chains. With fewer than three vertices there is nothing to simplify.
fn simplify_ring(pts: &[Coord2d], tolerance: f64) -> Vec<Coord2d> {
    let n = pts.len();
    if n < 3 {
        return pts.to_vec();
    }
    let anchor = pts[0];
    let mut far = 1;
    let mut far_d = -1.0;
    for (i, &p) in pts.iter().enumerate().skip(1) {
        // A non-finite distance means a vertex too far to measure: it is the
        // farthest, and the first such wins.
        match segment_distance(p, anchor, anchor) {
            None => {
                far = i;
                break;
            }
            Some(d) if d > far_d => {
                far = i;
                far_d = d;
            }
            Some(_) => {}
        }
    }
    // If the antipode is within tolerance of the anchor, every vertex is,
    // and the ring collapses to its anchor and antipode; the vouch refuses
    // that. Marking continues regardless so that the result is the RDP one.
    let mut keep = vec![false; n];
    mark_chain(pts, 0, far, tolerance, &mut keep);
    // The second arc runs from `far` through the end back to vertex 0: lay it
    // out as a contiguous chain by rotating.
    let mut second: Vec<Coord2d> = pts[far..].to_vec();
    second.push(pts[0]);
    let mut keep2 = vec![false; second.len()];
    mark_chain(&second, 0, second.len() - 1, tolerance, &mut keep2);
    for (k, &kept) in keep2.iter().enumerate().take(second.len() - 1) {
        if kept {
            keep[far + k] = true;
        }
    }
    pts.iter()
        .zip(&keep)
        .filter(|(_, k)| **k)
        .map(|(p, _)| *p)
        .collect()
}

pub(super) fn check_tolerance(tolerance: f64) -> Result<(), GeometryError> {
    if tolerance.is_finite() && tolerance >= 0.0 {
        Ok(())
    } else {
        Err(GeometryError::InvalidTolerance { tolerance })
    }
}

// ===========================================================================
// Polygons
// ===========================================================================

/// The rings of `polygon` after RDP, built as `LinearRing`s: exterior first.
///
/// # Errors
///
/// A ring that keeps fewer than three vertices cannot be a `LinearRing`;
/// its constructor's error is returned. Nothing else can fail here: a
/// subsequence of a ring without repeated consecutive vertices has none, and
/// coordinates are the input's own.
fn simplified_rings(polygon: &Polygon, tolerance: f64) -> Result<Vec<LinearRing>, GeometryError> {
    core::iter::once(polygon.exterior())
        .chain(polygon.holes())
        .map(|r| {
            LinearRing::new(CoordSeq::new(simplify_ring(
                r.coords().as_slice(),
                tolerance,
            )))
        })
        .collect()
}

/// The first ring of `rings` that touches or crosses itself: its index
/// (0 the exterior, 1.. the holes) and the two edges in contact.
fn self_contact(rings: &[LinearRing]) -> Option<(usize, (usize, usize))> {
    rings
        .iter()
        .enumerate()
        .find_map(|(ring, r)| first_self_contact(r.coords().as_slice()).map(|edges| (ring, edges)))
}

/// Construct and vouch one polygon.
///
/// `Ok(Vouched)` when every check passes; `Ok(Unvouched)` when the polygon
/// exists but a ring touches or crosses itself, with the first such ring as
/// the reason; `Err` when no polygon can be built from the result.
fn build_polygon(
    polygon: &Polygon,
    tolerance: f64,
) -> Result<MaybeVouched<Polygon>, GeometryError> {
    let mut rings = simplified_rings(polygon, tolerance)?;
    let contact = self_contact(&rings);
    let exterior = rings.remove(0);
    // `Polygon::new` is the rest of the vouch: certified winding, the
    // exact degeneracy test, hole placement and hole disjointness, all exact.
    let built = Polygon::new(exterior, rings)?;
    Ok(match contact {
        None => MaybeVouched::Vouched(built),
        Some((0, edges)) => MaybeVouched::Unvouched(Unvouched::new(
            built,
            UnvouchedReason::SelfIntersecting { edges },
        )),
        Some((ring, edges)) => MaybeVouched::Unvouched(Unvouched::new(
            built,
            UnvouchedReason::RingSelfIntersecting { ring, edges },
        )),
    })
}

impl Polygon {
    /// Ramer-Douglas-Peucker simplification, vouched: the result is a valid
    /// polygon or an error. See the [module documentation](self) for the
    /// contract, the vouch, and the tolerance.
    ///
    /// Each ring is simplified independently; every retained vertex is one
    /// of the input's own coordinates, unchanged.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::InvalidTolerance`] for a negative or non-finite
    ///   tolerance.
    /// - [`GeometryError::Unvouched`] with a `SelfIntersecting` reason
    ///   (exterior) or `RingSelfIntersecting` reason (a hole, by ring index)
    ///   when a simplified ring touches or crosses itself. The result exists
    ///   and is available from [`Polygon::simplify_unvouched`].
    /// - The constructor's error when no polygon can be built from the
    ///   result: a ring collapsed below three vertices
    ///   ([`GeometryError::InsufficientVertices`]), lost its area
    ///   ([`GeometryError::DegenerateRing`]) or inverted
    ///   ([`GeometryError::IncorrectWindingOrder`],
    ///   [`GeometryError::IncorrectHoleWindingOrder`]); a hole pushed outside
    ///   its exterior ([`GeometryError::HoleNotInsideExterior`]); holes that
    ///   now overlap ([`GeometryError::OverlappingHoles`]). These are errors
    ///   from both entry points.
    pub fn simplify(&self, tolerance: f64) -> Result<Polygon, GeometryError> {
        check_tolerance(tolerance)?;
        if tolerance == 0.0 {
            return Ok(self.clone());
        }
        build_polygon(self, tolerance)?.vouched()
    }

    /// [`Polygon::simplify`], returning a self-intersecting result wrapped in
    /// [`Unvouched`] instead of refusing it.
    ///
    /// The vouch still runs: a result it passes comes back as
    /// [`MaybeVouched::Vouched`], one that fails only on simplicity as
    /// [`MaybeVouched::Unvouched`] with the first self-contact as the reason
    /// (exterior first, then holes in order; the reason carries the two
    /// edges in contact).
    ///
    /// # Errors
    ///
    /// [`GeometryError::InvalidTolerance`], and every constructor error
    /// listed on [`Polygon::simplify`]: a result that is not a polygon at all
    /// cannot be returned as one, vouched or not.
    pub fn simplify_unvouched(
        &self,
        tolerance: f64,
    ) -> Result<MaybeVouched<Polygon>, GeometryError> {
        check_tolerance(tolerance)?;
        if tolerance == 0.0 {
            return Ok(MaybeVouched::Vouched(self.clone()));
        }
        build_polygon(self, tolerance)
    }
}

// ===========================================================================
// Multipolygons
// ===========================================================================

impl MultiPolygon {
    /// Ramer-Douglas-Peucker simplification of every part, vouched.
    ///
    /// Each part is simplified and vouched as a [`Polygon`], and the parts
    /// are then checked against one another by [`MultiPolygon::new`]: parts
    /// that were disjoint can overlap, or share an edge, after their
    /// boundaries move, so part validity is part of this vouch.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::InvalidTolerance`].
    /// - [`GeometryError::Unvouched`] with a
    ///   [`UnvouchedReason::PartSelfIntersecting`] reason naming the first
    ///   part and ring that touch or cross themselves; the result is
    ///   available from [`MultiPolygon::simplify_unvouched`].
    /// - Every per-part constructor error of [`Polygon::simplify`], and
    ///   [`GeometryError::OverlappingParts`] when two simplified parts
    ///   overlap or share a stretch of boundary. Errors from both entry
    ///   points.
    pub fn simplify(&self, tolerance: f64) -> Result<MultiPolygon, GeometryError> {
        check_tolerance(tolerance)?;
        if tolerance == 0.0 {
            return Ok(self.clone());
        }
        build_multipolygon(self, tolerance)?.vouched()
    }

    /// [`MultiPolygon::simplify`], returning a result whose only defect is a
    /// self-intersecting part wrapped in [`Unvouched`] instead of refusing
    /// it. See [`Polygon::simplify_unvouched`] for the shape of the answer.
    ///
    /// # Errors
    ///
    /// As [`MultiPolygon::simplify`], minus the `Unvouched` refusal.
    pub fn simplify_unvouched(
        &self,
        tolerance: f64,
    ) -> Result<MaybeVouched<MultiPolygon>, GeometryError> {
        check_tolerance(tolerance)?;
        if tolerance == 0.0 {
            return Ok(MaybeVouched::Vouched(self.clone()));
        }
        build_multipolygon(self, tolerance)
    }
}

fn build_multipolygon(
    mp: &MultiPolygon,
    tolerance: f64,
) -> Result<MaybeVouched<MultiPolygon>, GeometryError> {
    let mut parts = Vec::with_capacity(mp.num_polygons());
    let mut first_defect: Option<UnvouchedReason> = None;
    for (part, p) in mp.polygons().iter().enumerate() {
        match build_polygon(p, tolerance)? {
            MaybeVouched::Vouched(q) => parts.push(q),
            MaybeVouched::Unvouched(u) => {
                let (q, reason) = u.into_parts();
                if first_defect.is_none() {
                    let (ring, edges) = match reason {
                        UnvouchedReason::SelfIntersecting { edges } => (0, edges),
                        UnvouchedReason::RingSelfIntersecting { ring, edges } => (ring, edges),
                        other => {
                            return Err(GeometryError::Unvouched {
                                reason: other,
                                inputs: vec![part],
                            });
                        }
                    };
                    first_defect =
                        Some(UnvouchedReason::PartSelfIntersecting { part, ring, edges });
                }
                parts.push(q);
            }
        }
    }
    // Part validity, exact. Run on the self-intersecting result too: two
    // parts that have come to overlap are not a multipolygon whatever else
    // is wrong with them.
    let built = MultiPolygon::new(parts)?;
    Ok(match first_defect {
        None => MaybeVouched::Vouched(built),
        Some(reason) => MaybeVouched::Unvouched(Unvouched::new(built, reason)),
    })
}

// ===========================================================================
// Lines
// ===========================================================================

/// RDP over the vertices of a line string: an open chain, or a ring when
/// the line is closed.
fn simplify_line_coords(pts: &[Coord2d], tolerance: f64) -> Vec<Coord2d> {
    let n = pts.len();
    if n >= 4 && pts[0] == pts[n - 1] {
        // Closed: the closing vertex is the anchor's own repeat, so the ring
        // is the vertices without it, and the result is closed again.
        let mut out = simplify_ring(&pts[..n - 1], tolerance);
        out.push(pts[0]);
        return out;
    }
    simplify_chain(pts, tolerance)
}

impl LineString {
    /// Ramer-Douglas-Peucker simplification of the line.
    ///
    /// The endpoints are kept. A closed line (first vertex equal to the last)
    /// is simplified as a ring, so its closing point is not privileged over
    /// any other vertex; the result is closed again at the same point. There
    /// is no vouch beyond the constructor's: a line string may cross itself,
    /// so simplicity is not a validity condition for it.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::InvalidTolerance`] for a negative or non-finite
    ///   tolerance.
    /// - [`GeometryError::ConsecutiveDuplicateVertices`] if a closed line
    ///   collapses to its anchor alone, so that the result would be the two
    ///   equal vertices `a, a`. An open line always keeps its two distinct
    ///   endpoints.
    pub fn simplify(&self, tolerance: f64) -> Result<LineString, GeometryError> {
        check_tolerance(tolerance)?;
        if tolerance == 0.0 {
            return Ok(self.clone());
        }
        LineString::new(CoordSeq::new(simplify_line_coords(
            self.coords().as_slice(),
            tolerance,
        )))
    }
}

impl MultiLineString {
    /// [`LineString::simplify`] applied to every component, in order.
    ///
    /// # Errors
    ///
    /// As [`LineString::simplify`], for the first component that fails.
    pub fn simplify(&self, tolerance: f64) -> Result<MultiLineString, GeometryError> {
        check_tolerance(tolerance)?;
        if tolerance == 0.0 {
            return Ok(self.clone());
        }
        self.line_strings()
            .iter()
            .map(|l| l.simplify(tolerance))
            .collect::<Result<Vec<_>, _>>()
            .map(MultiLineString::new)
    }
}

#[cfg(test)]
mod tests;
