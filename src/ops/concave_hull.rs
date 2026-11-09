// src/ops/concave_hull.rs

//! [`concave_hull`]: a polygon containing every input point, tighter than
//! the convex hull, built by eroding long-edged border triangles from the
//! Delaunay triangulation.
//!
//! The method is JTS's `ConcaveHull` (Martin Davis; JTS 1.19,
//! `org.locationtech.jts.algorithm.hull.ConcaveHull`, EPL/EDL):
//! starting from the Delaunay triangulation, whose union is the
//! convex hull, repeatedly remove the border triangle with the longest
//! boundary edge, while that edge is longer than the threshold and the
//! removal keeps the region whole.
//!
//! # The threshold is a selection rule; the claim is exact
//!
//! Which triangles are eroded is decided by comparing edge lengths in
//! `f64`. That comparison selects among valid results and is never part of
//! the returned claim. The claim, that the result is a valid polygon whose
//! vertices are input points and which contains every input point, is
//! certified exactly: the boundary is traced as one cycle or refused, the ring goes
//! through the validating polygon constructor, and every input point is
//! located against the result by exact point location. A certification
//! failure, which the argument below says cannot happen, is
//! [`GeometryError::Unvouched`] with [`UnvouchedReason::NotAHull`], never a
//! plausible wrong polygon.
//!
//! # The removal rule, and why it keeps the contract
//!
//! A border triangle is removable exactly when:
//!
//! 1. **it has exactly one boundary edge**, and
//! 2. **the vertex opposite that edge is not on the boundary.**
//!
//! The rule follows JTS's `ConcaveHull`, which requires a removable triangle to
//! have exactly two adjacent triangles and refuses a *connecting* triangle, one
//! whose vertex opposite the boundary edge also lies on the boundary, as this would
//! create more than one result polygon."
//!
//! The rule preserves the invariant that the region is a triangulated
//! simple polygon containing every input point:
//!
//! - **Start.** The full triangulation's union is the convex hull, whose
//!   boundary is one simple cycle, and every input point is a vertex.
//! - **Removing an allowed triangle keeps a simple polygon.** The triangle
//!   meets the boundary in exactly its one boundary edge `(a, b)`: its
//!   apex `v` is off the boundary, so neither `v` nor the interior edges
//!   `(a, v)`, `(v, b)` touch it. The boundary cycle replaces `(a, b)` by
//!   `(a, v), (v, b)`: still one cycle, and still simple because `v` was
//!   not on it. The region stays connected: removing a triangle whose
//!   intersection with the boundary is one closed edge is an elementary
//!   collapse of the disk through that free edge.
//! - **No input point is lost.** Every input point is a vertex, no vertex
//!   lies inside a triangle, and a removal keeps all three of its vertices
//!   on the region: `a`, `b` keep their other boundary edges, and `v`
//!   joins the boundary.
//! - **The forbidden removals are forbidden for cause.** A triangle with
//!   two boundary edges is an ear; removing it would strand its shared
//!   vertex strictly outside the region, breaking containment. A triangle
//!   whose apex is already on the boundary would pinch the cycle into a
//!   figure eight. A triangle with three boundary edges is the whole region.
//!
//! The argument is combinatorial. It reads which edges and vertices are on
//! the boundary, never a coordinate, so the `f64` edge lengths choose
//! *which* allowed removal happens, and can never make an unallowed one.
//!
//! # Determinism
//!
//! Erosion processes boundary edges longest first, by squared length in
//! `f64`, with ties broken by the triangulation's canonical triangle order.
//! The triangulation itself depends only on the point set, so the hull
//! does too: permutations of the input give byte-identical results.
//!
//! # The threshold forms
//!
//! [`ConcaveHullParameter::MaxEdgeLength`] erodes boundary edges strictly
//! longer than the given length. [`ConcaveHullParameter::MaxEdgeLengthRatio`]
//! is the scale-free form, with the target length
//! `shortest + ratio · (longest − shortest)`; 0 "produces a concave hull of
//! minimum area that is still connected", 1 "produces the convex hull". One
//! difference here is that the comparison is strict, so the threshold at exactly
//! the longest boundary edge keeps the convex hull.

use std::collections::{BinaryHeap, HashMap};

use crate::GeometryError;
use crate::ops::basic::lex;
use crate::ops::delaunay::{Triangulation, delaunay};
use crate::ops::hull::convex_hull;
use crate::predicates::Relatable;
use crate::prepared::PreparedPolygonal;
use crate::topology::label::Location;
use crate::types::{Coord2d, CoordSeq, Geometry, LinearRing, Polygon, UnvouchedReason};

/// The erosion threshold of [`concave_hull`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ConcaveHullParameter {
    /// Erode boundary edges strictly longer than this length (map units).
    MaxEdgeLength(f64),
    /// The target length is `shortest + ratio · (longest − shortest)` over
    /// the Delaunay triangulation's edge lengths. In `[0, 1]`; 0 erodes
    /// everything the removal rule allows, 1 keeps the convex hull.
    MaxEdgeLengthRatio(f64),
}

/// A certification failure: refused, never returned (module documentation).
struct Broken;

/// The concave hull of the distinct vertices of `g`.
///
/// A valid polygon whose vertices are input points, containing every input
/// point (certified exactly), eroded from the convex hull by the
/// threshold's selection rule (see the documentation of src/ops/concave_hull.rs).
/// The result is canonical: it depends only on the point set, with its ring started
/// at its lexicographically smallest vertex. Collinear boundary points stay
/// on the ring (they are vertices of the region's boundary), so the
/// ratio-1 result equals the convex hull as a region, with every boundary
/// point as a vertex.
///
/// Degenerate input (fewer than three distinct points, or all collinear)
/// follows [`convex_hull`]'s convention: an empty `MultiPoint`, a `Point`,
/// or the `LineString` of the two extremes.
///
/// # Errors
///
/// - [`GeometryError::InvalidTolerance`] for a non-finite or negative
///   length, or a ratio outside `[0, 1]`.
/// - [`GeometryError::TooManyEntries`] as [`delaunay`] refuses it.
/// - [`GeometryError::Unvouched`] with [`UnvouchedReason::NotAHull`] if the
///   certification re-check fails, which is not expected to occur.
pub fn concave_hull(
    g: &impl Relatable,
    parameter: ConcaveHullParameter,
) -> Result<Geometry, GeometryError> {
    match parameter {
        ConcaveHullParameter::MaxEdgeLength(l) if !(l.is_finite() && l >= 0.0) => {
            return Err(GeometryError::InvalidTolerance { tolerance: l });
        }
        ConcaveHullParameter::MaxEdgeLengthRatio(r)
            if !(r.is_finite() && (0.0..=1.0).contains(&r)) =>
        {
            return Err(GeometryError::InvalidTolerance { tolerance: r });
        }
        _ => {}
    }
    let t = delaunay(g)?;
    if t.triangles().is_empty() {
        // Degenerate by the same classification convex_hull uses, so the
        // two agree byte for byte here.
        return convex_hull(g);
    }
    let threshold_sq = threshold_sq(&t, parameter);
    erode(&t, threshold_sq).map_err(|Broken| GeometryError::Unvouched {
        reason: UnvouchedReason::NotAHull,
        inputs: Vec::new(),
    })
}

/// The squared erosion threshold. A selection rule evaluated in `f64`,
/// with squared comparisons.
#[expect(
    clippy::float_cmp,
    reason = "the ratio endpoints 0 and 1 are compared exactly"
)]
fn threshold_sq(t: &Triangulation, parameter: ConcaveHullParameter) -> f64 {
    match parameter {
        ConcaveHullParameter::MaxEdgeLength(l) => l * l,
        ConcaveHullParameter::MaxEdgeLengthRatio(r) => {
            if r == 0.0 {
                return 0.0;
            }
            let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
            for e in edge_list(t) {
                let len = sq_len(t.points(), e[0], e[1]).sqrt();
                min = min.min(len);
                max = max.max(len);
            }
            // Ratio 1 targets twice the longest edge to ensure all
            // edges are included; otherwise min + r·(max − min).
            // Computed exactly so an independent recomputation of the target
            // agrees bit for bit.
            let target = if r == 1.0 {
                2.0 * max
            } else {
                r * (max - min) + min
            };
            target * target
        }
    }
}

/// Every distinct edge of the triangulation, as index pairs.
fn edge_list(t: &Triangulation) -> Vec<[u32; 2]> {
    let mut es: Vec<[u32; 2]> = t
        .triangles()
        .iter()
        .flat_map(|tr| (0..3).map(|k| sorted_pair(tr[k], tr[(k + 1) % 3])))
        .collect();
    es.sort_unstable();
    es.dedup();
    es
}

fn sorted_pair(a: u32, b: u32) -> [u32; 2] {
    if a < b { [a, b] } else { [b, a] }
}

fn sq_len(pts: &[Coord2d], a: u32, b: u32) -> f64 {
    let (p, q) = (pts[a as usize], pts[b as usize]);
    let (dx, dy) = (q.x - p.x, q.y - p.y);
    dx.mul_add(dx, dy * dy)
}

/// The erosion state over the canonical triangulation.
struct Erosion<'t> {
    t: &'t Triangulation,
    /// `adj[i][k]`: the triangle sharing edge `k` of triangle `i` (the edge
    /// from vertex `k` to vertex `k + 1`), with its own slot for that edge.
    adj: Vec<[Option<(u32, u8)>; 3]>,
    removed: Vec<bool>,
    /// Boundary edges incident to each point. A point is on the boundary
    /// exactly when its count is positive.
    border_degree: Vec<u32>,
    /// Pending boundary edges: squared length (as monotone bits, so the
    /// heap pops the longest), then the canonical tie-break: the smallest
    /// triangle index, then the smallest edge slot.
    heap: BinaryHeap<(u64, core::cmp::Reverse<u32>, core::cmp::Reverse<u8>)>,
}

/// Runs the erosion and certifies the result (module documentation).
fn erode(t: &Triangulation, threshold_sq: f64) -> Result<Geometry, Broken> {
    let mut e = Erosion::new(t)?;
    while let Some((bits, core::cmp::Reverse(tri), core::cmp::Reverse(slot))) = e.heap.pop() {
        // Strictly longer than the threshold erodes; the heap pops longest
        // first, and later pushes only follow removals, so the first edge
        // at or below the threshold ends the erosion.
        if f64::from_bits(bits) <= threshold_sq {
            break;
        }
        e.try_remove(tri, slot);
    }
    let ring = e.trace()?;
    let polygon = Polygon::new(
        LinearRing::new(CoordSeq::new(ring)).map_err(|_| Broken)?,
        Vec::new(),
    )
    .map_err(|_| Broken)?;
    // The exact containment certification: every input point in the closed
    // result. The prepared form answers each location through its edge
    // index; the answers are exact.
    let prepared = PreparedPolygonal::new(&polygon).map_err(|_| Broken)?;
    for p in t.points() {
        if prepared.locate(*p) == Location::Exterior {
            return Err(Broken);
        }
    }
    Ok(Geometry::Polygon(polygon))
}

impl<'t> Erosion<'t> {
    fn new(t: &'t Triangulation) -> Result<Erosion<'t>, Broken> {
        let n = t.triangles().len();
        let mut open: HashMap<[u32; 2], (u32, u8)> = HashMap::new();
        let mut adj: Vec<[Option<(u32, u8)>; 3]> = vec![[None; 3]; n];
        for (i, tr) in t.triangles().iter().enumerate() {
            let ti = u32::try_from(i).map_err(|_| Broken)?;
            for k in 0..3u8 {
                let key = sorted_pair(tr[usize::from(k)], tr[usize::from(k + 1) % 3]);
                match open.remove(&key) {
                    Some((o, j)) => {
                        adj[i][usize::from(k)] = Some((o, j));
                        adj[o as usize][usize::from(j)] = Some((ti, k));
                    }
                    None => {
                        if open.insert(key, (ti, k)).is_some() {
                            return Err(Broken);
                        }
                    }
                }
            }
        }
        let mut border_degree = vec![0u32; t.points().len()];
        let mut heap = BinaryHeap::new();
        for (i, tr) in t.triangles().iter().enumerate() {
            for k in 0..3u8 {
                if adj[i][usize::from(k)].is_none() {
                    let (a, b) = (tr[usize::from(k)], tr[usize::from(k + 1) % 3]);
                    border_degree[a as usize] += 1;
                    border_degree[b as usize] += 1;
                    let bits = sq_len(t.points(), a, b).to_bits();
                    heap.push((
                        bits,
                        core::cmp::Reverse(u32::try_from(i).map_err(|_| Broken)?),
                        core::cmp::Reverse(k),
                    ));
                }
            }
        }
        Ok(Erosion {
            t,
            adj,
            removed: vec![false; n],
            border_degree,
            heap,
        })
    }

    /// Is edge `k` of living triangle `i` on the current boundary?
    fn is_border(&self, i: u32, k: u8) -> bool {
        match self.adj[i as usize][usize::from(k)] {
            None => true,
            Some((o, _)) => self.removed[o as usize],
        }
    }

    /// Applies the removal rule to the popped boundary edge (module
    /// documentation). A triangle that fails the rule is left; its entry is
    /// never retried, because a vertex never leaves the boundary and a
    /// boundary edge never stops being one, so a refusal is final.
    fn try_remove(&mut self, tri: u32, slot: u8) {
        if self.removed[tri as usize] {
            return;
        }
        let border_count = (0..3u8).filter(|&k| self.is_border(tri, k)).count();
        if border_count != 1 {
            return;
        }
        let tr = self.t.triangles()[tri as usize];
        let apex = tr[usize::from(slot + 2) % 3];
        if self.border_degree[apex as usize] > 0 {
            return;
        }
        self.removed[tri as usize] = true;
        for k in 0..3u8 {
            let (a, b) = (tr[usize::from(k)], tr[usize::from(k + 1) % 3]);
            if k == slot {
                // The old boundary edge leaves the boundary with its
                // triangle.
                self.border_degree[a as usize] -= 1;
                self.border_degree[b as usize] -= 1;
            } else if let Some((o, j)) = self.adj[tri as usize][usize::from(k)] {
                // An interior edge becomes a boundary edge of the living
                // neighbor.
                self.border_degree[a as usize] += 1;
                self.border_degree[b as usize] += 1;
                let bits = sq_len(self.t.points(), a, b).to_bits();
                self.heap
                    .push((bits, core::cmp::Reverse(o), core::cmp::Reverse(j)));
            }
        }
    }

    /// Traces the region's boundary as one simple ring, counterclockwise,
    /// starting at its lexicographically smallest vertex. Living triangles
    /// are counterclockwise and traverse their edges with the region on the
    /// left, so the directed boundary edges chain region-left. A vertex
    /// with two outgoing boundary edges (a pinch), or more than one cycle,
    /// is a broken invariant and is refused.
    fn trace(&self) -> Result<Vec<Coord2d>, Broken> {
        let mut next: HashMap<u32, u32> = HashMap::new();
        for (i, tr) in self.t.triangles().iter().enumerate() {
            let ti = u32::try_from(i).map_err(|_| Broken)?;
            if self.removed[i] {
                continue;
            }
            for k in 0..3u8 {
                if self.is_border(ti, k)
                    && next
                        .insert(tr[usize::from(k)], tr[usize::from(k + 1) % 3])
                        .is_some()
                {
                    return Err(Broken);
                }
            }
        }
        let pts = self.t.points();
        let start = next
            .keys()
            .copied()
            .min_by(|&a, &b| lex(pts[a as usize], pts[b as usize]))
            .ok_or(Broken)?;
        let mut ring = Vec::with_capacity(next.len());
        let mut cur = start;
        for _ in 0..next.len() {
            ring.push(pts[cur as usize]);
            cur = *next.get(&cur).ok_or(Broken)?;
        }
        if cur != start || ring.len() != next.len() {
            return Err(Broken);
        }
        Ok(ring)
    }
}
