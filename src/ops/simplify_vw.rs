// src/ops/simplify_vw.rs

//! Visvalingam–Whyatt simplification, vouched.
//!
//! # The method
//!
//! Visvalingam and Whyatt (*Line generalisation by repeated elimination of
//! points*, The Cartographic Journal 30(1), 1993) rank every vertex by its
//! **effective area**: the area of the triangle it forms with its two
//! current neighbors. The vertex with the smallest effective area is
//! removed, its two neighbors' areas are recomputed against their new
//! neighbors, and the elimination repeats. [`simplify_vw`] stops when the
//! smallest effective area is at or above the threshold.
//!
//! The elimination, exactly:
//!
//! 1. **Fixed vertices.** A line's two endpoints are never removed. A closed
//!    line (first vertex equal to the last) keeps at least one vertex
//!    between them, so it stays two distinct points; unlike the
//!    Douglas–Peucker `simplify`, which treats a closed line as a ring, a
//!    line's endpoints are fixed here whether or not they coincide. A ring
//!    keeps at least three vertices: once it has three, none is removed.
//! 2. **The effective area** of a vertex `b` with current neighbors `a`
//!    (before) and `c` (after) is `|(b − a) × (c − a)| / 2`, computed in
//!    `f64` exactly as written: `((b.x − a.x)·(c.y − a.y) − (b.y − a.y)·(c.x −
//!    a.x)).abs() · 0.5`. A non-finite value (products overflowing at
//!    coordinates above about `1e154`) counts as `+∞`, so the vertex is kept.
//!    A ring's vertices all have neighbors, cyclically.
//! 3. **A vertex whose two neighbors are the same point is not removable**
//!    while they are: removing it would leave two equal consecutive
//!    vertices, which no line or ring may have. (Its area is then 0, so
//!    without this rule it would be the first to go.) It becomes removable
//!    again if a neighbor changes.
//! 4. **Selection.** Among the removable vertices, the one with the smallest
//!    effective area is taken; ties go to the vertex earlier in the input,
//!    so the result is deterministic.
//! 5. **Stop or remove.** The selected vertex is removed if its area is
//!    below the threshold, or is exactly `0.0`; otherwise the elimination
//!    stops, since every other removable vertex has at least that area.
//! 6. **After a removal,** the two neighbors' areas are recomputed against
//!    their new neighbors.
//!
//! **Zero areas.** The paper's step 2 deletes every point of zero area
//! first. Here zero-area vertices go through the same queue, which takes
//! them first because 0 is the smallest area; a vertex made collinear by an
//! earlier removal is then removed too. So a threshold of 0 removes exactly
//! the vertices whose effective area computes as `0.0` in `f64`: collinear
//! ones, or ones within rounding of collinear.
//!
//! **The paper's monotone rule.** When a point is selected, "if its
//! calculated area is less than that of the last point to be eliminated,
//! use the latter's area instead. This ensures that the current point
//! cannot be eliminated without eliminating previously eliminated points"
//! (Visvalingam & Whyatt 1993, the algorithm's step 3). The rule fixes the
//! area *recorded* for the eliminated point, the rank the paper stores for
//! filtering at run time; the selection itself is on the calculated areas.
//! The internal routine records exactly that. For a threshold run it never
//! changes which vertices are removed: by induction, every recorded area so
//! far is below the threshold (or 0 at a threshold of 0), so the *recorded*
//! area `max(calculated, last)` is below the threshold exactly when the
//! calculated one is.
//!
//! # Selection rule, and the vouch
//!
//! **The effective areas are a selection rule, computed in `f64`.** They decide
//! which vertices are removed. Comparing triangle areas exactly needs exact
//! triangle-area comparisons that `arithmetica` does not yet offer. An area
//! within rounding of the threshold, or tied with another, may go either way.
//! That choice is never a correctness claim, because of the vouch.
//!
//! **The result is vouched, or refused.** Each line and ring is simplified
//! independently, as [`simplify`](crate::types::Polygon::simplify) does,
//! and the result is built through the checked constructors and then judged
//! by [`validity_issues`]: the same exact checks
//! [`validate`](super::validate()) makes. A result with any issue is
//! refused with [`UnvouchedReason::SimplifiedInvalid`]. A result that
//! cannot be built at all (a ring left with three collinear vertices) is
//! refused with its constructor's error. Visvalingam–Whyatt does **not**
//! preserve topology: rings can come to cross, holes to escape, parts to
//! overlap. For topology preserved by construction, use
//! [`simplify_preserve_topology`](super::simplify_preserve_topology), or
//! `coverage_simplify` for a coverage.
//!
//! # The threshold
//!
//! The threshold is an **area**, in the coordinates' units squared: the
//! quantity the method compares. It must be finite and at least 0.
//!
//! JTS's `VWSimplifier` (`org.locationtech.jts.simplify`) takes a distance
//! tolerance instead, and its documentation states that "this is converted
//! to an area tolerance by squaring it". To match a JTS call with distance
//! tolerance `t`, pass `t * t`.
//!
//! Every vertex of the result is a vertex of the input; nothing is
//! constructed.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use super::simplify::check_tolerance;
use super::simplify_topology::sealed::{self, Flat, Kept};
use super::validity_issues;
use crate::GeometryError;
use crate::predicates::Relatable;
use crate::types::{
    Coord2d, GeometryCollection, LineString, MaybeVouched, MultiLineString, MultiPolygon, Polygon,
    Unvouched, UnvouchedReason,
};

// ===========================================================================
// The public entry point
// ===========================================================================

/// The geometry types [`simplify_vw`] accepts: those with vertices to
/// remove. Sealed.
///
/// Implemented for [`LineString`], [`MultiLineString`], [`Polygon`],
/// [`MultiPolygon`] and [`GeometryCollection`]. A collection's points and
/// multipoints are kept unchanged.
pub trait Simplifiable: sealed::Sealed + Relatable + Clone {}

impl Simplifiable for LineString {}
impl Simplifiable for MultiLineString {}
impl Simplifiable for Polygon {}
impl Simplifiable for MultiPolygon {}
impl Simplifiable for GeometryCollection {}

/// Visvalingam–Whyatt simplification, vouched: repeatedly remove the vertex
/// whose effective area (the triangle it forms with its two current
/// neighbors) is smallest, while that area is below `threshold`. See the
/// [module documentation](self) for the elimination, the selection rule and
/// the vouch.
///
/// The result has the same type as `g`, and every one of its vertices is a
/// vertex of `g`. Each line and ring is simplified independently. Line
/// endpoints are kept; rings keep at least three vertices.
///
/// `threshold` is an area (the coordinates' units squared).
///
/// ```
/// use geometres::ops::simplify_vw;
/// use geometres::types::{Coord2d, CoordSeq, LineString};
///
/// let line = LineString::new(CoordSeq::new(vec![
///     Coord2d::new(0.0, 0.0),
///     Coord2d::new(5.0, 0.1), // effective area 0.5
///     Coord2d::new(10.0, 0.0),
///     Coord2d::new(10.0, 10.0),
/// ]))?;
/// let out = simplify_vw(&line, 1.0)?;
/// assert_eq!(out.num_points(), 3);
/// # Ok::<(), geometres::GeometryError>(())
/// ```
///
/// # Errors
///
/// - [`GeometryError::InvalidTolerance`] for a negative or non-finite
///   threshold.
/// - [`GeometryError::Unvouched`] with
///   [`UnvouchedReason::SimplifiedInvalid`] when the simplified geometry
///   has a validity issue (a ring that crosses or touches itself, a hole
///   outside its exterior, overlapping parts, and so on).
/// - A constructor's error when the simplified geometry cannot be built: a
///   ring left with three collinear vertices
///   ([`GeometryError::DegenerateRing`]), or an inverted one.
pub fn simplify_vw<G: Simplifiable>(g: &G, threshold: f64) -> Result<G, GeometryError> {
    check_tolerance(threshold)?;
    let out = simplify_vw_raw(g, threshold)?;
    match validity_issues(&out).first().map(|i| i.issue) {
        None => Ok(out),
        Some(issue) => MaybeVouched::Unvouched(Unvouched::new(
            out,
            UnvouchedReason::SimplifiedInvalid { issue },
        ))
        .vouched(),
    }
}

/// Flatten, eliminate on every chain, rebuild through the constructors. No
/// validity check beyond the constructors'.
pub(super) fn simplify_vw_raw<G: Simplifiable>(g: &G, threshold: f64) -> Result<G, GeometryError> {
    let mut flat = Flat::default();
    g.flatten(&mut flat);
    let kept = kept_chains(&flat, threshold);
    g.rebuild(&mut Kept::new(kept))
}

/// The kept vertices of every chain of `flat`, in order.
pub(super) fn kept_chains(flat: &Flat, threshold: f64) -> Vec<Vec<Coord2d>> {
    flat.chains
        .iter()
        .map(|ch| {
            let pts = &flat.verts[ch.start..ch.start + ch.len];
            let kind = if ch.closed {
                ChainKind::Ring
            } else {
                ChainKind::Line
            };
            let e = eliminate(pts, kind, threshold, &mut |_, _, _| true);
            pts.iter()
                .zip(&e.keep)
                .filter_map(|(&p, &k)| k.then_some(p))
                .collect()
        })
        .collect()
}

// ===========================================================================
// The elimination core (crate-private; `coverage_simplify` adds a veto)
// ===========================================================================

/// How a chain's ends behave.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ChainKind {
    /// An open or closed line: endpoints fixed.
    Line,
    /// A ring, stored without its closing vertex: cyclic, at least three kept.
    Ring,
}

/// What the elimination did to one chain.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Elimination {
    /// `keep[i]` for each input vertex.
    pub(crate) keep: Vec<bool>,
    /// The removed vertices in removal order, each with the area the paper
    /// records for it: its calculated area, raised to the previous recorded
    /// area when smaller (the monotone rule of the module documentation).
    pub(crate) removed: Vec<(usize, f64)>,
}

/// The effective area of `b` between `a` and `c`: `|(b − a) × (c − a)| / 2`
/// in `f64`, in exactly this order of operations (the Python reference
/// recomputes it bit for bit), with a non-finite value read as `+∞`.
pub(crate) fn effective_area(a: Coord2d, b: Coord2d, c: Coord2d) -> f64 {
    let cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    let area = cross.abs() * 0.5;
    if area.is_finite() {
        area
    } else {
        f64::INFINITY
    }
}

/// Visvalingam–Whyatt elimination of one chain (module documentation,
/// steps 1–6), with a veto.
///
/// Before each removal, `allow(i, before, after)` is asked, with the
/// vertex's index and its current neighbors' indices. If it answers
/// `false`, the vertex is fixed for the rest of the elimination and the next
/// one is tried. [`simplify_vw`] passes a veto that always allows.
///
/// A priority queue keyed by `(area, index)`, with stale entries skipped by
/// a per-vertex stamp: `O(n log n)`. It selects exactly what a scan over
/// every removable vertex would, in the same order.
pub(crate) fn eliminate(
    pts: &[Coord2d],
    kind: ChainKind,
    threshold: f64,
    allow: &mut dyn FnMut(usize, usize, usize) -> bool,
) -> Elimination {
    let n = pts.len();
    let mut keep = vec![true; n];
    let mut removed = Vec::new();
    let min_alive = match kind {
        ChainKind::Ring => 3,
        ChainKind::Line if n >= 2 && pts[0] == pts[n - 1] => 3,
        ChainKind::Line => 2,
    };
    if n <= min_alive {
        return Elimination { keep, removed };
    }
    let mut prev: Vec<usize> = (0..n).map(|i| if i == 0 { n - 1 } else { i - 1 }).collect();
    let mut next: Vec<usize> = (0..n).map(|i| if i + 1 == n { 0 } else { i + 1 }).collect();
    let mut fixed = vec![false; n];
    if kind == ChainKind::Line {
        fixed[0] = true;
        fixed[n - 1] = true;
    }
    let mut stamp = vec![0_u32; n];
    let mut heap: BinaryHeap<Reverse<(u64, usize, u32)>> = BinaryHeap::with_capacity(n);
    let entry = |i: usize, prev: &[usize], next: &[usize], stamp: &[u32]| {
        let (a, c) = (pts[prev[i]], pts[next[i]]);
        // Rule 3: equal neighbors make the vertex unremovable for now.
        (a != c).then(|| Reverse((effective_area(a, pts[i], c).to_bits(), i, stamp[i])))
    };
    for i in (0..n).filter(|&i| !fixed[i]) {
        if let Some(e) = entry(i, &prev, &next, &stamp) {
            heap.push(e);
        }
    }
    let mut alive = n;
    let mut last = 0.0_f64;
    while alive > min_alive {
        let Some(Reverse((bits, i, s))) = heap.pop() else {
            break;
        };
        if !keep[i] || fixed[i] || s != stamp[i] {
            continue; // stale
        }
        // Areas are non-negative, so their bit patterns order as the values.
        let area = f64::from_bits(bits);
        if !(area < threshold || area == 0.0) {
            break; // step 5: every other removable vertex has at least this area
        }
        let (b, a) = (prev[i], next[i]);
        if !allow(i, b, a) {
            fixed[i] = true;
            continue;
        }
        keep[i] = false;
        alive -= 1;
        next[b] = a;
        prev[a] = b;
        last = area.max(last);
        removed.push((i, last));
        for j in [b, a] {
            stamp[j] = stamp[j].wrapping_add(1);
            if !fixed[j]
                && let Some(e) = entry(j, &prev, &next, &stamp)
            {
                heap.push(e);
            }
        }
    }
    Elimination { keep, removed }
}
