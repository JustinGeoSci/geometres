// src/exact/kernel.rs

//! Exact ring-level geometric primitives, built directly on
//! `arithmetica::orient2d`.
//!
//! # Why this sits below `topology`
//!
//! `Polygon::new` must decide whether a hole lies inside its exterior, and
//! `LinearRing::is_simple` must decide whether two edges cross. Both are
//! questions about rings, and both need an exact sign test. The type layer
//! sits below `topology`, so the exact primitives it needs live here, and
//! `topology` imports them from here rather than keeping its own copies.
//!
//! Every sign decision is `orient2d` on stored coordinates, and every
//! ordering decision compares stored coordinates. No intersection point is
//! ever constructed.
//!
//! # Finiteness
//!
//! Every coordinate reaching this module comes from a validated `LinearRing`,
//! whose constructor runs `check_coords`, so all values are finite and
//! `orient2d` cannot return `Err`. The `unwrap_or(Sign::Zero)` in [`orient`]
//! is therefore unreachable. The crate follows the same convention wherever
//! the inputs are finite by construction.
//!
//! That argument holds only for callers inside this crate. The items other
//! crates need are public through the parent module, [`crate::exact`]: those
//! taking validated types are re-exported unchanged,
//! and those taking bare coordinates are wrapped in a finiteness check that
//! returns an error instead of reaching the fallback.

use arithmetica::{Sign, orient2d};

use crate::types::Coord2d;
use crate::types::LinearRing;

// ===========================================================================
// Point and segment primitives
// ===========================================================================

/// Exact orientation of `c` relative to the directed line `a → b`.
///
/// `Positive` is left of (counterclockwise from) `a → b`, `Negative` is
/// right, `Zero` is on the line. See the module note on why the `Err`
/// fallback cannot be reached.
#[inline]
pub(crate) fn orient(a: Coord2d, b: Coord2d, c: Coord2d) -> Sign {
    orient2d([a.x, a.y], [b.x, b.y], [c.x, c.y]).unwrap_or(Sign::Zero)
}

/// `true` if `p` lies on the closed segment `[a, b]`.
///
/// The bounding-box test is exact (it compares stored coordinates) and is
/// applied first only because it is cheaper than `orient2d`.
pub(crate) fn on_segment(p: Coord2d, a: Coord2d, b: Coord2d) -> bool {
    p.x >= a.x.min(b.x)
        && p.x <= a.x.max(b.x)
        && p.y >= a.y.min(b.y)
        && p.y <= a.y.max(b.y)
        && orient(a, b, p) == Sign::Zero
}

/// `true` if segments `ab` and `cd` cross at a single point interior to
/// both.
///
/// Each endpoint of one segment must lie strictly on opposite sides of the
/// other segment's line. A `Zero` from either side test means an endpoint is
/// on the other line (a T-junction, a shared endpoint, or collinearity), and
/// none of those is a proper crossing. Collinear overlap is therefore reported
/// as *not* crossing, the semantics `LinearRing::is_simple` documents and the
/// ring-placement tests below rely on. All four sign tests are exact.
pub(crate) fn segments_properly_cross(a: Coord2d, b: Coord2d, c: Coord2d, d: Coord2d) -> bool {
    let d1 = orient(c, d, a);
    let d2 = orient(c, d, b);
    if d1 == Sign::Zero || d2 == Sign::Zero || d1 == d2 {
        return false;
    }
    let d3 = orient(a, b, c);
    let d4 = orient(a, b, d);
    d3 != Sign::Zero && d4 != Sign::Zero && d3 != d4
}

// ===========================================================================
// Point-in-ring
// ===========================================================================

/// Where a point lies relative to a single ring, decided exactly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RingLocation {
    /// Strictly inside the region the ring bounds.
    Inside,
    /// On a ring edge, endpoints included.
    OnBoundary,
    /// Strictly outside.
    Outside,
}

/// Crossing-number test of a point known *not* to be on `ring`.
///
/// Casts a ray toward +x. An edge counts when it straddles the horizontal
/// line through `p` under the half-open rule `(a.y > p.y) != (b.y > p.y)`,
/// which counts a ray passing exactly through a vertex exactly once. For a
/// straddling edge, the sign of `orient(a, b, p)` decides, exactly, whether
/// the edge is to the right of `p`: `p` left of an upward edge, or right of a
/// downward edge, means the edge crosses the ray.
///
/// Winding-agnostic: the parity of crossings does not depend on the ring's
/// direction, so this is correct for exteriors and holes alike.
pub(crate) fn ring_contains_off_boundary(p: Coord2d, ring: &LinearRing) -> bool {
    let mut inside = false;
    for (a, b) in ring.edges() {
        if edge_crosses_ray(a, b, p) {
            inside = !inside;
        }
    }
    inside
}

/// Whether the edge `a → b` crosses the ray from `p` toward +x, under the
/// rule [`ring_contains_off_boundary`] states: the edge straddles the
/// horizontal line through `p` (half-open, `(a.y > p.y) != (b.y > p.y)`) and
/// lies to the right of `p` there, decided by the sign of `orient(a, b, p)`.
///
/// One function, so that the whole-ring test above and the edge-filtered
/// test of a prepared mask (`topology::areal::IndexedAreal`) toggle on
/// exactly the same edges. An edge whose closed box does not meet the ray's
/// box (a `y` range excluding `p.y`, or both endpoints left of `p`) never
/// crosses: the straddle test fails, or `p` is to the right of the edge's
/// line within its `y` range.
pub(crate) fn edge_crosses_ray(a: Coord2d, b: Coord2d, p: Coord2d) -> bool {
    if (a.y > p.y) == (b.y > p.y) {
        return false;
    }
    let s = orient(a, b, p);
    if b.y > a.y {
        s == Sign::Positive
    } else {
        s == Sign::Negative
    }
}

/// Exact location of `p` relative to `ring`.
pub(crate) fn locate_in_ring(p: Coord2d, ring: &LinearRing) -> RingLocation {
    if ring.edges().any(|(a, b)| on_segment(p, a, b)) {
        return RingLocation::OnBoundary;
    }
    if ring_contains_off_boundary(p, ring) {
        RingLocation::Inside
    } else {
        RingLocation::Outside
    }
}

// ===========================================================================
// Ring-against-ring placement
// ===========================================================================

/// Where one simple ring lies relative to another, decided exactly.
///
/// For two simple closed curves whose boundaries meet in at most finitely
/// many points there are exactly three possibilities: one is inside the
/// other, they are separate, or one passes from inside the other to outside
/// it. A fourth case is separated out because no valid placement allows it:
/// the boundaries share a stretch of positive length.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RingPlacement {
    /// `inner`'s interior is entirely within the region `outer` bounds.
    /// `inner` may touch `outer`'s boundary at isolated points.
    Inside,
    /// `inner`'s interior is entirely outside the region `outer` bounds
    /// (which includes `inner` surrounding `outer`). It may touch at
    /// isolated points.
    Outside,
    /// `inner` is partly inside and partly outside: an edge properly crosses
    /// `outer`'s boundary, or `inner` passes through it at a vertex.
    Crossing,
    /// The boundaries share a stretch of positive length, as two identical
    /// rings do. OGC allows rings to touch at points only, so this is
    /// never a valid hole placement, hole pair, or pair of disjoint parts.
    SharedEdge,
}

/// Which side of `outer` one sub-edge of `inner` lies on.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SubEdgeSide {
    Inside,
    Outside,
    Along,
}

/// The side of `outer` the open segment `s → e` lies on, given that no
/// boundary point of `outer` lies strictly inside it except, possibly, along
/// its whole length.
fn sub_edge_side(s: Coord2d, e: Coord2d, outer: &LinearRing) -> SubEdgeSide {
    for p in [s, e] {
        match locate_in_ring(p, outer) {
            RingLocation::Inside => return SubEdgeSide::Inside,
            RingLocation::Outside => return SubEdgeSide::Outside,
            RingLocation::OnBoundary => {}
        }
    }
    if edge_along_boundary(s, e, outer) {
        return SubEdgeSide::Along;
    }
    if chord_departs_inward(s, e, outer) {
        SubEdgeSide::Inside
    } else {
        SubEdgeSide::Outside
    }
}

/// `outer`'s vertices lying strictly inside the segment `a → b`, in order
/// from `a`. Every comparison is of stored coordinates.
fn cuts_along(a: Coord2d, b: Coord2d, outer: &LinearRing) -> Vec<Coord2d> {
    let mut cuts: Vec<Coord2d> = outer
        .coords()
        .as_slice()
        .iter()
        .copied()
        .filter(|&v| v != a && v != b && on_segment(v, a, b))
        .collect();
    // The points are collinear with `a → b`, so one coordinate orders them:
    // `x` unless the segment is vertical. Ascending or descending to follow
    // the segment's direction.
    #[allow(clippy::float_cmp)] // stored coordinates, compared exactly
    let key = |c: &Coord2d| if a.x == b.x { c.y } else { c.x };
    let forward = key(&b) > key(&a);
    cuts.sort_by(|p, q| {
        let o = key(p).total_cmp(&key(q));
        if forward { o } else { o.reverse() }
    });
    cuts.dedup();
    cuts
}

/// Classify `inner` against `outer`, exactly, for every way two simple rings
/// can meet.
///
/// # Method
///
/// A proper crossing decides `Crossing` outright. Otherwise every edge of
/// `inner` is cut at each vertex of `outer` lying strictly inside it, and
/// each resulting sub-edge is classified on its own:
///
/// - if either end is off `outer`'s boundary, that end's location is the
///   whole sub-edge's. No boundary point lies strictly inside the
///   sub-edge, since such a point would be a vertex of `outer` (a cut), a
///   proper crossing (already excluded), or the start of a shared stretch,
///   which begins at a vertex of `outer` or of `inner`;
/// - if both ends are on the boundary and on one edge of `outer`, the
///   sub-edge lies *along* that edge;
/// - otherwise it is a chord, and the side it departs into from its start
///   is decided by an exact wedge test there (`chord_departs_inward`).
///
/// The answer is `SharedEdge` if any sub-edge is **along** the boundary, and
/// otherwise `Inside`, `Outside` or `Crossing` by whether the sub-edges are
/// all inside, all outside, or mixed.
///
/// # Why one vertex does not decide
///
/// Without a proper crossing, it is tempting to locate one vertex of `inner`
/// off `outer`'s boundary and take that side for the whole ring. The
/// argument fails: `inner` can pass from one side to the other *through a
/// vertex*, either its own vertex on `outer`'s boundary or `outer`'s vertex
/// on its edge, and neither is a proper crossing. A hole leaving an L-shaped
/// exterior through two boundary vertices is one such case. Classifying
/// every sub-edge also detects a shared stretch.
///
/// Every decision is `orient` on stored coordinates or a comparison of
/// stored coordinates. No point is constructed. `O(n·m)`.
pub fn place_ring(inner: &LinearRing, outer: &LinearRing) -> RingPlacement {
    if rings_properly_cross(inner, outer) {
        return RingPlacement::Crossing;
    }
    let (mut inside, mut outside, mut along) = (false, false, false);
    for (a, b) in inner.edges() {
        let mut s = a;
        for e in cuts_along(a, b, outer)
            .into_iter()
            .chain(core::iter::once(b))
        {
            match sub_edge_side(s, e, outer) {
                SubEdgeSide::Inside => inside = true,
                SubEdgeSide::Outside => outside = true,
                SubEdgeSide::Along => along = true,
            }
            s = e;
        }
    }
    if along {
        RingPlacement::SharedEdge
    } else if inside && outside {
        RingPlacement::Crossing
    } else if inside {
        RingPlacement::Inside
    } else {
        RingPlacement::Outside
    }
}

/// The first pair of edges of a closed ring that touch anywhere other than at
/// the vertex two consecutive edges share, by edge index; `None` if the ring
/// is simple in the OGC sense.
///
/// # Why proper crossing is not enough
///
/// A ring can fail to be simple without any two edges crossing: it can pass
/// through one point twice (a pinch), a vertex can sit on a non-adjacent
/// edge (a T-touch), or two consecutive edges can fold back over each other
/// (a spike). [`segments_properly_cross`] reports none of these, by design,
/// for the placement tests that rely on its semantics. A simplicity test
/// built on proper crossings alone would accept a union whose hole touches
/// its exterior at a point, written as one self-touching ring. Rounding can
/// produce all three contacts from an exact ring that had none: two distinct
/// exact vertices rounding to one coordinate is a pinch, and a vertex
/// rounding onto a nearby edge is a T-touch. The vouches therefore use this
/// test.
///
/// # Exactness
///
/// Every decision is [`orient`] on stored coordinates or a comparison of
/// stored coordinates: [`segments_properly_cross`] for crossings,
/// [`on_segment`] on the closed segment for every touching configuration.
/// No intersection point is constructed. The bounding-box prefilter compares
/// stored coordinates and is exact; it only skips pairs that cannot touch.
///
/// `O(n²)` in the vertex count, as the proper-crossing test it replaces.
pub(crate) fn first_self_contact(v: &[Coord2d]) -> Option<(usize, usize)> {
    if v.len() < 3 {
        return Some((0, 0));
    }
    self_contacts(v, true).first().copied()
}

/// Every pair of edges of a closed ring that touch anywhere other than at
/// the vertex two consecutive edges share, by edge index, in order of the
/// first edge; the enumerating form of [`first_self_contact`], for
/// diagnostics that report every issue (`ops::validate`).
///
/// With `first_only` the scan stops at the first contact, which is what the
/// vouches need. The tests are the same exact ones in either mode. A ring
/// with fewer than three vertices has no edges to compare and yields nothing;
/// callers report that condition themselves.
pub(crate) fn self_contacts(v: &[Coord2d], first_only: bool) -> Vec<(usize, usize)> {
    if v.len() > SWEEP_MIN_VERTICES {
        return self_contacts_sweep(v, first_only);
    }
    self_contacts_scan(v, first_only)
}

/// Vertices above which [`self_contacts`] finds its candidate pairs by a
/// sweep instead of the double loop. A performance choice only: both give
/// the same pairs in the same order.
const SWEEP_MIN_VERTICES: usize = 64;

/// [`self_contacts`] by the double loop: every pair, box prefilter first.
fn self_contacts_scan(v: &[Coord2d], first_only: bool) -> Vec<(usize, usize)> {
    let n = v.len();
    let mut found = Vec::new();
    if n < 3 {
        return found;
    }
    let edge = |i: usize| (v[i], v[(i + 1) % n]);
    let boxes_meet = |(a, b): (Coord2d, Coord2d), (c, d): (Coord2d, Coord2d)| {
        a.x.min(b.x) <= c.x.max(d.x)
            && c.x.min(d.x) <= a.x.max(b.x)
            && a.y.min(b.y) <= c.y.max(d.y)
            && c.y.min(d.y) <= a.y.max(b.y)
    };

    for i in 0..n {
        let (a, b) = edge(i);

        // Consecutive pair (i, i+1), sharing b: a fold-back puts the far end
        // of one edge on the other. A collinear continuation does not.
        let (_, c) = edge((i + 1) % n);
        if on_segment(c, a, b) || on_segment(a, b, c) {
            found.push((i, (i + 1) % n));
            if first_only {
                return found;
            }
        }

        // Non-consecutive pairs: any common point at all is a contact.
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue; // consecutive through the closing edge, handled above
            }
            let (c, d) = edge(j);
            if !boxes_meet((a, b), (c, d)) {
                continue;
            }
            if segments_properly_cross(a, b, c, d)
                || on_segment(a, c, d)
                || on_segment(b, c, d)
                || on_segment(c, a, b)
                || on_segment(d, a, b)
            {
                found.push((i, j));
                if first_only {
                    return found;
                }
            }
        }
    }
    found
}

/// [`self_contacts`] with the non-consecutive candidates found by a sweep
/// over the edges sorted by their least `x`.
///
/// The double loop tests a non-consecutive pair only when the two edges'
/// closed boxes meet, so the pairs it can report are among the pairs whose
/// boxes meet. The sweep enumerates exactly those pairs: an edge stays
/// active while its greatest `x` is at least the least `x` of the edge
/// being added (closed comparison, so boxes touching at one `x` meet), and
/// each active edge is then checked for overlap in `y`. Every candidate is
/// tested by the same exact contact test, the consecutive pairs by the same
/// fold-back test, and the contacts are put in the double loop's order:
/// by first edge, the consecutive pair before the others, then by second
/// edge. So the result is the double loop's, pair for pair; with
/// `first_only` it is its first pair.
fn self_contacts_sweep(v: &[Coord2d], first_only: bool) -> Vec<(usize, usize)> {
    let n = v.len();
    let edge = |i: usize| (v[i], v[(i + 1) % n]);
    let lo: Vec<f64> = (0..n).map(|i| v[i].x.min(v[(i + 1) % n].x)).collect();
    let hi: Vec<f64> = (0..n).map(|i| v[i].x.max(v[(i + 1) % n].x)).collect();
    let lo_x = |i: usize| lo[i];
    let hi_x = |i: usize| hi[i];
    // (first edge, 0 for the consecutive pair or 1 for the others, second).
    let mut found: Vec<(usize, u8, usize)> = Vec::new();
    for i in 0..n {
        let (a, b) = edge(i);
        let (_, c) = edge((i + 1) % n);
        if on_segment(c, a, b) || on_segment(a, b, c) {
            found.push((i, 0, (i + 1) % n));
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| lo_x(i).total_cmp(&lo_x(j)));
    let mut active: Vec<usize> = Vec::new();
    for &k in &order {
        let x = lo_x(k);
        active.retain(|&i| hi_x(i) >= x);
        let (a, b) = edge(k);
        for &m in &active {
            let (i, j) = (k.min(m), k.max(m));
            if j < i + 2 || (i == 0 && j == n - 1) {
                continue; // consecutive: the fold-back test above
            }
            let (c, d) = edge(m);
            if a.y.min(b.y) > c.y.max(d.y) || c.y.min(d.y) > a.y.max(b.y) {
                continue;
            }
            let ((p, q), (r, t)) = (edge(i), edge(j));
            if segments_properly_cross(p, q, r, t)
                || on_segment(p, r, t)
                || on_segment(q, r, t)
                || on_segment(r, p, q)
                || on_segment(t, p, q)
            {
                found.push((i, 1, j));
            }
        }
        active.push(k);
    }
    found.sort_unstable();
    if first_only {
        found.truncate(1);
    }
    found.into_iter().map(|(i, _, j)| (i, j)).collect()
}

/// `true` if any edge of `a` properly crosses any edge of `b`. `O(n·m)`.
pub(crate) fn rings_properly_cross(a: &LinearRing, b: &LinearRing) -> bool {
    let b_edges: Vec<(Coord2d, Coord2d)> = b.edges().collect();
    a.edges().any(|(p, q)| {
        b_edges
            .iter()
            .any(|&(r, s)| segments_properly_cross(p, q, r, s))
    })
}

/// `true` if the regions bounded by `a` and `b` have disjoint interiors and
/// their boundaries share no stretch.
///
/// Both directions are required: `a` outside `b` does not imply `b` outside
/// `a`. If `b` lies *inside* `a`, then `a` is outside `b`, yet the interiors
/// overlap. Two identical rings, or two rings sharing an edge
/// stretch, are `SharedEdge` and are reported as *not* disjoint: OGC allows two
/// holes of a polygon to touch at points only. Rings touching at isolated
/// points are disjoint.
pub(crate) fn rings_interior_disjoint(a: &LinearRing, b: &LinearRing) -> bool {
    place_ring(a, b) == RingPlacement::Outside && place_ring(b, a) == RingPlacement::Outside
}

// ===========================================================================
// Between-ring diagnostics (`ops::validate`)
// ===========================================================================

/// The sub-edges of `inner` against `outer`, in ring order: each edge of
/// `inner` (by index) cut at the vertices of `outer` strictly inside it, with
/// the side of `outer` each piece lies on. This is the classification
/// [`place_ring`] makes, kept rather than folded into one verdict.
///
/// Only for the diagnostics below, which run on rings already known to be
/// misplaced; `place_ring` itself stays allocation-light.
fn sub_edges(
    inner: &LinearRing,
    outer: &LinearRing,
) -> Vec<(usize, Coord2d, Coord2d, SubEdgeSide)> {
    let mut out = Vec::new();
    for (i, (a, b)) in inner.edges().enumerate() {
        let mut s = a;
        for e in cuts_along(a, b, outer)
            .into_iter()
            .chain(core::iter::once(b))
        {
            out.push((i, s, e, sub_edge_side(s, e, outer)));
            s = e;
        }
    }
    out
}

/// The first edge of `ring`, by index, whose closed segment holds `p`.
fn edge_holding(p: Coord2d, ring: &LinearRing) -> Option<usize> {
    ring.edges().position(|(a, b)| on_segment(p, a, b))
}

/// A pair of edges, `(edge of inner, edge of outer)`, that witnesses why
/// `inner` is not cleanly inside or outside `outer`. For rings that
/// [`place_ring`] classifies as `Crossing` or `SharedEdge`; exact, and
/// nothing is constructed.
///
/// In order of preference:
///
/// 1. **A shared stretch:** the first sub-edge of `inner` lying along
///    `outer`, with the edge of `outer` it lies on.
/// 2. **A proper crossing:** the first pair of edges that properly cross.
/// 3. **A crossing through a point:** where `inner`'s sub-edges change from
///    inside to outside (or back), that change happens at a point of
///    `outer`'s boundary, since no sub-edge crosses it. The edge of `inner`
///    leaving that point, and the first edge of `outer` holding it.
/// 4. Otherwise, the first pair of edges with any common point. Only rings
///    that touch without crossing reach this case, and callers do not pass
///    them.
///    `None` if the rings have no common point at all.
pub(crate) fn contact_witness(inner: &LinearRing, outer: &LinearRing) -> Option<(usize, usize)> {
    let pieces = sub_edges(inner, outer);
    for &(i, s, e, side) in &pieces {
        if side == SubEdgeSide::Along
            && let Some(j) = outer
                .edges()
                .position(|(p, q)| on_segment(s, p, q) && on_segment(e, p, q))
        {
            return Some((i, j));
        }
    }
    let outer_edges: Vec<(Coord2d, Coord2d)> = outer.edges().collect();
    for (i, (a, b)) in inner.edges().enumerate() {
        if let Some(j) = outer_edges
            .iter()
            .position(|&(c, d)| segments_properly_cross(a, b, c, d))
        {
            return Some((i, j));
        }
    }
    let n = pieces.len();
    for k in 0..n {
        let prev = pieces[(k + n - 1) % n].3;
        let (i, s, _, side) = pieces[k];
        if prev != side
            && let Some(j) = edge_holding(s, outer)
        {
            return Some((i, j));
        }
    }
    for (i, (a, b)) in inner.edges().enumerate() {
        if let Some(j) = outer_edges.iter().position(|&(c, d)| {
            segments_properly_cross(a, b, c, d)
                || on_segment(a, c, d)
                || on_segment(b, c, d)
                || on_segment(c, a, b)
                || on_segment(d, a, b)
        }) {
            return Some((i, j));
        }
    }
    None
}

/// `true` if `a` and `b` bound the same point set: every piece of each lies
/// along the other. Identical rings, rings listed from a different start or
/// in the other direction, and rings differing only by vertices collinear
/// within an edge all qualify. Exact.
pub(crate) fn same_point_set(a: &LinearRing, b: &LinearRing) -> bool {
    sub_edges(a, b)
        .iter()
        .all(|&(_, _, _, side)| side == SubEdgeSide::Along)
        && sub_edges(b, a)
            .iter()
            .all(|&(_, _, _, side)| side == SubEdgeSide::Along)
}

/// The points where two rings meet, for rings that neither cross nor share
/// a stretch: each is a vertex of `a` on `b`'s boundary or a vertex of `b`
/// on `a`'s, in that order of discovery, without repeats.
///
/// That these are all the common points is exact geometry, not a sampling:
/// two segments meeting at a single point interior to both, without being
/// collinear, cross properly, which the precondition excludes; and collinear
/// segments meeting in one point meet at an endpoint. So every common point
/// is a vertex of one ring or the other, an input coordinate. Points are
/// compared by value (`==`), so `-0.0` and `0.0` are one point.
pub(crate) fn ring_touch_points(a: &LinearRing, b: &LinearRing) -> Vec<Coord2d> {
    let mut points: Vec<Coord2d> = Vec::new();
    let candidates = a
        .coords()
        .as_slice()
        .iter()
        .filter(|&&v| locate_in_ring(v, b) == RingLocation::OnBoundary)
        .chain(
            b.coords()
                .as_slice()
                .iter()
                .filter(|&&v| locate_in_ring(v, a) == RingLocation::OnBoundary),
        );
    for &p in candidates {
        if !points.contains(&p) {
            points.push(p);
        }
    }
    points
}

/// `true` if segment `ab` lies along a single edge of `ring` (both endpoints
/// on the same edge).
fn edge_along_boundary(a: Coord2d, b: Coord2d, ring: &LinearRing) -> bool {
    ring.edges()
        .any(|(p, q)| on_segment(a, p, q) && on_segment(b, p, q))
}

/// For `a` on `ring`'s boundary and `b` not along the same edge: does the
/// direction `a → b` depart into `ring`'s interior?
///
/// The interior side is given by the ring's orientation: counterclockwise
/// puts it on the left of every edge, clockwise on the right. At a vertex of
/// the ring the departure is inward if it falls *within* the interior angle,
/// which for a convex corner is "left of both incident edges" and for a
/// reflex corner is "left of either" (mirrored for clockwise).
fn chord_departs_inward(a: Coord2d, b: Coord2d, ring: &LinearRing) -> bool {
    let ccw = ring.signed_area() > 0.0;
    let inward = |s: Sign| -> bool {
        match s {
            Sign::Positive => ccw,
            Sign::Negative => !ccw,
            Sign::Zero => false,
        }
    };
    let pts = ring.coords().as_slice();
    let n = pts.len();

    // Is `a` a vertex of the ring?
    if let Some(k) = pts.iter().position(|&v| v == a) {
        let p = pts[(k + n - 1) % n];
        let q = pts[(k + 1) % n];
        let convex = inward(orient(p, a, q));
        let left_of_in = inward(orient(p, a, b));
        let left_of_out = inward(orient(a, q, b));
        return if convex {
            left_of_in && left_of_out
        } else {
            left_of_in || left_of_out
        };
    }
    // Otherwise `a` lies in the interior of exactly one edge.
    for (p, q) in ring.edges() {
        if on_segment(a, p, q) {
            return inward(orient(p, q, b));
        }
    }
    false
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used)]
// test-lint-posture: see the workspace lint comment in Cargo.toml
// `unreadable_literal` is allowed here as in the integration suites: the
// regression fixture's coordinates are transcribed exactly as proptest
// reported them, so they can be checked against the shrink output.
#[allow(
    clippy::float_cmp,
    clippy::unreadable_literal,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::range_plus_one,
    clippy::no_effect_underscore_binding
)]
mod tests {
    use super::*;
    use crate::types::CoordSeq;

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    fn ring(pts: &[(f64, f64)]) -> LinearRing {
        LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
    }

    fn unit_square() -> LinearRing {
        ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    }

    // --- segments ---

    #[test]
    fn proper_crossing_is_detected_and_touching_is_not() {
        assert!(segments_properly_cross(
            c(0.0, 0.0),
            c(4.0, 4.0),
            c(0.0, 4.0),
            c(4.0, 0.0)
        ));
        // T-junction: endpoint on the other segment's interior.
        assert!(!segments_properly_cross(
            c(0.0, 0.0),
            c(4.0, 0.0),
            c(2.0, 0.0),
            c(2.0, 3.0)
        ));
        // Shared endpoint.
        assert!(!segments_properly_cross(
            c(0.0, 0.0),
            c(3.0, 3.0),
            c(3.0, 3.0),
            c(6.0, 0.0)
        ));
        // Collinear overlap, documented as not a crossing.
        assert!(!segments_properly_cross(
            c(0.0, 0.0),
            c(4.0, 0.0),
            c(2.0, 0.0),
            c(6.0, 0.0)
        ));
        // Parallel.
        assert!(!segments_properly_cross(
            c(0.0, 0.0),
            c(4.0, 0.0),
            c(0.0, 1.0),
            c(4.0, 1.0)
        ));
    }

    // --- point in ring ---

    #[test]
    fn locate_in_ring_distinguishes_all_three() {
        let r = unit_square();
        assert_eq!(locate_in_ring(c(5.0, 5.0), &r), RingLocation::Inside);
        assert_eq!(locate_in_ring(c(15.0, 5.0), &r), RingLocation::Outside);
        assert_eq!(locate_in_ring(c(10.0, 5.0), &r), RingLocation::OnBoundary);
        assert_eq!(
            locate_in_ring(c(0.0, 0.0), &r),
            RingLocation::OnBoundary,
            "vertex"
        );
        // Ray through a vertex must count once, not twice.
        assert_eq!(
            locate_in_ring(c(5.0, 10.0 - 1e-9), &r),
            RingLocation::Inside
        );
        assert_eq!(locate_in_ring(c(-1.0, 10.0), &r), RingLocation::Outside);
    }

    #[test]
    fn locate_is_winding_agnostic() {
        let cw = ring(&[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)]);
        assert_eq!(locate_in_ring(c(5.0, 5.0), &cw), RingLocation::Inside);
        assert_eq!(locate_in_ring(c(15.0, 5.0), &cw), RingLocation::Outside);
    }

    // --- ring placement ---

    #[test]
    fn place_ring_covers_the_four_pictures() {
        let outer = unit_square();
        let inside = ring(&[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)]);
        let disjoint = ring(&[(20.0, 20.0), (24.0, 20.0), (24.0, 24.0), (20.0, 24.0)]);
        let overlapping = ring(&[(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0)]);
        let containing = ring(&[(-5.0, -5.0), (15.0, -5.0), (15.0, 15.0), (-5.0, 15.0)]);

        assert_eq!(place_ring(&inside, &outer), RingPlacement::Inside);
        assert_eq!(place_ring(&disjoint, &outer), RingPlacement::Outside);
        assert_eq!(place_ring(&overlapping, &outer), RingPlacement::Crossing);
        assert_eq!(
            place_ring(&containing, &outer),
            RingPlacement::Outside,
            "a ring that surrounds `outer` has every vertex outside it"
        );
        assert_eq!(
            place_ring(&outer, &containing),
            RingPlacement::Inside,
            "and from the other side, `outer` is inside it"
        );
    }

    #[test]
    fn touching_at_a_vertex_is_still_inside() {
        let outer = unit_square();
        // Hole touches the exterior at (0,0) and is otherwise inside.
        let touching = ring(&[(0.0, 0.0), (3.0, 1.0), (1.0, 3.0)]);
        assert_eq!(place_ring(&touching, &outer), RingPlacement::Inside);
    }

    #[test]
    fn identical_rings_share_their_whole_boundary_and_are_not_disjoint() {
        let a = unit_square();
        let b = unit_square();
        assert_eq!(place_ring(&a, &b), RingPlacement::SharedEdge);
        assert!(!rings_interior_disjoint(&a, &b));
    }

    #[test]
    fn disjointness_requires_both_directions() {
        let big = unit_square();
        let small = ring(&[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)]);
        // Every vertex of `big` is outside `small`, but `small` is inside
        // `big`, so they are not disjoint. One direction alone gets this
        // wrong.
        assert_eq!(place_ring(&big, &small), RingPlacement::Outside);
        assert!(!rings_interior_disjoint(&big, &small));

        let far = ring(&[(20.0, 20.0), (24.0, 20.0), (24.0, 24.0), (20.0, 24.0)]);
        assert!(rings_interior_disjoint(&big, &far));
    }

    #[test]
    fn rings_touching_at_one_vertex_are_disjoint() {
        let a = ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]);
        let b = ring(&[(4.0, 4.0), (8.0, 4.0), (8.0, 8.0), (4.0, 8.0)]);
        assert!(rings_interior_disjoint(&a, &b));
    }

    /// The hole from the failing `offset_properties` case 2: twenty times the
    /// exterior's area and straddling its boundary.
    #[test]
    fn regression_offset_properties_case_2_hole_crosses_exterior() {
        let ext = ring(&[
            (-4.003133454020272, 16.42962109085804),
            (-9.590873621482961, -8.093162064591926),
            (-15.530846932663101, -34.6057195833789),
        ]);
        let hole = ring(&[
            (-3.1430959461392274e-17, -0.17110217402152594),
            (-5.647170630321792, 6.915789436662523e-16),
            (3.41673958700421e-16, 5.579959200290364),
            (2.949389742165336, 0.0),
        ]);
        assert_ne!(place_ring(&hole, &ext), RingPlacement::Inside);
    }

    // ==================================================================
    // Rings touching at vertices: chords, and rings passing through a
    // vertex. A placement decided by one vertex answers `Boundary` for the
    // first three cases and the wrong side for the passing-through ones.
    // ==================================================================

    /// A diamond inscribed in a square hole, touching it at every vertex:
    /// every edge is a chord, and every chord departs inward.
    #[test]
    fn an_inscribed_diamond_is_inside_the_square() {
        let square_cw = ring(&[(2.0, 2.0), (2.0, 6.0), (6.0, 6.0), (6.0, 2.0)]);
        let diamond_ccw = ring(&[(6.0, 4.0), (4.0, 6.0), (2.0, 4.0), (4.0, 2.0)]);
        assert_eq!(place_ring(&diamond_ccw, &square_cw), RingPlacement::Inside);
        // The square's vertices are strictly outside the diamond.
        assert_eq!(place_ring(&square_cw, &diamond_ccw), RingPlacement::Outside);
        // Orientation of `outer` does not change the answer.
        let square_ccw = ring(&[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]);
        assert_eq!(place_ring(&diamond_ccw, &square_ccw), RingPlacement::Inside);
    }

    /// Rings with every vertex on a concave ring, around its reflex vertex.
    ///
    /// The second case is the one a test of the first chord alone gets
    /// wrong: that chord departs inward, so the triangle appears to lie
    /// *within* the concave ring. Its top edge `(6,6) → (2,6)` spans the
    /// notch, and the point `(4, 5)` is inside the triangle and outside the
    /// ring, so it crosses.
    #[test]
    fn rings_touching_a_concave_ring_at_every_vertex_are_classified_by_every_edge() {
        // (0,0)→(8,0)→(8,8)→(4,4)→(0,8): a square with a V notch from the
        // top, reflex at (4,4).
        let concave = ring(&[(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (4.0, 4.0), (0.0, 8.0)]);
        // Apex at the reflex vertex, two edges along the notch's walls.
        let notch = ring(&[(4.0, 4.0), (6.0, 6.0), (2.0, 6.0)]);
        assert_eq!(place_ring(&notch, &concave), RingPlacement::SharedEdge);
        // Apex on the bottom edge, top edge across the notch.
        let spans = ring(&[(4.0, 0.0), (6.0, 6.0), (2.0, 6.0)]);
        assert_eq!(place_ring(&spans, &concave), RingPlacement::Crossing);
        // Apex at the reflex vertex, departing downward into the body; the
        // other two vertices on the side edges.
        let below = ring(&[(4.0, 4.0), (0.0, 2.0), (8.0, 2.0)]);
        assert_eq!(place_ring(&below, &concave), RingPlacement::Inside);
    }

    /// An L-shaped exterior and a triangle with one vertex inside and two on
    /// the notch's edges. Its third edge runs through the notch, outside,
    /// with no proper crossing anywhere.
    #[test]
    fn a_ring_leaving_through_two_boundary_vertices_is_crossing() {
        let l = ring(&[
            (0.0, 0.0),
            (8.0, 0.0),
            (8.0, 4.0),
            (4.0, 4.0),
            (4.0, 8.0),
            (0.0, 8.0),
        ]);
        let tri = ring(&[(3.0, 3.0), (4.0, 6.0), (6.0, 4.0)]);
        assert_eq!(place_ring(&tri, &l), RingPlacement::Crossing);
    }

    /// `inner` passing straight through `outer`'s reflex vertex: the vertex
    /// lies in the middle of an edge of `inner`, which goes from inside to
    /// outside there.
    #[test]
    fn a_ring_passing_through_a_reflex_vertex_is_crossing() {
        let l = ring(&[
            (0.0, 0.0),
            (8.0, 0.0),
            (8.0, 4.0),
            (4.0, 4.0),
            (4.0, 8.0),
            (0.0, 8.0),
        ]);
        // Edge (2,2)→(6,6) passes through the reflex vertex (4,4) into the
        // notch.
        let tri = ring(&[(2.0, 2.0), (6.0, 6.0), (2.0, 6.0)]);
        assert_eq!(place_ring(&tri, &l), RingPlacement::Crossing);
        // Touching the reflex vertex from inside, without passing through,
        // is still inside.
        let touch = ring(&[(2.0, 2.0), (4.0, 4.0), (2.0, 6.0)]);
        assert_eq!(place_ring(&touch, &l), RingPlacement::Inside);
    }

    /// A chord whose two ends are on the boundary and whose middle passes
    /// through a reflex vertex of `outer`, from inside to outside. Neither end
    /// and no proper crossing reveals it; only cutting the edge at the vertex
    /// does.
    #[test]
    fn a_chord_through_a_reflex_vertex_is_cut_there() {
        // A U: notch x ∈ (4, 8), y > 2, reflex vertices (4,2) and (8,2).
        let u = ring(&[
            (0., 0.),
            (12., 0.),
            (12., 8.),
            (8., 8.),
            (8., 2.),
            (4., 2.),
            (4., 8.),
            (0., 8.),
        ]);
        // (0,0) → (8,4) runs along y = x/2: inside up to (4,2), then across
        // the notch to its right wall. (8,4) → (12,1) → (0,0) stay inside.
        let tri = ring(&[(0., 0.), (8., 4.), (12., 1.)]);
        assert_eq!(place_ring(&tri, &u), RingPlacement::Crossing);
    }

    /// Rings sharing a stretch of an edge: a collinear overlap, which no
    /// proper-crossing test detects.
    #[test]
    fn rings_sharing_an_edge_stretch_are_shared_edge_not_disjoint() {
        let a = ring(&[(1.0, 2.0), (1.0, 6.0), (3.0, 6.0), (3.0, 2.0)]);
        let b = ring(&[(3.0, 4.0), (3.0, 6.0), (6.0, 6.0), (6.0, 4.0)]);
        assert_eq!(place_ring(&a, &b), RingPlacement::SharedEdge);
        assert!(!rings_interior_disjoint(&a, &b));
        // A hole lying along part of the exterior's edge.
        let hole = ring(&[(0.0, 2.0), (0.0, 4.0), (3.0, 4.0), (3.0, 2.0)]);
        assert_eq!(place_ring(&hole, &unit_square()), RingPlacement::SharedEdge);
    }

    // ==================================================================
    // first_self_contact
    // ==================================================================

    fn coords(pts: &[(f64, f64)]) -> Vec<Coord2d> {
        pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect()
    }

    /// Collinear consecutive vertices are simple: the arrangement splits
    /// edges at every crossing, so output rings carry them routinely, and a
    /// simplicity test that refused them would refuse most overlays.
    #[test]
    fn collinear_consecutive_vertices_are_simple() {
        let v = coords(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (4.0, 6.0),
            (0.0, 6.0),
        ]);
        assert_eq!(first_self_contact(&v), None);
    }

    #[test]
    fn a_pinched_ring_is_not_simple() {
        // an inverted hole touching its exterior, reduced
        let v = coords(&[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 6.0),
            (3.0, 3.0),
            (4.0, 2.0),
            (2.0, 2.0),
            (3.0, 3.0),
            (0.0, 6.0),
        ]);
        assert!(first_self_contact(&v).is_some());
    }

    #[test]
    fn a_proper_crossing_is_still_caught() {
        let v = coords(&[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0)]);
        assert!(first_self_contact(&v).is_some());
    }

    /// The sweep and the double loop give the same contacts in the same
    /// order, on rings of many vertices on a small integer grid (so edges
    /// touch, overlap and cross, many at the edges of each other's boxes),
    /// in both modes.
    #[test]
    fn the_sweep_reports_what_the_double_loop_reports() {
        let mut state = 0x5EED_5EEDu64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        let mut with_contacts = 0;
        for case in 0..300 {
            let n = 65 + (next() % 200) as usize;
            let span = if case % 3 == 0 { 6 } else { 40 };
            let mut v: Vec<Coord2d> = (0..n)
                .map(|_| Coord2d::new((next() % span) as f64, (next() % span) as f64))
                .collect();
            v.dedup();
            if v.len() <= SWEEP_MIN_VERTICES {
                continue;
            }
            let all = self_contacts_sweep(&v, false);
            assert_eq!(all, self_contacts_scan(&v, false), "case {case}");
            assert_eq!(
                self_contacts_sweep(&v, true),
                self_contacts_scan(&v, true),
                "case {case}"
            );
            with_contacts += usize::from(!all.is_empty());
        }
        // A convex polygon of many vertices: no contact at all.
        let ring: Vec<Coord2d> = (0..400)
            .map(|k| {
                let t = std::f64::consts::TAU * f64::from(k) / 400.0;
                Coord2d::new(1000.0 * t.cos(), 1000.0 * t.sin())
            })
            .collect();
        assert!(self_contacts_sweep(&ring, false).is_empty());
        assert!(self_contacts_scan(&ring, false).is_empty());
        assert!(with_contacts > 200);
    }
}
