// src/ops/simplify_topology.rs

//! Topology-preserving simplification: Douglas–Peucker vertex reduction
//! that removes a vertex only where doing so changes no topological
//! relationship, decided exactly.
//!
//! # Why it exists
//!
//! [`Polygon::simplify`] is plain Ramer–Douglas–Peucker, vouched. RDP is not
//! topology-preserving, and at the tolerances simplification is used at it
//! routinely collapses small rings, pushes holes across exteriors and makes
//! parts overlap; its vouch rightly refuses those results. This operation
//! never produces them.
//!
//! # The promise
//!
//! For `simplify_preserve_topology(g, tolerance)`:
//!
//! 1. **Every output vertex is an input vertex,** unchanged, bit for bit.
//!    Nothing is constructed: vertices are only removed.
//! 2. **Every removed vertex lies within `tolerance` of the chord that
//!    replaced it,** as evaluated in `f64`. Which vertices to *try* removing
//!    is a selection rule, the same `f64` distance test as
//!    [`Polygon::simplify`]. It cannot make a result wrong, only different,
//!    because every removal is then decided exactly.
//! 3. **No two elements gain or lose a contact.** The elements are the
//!    rings of every polygon, the lines, and the points, all treated as
//!    point sets. Elements that were disjoint stay disjoint; where two
//!    touched, they still touch at the same points; nothing crosses that
//!    did not cross, and a line's own self-intersections are unchanged.
//! 4. **No element moves across another.** No ring or line passes over any
//!    part of another element. This keeps every hole inside its exterior,
//!    every part outside the others, and every point and vertex on the same
//!    side of every ring it was on.
//! 5. **The output is valid whenever the input is:** it is built by the
//!    checked constructors, and [`validity_issues`] of it is empty.
//!
//! For an invalid input the result is still computed, since operations do
//! not refuse invalid input, but it carries no vouch:
//! [`simplify_preserve_topology`] refuses with
//! [`UnvouchedReason::InvalidInput`], and
//! [`simplify_preserve_topology_unvouched`] returns it as
//! [`MaybeVouched::Unvouched`]. Promises 1 to 4 hold for it all the same, so
//! the result has exactly the input's defects and no new ones.
//!
//! # The algorithm
//!
//! **Fixed vertices** are never removed:
//!
//! - the endpoints of every line;
//! - every vertex whose coordinate occurs more than once in the input
//!   (vertices of any element, and points), compared by exact equality: a
//!   vertex shared by two elements, a touch between parts, a hole touching
//!   its exterior at a common vertex, a pinch, a point lying on a vertex;
//! - every vertex lying on another segment, other than at that segment's
//!   endpoints, decided by `orient2d`: a hole vertex resting on the inside
//!   of an exterior edge, or a line ending on another line.
//!
//! **Candidates** come from Douglas–Peucker on each element, between
//! consecutive fixed vertices. A ring with one fixed vertex, or none (then
//! vertex 0 is its anchor), is one run from the anchor all the way round
//! back to it; its chord has zero length, so it is never accepted, and the
//! run splits at its farthest vertex, which is how [`Polygon::simplify`]
//! anchors a ring too.
//!
//! ## Accepting a chord
//!
//! Let a run `v_i … v_j` be a candidate, its chord
//! `c = v_i v_j`, and `P` the region bounded by the run and the chord, by
//! the even-odd rule (the closed curve may cross itself). Let `E` be every
//! other current segment of every element, including the same element's
//! segments outside the run, and every point. The chord is accepted exactly
//! when
//!
//! > `E ∩ cl(P) ⊆ {v_i, v_j}`, the run does not touch itself, and
//! > `v_i ≠ v_j`.
//!
//! It is decided exactly, with `orient2d` and comparisons of stored
//! coordinates, on the candidates the index returns for the run's box,
//! which contains `cl(P)`. Four checks make up the condition:
//!
//! - **(a) the chord:** no segment of `E` meets `c` other than at `v_i` or
//!   `v_j`, and none overlaps it along a stretch; no point lies on it;
//! - **(b) the region:** no endpoint of a segment of `E`, and no point,
//!   lies inside `P`, other than `v_i` and `v_j`;
//! - **(c) the run:** no segment of `E` meets the run other than at `v_i`
//!   or `v_j`; no point lies on it;
//! - **(d) the run itself:** no two segments of the run touch, except
//!   consecutive ones at their common vertex.
//!
//! These give the condition: a segment of `E` that meets `cl(P)` at a point
//! other than `v_i`, `v_j` either has an endpoint in `P` (b), or meets
//! `P`'s boundary, which is the chord and the run, at such a point (a, c).
//!
//! A rejected chord is split at the run's farthest vertex, and both halves
//! are tried in turn, as Douglas–Peucker recurses. Elements are processed
//! in a fixed order (each polygon's exterior, then its holes, polygon by
//! polygon, in the order of the input; lines in order), and runs within an
//! element by the recursion, so the result is deterministic.
//!
//! # Why the condition preserves topology
//!
//! Accepting one chord replaces the run by `c` in one element and leaves
//! every other segment where it is. Promises 3 and 4 then follow, one chord
//! at a time.
//!
//! ## Contacts (promise 3)
//!
//! By (c), the run meets the rest of the input
//! only at `v_i` and `v_j`; by (a), so does the chord. Both endpoints stay.
//! So the element's intersection with the rest of the input, as a point
//! set, is the same before and after: the unchanged part, plus whichever of
//! `v_i`, `v_j` the rest passes through, which is the same on both sides.
//! (d) and the same-element segments in `E` make the same true of the
//! element's intersection with itself.
//!
//! At `v_i` (or `v_j`), whether another element *crosses* or only *touches*
//! there depends on the order of the segments around the point, and the
//! element's outgoing segment turns from the run's first segment to the
//! chord. Near `v_i`, `P` is one of the two wedges between those two
//! directions. The element's own incoming segment lies in `E`, so not in
//! `P`, by (b) and (c): it is in the other wedge. So the turn sweeps
//! exactly `P`'s wedge, which holds nothing of `E`. No segment at `v_i`
//! changes side.
//!
//! ## No move across (promise 4)
//!
//! For a ring, the interiors before and after differ by exactly the interior of
//! `P` under the even-odd rule: the symmetric difference of the regions bounded
//! by two closed curves that share everything but the run and the chord is the
//! region bounded by the run and the chord. Nothing of any other element is in
//! `cl(P)` except `v_i`, `v_j`, which lie on both rings. So every point of
//! every other element is inside, on, or outside the ring after exactly as
//! before. That keeps holes inside their exteriors (and outside one another)
//! and parts outside one another, and it is how a ring is kept from sweeping
//! over a whole small element, which has no segment crossing the ring at all:
//! such an element would have a vertex inside `P`, which (b) refuses.
//!
//! ## Degeneration
//!
//! A ring cannot collapse to two vertices or to zero
//! area: either would need a chord overlapping a remaining segment along a
//! stretch, which (a) refuses. A chord with `v_i = v_j` is never accepted.
//!
//! ## Validity (promise 5)
//!
//! OGC validity is a statement about contacts:
//! simple rings, holes inside the exterior and interior-disjoint, parts
//! meeting at points, a connected interior. Promises 3 and 4 preserve every
//! contact and every side relation, so a valid input gives a valid output.
//!
//! # The vouch
//!
//! The output is built through the checked constructors, and
//! [`validity_issues`] is run on it. By the argument above it is always
//! clean for a valid input, and the check is not expected ever to fire. It
//! is kept as a defense: if it fires, the result is refused with
//! [`UnvouchedReason::InvalidResult`] rather than returned as valid. A
//! constructor that refuses the output (which the same argument rules out
//! for a valid input) returns its error from both entry points.
//!
//! # The candidate index
//!
//! The current segments change as runs are replaced, so they are held in
//! the crate's dynamic R\*-tree: an accepted chord removes the run's
//! segments and inserts itself. The index only filters candidates, with no
//! false negatives: every segment whose box meets the run's box is
//! returned, and every decision is then made exactly on the candidate.
//!
//! # Tolerance
//!
//! As [`Polygon::simplify`]: a vertex may go when its distance from the
//! chord is **not greater than** `tolerance`; zero is the identity; a
//! negative or non-finite tolerance is [`GeometryError::InvalidTolerance`].

use crate::GeometryError;
use crate::exact::kernel::{edge_crosses_ray, on_segment, orient, segments_properly_cross};
use crate::index::SegmentIndex;
use crate::ops::simplify::{check_tolerance, segment_distance};
use crate::ops::validate::validity_issues;
use crate::predicates::Relatable;
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MaybeVouched,
    MultiLineString, MultiPolygon, Polygon, Unvouched, UnvouchedReason,
};
use arithmetica::Sign;

#[cfg(test)]
mod tests;

// ===========================================================================
// The public entry points
// ===========================================================================

/// The geometry types [`simplify_preserve_topology`] accepts: those with
/// vertices to remove. Sealed.
///
/// Implemented for [`LineString`], [`MultiLineString`], [`Polygon`],
/// [`MultiPolygon`] and [`GeometryCollection`]. A collection's points and
/// multipoints are kept unchanged, as obstacles.
pub trait TopologySimplifiable: sealed::Sealed + Relatable + Clone {}

impl TopologySimplifiable for LineString {}
impl TopologySimplifiable for MultiLineString {}
impl TopologySimplifiable for Polygon {}
impl TopologySimplifiable for MultiPolygon {}
impl TopologySimplifiable for GeometryCollection {}

/// Douglas–Peucker simplification that preserves topology, vouched: every
/// vertex removed is removed only where doing so changes no topological
/// relationship, decided exactly. See the [module documentation](self) for
/// the promise, the algorithm and the proof.
///
/// The result has the same type as `g`, and every one of its vertices is a
/// vertex of `g`.
///
/// # Errors
///
/// - [`GeometryError::InvalidTolerance`] for a negative or non-finite
///   tolerance.
/// - [`GeometryError::Unvouched`] with [`UnvouchedReason::InvalidInput`]
///   when `g` has a validity issue; the result exists and is available from
///   [`simplify_preserve_topology_unvouched`]. With
///   [`UnvouchedReason::InvalidResult`] if the defensive check of the output
///   ever fires (not expected; see the module documentation).
/// - [`GeometryError::TooManyEntries`] for an input with more than
///   `u32::MAX` vertices and points.
/// - A constructor's error if the output cannot be built, which the
///   argument rules out for a valid input.
pub fn simplify_preserve_topology<G: TopologySimplifiable>(
    g: &G,
    tolerance: f64,
) -> Result<G, GeometryError> {
    simplify_preserve_topology_unvouched(g, tolerance)?.vouched()
}

/// [`simplify_preserve_topology`], returning a result that is not vouched
/// for as [`MaybeVouched::Unvouched`] instead of refusing it.
///
/// A valid input gives [`MaybeVouched::Vouched`]. An invalid input gives
/// the computed result with [`UnvouchedReason::InvalidInput`] and the
/// input's first validity issue: promises 1 to 4 hold for it, so its
/// defects are exactly the input's.
///
/// # Errors
///
/// As [`simplify_preserve_topology`], minus the `Unvouched` refusal.
pub fn simplify_preserve_topology_unvouched<G: TopologySimplifiable>(
    g: &G,
    tolerance: f64,
) -> Result<MaybeVouched<G>, GeometryError> {
    check_tolerance(tolerance)?;
    let input_issue = validity_issues(g).first().map(|i| i.issue);
    let out = if tolerance == 0.0 {
        g.clone()
    } else {
        simplify_raw(g, tolerance)?
    };
    Ok(match input_issue {
        Some(issue) => {
            MaybeVouched::Unvouched(Unvouched::new(out, UnvouchedReason::InvalidInput { issue }))
        }
        None => vouch(out),
    })
}

/// The defensive check: `out`, vouched when [`validity_issues`] finds
/// nothing, refused with its first issue otherwise.
fn vouch<G: Relatable>(out: G) -> MaybeVouched<G> {
    match validity_issues(&out).first().map(|i| i.issue) {
        None => MaybeVouched::Vouched(out),
        Some(issue) => MaybeVouched::Unvouched(Unvouched::new(
            out,
            UnvouchedReason::InvalidResult { issue },
        )),
    }
}

/// Flatten, simplify every element, rebuild through the constructors.
fn simplify_raw<G: TopologySimplifiable>(g: &G, tolerance: f64) -> Result<G, GeometryError> {
    let mut flat = sealed::Flat::default();
    g.flatten(&mut flat);
    let kept = Simplifier::new(&flat, tolerance)?.run();
    g.rebuild(&mut sealed::Kept::new(kept))
}

// ===========================================================================
// Flattening and rebuilding, per type (sealed)
// ===========================================================================

pub(super) mod sealed {
    use super::{
        Coord2d, CoordSeq, Geometry, GeometryCollection, GeometryError, LineString, LinearRing,
        MultiLineString, MultiPolygon, Polygon,
    };

    /// One chain: a ring (closed, stored without a closing vertex) or a line.
    #[derive(Clone, Copy, Debug)]
    pub struct Chain {
        pub(in crate::ops) start: usize,
        pub(in crate::ops) len: usize,
        pub(in crate::ops) closed: bool,
    }

    /// The input as chains of vertices, in processing order, and points.
    #[derive(Debug, Default)]
    pub struct Flat {
        pub(in crate::ops) verts: Vec<Coord2d>,
        pub(in crate::ops) chains: Vec<Chain>,
        pub(in crate::ops) points: Vec<Coord2d>,
    }

    impl Flat {
        fn push(&mut self, coords: &[Coord2d], closed: bool) {
            self.chains.push(Chain {
                start: self.verts.len(),
                len: coords.len(),
                closed,
            });
            self.verts.extend_from_slice(coords);
        }
        fn ring(&mut self, r: &LinearRing) {
            self.push(r.coords().as_slice(), true);
        }
        fn line(&mut self, l: &LineString) {
            self.push(l.coords().as_slice(), false);
        }
        fn polygon(&mut self, p: &Polygon) {
            self.ring(p.exterior());
            for h in p.holes() {
                self.ring(h);
            }
        }
        /// A collection element: chains for rings and lines, obstacles for
        /// points. Nested collections are walked in place.
        fn geometry(&mut self, g: &Geometry) {
            match g {
                Geometry::Point(p) => self.points.push(p.coord()),
                Geometry::MultiPoint(m) => self.points.extend(m.points().iter().map(|p| p.coord())),
                Geometry::LineString(l) => self.line(l),
                Geometry::MultiLineString(m) => m.line_strings().iter().for_each(|l| self.line(l)),
                Geometry::Polygon(p) => self.polygon(p),
                Geometry::MultiPolygon(m) => m.polygons().iter().for_each(|p| self.polygon(p)),
                Geometry::GeometryCollection(c) => {
                    c.geometries().iter().for_each(|e| self.geometry(e));
                }
            }
        }
    }

    /// The kept vertices of each chain, consumed in the order they were
    /// flattened.
    pub struct Kept {
        chains: std::vec::IntoIter<Vec<Coord2d>>,
    }

    impl Kept {
        pub(in crate::ops) fn new(chains: Vec<Vec<Coord2d>>) -> Self {
            Self {
                chains: chains.into_iter(),
            }
        }
        /// The next chain's vertices. Rebuilding consumes exactly the chains
        /// flattening produced, in the same order, so the list cannot run
        /// out; were it to, the empty chain would fail its constructor, an
        /// error rather than a panic.
        fn next(&mut self) -> CoordSeq {
            CoordSeq::new(self.chains.next().unwrap_or_default())
        }
        fn ring(&mut self) -> Result<LinearRing, GeometryError> {
            LinearRing::new(self.next())
        }
        fn line(&mut self) -> Result<LineString, GeometryError> {
            LineString::new(self.next())
        }
        fn polygon(&mut self, p: &Polygon) -> Result<Polygon, GeometryError> {
            let exterior = self.ring()?;
            let holes = p
                .holes()
                .iter()
                .map(|_| self.ring())
                .collect::<Result<Vec<_>, _>>()?;
            Polygon::new(exterior, holes)
        }
        fn geometry(&mut self, g: &Geometry) -> Result<Geometry, GeometryError> {
            Ok(match g {
                Geometry::Point(_) | Geometry::MultiPoint(_) => g.clone(),
                Geometry::LineString(_) => Geometry::LineString(self.line()?),
                Geometry::MultiLineString(m) => Geometry::MultiLineString(self.lines(m)?),
                Geometry::Polygon(p) => Geometry::Polygon(self.polygon(p)?),
                Geometry::MultiPolygon(m) => Geometry::MultiPolygon(self.polygons(m)?),
                Geometry::GeometryCollection(c) => {
                    Geometry::GeometryCollection(self.collection(c)?)
                }
            })
        }
        fn lines(&mut self, m: &MultiLineString) -> Result<MultiLineString, GeometryError> {
            let lines = m
                .line_strings()
                .iter()
                .map(|_| self.line())
                .collect::<Result<Vec<_>, _>>()?;
            Ok(MultiLineString::new(lines))
        }
        fn polygons(&mut self, m: &MultiPolygon) -> Result<MultiPolygon, GeometryError> {
            let parts = m
                .polygons()
                .iter()
                .map(|p| self.polygon(p))
                .collect::<Result<Vec<_>, _>>()?;
            MultiPolygon::new(parts)
        }
        fn collection(
            &mut self,
            c: &GeometryCollection,
        ) -> Result<GeometryCollection, GeometryError> {
            let elements = c
                .geometries()
                .iter()
                .map(|e| self.geometry(e))
                .collect::<Result<Vec<_>, _>>()?;
            GeometryCollection::new(elements)
        }
    }

    /// Flattening to chains and rebuilding from the kept vertices.
    pub trait Sealed: Sized {
        /// Appends the chains and points of `self`, in processing order.
        fn flatten(&self, flat: &mut Flat);
        /// The value of the same shape as `self`, built from the kept
        /// vertices through the checked constructors.
        fn rebuild(&self, kept: &mut Kept) -> Result<Self, GeometryError>;
    }

    impl Sealed for LineString {
        fn flatten(&self, flat: &mut Flat) {
            flat.line(self);
        }
        fn rebuild(&self, kept: &mut Kept) -> Result<Self, GeometryError> {
            kept.line()
        }
    }

    impl Sealed for MultiLineString {
        fn flatten(&self, flat: &mut Flat) {
            self.line_strings().iter().for_each(|l| flat.line(l));
        }
        fn rebuild(&self, kept: &mut Kept) -> Result<Self, GeometryError> {
            kept.lines(self)
        }
    }

    impl Sealed for Polygon {
        fn flatten(&self, flat: &mut Flat) {
            flat.polygon(self);
        }
        fn rebuild(&self, kept: &mut Kept) -> Result<Self, GeometryError> {
            kept.polygon(self)
        }
    }

    impl Sealed for MultiPolygon {
        fn flatten(&self, flat: &mut Flat) {
            self.polygons().iter().for_each(|p| flat.polygon(p));
        }
        fn rebuild(&self, kept: &mut Kept) -> Result<Self, GeometryError> {
            kept.polygons(self)
        }
    }

    impl Sealed for GeometryCollection {
        fn flatten(&self, flat: &mut Flat) {
            self.geometries().iter().for_each(|g| flat.geometry(g));
        }
        fn rebuild(&self, kept: &mut Kept) -> Result<Self, GeometryError> {
            kept.collection(self)
        }
    }
}

// ===========================================================================
// The simplifier
// ===========================================================================

/// "No next vertex": the end of a line.
const END: u32 = u32::MAX;

/// What an index id names: a current segment, or a point.
#[derive(Clone, Copy, Debug)]
pub(in crate::ops) enum Obstacle {
    Segment(Coord2d, Coord2d),
    Point(Coord2d),
}

/// The state of one simplification.
///
/// Vertices are numbered globally, chain after chain (`gid`). A segment is
/// named by the id of its start vertex, and ends at that vertex's `next`;
/// a point is named by `verts.len()` plus its position. Every id fits in
/// `u32`, checked once at the start, and [`END`] is never an id.
struct Simplifier<'a> {
    verts: &'a [Coord2d],
    chains: &'a [sealed::Chain],
    points: &'a [Coord2d],
    tolerance: f64,
    /// The next current vertex of the same chain, or [`END`].
    next: Vec<u32>,
    alive: Vec<bool>,
    fixed: Vec<bool>,
    index: SegmentIndex,
    /// Scratch: the start vertices of the segments of the run under test.
    in_run: Vec<bool>,
}

/// A candidate run: its vertices' ids, `v_i` first and `v_j` last.
pub(in crate::ops) struct Run {
    pub(in crate::ops) ids: Vec<u32>,
    pub(in crate::ops) pts: Vec<Coord2d>,
    pub(in crate::ops) lo: Coord2d,
    pub(in crate::ops) hi: Coord2d,
}

impl Run {
    fn vi(&self) -> Coord2d {
        self.pts[0]
    }
    fn vj(&self) -> Coord2d {
        self.pts[self.pts.len() - 1]
    }
    /// `v_i` or `v_j`: the points the run and the chord may share with the
    /// rest of the input.
    fn allowed(&self, p: Coord2d) -> bool {
        p == self.vi() || p == self.vj()
    }
    fn segments(&self) -> impl Iterator<Item = (Coord2d, Coord2d)> + '_ {
        self.pts.windows(2).map(|w| (w[0], w[1]))
    }
    /// The run's box, which contains the chord and so `cl(P)`; a point
    /// outside it is outside `P`.
    fn contains_box(&self, p: Coord2d) -> bool {
        p.x >= self.lo.x && p.x <= self.hi.x && p.y >= self.lo.y && p.y <= self.hi.y
    }
}

/// The closed boxes of two segments meet. Exact: comparisons of stored
/// coordinates; it only skips pairs that cannot touch.
fn boxes_meet((a, b): (Coord2d, Coord2d), (c, d): (Coord2d, Coord2d)) -> bool {
    a.x.min(b.x) <= c.x.max(d.x)
        && c.x.min(d.x) <= a.x.max(b.x)
        && a.y.min(b.y) <= c.y.max(d.y)
        && c.y.min(d.y) <= a.y.max(b.y)
}

/// Segments `ab` and `pq` are collinear and overlap along a stretch of
/// positive length. Exact: two `orient2d` signs, then interval comparison
/// of stored coordinates along an axis on which `p ≠ q` (on a line, the
/// order along that axis is the order along the line).
#[allow(clippy::float_cmp)] // exact comparison of stored coordinates is intended
fn collinear_overlap(a: Coord2d, b: Coord2d, p: Coord2d, q: Coord2d) -> bool {
    if orient(p, q, a) != Sign::Zero || orient(p, q, b) != Sign::Zero {
        return false;
    }
    let along = |c: Coord2d| if p.x == q.x { c.y } else { c.x };
    let (s0, s1) = (along(a).min(along(b)), along(a).max(along(b)));
    let (t0, t1) = (along(p).min(along(q)), along(p).max(along(q)));
    s0.max(t0) < s1.min(t1)
}

/// The segment `ab` meets the segment `pq` at a point that is not allowed,
/// or overlaps it along a stretch. Every contact configuration is covered:
/// a proper crossing, an endpoint of either on the other (a T-junction, a
/// shared vertex), and collinear overlap. Exact: `orient2d` on stored
/// coordinates and comparisons; nothing is constructed.
fn meets_outside(
    (a, b): (Coord2d, Coord2d),
    (p, q): (Coord2d, Coord2d),
    allowed: impl Fn(Coord2d) -> bool,
) -> bool {
    segments_properly_cross(a, b, p, q)
        || (on_segment(a, p, q) && !allowed(a))
        || (on_segment(b, p, q) && !allowed(b))
        || (on_segment(p, a, b) && !allowed(p))
        || (on_segment(q, a, b) && !allowed(q))
        || collinear_overlap(a, b, p, q)
}

impl<'a> Simplifier<'a> {
    fn new(flat: &'a sealed::Flat, tolerance: f64) -> Result<Self, GeometryError> {
        let count = flat.verts.len() + flat.points.len();
        // Ids run to `count - 1`; `END` (u32::MAX) must stay free.
        if u32::try_from(count).map_or(true, |n| n == u32::MAX) {
            return Err(GeometryError::TooManyEntries { count });
        }
        let n = flat.verts.len();
        let mut s = Self {
            verts: &flat.verts,
            chains: &flat.chains,
            points: &flat.points,
            tolerance,
            next: vec![END; n],
            alive: vec![true; n],
            fixed: vec![false; n],
            index: SegmentIndex::new(),
            in_run: vec![false; n],
        };
        s.link_and_index();
        s.fix_vertices();
        Ok(s)
    }

    /// An id as `u32`. Every id is below the count checked in `new`.
    fn id(i: usize) -> u32 {
        u32::try_from(i).unwrap_or(END)
    }

    fn at(&self, id: u32) -> Coord2d {
        self.verts[id as usize]
    }

    /// Links each chain's vertices and indexes every segment and point.
    fn link_and_index(&mut self) {
        for ch in self.chains {
            for k in 0..ch.len {
                let here = ch.start + k;
                let succ = if k + 1 < ch.len {
                    Some(here + 1)
                } else if ch.closed {
                    Some(ch.start)
                } else {
                    None
                };
                if let Some(t) = succ {
                    self.next[here] = Self::id(t);
                    self.index
                        .insert(self.verts[here], self.verts[t], Self::id(here));
                }
            }
        }
        let n = self.verts.len();
        for (k, &p) in self.points.iter().enumerate() {
            self.index.insert(p, p, Self::id(n + k));
        }
    }

    /// What an index id names, in the current state.
    fn obstacle(&self, id: u32) -> Obstacle {
        let n = self.verts.len();
        let i = id as usize;
        if i < n {
            Obstacle::Segment(self.verts[i], self.at(self.next[i]))
        } else {
            Obstacle::Point(self.points[i - n])
        }
    }

    /// Marks the fixed vertices: line endpoints, repeated coordinates, and
    /// vertices on another segment (see the module documentation).
    fn fix_vertices(&mut self) {
        for ch in self.chains {
            if !ch.closed {
                self.fixed[ch.start] = true;
                self.fixed[ch.start + ch.len - 1] = true;
            }
        }
        self.fix_repeated();
        self.fix_on_segment();
    }

    /// Every vertex whose coordinate occurs more than once among all
    /// vertices and points. Sorted by the exact lexicographic order, in
    /// which `-0.0` and `0.0` are equal, as `==` has them; equal coordinates
    /// are then adjacent.
    fn fix_repeated(&mut self) {
        let n = self.verts.len();
        let coord = |k: usize| {
            if k < n {
                self.verts[k]
            } else {
                self.points[k - n]
            }
        };
        let mut order: Vec<usize> = (0..n + self.points.len()).collect();
        order.sort_by(|&i, &j| {
            let (p, q) = (coord(i), coord(j));
            p.x.partial_cmp(&q.x)
                .unwrap_or(core::cmp::Ordering::Equal)
                .then(p.y.partial_cmp(&q.y).unwrap_or(core::cmp::Ordering::Equal))
        });
        let mut start = 0;
        while start < order.len() {
            let mut end = start + 1;
            while end < order.len() && coord(order[end]) == coord(order[start]) {
                end += 1;
            }
            if end - start > 1 {
                for &k in &order[start..end] {
                    if k < n {
                        self.fixed[k] = true;
                    }
                }
            }
            start = end;
        }
    }

    /// Every vertex lying on a segment of which it is not an endpoint.
    fn fix_on_segment(&mut self) {
        for i in 0..self.verts.len() {
            let p = self.verts[i];
            let mut hit = false;
            self.index.candidates(p, p, |id| {
                if let Obstacle::Segment(a, b) = self.obstacle(id) {
                    hit |= p != a && p != b && on_segment(p, a, b);
                }
            });
            self.fixed[i] |= hit;
        }
    }

    /// Simplifies every chain, in order, and returns each chain's kept
    /// vertices.
    fn run(mut self) -> Vec<Vec<Coord2d>> {
        for c in 0..self.chains.len() {
            let ch = self.chains[c];
            let mut anchors: Vec<usize> =
                (0..ch.len).filter(|&k| self.fixed[ch.start + k]).collect();
            if ch.closed {
                if anchors.is_empty() {
                    anchors.push(0);
                }
                for (k, &lo) in anchors.iter().enumerate() {
                    let hi = anchors.get(k + 1).copied().unwrap_or(anchors[0] + ch.len);
                    self.simplify_run(ch, lo, hi);
                }
            } else {
                for w in anchors.windows(2) {
                    self.simplify_run(ch, w[0], w[1]);
                }
            }
        }
        self.chains
            .iter()
            .map(|ch| {
                (ch.start..ch.start + ch.len)
                    .filter(|&g| self.alive[g])
                    .map(|g| self.verts[g])
                    .collect()
            })
            .collect()
    }

    /// The global id of position `k` of a chain, positions running past the
    /// end of a ring wrapping round.
    fn gid(ch: sealed::Chain, k: usize) -> usize {
        ch.start + k % ch.len
    }

    /// Douglas–Peucker on the run between positions `lo` and `hi` of `ch`,
    /// with an explicit stack (a parcel ring can have thousands of
    /// vertices). Every vertex strictly between `lo` and `hi` is alive when
    /// its run is tried: runs are disjoint, and a run is tried before any
    /// of its parts.
    fn simplify_run(&mut self, ch: sealed::Chain, lo: usize, hi: usize) {
        let mut stack = vec![(lo, hi)];
        while let Some((lo, hi)) = stack.pop() {
            if hi < lo + 2 {
                continue;
            }
            let (a, b) = (self.verts[Self::gid(ch, lo)], self.verts[Self::gid(ch, hi)]);
            let (split, within) = self.farthest(ch, lo, hi, a, b);
            if within && a != b {
                let run = self.make_run(ch, lo, hi);
                if self.accept(&run) {
                    self.commit(&run);
                    continue;
                }
            }
            // Halves pushed so that the first is tried first.
            stack.push((split, hi));
            stack.push((lo, split));
        }
    }

    /// The farthest interior vertex from the chord `ab` (the first one
    /// attaining the maximum, or the first whose distance could not be
    /// computed, which is kept conservatively), and whether every interior
    /// vertex is within tolerance. The `f64` selection rule.
    fn farthest(
        &self,
        ch: sealed::Chain,
        lo: usize,
        hi: usize,
        a: Coord2d,
        b: Coord2d,
    ) -> (usize, bool) {
        let mut best: Option<(usize, f64)> = None;
        for k in lo + 1..hi {
            match segment_distance(self.verts[Self::gid(ch, k)], a, b) {
                None => return (k, false),
                Some(d) => {
                    if best.is_none_or(|(_, bd)| d > bd) {
                        best = Some((k, d));
                    }
                }
            }
        }
        // `hi >= lo + 2`, so the loop ran and `best` is set.
        best.map_or((lo + 1, false), |(k, d)| (k, d <= self.tolerance))
    }

    fn make_run(&self, ch: sealed::Chain, lo: usize, hi: usize) -> Run {
        let ids: Vec<u32> = (lo..=hi).map(|k| Self::id(Self::gid(ch, k))).collect();
        let pts: Vec<Coord2d> = ids.iter().map(|&i| self.at(i)).collect();
        let mut lo_c = pts[0];
        let mut hi_c = pts[0];
        for p in &pts {
            lo_c = Coord2d::new(lo_c.x.min(p.x), lo_c.y.min(p.y));
            hi_c = Coord2d::new(hi_c.x.max(p.x), hi_c.y.max(p.y));
        }
        Run {
            ids,
            pts,
            lo: lo_c,
            hi: hi_c,
        }
    }

    /// The acceptance test: `E ∩ cl(P) ⊆ {v_i, v_j}` and (d). Exact.
    fn accept(&mut self, run: &Run) -> bool {
        if run_touches_itself(run) {
            return false;
        }
        let m = run.ids.len() - 1;
        for &id in &run.ids[..m] {
            self.in_run[id as usize] = true;
        }
        let mut candidates = Vec::new();
        self.index
            .candidates(run.lo, run.hi, |id| candidates.push(id));
        let ok = candidates.iter().all(|&id| {
            (id as usize) < self.in_run.len() && self.in_run[id as usize]
                || clear_of(run, self.obstacle(id))
        });
        for &id in &run.ids[..m] {
            self.in_run[id as usize] = false;
        }
        ok
    }

    /// Replaces the run by its chord: the run's segments leave the index,
    /// the chord enters it, and the interior vertices are dropped.
    fn commit(&mut self, run: &Run) {
        let m = run.ids.len() - 1;
        for k in 0..m {
            let removed = self.index.remove(run.pts[k], run.pts[k + 1], run.ids[k]);
            debug_assert!(removed, "a current segment is in the index");
        }
        for &id in &run.ids[1..m] {
            self.alive[id as usize] = false;
        }
        self.next[run.ids[0] as usize] = run.ids[m];
        self.index.insert(run.vi(), run.vj(), run.ids[0]);
    }
}

// ===========================================================================
// The four checks
// ===========================================================================

/// One obstacle of `E` meets `cl(P)` at most at `v_i`, `v_j`.
pub(in crate::ops) fn clear_of(run: &Run, o: Obstacle) -> bool {
    match o {
        Obstacle::Point(p) => run.allowed(p) || point_clear(run, p),
        Obstacle::Segment(a, b) => {
            !chord_meets(run, (a, b))
                && !run_meets(run, (a, b))
                && [a, b]
                    .iter()
                    .all(|&e| run.allowed(e) || !inside_region(run, e))
        }
    }
}

/// (a) The segment meets the chord other than at `v_i`, `v_j`, or overlaps
/// it.
fn chord_meets(run: &Run, s: (Coord2d, Coord2d)) -> bool {
    let c = (run.vi(), run.vj());
    boxes_meet(s, c) && meets_outside(s, c, |p| run.allowed(p))
}

/// (c) The segment meets the run other than at `v_i`, `v_j`.
fn run_meets(run: &Run, s: (Coord2d, Coord2d)) -> bool {
    run.segments()
        .any(|r| boxes_meet(s, r) && meets_outside(s, r, |p| run.allowed(p)))
}

/// (b) The point, known not to be `v_i` or `v_j`, lies inside `P`.
///
/// Called on an endpoint that (a) and (c) have already shown is not on
/// `P`'s boundary (a point on the chord or the run is a contact they
/// report), so the even-odd crossing count decides interior or exterior.
/// The count is the crossing-number test of point-in-ring location, with
/// its half-open rule, over the run's segments and the chord.
fn inside_region(run: &Run, p: Coord2d) -> bool {
    if !run.contains_box(p) {
        return false;
    }
    let mut inside = edge_crosses_ray(run.vj(), run.vi(), p);
    for (a, b) in run.segments() {
        if edge_crosses_ray(a, b, p) {
            inside = !inside;
        }
    }
    inside
}

/// A point obstacle, other than `v_i`, `v_j`, is clear of `cl(P)`: not on
/// the chord (a), not on the run (c), not inside (b).
fn point_clear(run: &Run, p: Coord2d) -> bool {
    !on_segment(p, run.vi(), run.vj())
        && !run.segments().any(|(a, b)| on_segment(p, a, b))
        && !inside_region(run, p)
}

/// (d) Two segments of the run touch other than consecutive ones at their
/// common vertex: a fold-back of consecutive segments, or any contact of
/// others. Exact, as the ring self-contact test.
pub(in crate::ops) fn run_touches_itself(run: &Run) -> bool {
    let s: Vec<(Coord2d, Coord2d)> = run.segments().collect();
    for i in 0..s.len() {
        if i + 1 < s.len() {
            let ((a, b), (_, c)) = (s[i], s[i + 1]);
            if on_segment(c, a, b) || on_segment(a, b, c) {
                return true;
            }
        }
        for j in i + 2..s.len() {
            if boxes_meet(s[i], s[j]) && meets_outside(s[i], s[j], |_| false) {
                return true;
            }
        }
    }
    false
}
