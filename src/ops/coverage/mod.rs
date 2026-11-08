// src/ops/coverage/mod.rs

//! Polygonal coverages and noded linework: [`coverage_validate`] (with an
//! optional narrow-gap check, [`coverage_validate_with_gap`]),
//! [`coverage_union`], and [`polygonize_noded`].
//!
//! # Why these belong here
//!
//! Every operation in this module works on input that is **already exactly
//! noded**: adjacent polygons share their vertices bit for bit, and lines meet
//! only at shared endpoints. So every output coordinate is an input
//! coordinate. Nothing is intersected, rounded or moved. The versions that
//! node general input first, and therefore construct coordinates, belong to
//! a constructive crate rather than this one.
//!
//! # The two exact building blocks
//!
//! ## Edge identity
//!
//! Two directed edges are the same undirected edge when
//! their endpoint pairs are equal, in either order, by exact coordinate
//! equality (`-0.0` equals `0.0`, as coordinates do everywhere in this
//! crate).
//!
//! ## **Angular order at a vertex**
//!
//! This is the one numerical primitive the module adds. The edges leaving a
//! vertex `v` are sorted *counterclockwise* from the `+x` direction:
//!
//! 1. Each direction `a − v` is classified into a half-turn: the upper
//!    half `[0, π)` (`dy > 0`, or `dy = 0` and `dx > 0`) or the lower half
//!    `[π, 2π)`. The signs of `dx = fl(a.x − v.x)` and `dy` are the signs of
//!    the exact differences: subtraction is correctly rounded, rounding
//!    never changes a sign, and a nonzero difference never rounds to zero
//!    (a difference too small to be normal is exact, by gradual underflow).
//!    Opposite directions always fall in different halves.
//! 2. Within one half, every pair of directions spans less than a half
//!    turn, so `b` comes after `a` exactly when the exact predicate
//!    `orient2d(v, a, b) > 0` holds.
//!
//! Nothing is constructed, no angle is computed, and every comparison is
//! exact. Two distinct edges in the same direction from `v` would overlap
//! collinearly; the validations exclude that before the order is used.
//!
//! # Tracing faces
//!
//! Faces are traced with the face on the left of each directed edge: from
//! the edge `u → v`, the next edge is the one **first clockwise** from the
//! reverse direction `v → u` among the edges leaving `v`, which is the
//! tightest left turn. (That is the standard half-edge successor; the next edge
//! *counterclockwise* from the reverse would trace faces on the right.)
//!
//! A face whose boundary touches itself (a hole touching its exterior at a
//! vertex, two parts of a union meeting at a corner) has **one** boundary
//! walk that visits the touching vertex twice; no turn rule can split it,
//! because the walk is a property of the embedding. So each walk is split at
//! repeated vertices into simple closed loops, and each loop is classified
//! by the certified sign of its area: positive (counterclockwise) is an
//! exterior, negative (clockwise) a hole. Rings that touch at a vertex are
//! valid, so the split loops are the valid rings OGC asks for.
//!
//! Holes are assigned to the innermost exterior that contains them, by exact
//! ring placement ([`place_ring`](crate::exact::place_ring)); "innermost" is
//! decided by placement too, never by comparing approximate areas.

mod angular;
pub(in crate::ops) mod graph;
mod polygonize;
mod union;
mod validate;

#[cfg(test)]
mod tests;

use crate::error::GeometryError;
use crate::types::{MultiLineString, MultiPolygon, Polygon};

// ===========================================================================
// Coverage validation
// ===========================================================================

/// One reason a set of polygons is not a valid coverage.
///
/// Indices refer to the input slice: `polygon` and `other` are positions in
/// it, `ring` is `0` for the exterior and `1..` for the holes in order, and
/// `edge` is the edge starting at that ring's vertex `edge`.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CoverageIssue {
    /// A ring of the polygon touches or crosses itself, decided by the exact
    /// self-contact test. Polygon constructors do not check simplicity, and
    /// a coverage of non-simple polygons has no well-defined shared
    /// boundary, so it is not a valid coverage.
    NonSimpleRing {
        /// The polygon.
        polygon: usize,
        /// The ring (`0` exterior, `1..` holes).
        ring: usize,
    },
    /// The interiors of polygons `a` and `b` meet (`a < b`), decided by the
    /// exact `relate` (`II ≠ F`). Duplicate polygons, a polygon crossing
    /// into another, and one inside another's interior all show here.
    Overlap {
        /// The first polygon.
        a: usize,
        /// The second polygon.
        b: usize,
    },
    /// An edge meets polygon `other`'s boundary somewhere other than along
    /// an identical edge or at its own endpoints: it overlaps one of
    /// `other`'s edges collinearly without having the same two endpoints,
    /// or a vertex of `other` lies in its interior. The coverage is not
    /// noded there, so the shared boundary does not match edge for edge.
    UnmatchedEdge {
        /// The polygon the edge belongs to.
        polygon: usize,
        /// The ring of the edge (`0` exterior, `1..` holes).
        ring: usize,
        /// The edge's index in its ring.
        edge: usize,
        /// The polygon whose boundary it meets.
        other: usize,
    },
    /// Two boundary edges of different polygons lie within the gap width of
    /// each other, are not the same edge, and do not touch: a gap (or a
    /// gore) narrower than the width. Reported only by
    /// [`coverage_validate_with_gap`]. The first edge is the one earlier in
    /// `(polygon, ring, edge)` order.
    NarrowGap {
        /// The first edge's polygon.
        polygon: usize,
        /// The first edge's ring.
        ring: usize,
        /// The first edge's index in its ring.
        edge: usize,
        /// The second edge's polygon.
        other: usize,
        /// The second edge's ring.
        other_ring: usize,
        /// The second edge's index in its ring.
        other_edge: usize,
    },
}

/// Every reason `polys` is not a valid polygonal coverage, sorted; empty
/// when it is one.
///
/// A coverage is valid when every ring is simple
/// ([`CoverageIssue::NonSimpleRing`]), no two interiors meet
/// ([`CoverageIssue::Overlap`]), and every shared stretch of boundary is made
/// of identical edges ([`CoverageIssue::UnmatchedEdge`]), so adjacent
/// polygons share their vertices bit for bit. Gaps between polygons are allowed; use
/// [`coverage_validate_with_gap`] to find narrow ones.
///
/// Exact: overlaps by `relate`, edge identity by coordinate equality,
/// collinearity and vertex-on-edge by `orient2d` on stored coordinates.
/// Candidate pairs come from a box index (no false negatives).
///
/// ```
/// use geometres::ops::coverage_validate;
/// use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon};
///
/// let square = |x: f64| -> Result<Polygon, geometres::GeometryError> {
///     Polygon::new_simple(LinearRing::new(CoordSeq::new(vec![
///         Coord2d::new(x, 0.0), Coord2d::new(x + 1.0, 0.0),
///         Coord2d::new(x + 1.0, 1.0), Coord2d::new(x, 1.0),
///     ]))?)
/// };
/// // Two unit squares sharing the edge x = 1: a valid coverage.
/// assert!(coverage_validate(&[square(0.0)?, square(1.0)?]).is_empty());
/// # Ok::<(), geometres::GeometryError>(())
/// ```
#[must_use]
pub fn coverage_validate(polys: &[Polygon]) -> Vec<CoverageIssue> {
    validate::issues(polys, None)
}

/// [`coverage_validate`], plus [`CoverageIssue::NarrowGap`] for every pair
/// of boundary edges of different polygons that lie within `width` of each
/// other (distance `≤ width`, decided exactly), are not the same edge, and
/// do not touch.
///
/// The test is `ops::distance`'s exact segment-to-segment composition: the
/// edges do not touch, so their distance is attained at an endpoint of one
/// of them, and each endpoint-to-segment comparison is `arithmetica`'s exact
/// predicate. "This gap is narrower than `width`" is therefore never decided
/// in `f64`. Candidates come from each edge's box grown by `width` and
/// rounded outward: a box rounded inward could miss a pair, while one
/// rounded outward only adds candidates, which the exact test rejects.
///
/// Gap detection is finding narrow gaps of at most a maximum gap width, including
/// gaps separating disjoint coverage regions and narrow gores. `width = 0` reports
/// no gap (a distance of zero means the edges touch).
///
/// # Errors
///
/// [`GeometryError::InvalidTolerance`] if `width` is negative or not
/// finite.
pub fn coverage_validate_with_gap(
    polys: &[Polygon],
    width: f64,
) -> Result<Vec<CoverageIssue>, GeometryError> {
    crate::ops::distance::check_distance(width)?;
    Ok(validate::issues(polys, Some(width)))
}

// ===========================================================================
// Coverage union
// ===========================================================================

/// The union of a valid coverage, by removing shared edges. It takes no
/// overlay and no intersection, and adds no vertex.
///
/// Every ring edge is collected with its direction. In a valid coverage an
/// edge shared by two polygons appears exactly twice, in opposite
/// directions (exteriors run counterclockwise and holes clockwise, so the
/// interior is always on the left), and lies inside the union: it is
/// dropped. An edge appearing once is on the union's boundary. The boundary
/// edges are traced into rings with the union's interior on the left (the
/// module documentation), classified by the certified sign of their area,
/// and holes are assigned to exteriors by exact placement. The result is
/// built through `Polygon::new` and `MultiPolygon::new`.
///
/// Collinear vertices left on a dissolved boundary are **kept**: removing
/// them would be simplification, and every output vertex is an input
/// vertex, bit for bit.
///
/// # Errors
///
/// - [`GeometryError::InvalidCoverage`] carrying the first issue
///   [`coverage_validate`] finds, if the input is not a valid coverage. For
///   general input, use a constructive union.
/// - [`GeometryError::NonFiniteResult`] if a ring of the union has a signed
///   area beyond the range of `f64`. Every input polygon has a finite area,
///   but their union can exceed it (coordinates near 1e154 and beyond), and
///   no `Polygon` can carry such a ring.
pub fn coverage_union(polys: &[Polygon]) -> Result<MultiPolygon, GeometryError> {
    if let Some(&first) = validate::issues(polys, None).first() {
        return Err(GeometryError::InvalidCoverage { first });
    }
    union::dissolve(polys)
}

// ===========================================================================
// Polygonize noded linework
// ===========================================================================

/// What [`polygonize_noded`] returns: the faces of a noded line arrangement,
/// and the linework that bounds no face.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Polygonized {
    /// The bounded faces, each with the rings of the components nested in
    /// it as holes. Faces are interior-disjoint but may share edges, so they
    /// are not a `MultiPolygon`. Sorted by exterior vertex sequence.
    pub polygons: Vec<Polygon>,
    /// Edges with a free end, removed repeatedly until none is left, merged
    /// into maximal chains.
    pub dangles: MultiLineString,
    /// Edges with the same face on both sides (bridges), merged into
    /// maximal chains.
    pub cut_edges: MultiLineString,
    /// Closed rings that cannot form a valid polygon, as closed lines.
    /// Every traced loop of noded input is a simple ring of at least three
    /// non-collinear vertices whose holes are placed exactly, so
    /// `Polygon::new` accepts it, with one exception: a loop whose signed
    /// area is beyond the range of `f64` (coordinates near 1e154 and
    /// beyond), which no `Polygon` can carry. Such a loop is reported here
    /// rather than dropped.
    pub invalid_rings: MultiLineString,
}

/// The faces of the noded linework of `lines`, as polygons, with their
/// by-products: dangles, cut edges and invalid rings.
///
/// The linework is every segment of every line and ring in `lines`
/// (points are ignored); identical segments, in either direction, are one
/// edge, and **every vertex is a node**.
///
/// **The linework must be noded:** two segments may meet only at an endpoint of
/// both, with no crossing, no collinear overlap, and no endpoint in the
/// interior of another segment. That is checked exactly before anything is
/// traced.
///
/// Dangles are removed repeatedly; then cut edges (edges whose two sides
/// trace the same face); the remaining faces are traced with the face on
/// the left. A face's outer loop is its exterior; the loops of components
/// nested inside it are its holes, assigned by exact placement. The faces
/// of a nested component are output as well: a hole and the island in it
/// are both faces. Every output coordinate is an input coordinate.
///
/// # Errors
///
/// [`GeometryError::NotNoded`], naming the first offending pair of
/// segments, if the linework is not noded. Node the input first with a
/// constructive operation.
pub fn polygonize_noded(
    lines: &impl crate::predicates::Relatable,
) -> Result<Polygonized, GeometryError> {
    polygonize::run(crate::predicates::borrowed(lines))
}
