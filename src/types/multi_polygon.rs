// src/types/multi_polygon.rs

//! Multi-polygon geometry (`MultiPolygon`).

use crate::error::GeometryError;
use crate::exact::kernel::{RingPlacement, place_ring};
use crate::types::polygon::Polygon;

/// A collection of zero or more `Polygon` values whose parts meet at most at
/// finitely many points.
///
/// # Invariants
///
/// - Every element is a valid `Polygon`, enforced at `Polygon` construction.
/// - **No two parts overlap**: their interiors are disjoint and their
///   boundaries share no stretch of positive length, so they meet at most at
///   finitely many isolated points (OGC). Enforced by [`MultiPolygon::new`].
/// - The collection may be empty.
///
/// As everywhere in `types`, **ring simplicity is not checked** (the
/// construction contract in [`types`](crate::types)):
/// [`ops::validity_issues`](crate::ops::validity_issues) reports a
/// self-touching ring in any part.
///
/// Every multipart operation rests on the second invariant: the DE-9IM
/// matrix (which analyzes all parts' rings as one ring set) and overlay's
/// winding argument both assume the parts do not overlap. A value that broke
/// it could not be analyzed correctly by any of them, so it cannot be built.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiPolygon {
    polygons: Vec<Polygon>,
}

/// `(min_x, min_y, max_x, max_y)` of a polygon's exterior.
fn envelope(p: &Polygon) -> (f64, f64, f64, f64) {
    p.exterior().bbox()
}

/// `true` if part `q` lies inside one of `p`'s holes (touching it at points
/// at most).
pub(crate) fn in_a_hole_of(q: &Polygon, p: &Polygon) -> bool {
    p.holes()
        .iter()
        .any(|h| place_ring(q.exterior(), h) == RingPlacement::Inside)
}

/// `true` if two valid polygons meet at most at finitely many points.
///
/// Decided by exact ring placement (`exact::place_ring`, exact for every
/// way two rings can meet):
///
/// - the exteriors are mutually `Outside`: the parts are side by side,
///   touching at points at most; or
/// - one exterior is `Inside` the other and lies `Inside` one of the other's
///   holes: an island in a hole, touching it at points at most.
///
/// Anything else (crossing rings, a shared stretch, one part's exterior
/// inside the other's body) is an overlap. The predicate is equivalent to
/// the DE-9IM criterion `II = F ∧ BB ≠ 1` (interiors disjoint, no shared
/// boundary stretch); `parts_valid_agrees_with_the_de9im_criterion` checks
/// the equivalence on every pair of the polygon corpora. It is stated in
/// ring terms because `types` sits below `topology`.
pub(crate) fn parts_meet_only_at_points(p: &Polygon, q: &Polygon) -> bool {
    match (
        place_ring(q.exterior(), p.exterior()),
        place_ring(p.exterior(), q.exterior()),
    ) {
        (RingPlacement::Outside, RingPlacement::Outside) => true,
        (RingPlacement::Inside, RingPlacement::Outside) => in_a_hole_of(q, p),
        (RingPlacement::Outside, RingPlacement::Inside) => in_a_hole_of(p, q),
        _ => false,
    }
}

/// The first pair of parts that overlap, or `None` if every pair meets at
/// most at finitely many points. This is *the* check behind
/// [`MultiPolygon::new`].
///
/// A closed bounding-box prefilter first (exact: it compares stored
/// coordinates, and envelopes that only touch still go to the exact test,
/// since touching parts may share a point), then
/// [`parts_meet_only_at_points`].
pub(crate) fn first_overlapping_parts(polygons: &[Polygon]) -> Option<(usize, usize)> {
    let refs: Vec<&Polygon> = polygons.iter().collect();
    first_overlapping_parts_ref(&refs)
}

/// As [`first_overlapping_parts`], over borrowed parts: the form the
/// collection decomposition uses, whose polygon leaves are
/// scattered across a `GeometryCollection` rather than held in one vector.
/// One body serves both, so the restriction a collection's polygonal leaves
/// must satisfy is *the* multipolygon part rule, not a reimplementation.
pub(crate) fn first_overlapping_parts_ref(polygons: &[&Polygon]) -> Option<(usize, usize)> {
    let envs: Vec<(f64, f64, f64, f64)> = polygons.iter().map(|p| envelope(p)).collect();
    for i in 0..polygons.len() {
        for j in (i + 1)..polygons.len() {
            let (a, b) = (envs[i], envs[j]);
            let boxes_meet = a.0 <= b.2 && b.0 <= a.2 && a.1 <= b.3 && b.1 <= a.3;
            if boxes_meet && !parts_meet_only_at_points(polygons[i], polygons[j]) {
                return Some((i, j));
            }
        }
    }
    None
}

impl MultiPolygon {
    /// A multipolygon from parts already known valid together: a
    /// validated multipolygon's parts, permuted, each with its rings
    /// rotated. The part conditions are pairwise and on point sets, so
    /// neither permutation nor rotation affects them. Crate-private.
    pub(crate) fn from_valid(polygons: Vec<Polygon>) -> Self {
        Self { polygons }
    }

    /// As [`MultiPolygon::new`], with the offending pair of part indices as
    /// the error itself.
    ///
    /// The check is the one `new` runs, exactly: a closed bounding-box
    /// prefilter, then exact ring placement. This entry point exists for a
    /// caller that reports the pair in its own terms, such as indices of its
    /// own, and so needs the pair without a `GeometryError` to unpack. No
    /// multipolygon is ever built from unchecked parts.
    ///
    /// # Errors
    ///
    /// `(first, second)`, the first pair of parts whose interiors overlap or
    /// whose boundaries share a stretch, in the order `new` reports it as
    /// [`GeometryError::OverlappingParts`].
    pub fn from_parts(polygons: Vec<Polygon>) -> Result<Self, (usize, usize)> {
        match first_overlapping_parts(&polygons) {
            Some(pair) => Err(pair),
            None => Ok(Self { polygons }),
        }
    }

    /// Construct a `MultiPolygon`, checking that no two parts overlap.
    ///
    /// Every pair of parts is checked: a closed bounding-box prefilter, then
    /// the exact ring-placement predicate (equivalent to the DE-9IM criterion
    /// `II = F ∧ BB ≠ 1`). `O(n²)` pairs in the part count, each `O(e·f)`
    /// in the two parts' edge counts when the boxes meet.
    ///
    /// An empty vector is accepted; it represents an empty geometry.
    ///
    /// # Errors
    ///
    /// [`GeometryError::OverlappingParts`] naming the first pair of parts
    /// whose interiors overlap or whose boundaries share a stretch.
    pub fn new(polygons: Vec<Polygon>) -> Result<Self, GeometryError> {
        match first_overlapping_parts(&polygons) {
            Some((first, second)) => Err(GeometryError::OverlappingParts { first, second }),
            None => Ok(Self { polygons }),
        }
    }

    /// The empty multipolygon.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            polygons: Vec::new(),
        }
    }

    /// The parts of `operand` that `keep` accepts, in their original order,
    /// as a multipolygon.
    ///
    /// No validation runs, and none is needed: a valid operand's parts meet
    /// at most at finitely many points pairwise (a `Polygon` has one part; a
    /// `MultiPolygon` passed its constructor), and every pair of a subset is
    /// a pair of the whole. So the result is valid by construction, without
    /// the `O(n²)` check of [`MultiPolygon::new`]. The order is kept, so a
    /// computation that depends on part order sees the retained parts in
    /// the same relative order as in `operand`.
    ///
    /// `keep` is called once per part, with the part's index in `operand`
    /// and the part itself. The retained parts are cloned.
    ///
    /// ```
    /// use geometres::types::{Coord2d, CoordSeq, LinearRing, MultiPolygon, Polygon};
    ///
    /// let square = |x: f64| {
    ///     let ring = LinearRing::new(CoordSeq::new(vec![
    ///         Coord2d::new(x, 0.0), Coord2d::new(x + 1.0, 0.0),
    ///         Coord2d::new(x + 1.0, 1.0), Coord2d::new(x, 1.0),
    ///     ]))?;
    ///     Polygon::new_simple(ring)
    /// };
    /// let mp = MultiPolygon::new(vec![square(0.0)?, square(2.0)?, square(4.0)?])?;
    /// let odd = MultiPolygon::from_retained_parts(&mp, |i, _| i % 2 == 0);
    /// assert_eq!(odd.polygons(), [mp.polygons()[0].clone(), mp.polygons()[2].clone()]);
    /// # Ok::<(), geometres::GeometryError>(())
    /// ```
    #[must_use]
    pub fn from_retained_parts(
        operand: &impl crate::types::Polygonal,
        mut keep: impl FnMut(usize, &Polygon) -> bool,
    ) -> Self {
        let polygons = operand
            .parts()
            .iter()
            .enumerate()
            .filter(|(i, p)| keep(*i, p))
            .map(|(_, p)| p.clone())
            .collect();
        Self { polygons }
    }

    /// Return a slice over the constituent polygons.
    pub fn polygons(&self) -> &[Polygon] {
        &self.polygons
    }

    /// Return the number of polygons in the collection.
    pub fn num_polygons(&self) -> usize {
        self.polygons.len()
    }

    /// Return `true` if the collection contains no polygons.
    pub fn is_empty(&self) -> bool {
        self.polygons.is_empty()
    }

    /// Return the single `Polygon` if this collection holds exactly one, else
    /// `None` (including the empty case).
    ///
    /// This is the ergonomic path for operations such as outward offset that
    /// return a `MultiPolygon` for generality but produce a single polygon in
    /// the common case; the caller unwraps the simple result without matching
    /// on a vector length.
    #[must_use]
    pub fn try_into_polygon(self) -> Option<Polygon> {
        if self.polygons.len() == 1 {
            self.polygons.into_iter().next()
        } else {
            None
        }
    }
}

/// One polygon is always a valid multipolygon.
impl From<Polygon> for MultiPolygon {
    fn from(p: Polygon) -> Self {
        Self { polygons: vec![p] }
    }
}

impl TryFrom<MultiPolygon> for Polygon {
    type Error = crate::error::GeometryError;

    /// Succeeds only when the collection contains exactly one polygon;
    /// otherwise returns [`GeometryError::NotASinglePolygon`](crate::error::GeometryError::NotASinglePolygon).
    fn try_from(mp: MultiPolygon) -> Result<Polygon, Self::Error> {
        mp.try_into_polygon()
            .ok_or(crate::error::GeometryError::NotASinglePolygon)
    }
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used)]
// test-lint-posture: see the workspace lint comment in Cargo.toml
#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::range_plus_one,
    clippy::no_effect_underscore_binding
)]
mod tests {
    use super::*;
    use crate::types::coord::{Coord2d, CoordSeq};
    use crate::types::linear_ring::LinearRing;

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    fn seq(pts: &[(f64, f64)]) -> CoordSeq {
        CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())
    }

    fn ccw_square(ox: f64, oy: f64, side: f64) -> Polygon {
        let ring = LinearRing::new(seq(&[
            (ox, oy),
            (ox + side, oy),
            (ox + side, oy + side),
            (ox, oy + side),
        ]))
        .unwrap();
        Polygon::new_simple(ring).unwrap()
    }

    #[test]
    fn empty_multi_polygon_is_valid() {
        let mp = MultiPolygon::new(vec![]).unwrap();
        assert!(mp.is_empty());
        assert_eq!(MultiPolygon::empty(), mp);
        assert_eq!(mp.num_polygons(), 0);
        assert_eq!(mp.polygons(), &[]);
    }

    #[test]
    fn single_polygon_accessors() {
        let p = ccw_square(0.0, 0.0, 1.0);
        let mp = MultiPolygon::new(vec![p.clone()]).unwrap();
        assert_eq!(MultiPolygon::from(p.clone()), mp);
        assert!(!mp.is_empty());
        assert_eq!(mp.num_polygons(), 1);
        assert_eq!(mp.polygons(), &[p]);
    }

    #[test]
    fn multiple_polygons_round_trip() {
        let p1 = ccw_square(0.0, 0.0, 1.0);
        let p2 = ccw_square(2.0, 0.0, 1.0);
        let mp = MultiPolygon::new(vec![p1.clone(), p2.clone()]).unwrap();
        assert_eq!(mp.num_polygons(), 2);
        assert_eq!(mp.polygons()[0], p1);
        assert_eq!(mp.polygons()[1], p2);
    }

    // ===================================================================
    // from_retained_parts
    // ===================================================================

    /// Every subset of a valid multipolygon's parts, kept in order, equals
    /// what `MultiPolygon::new` builds from the same parts: same parts,
    /// same order, and `new` accepts it.
    #[test]
    fn retained_parts_are_every_subset_in_order() {
        // Four parts, two of which touch at a corner (a point contact is
        // valid), so subsets include touching and separated pairs.
        let parts = vec![
            ccw_square(0.0, 0.0, 1.0),
            ccw_square(1.0, 1.0, 1.0),
            ccw_square(3.0, 0.0, 1.0),
            ccw_square(5.0, 5.0, 2.0),
        ];
        let mp = MultiPolygon::new(parts.clone()).unwrap();
        for mask in 0u32..16 {
            let got = MultiPolygon::from_retained_parts(&mp, |i, _| mask & (1 << i) != 0);
            let want: Vec<Polygon> = (0..4)
                .filter(|i| mask & (1 << i) != 0)
                .map(|i| parts[i].clone())
                .collect();
            assert_eq!(got.polygons(), &want[..], "subset {mask:04b}");
            assert_eq!(MultiPolygon::new(want).unwrap(), got);
        }
    }

    #[test]
    fn retained_parts_of_a_polygon_are_it_or_nothing() {
        let p = ccw_square(0.0, 0.0, 1.0);
        assert_eq!(
            MultiPolygon::from_retained_parts(&p, |_, _| true),
            MultiPolygon::from(p.clone())
        );
        assert!(MultiPolygon::from_retained_parts(&p, |_, _| false).is_empty());
    }

    #[test]
    fn retained_parts_see_each_index_once_in_order() {
        let mp = MultiPolygon::new(vec![
            ccw_square(0.0, 0.0, 1.0),
            ccw_square(2.0, 0.0, 1.0),
            ccw_square(4.0, 0.0, 1.0),
        ])
        .unwrap();
        let mut seen = Vec::new();
        let _ = MultiPolygon::from_retained_parts(&mp, |i, p| {
            seen.push((i, p.clone()));
            false
        });
        assert_eq!(seen.len(), 3);
        for (k, (i, p)) in seen.iter().enumerate() {
            assert_eq!((*i, p), (k, &mp.polygons()[k]));
        }
    }
}
