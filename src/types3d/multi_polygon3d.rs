// src/types3d/multi_polygon3d.rs

//! Multi-polygon geometry in 3D (`MultiPolygon3d`).
//!
//! The type extends the 2D multipart areal rule, and its topological
//! boundary, to the XY projection rather than inventing a rule of its own.

use crate::error::GeometryError;
use crate::types::{Polygon, first_overlapping_parts};
use crate::types3d::Polygon3d;

/// A collection of zero or more `Polygon3d` values whose XY projections meet
/// at most at finitely many points.
///
/// # The rule is the 2D rule on the projection
///
/// Two parts may not overlap **in plan**: their XY projections must have
/// disjoint interiors and share no boundary stretch of positive length, so
/// they meet at most at isolated points. This is exactly
/// [`MultiPolygon`](crate::types::MultiPolygon)'s invariant applied to each
/// part's projection, and it is what every topology operation in this crate
/// requires, since all of them work on the XY projection alone (see
/// [`types3d`](crate::types3d)).
///
/// It is **not a rule about 3D solids or surfaces**. Two parts at different
/// heights whose projections overlap are refused, even though as surfaces in
/// space they are disjoint. A caller with stacked geometry keeps the layers
/// as separate values.
///
/// # Invariants
///
/// - Every element is a valid `Polygon3d`, enforced at `Polygon3d`
///   construction, which also makes every ring's XY projection a valid
///   `LinearRing`.
/// - **No two parts overlap in plan**, decided exactly on the projections
///   by the same ring-placement predicate `MultiPolygon::new` uses. Enforced
///   by [`MultiPolygon3d::new`].
/// - The collection may be empty.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiPolygon3d {
    polygons: Vec<Polygon3d>,
}

impl MultiPolygon3d {
    /// Construct a `MultiPolygon3d`, checking that no two parts overlap in
    /// plan.
    ///
    /// Each part is projected to a 2D `Polygon` (every coordinate of the
    /// projection is an input `x`, `y` pair; nothing is computed), and the
    /// projections are checked pairwise with the exact ring-placement
    /// predicate behind `MultiPolygon::new`: a closed bounding-box prefilter,
    /// then exact placement. `O(n²)` pairs in the part count.
    ///
    /// An empty vector is accepted; it represents an empty geometry.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::OverlappingParts`] naming the first pair of parts
    ///   whose projections have overlapping interiors or share a boundary
    ///   stretch.
    /// - Any error `Polygon::new` raises on a projection. For a valid
    ///   `Polygon3d` there is none (`Polygon3d::new` decides winding by the
    ///   same certified sum, degeneracy by the same exact test, and hole
    ///   placement on the same projection), but the error is
    ///   propagated rather than assumed away.
    pub fn new(polygons: Vec<Polygon3d>) -> Result<Self, GeometryError> {
        let projected = polygons
            .iter()
            .map(Polygon3d::xy_projection)
            .collect::<Result<Vec<Polygon>, _>>()?;
        match first_overlapping_parts(&projected) {
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

    /// Return a slice over the constituent polygons.
    pub fn polygons(&self) -> &[Polygon3d] {
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
}

/// One polygon is always a valid multipolygon.
impl From<Polygon3d> for MultiPolygon3d {
    fn from(p: Polygon3d) -> Self {
        Self { polygons: vec![p] }
    }
}

impl TryFrom<MultiPolygon3d> for Polygon3d {
    type Error = GeometryError;

    /// Succeeds only when the collection contains exactly one polygon;
    /// otherwise returns [`GeometryError::NotASinglePolygon`].
    fn try_from(mp: MultiPolygon3d) -> Result<Polygon3d, Self::Error> {
        if mp.polygons.len() == 1 {
            mp.polygons
                .into_iter()
                .next()
                .ok_or(GeometryError::NotASinglePolygon)
        } else {
            Err(GeometryError::NotASinglePolygon)
        }
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
    use crate::types3d::{Coord3d, CoordSeq3d, LinearRing3d};

    fn ring(pts: &[(f64, f64, f64)]) -> LinearRing3d {
        LinearRing3d::new(CoordSeq3d::new(
            pts.iter().map(|&(x, y, z)| Coord3d::new(x, y, z)).collect(),
        ))
        .unwrap()
    }

    /// Counterclockwise square of side `side` at `(ox, oy)`, at height `z`.
    fn square(ox: f64, oy: f64, side: f64, z: f64) -> Polygon3d {
        Polygon3d::new_simple(ring(&[
            (ox, oy, z),
            (ox + side, oy, z),
            (ox + side, oy + side, z),
            (ox, oy + side, z),
        ]))
        .unwrap()
    }

    #[test]
    fn empty_is_valid() {
        let mp = MultiPolygon3d::new(vec![]).unwrap();
        assert!(mp.is_empty());
        assert_eq!(MultiPolygon3d::empty(), mp);
        assert_eq!(mp.num_polygons(), 0);
        assert_eq!(mp.polygons(), &[]);
    }

    #[test]
    fn disjoint_parts_are_accepted() {
        let a = square(0.0, 0.0, 1.0, 0.0);
        let b = square(2.0, 0.0, 1.0, 5.0);
        let mp = MultiPolygon3d::new(vec![a.clone(), b.clone()]).unwrap();
        assert_eq!(mp.num_polygons(), 2);
        assert_eq!(mp.polygons(), &[a, b]);
    }

    #[test]
    fn parts_touching_at_a_vertex_are_accepted() {
        let a = square(0.0, 0.0, 1.0, 0.0);
        let b = square(1.0, 1.0, 1.0, 0.0);
        assert!(MultiPolygon3d::new(vec![a, b]).is_ok());
    }

    #[test]
    fn parts_sharing_an_edge_are_refused() {
        let a = square(0.0, 0.0, 1.0, 0.0);
        let b = square(1.0, 0.0, 1.0, 0.0);
        assert_eq!(
            MultiPolygon3d::new(vec![a, b]).unwrap_err(),
            GeometryError::OverlappingParts {
                first: 0,
                second: 1
            }
        );
    }

    #[test]
    fn overlapping_parts_are_refused_with_the_first_pair() {
        let a = square(0.0, 0.0, 1.0, 0.0);
        let far = square(5.0, 5.0, 1.0, 0.0);
        let b = square(0.5, 0.5, 1.0, 0.0);
        assert_eq!(
            MultiPolygon3d::new(vec![a, far, b]).unwrap_err(),
            GeometryError::OverlappingParts {
                first: 0,
                second: 2
            }
        );
    }

    /// The rule is on the projection: stacked parts at different heights
    /// overlap in plan and are refused.
    #[test]
    fn parts_overlapping_in_plan_at_different_heights_are_refused() {
        let ground = square(0.0, 0.0, 1.0, 0.0);
        let roof = square(0.0, 0.0, 1.0, 10.0);
        assert_eq!(
            MultiPolygon3d::new(vec![ground, roof]).unwrap_err(),
            GeometryError::OverlappingParts {
                first: 0,
                second: 1
            }
        );
    }

    #[test]
    fn an_island_in_a_hole_is_accepted() {
        let hole = ring(&[
            (0.2, 0.2, 0.0),
            (0.2, 0.8, 0.0),
            (0.8, 0.8, 0.0),
            (0.8, 0.2, 0.0),
        ]);
        let outer = Polygon3d::new(
            ring(&[
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 1.0, 0.0),
            ]),
            vec![hole],
        )
        .unwrap();
        let island = square(0.4, 0.4, 0.2, 3.0);
        assert!(MultiPolygon3d::new(vec![outer, island]).is_ok());
    }

    #[test]
    fn from_and_try_from_round_trip_one_part() {
        let p = square(0.0, 0.0, 1.0, 1.0);
        let mp = MultiPolygon3d::from(p.clone());
        assert_eq!(mp.num_polygons(), 1);
        assert_eq!(Polygon3d::try_from(mp).unwrap(), p);
    }

    #[test]
    fn try_from_refuses_zero_and_two_parts() {
        assert_eq!(
            Polygon3d::try_from(MultiPolygon3d::empty()).unwrap_err(),
            GeometryError::NotASinglePolygon
        );
        let two = MultiPolygon3d::new(vec![square(0.0, 0.0, 1.0, 0.0), square(3.0, 0.0, 1.0, 0.0)])
            .unwrap();
        assert_eq!(
            Polygon3d::try_from(two).unwrap_err(),
            GeometryError::NotASinglePolygon
        );
    }
}
