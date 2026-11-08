// src/ops/mod.rs

//! Non-predicate geometric operations: area, length, validation, simplicity
//! and vouched simplification ([`simplify`]).
//!
//! Constructive operations (offset, centroid, interior point and
//! tessellation, among others) are not here, because each returns
//! coordinates that are not input coordinates.
//!
//! # Exactness contract
//!
//! These operations return `f64` measurements, not geometric verdicts, so they
//! are approximate by nature. What each guarantees:
//!
//! | Function | Guarantee |
//! |---|---|
//! | [`signed_area_ring`], [`signed_area`], [`area`] | Computed as if in twice the working precision (Ogita-Rump-Oishi `Dot2`); see `area.rs` for the bound |
//! | [`length`], [`ring_perimeter`], [`perimeter`], [`edge_length`], [`polygon_length`] | Compensated sum of `hypot` edge lengths; no proven bound, accuracy measured |
//! | [`collection_area`], [`collection_length`] | Sum over the collection's elements (overlaps count twice); certified zero exactly when there is nothing to measure; bound composed from the terms'; see `ops::measure` |
//! | [`simplify`] (methods on `Polygon`, `MultiPolygon`, `LineString`, `MultiLineString`) | Vertex choice is `f64`; the result is **vouched exactly** or refused; see `ops::simplify` |
//! | [`simplify_vw()`] | Vertex choice is `f64` (effective areas, a selection rule); the result is **vouched exactly** by the constructors and `validity_issues`, or refused; see `ops::simplify_vw` |
//! | [`simplify_preserve_topology`] (and its `_unvouched` form) | Vertex choice is `f64`; every topological check is **exact**, and the result is vouched exactly; see [`simplify_topology`] |
//! | [`validate()`], [`validate_line`], [`validate_multipolygon`], [`is_valid`], [`validity_issues`] | **Exact**: certified area sign for orientation, exact self-contact and ring placement by index, exact touch graph for a disconnected interior; every issue reported, none constructed; see `ops::validate` |
//! | [`is_simple`] | **Exact**: OGC simplicity per type by exact orientation and coordinate comparison; see `ops::simple` |
//! | [`split_pinches`] | **Exact**: vertices compared by stored bits, no coordinate constructed; the loops are runs of the ring's own edges; see `ops::pinch` |
//! | [`coverage_validate`], [`coverage_validate_with_gap`] | **Exact**: overlaps by `relate`, edge identity by coordinate equality, collinearity and vertex-on-edge by `orient2d`, gap width by `ops::distance`'s exact segment test; see `ops::coverage` |
//! | [`coverage_simplify()`], [`coverage_simplify_inner`] | Vertex choice is `f64` (effective areas); every removal is decided by **exact** topology checks, shared edges are simplified once, and the result is **checked exactly** by `coverage_validate`; see `ops::coverage_simplify` |
//! | [`coverage_union`], [`polygonize_noded`] | **Exact**, constructing nothing: every output vertex is an input vertex; faces are traced by an exact angular order and classified by certified area signs; see `ops::coverage` |
//! | [`dwithin`] | **Exact**: the intersection test and every element-pair comparison are exact predicates; see `ops::distance` |
//! | [`delaunay()`] | **Exact and canonical**: every decision is `orient2d` or `incircle` on input coordinates; cocircular ties are cut by a fixed fan rule, so the result depends only on the point set; see `ops::delaunay` |
//! | [`concave_hull()`] | The threshold is a **selection rule** in `f64`; the claim (a valid polygon of input vertices containing every input point) is **certified exactly** or refused; see `ops::concave_hull` |
//! | [`distance()`] | Exactly `0.0` if and only if the operands intersect; otherwise positive and within `max(64·u·M, 2⁻¹⁰⁷⁴)` of the true distance; see `ops::distance` |
//! | [`hausdorff_within`], [`frechet_within`] | **Exact**: every point-to-element and point-to-point comparison is an exact `arithmetica` predicate; see `ops::hausdorff`, `ops::frechet` |
//! | [`hausdorff_distance`], [`frechet_distance`] | Exactly `0.0` when, and only when, the exact zero test holds; otherwise within `64·u·M` of the true value (Hausdorff: or `2⁻¹⁰⁷⁴`, if larger); see `ops::hausdorff`, `ops::frechet` |
//!
//! None of these functions panics. Inputs that make an intermediate result
//! overflow produce a non-finite `f64`; a finite wrong value is not produced
//! by overflow.

mod area;
pub(crate) mod basic;
#[cfg(test)]
mod basic_tests;
mod concave_hull;
#[cfg(test)]
mod concave_hull_tests;
pub mod coverage;
pub mod coverage_simplify;
#[cfg(test)]
mod coverage_simplify_tests;
mod delaunay;
#[cfg(test)]
mod delaunay_tests;
pub mod distance;
pub mod frechet;
#[cfg(test)]
mod frechet_tests;
pub mod hausdorff;
#[cfg(test)]
mod hausdorff_tests;
mod hull;
mod length;
mod linework;
mod measure;
pub mod pinch;
mod shared;
mod simple;
pub mod simplify;
pub mod simplify_topology;
pub mod simplify_vw;
#[cfg(test)]
mod simplify_vw_tests;
pub mod validate;

pub use area::{area, signed_area, signed_area_ring};
pub use basic::{Reversible, boundary, normalize, reverse, unique_points};
pub use concave_hull::{ConcaveHullParameter, concave_hull};
pub use coverage::{
    CoverageIssue, Polygonized, coverage_union, coverage_validate, coverage_validate_with_gap,
    polygonize_noded,
};
pub use coverage_simplify::{coverage_simplify, coverage_simplify_inner};
pub use delaunay::{Triangulation, delaunay};
pub use distance::{distance, dwithin};
pub use frechet::{frechet_distance, frechet_within};
pub use hausdorff::{hausdorff_distance, hausdorff_within};
pub use hull::convex_hull;
pub use linework::{line_merge, line_merge_directed, line_sequence, split_by_point};
pub use shared::{SharedPaths, shared_paths};
// The one certified sum, for `Polygon3d`'s winding on projected edges.
#[cfg(feature = "3d")]
pub(crate) use area::certified_twice_area_edges;
pub use length::{edge_length, length, perimeter, ring_perimeter};
pub use measure::{collection_area, collection_length, polygon_length};
pub use pinch::split_pinches;
pub use simple::is_simple;
pub use simplify_topology::{
    TopologySimplifiable, simplify_preserve_topology, simplify_preserve_topology_unvouched,
};
pub use simplify_vw::{Simplifiable, simplify_vw};
pub use validate::{
    LocatedIssue, ValidityIssue, is_valid, validate, validate_line, validate_multipolygon,
    validity_issues,
};

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
    use crate::GeometryError;
    use crate::types::{Coord2d, CoordSeq, LineString, LinearRing, Polygon};

    fn seq(pts: &[(f64, f64)]) -> CoordSeq {
        CoordSeq::new(pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect())
    }

    fn ring(pts: &[(f64, f64)]) -> LinearRing {
        LinearRing::new(seq(pts)).unwrap()
    }

    fn poly(pts: &[(f64, f64)]) -> Polygon {
        Polygon::new_simple(ring(pts)).unwrap()
    }

    const UNIT_CCW: [(f64, f64); 4] = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    const UNIT_CW: [(f64, f64); 4] = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)];

    // signed_area, signed_area (CW), area.
    #[test]
    fn signed_area_of_a_ccw_unit_square_is_one() {
        assert_eq!(signed_area(&poly(&UNIT_CCW)), 1.0);
        assert_eq!(signed_area_ring(&ring(&UNIT_CCW)), 1.0);
    }

    #[test]
    fn signed_area_of_a_cw_unit_square_is_minus_one() {
        // A CW ring is not a valid exterior, so use the ring function.
        assert_eq!(signed_area_ring(&ring(&UNIT_CW)), -1.0);
    }

    #[test]
    fn area_is_one_regardless_of_winding() {
        assert_eq!(area(&poly(&UNIT_CCW)), 1.0);
    }

    #[test]
    fn length_of_a_unit_segment_is_one() {
        let l = LineString::new(seq(&[(0.0, 0.0), (1.0, 0.0)])).unwrap();
        assert_eq!(length(&l), 1.0);
    }

    #[test]
    fn length_and_perimeter_of_known_shapes() {
        let l = LineString::new(seq(&[(0.0, 0.0), (3.0, 4.0), (3.0, 10.0)])).unwrap();
        assert_eq!(length(&l), 11.0);
        assert_eq!(perimeter(&poly(&UNIT_CCW)), 4.0);
        assert_eq!(ring_perimeter(&ring(&UNIT_CW)), 4.0);
    }

    #[test]
    fn a_hole_that_exceeds_the_exterior_is_rejected_at_construction() {
        // The hole surrounding the exterior is `RingPlacement::Outside` from
        // the hole's point of view (every one of its vertices is outside the
        // unit square), and `Polygon::new` refuses it before any operation
        // can be asked.
        let big_hole = ring(&[(-5.0, -5.0), (-5.0, 5.0), (5.0, 5.0), (5.0, -5.0)]);
        assert_eq!(
            Polygon::new(ring(&UNIT_CCW), vec![big_hole]).unwrap_err(),
            GeometryError::HoleNotInsideExterior { hole_index: 0 }
        );
    }

    #[test]
    fn area_far_from_the_origin_is_exact_where_the_plain_sum_is_not() {
        // A unit square at (1e8, 1e8): products are ~1e16, where one ulp is 2,
        // so the plain shoelace sum cannot represent an area of 1.
        let o = 1e8;
        let pts = [(o, o), (o + 1.0, o), (o + 1.0, o + 1.0), (o, o + 1.0)];
        let r = ring(&pts);
        assert_eq!(signed_area_ring(&r), 1.0);
        let mut plain = 0.0_f64;
        for (a, b) in r.edges() {
            plain += a.x * b.y;
            plain -= b.x * a.y;
        }
        assert_ne!(
            plain * 0.5,
            1.0,
            "the plain shoelace is expected to lose this"
        );
        assert_eq!(r.signed_area(), 1.0, "LinearRing::signed_area is exact");
    }
}
