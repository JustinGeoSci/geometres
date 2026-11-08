// src/topology/point_analysis.rs

//! Exact DE-9IM computation for zero-dimensional geometry (`Point`,
//! `MultiPoint`, and empty geometry) against any other geometry.
//!
//! A finite set of points `P` has an interior equal to `P` itself, an empty
//! boundary (OGC: the boundary of a zero-dimensional geometry is empty), and
//! the whole plane minus `P` as its exterior. Against a target `T` that gives
//! the matrix directly from point locations:
//!
//! | Cell | Rule |
//! |---|---|
//! | `II`, `IB`, `IE` | `Point` if some `p ∈ P` is in the interior, boundary, or exterior of `T`; else `Empty` |
//! | `BI`, `BB`, `BE` | always `Empty` (no boundary) |
//! | `EI` | dimension of `interior(T)` unless `P` can cover it: `Line` for lines, `Area` for polygons, and for a point-set target `Point` if some target point is not in `P` |
//! | `EB` | `Line` for a polygon target (a finite set cannot cover a curve); for a line target `Point` if some line boundary point is not in `P`; `Empty` for a point-set target |
//! | `EE` | always `Area` (the complement of a closed set of dimension ≤ 2 is open and dense) |
//!
//! An empty geometry is the empty point set. The same table applies and gives
//! `FFFFFF…2` rows, which is the OGC result for an empty geometry of any
//! type: the dimension of an empty geometry never appears in a DE-9IM cell.
//!
//! # Exactness
//!
//! Locating a point against a polygon uses [`locate`] (exact crossing number
//! with `arithmetica::orient2d` for every edge decision); against lines it
//! uses [`LineSet::locate`] (exact `on_segment`). Equality of point
//! coordinates is `==` on the original `f64` values. Nothing is constructed
//! or rounded, so a point exactly on a boundary is `Boundary` and a point one
//! ulp off it is not.

use crate::topology::areal::AreaOperand;
#[cfg(test)]
use crate::topology::areal::Areal;
use crate::topology::de9im::IntersectionMatrix;
use crate::topology::label::{Dimension, Location};
use crate::topology::line_analysis::LineSet;
use crate::types::Coord2d;

/// The geometry a point set is compared against.
pub(crate) enum Target<'a> {
    /// Another finite set of points (possibly empty).
    Points(&'a [Coord2d]),
    /// A non-empty linear geometry.
    Lines(&'a LineSet<'a>),
    /// A polygon or a valid multipolygon, indexed or not.
    Area(&'a dyn AreaOperand),
}

impl Target<'_> {
    /// Exact location of `p` relative to the target.
    fn locate(&self, p: Coord2d) -> Location {
        match self {
            Target::Points(q) => {
                if q.contains(&p) {
                    Location::Interior
                } else {
                    Location::Exterior
                }
            }
            Target::Lines(l) => l.locate(p),
            Target::Area(area) => area.locate(p),
        }
    }
}

/// Computes the DE-9IM matrix of the point set `pts` (A) against `target` (B).
///
/// `pts` may be empty. Duplicate points are permitted and have no effect.
/// All coordinates must be finite; the caller (`compute_geometry`) rejects
/// non-finite points before reaching this function, because `orient2d` would
/// otherwise report a non-finite input as a zero sign and the result would be
/// a silent wrong answer.
pub(crate) fn points_matrix(pts: &[Coord2d], target: &Target<'_>) -> IntersectionMatrix {
    let (mut ii, mut ib, mut ie) = (false, false, false);
    for &p in pts {
        match target.locate(p) {
            Location::Interior => ii = true,
            Location::Boundary => ib = true,
            Location::Exterior => ie = true,
        }
    }

    // Target interior not covered by a finite point set: a curve or a region
    // always has points outside `pts`; a point-set target needs a check.
    let ei = match target {
        Target::Points(q) => q.iter().any(|x| !pts.contains(x)),
        Target::Lines(_) | Target::Area(_) => true,
    };
    let ei_dim = match target {
        Target::Points(_) => Dimension::Point,
        Target::Lines(_) => Dimension::Line,
        Target::Area(_) => Dimension::Area,
    };
    let eb_dim = match target {
        Target::Points(_) => Dimension::Empty,
        Target::Lines(l) => {
            if l.boundary().iter().any(|b| !pts.contains(b)) {
                Dimension::Point
            } else {
                Dimension::Empty
            }
        }
        Target::Area(_) => Dimension::Line,
    };

    let dim = |flag: bool, d: Dimension| if flag { d } else { Dimension::Empty };
    IntersectionMatrix {
        interior_interior: dim(ii, Dimension::Point),
        interior_boundary: dim(ib, Dimension::Point),
        interior_exterior: dim(ie, Dimension::Point),
        boundary_interior: Dimension::Empty,
        boundary_boundary: Dimension::Empty,
        boundary_exterior: Dimension::Empty,
        exterior_interior: dim(ei, ei_dim),
        exterior_boundary: eb_dim,
        exterior_exterior: Dimension::Area,
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
    use crate::types::{CoordSeq, LineString, LinearRing, Polygon};

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    fn line(pts: &[(f64, f64)]) -> LineString {
        LineString::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
    }

    fn square() -> Polygon {
        let ring = LinearRing::new(CoordSeq::new(vec![
            c(0.0, 0.0),
            c(4.0, 0.0),
            c(4.0, 4.0),
            c(0.0, 4.0),
        ]))
        .unwrap();
        Polygon::new_simple(ring).unwrap()
    }

    #[test]
    fn point_against_polygon_each_location() {
        let sq = square();
        let m = |x, y| points_matrix(&[c(x, y)], &Target::Area(&Areal::single(&sq))).to_string();
        assert_eq!(m(2.0, 2.0), "0FFFFF212", "interior");
        assert_eq!(m(4.0, 2.0), "F0FFFF212", "edge");
        assert_eq!(m(4.0, 4.0), "F0FFFF212", "vertex");
        assert_eq!(m(9.0, 9.0), "FF0FFF212", "exterior");
    }

    #[test]
    fn point_against_open_and_closed_line() {
        let open = line(&[(0.0, 0.0), (4.0, 4.0)]);
        let set = LineSet::single(&open);
        let m = |x, y| points_matrix(&[c(x, y)], &Target::Lines(&set)).to_string();
        assert_eq!(
            m(2.0, 2.0),
            "0FFFFF102",
            "interior point; both ends uncovered"
        );
        assert_eq!(m(0.0, 0.0), "F0FFFF102", "at the start, a boundary point");
        assert_eq!(m(1.0, 3.0), "FF0FFF102", "off the line");

        let closed = line(&[(0.0, 0.0), (4.0, 0.0), (2.0, 3.0), (0.0, 0.0)]);
        let cset = LineSet::single(&closed);
        assert_eq!(
            points_matrix(&[c(0.0, 0.0)], &Target::Lines(&cset)).to_string(),
            "0FFFFF1F2",
            "a closed line has no boundary, so its start is an interior point"
        );
    }

    #[test]
    fn point_set_against_point_set() {
        let a = [c(1.0, 1.0), c(1.0, 1.0)];
        let b = [c(1.0, 1.0)];
        assert_eq!(
            points_matrix(&a, &Target::Points(&b)).to_string(),
            "0FFFFFFF2",
            "duplicates are one point"
        );
        let b2 = [c(1.0, 1.0), c(2.0, 2.0)];
        assert_eq!(
            points_matrix(&a, &Target::Points(&b2)).to_string(),
            "0FFFFF0F2"
        );
        assert_eq!(
            points_matrix(&b2, &Target::Points(&a)).to_string(),
            "0F0FFFFF2"
        );
    }

    #[test]
    fn empty_point_set_matches_the_ogc_empty_rows() {
        let sq = square();
        assert_eq!(
            points_matrix(&[], &Target::Area(&Areal::single(&sq))).to_string(),
            "FFFFFF212"
        );
        assert_eq!(
            points_matrix(&[], &Target::Points(&[])).to_string(),
            "FFFFFFFF2"
        );
        let open = line(&[(0.0, 0.0), (1.0, 1.0)]);
        let set = LineSet::single(&open);
        assert_eq!(
            points_matrix(&[], &Target::Lines(&set)).to_string(),
            "FFFFFF102"
        );
    }

    #[test]
    fn boundary_point_of_a_line_that_is_covered_by_the_set_leaves_eb_empty() {
        // Both endpoints of the line are in the point set, so no line
        // boundary point is outside it: EB = F.
        let open = line(&[(0.0, 0.0), (4.0, 4.0)]);
        let set = LineSet::single(&open);
        let pts = [c(0.0, 0.0), c(4.0, 4.0)];
        // IB = 0 (both points sit on the line boundary), EI = 1, EB = F.
        assert_eq!(
            points_matrix(&pts, &Target::Lines(&set)).to_string(),
            "F0FFFF1F2"
        );
    }
}
