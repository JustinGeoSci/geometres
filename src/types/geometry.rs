// src/types/geometry.rs

//! Heterogeneous geometry value (`Geometry`).

use crate::types::{
    GeometryCollection, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};

/// A heterogeneous geometry value.
///
/// `Geometry` is `#[non_exhaustive]` so that new variants can be added in
/// minor versions without breaking downstream match expressions. Callers
/// must include a catch-all `_ =>` arm.
///
/// Each variant wraps an already-validated geometry type, so constructing a
/// `Geometry` is infallible and adds no invariants of its own.
///
/// # `GeometryCollection`
///
/// A collection of geometries of any types, nested at most
/// [`GeometryCollection::MAX_NESTING`] levels deep.
/// [`IntersectionMatrix::compute_geometry`] and the predicates accept it
/// when its polygonal leaves together form a valid `MultiPolygon` (line and
/// point leaves are unrestricted), and refuse it with
/// `Err(GeometryError::CollectionElementsOverlap)` otherwise.
///
/// [`IntersectionMatrix::compute_geometry`]: crate::topology::de9im::IntersectionMatrix::compute_geometry
#[non_exhaustive]
#[derive(Clone, Debug, PartialEq)]
pub enum Geometry {
    /// A single point.
    Point(Point),
    /// An open or closed polyline.
    LineString(LineString),
    /// A polygon with optional holes.
    Polygon(Polygon),
    /// A collection of points.
    MultiPoint(MultiPoint),
    /// A collection of polylines.
    MultiLineString(MultiLineString),
    /// A collection of polygons whose parts meet at most at finitely many
    /// points, as validated by `MultiPolygon::new`.
    MultiPolygon(MultiPolygon),
    /// A collection of geometries of any types, including other
    /// collections.
    GeometryCollection(GeometryCollection),
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
    use crate::types::Coord2d;

    #[test]
    fn variants_round_trip_by_pattern_match() {
        let pt = Point::new(Coord2d::new(1.0, 2.0)).unwrap();
        let g = Geometry::Point(pt);
        assert!(matches!(g, Geometry::Point(p) if p == pt));
        assert_eq!(g.clone(), g);
        assert_ne!(g, Geometry::MultiPoint(MultiPoint::empty()));
        assert_ne!(
            Geometry::MultiLineString(MultiLineString::empty()),
            Geometry::MultiPoint(MultiPoint::empty())
        );
    }
}
