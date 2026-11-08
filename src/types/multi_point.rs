// src/types/multi_point.rs

//! Multi-point geometry (`MultiPoint`).

use crate::types::point::Point;

/// A collection of zero or more [`Point`] values.
///
/// `MultiPoint::new` is infallible: the collection adds no constraint of its
/// own, and every element is a [`Point`], which is already guaranteed finite.
/// An empty `MultiPoint` is valid and represents an empty geometry.
///
/// # Why construction cannot fail
///
/// Unlike `LineString::new` and `Polygon::new`, which validate a raw
/// coordinate sequence (finiteness, minimum length, duplicates, winding),
/// `MultiPoint` is assembled from elements that are already validated
/// geometries. Finiteness is checked once, in [`Point::new`], so it cannot be
/// skipped here. The collection adds no constraint of its own: there is no
/// minimum count, no ordering requirement, and duplicate points are
/// permitted.
///
/// # Invariants
///
/// - Every element is a `Point`, so every coordinate is finite. This also
///   holds for a `MultiPoint` produced by decoding the native binary format,
///   because each element is decoded through [`Point::new`].
/// - The collection itself may be empty.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiPoint {
    points: Vec<Point>,
}

impl MultiPoint {
    /// Construct a `MultiPoint` from a `Vec<Point>`.
    ///
    /// Infallible. An empty vector is accepted and represents an empty
    /// geometry. Element order is preserved.
    pub fn new(points: Vec<Point>) -> Self {
        Self { points }
    }

    /// Construct an empty `MultiPoint` (zero points).
    ///
    /// Equivalent to `MultiPoint::new(vec![])`.
    pub fn empty() -> Self {
        Self { points: vec![] }
    }

    /// Return a slice over the constituent points, in construction order.
    pub fn points(&self) -> &[Point] {
        &self.points
    }

    /// Return the number of points in the collection.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Return `true` if the collection contains no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
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
    use crate::types::coord::Coord2d;

    fn p(x: f64, y: f64) -> Point {
        Point::new(Coord2d::new(x, y)).unwrap()
    }

    #[test]
    fn single_point_round_trips_through_points() {
        let pt = p(1.0, 2.0);
        let mp = MultiPoint::new(vec![pt]);
        assert_eq!(mp.points(), &[pt]);
        assert_eq!(mp.num_points(), 1);
        assert!(!mp.is_empty());
    }

    #[test]
    fn multiple_points_preserve_order() {
        let pts = vec![p(3.0, 3.0), p(1.0, 1.0), p(2.0, 2.0)];
        let mp = MultiPoint::new(pts.clone());
        assert_eq!(mp.points(), pts.as_slice());
        assert_eq!(mp.num_points(), 3);
    }

    #[test]
    fn duplicate_points_are_permitted() {
        let mp = MultiPoint::new(vec![p(1.0, 1.0), p(1.0, 1.0)]);
        assert_eq!(mp.num_points(), 2);
    }

    #[test]
    fn empty_is_empty() {
        let mp = MultiPoint::empty();
        assert!(mp.is_empty());
        assert_eq!(mp.num_points(), 0);
        assert_eq!(mp.points(), &[]);
    }

    #[test]
    fn empty_equals_new_with_empty_vec() {
        assert_eq!(MultiPoint::empty(), MultiPoint::new(vec![]));
        assert!(MultiPoint::new(vec![]).is_empty());
    }

    /// A non-finite element must be refused on decode.
    ///
    /// `MultiPoint::new` is infallible because its elements are
    /// already `Point`s and therefore already validated; the property worth
    /// testing is that a hostile stream cannot smuggle an invalid element in
    /// through the collection, which it cannot, because each element decodes
    /// through `Point::new`.
    #[test]
    fn decoding_a_non_finite_element_fails() {
        use crate::codec::{CodecError, FORMAT_VERSION, decode, encode, tag};
        use crate::error::GeometryError;

        let good = MultiPoint::new(vec![p(1.0, 2.0), p(3.0, 4.0)]);
        assert_eq!(decode::<MultiPoint>(&encode(&good)).unwrap(), good);

        // Hand-built: two points, the second with a NaN x.
        let mut bytes = vec![FORMAT_VERSION, tag::MULTI_POINT];
        bytes.extend_from_slice(&2u32.to_le_bytes());
        for (x, y) in [(1.0_f64, 2.0_f64), (f64::NAN, 4.0)] {
            bytes.extend_from_slice(&x.to_le_bytes());
            bytes.extend_from_slice(&y.to_le_bytes());
        }
        assert_eq!(
            decode::<MultiPoint>(&bytes).unwrap_err(),
            CodecError::Geometry(GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 0
            })
        );
    }
}
