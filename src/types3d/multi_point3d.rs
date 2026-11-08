// src/types3d/multi_point3d.rs

//! Collection of 3D points (`MultiPoint3d`).
//!
//! Follows [`MultiPoint`](crate::types::MultiPoint) exactly, with an
//! infallible constructor over already-validated elements.

use crate::types3d::Point3d;

/// A collection of zero or more 3D points.
///
/// Semantics match OGC `MultiPoint`: no minimum count, no ordering
/// requirement, duplicates permitted.
///
/// # Invariants
///
/// - The collection itself may be empty.
/// - Every element is a `Point3d`, so every coordinate is finite. This also
///   holds for a `MultiPoint3d` produced by decoding the native binary format,
///   because each element is decoded through [`Point3d::new`].
///
/// The finiteness guarantee matches 2D [`MultiPoint`](crate::types::MultiPoint).
/// It follows from the element type: [`Point3d::new`] rejects a non-finite
/// coordinate, so this type needs no check of its own.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiPoint3d {
    points: Vec<Point3d>,
}

impl MultiPoint3d {
    /// Construct a `MultiPoint3d` from a `Vec<Point3d>`.
    ///
    /// Infallible. An empty vector is accepted and represents an empty
    /// geometry. Element order is preserved.
    pub fn new(points: Vec<Point3d>) -> Self {
        Self { points }
    }

    /// Construct an empty `MultiPoint3d` (zero points).
    pub fn empty() -> Self {
        Self { points: vec![] }
    }

    /// Return a slice over the constituent points, in construction order.
    pub fn points(&self) -> &[Point3d] {
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
    use crate::types3d::Coord3d;

    fn p(x: f64, y: f64, z: f64) -> Point3d {
        Point3d::new(Coord3d::new(x, y, z)).unwrap()
    }

    #[test]
    fn round_trips_through_its_accessors() {
        let pts = vec![p(1.0, 2.0, 3.0), p(4.0, 5.0, 6.0)];
        let m = MultiPoint3d::new(pts.clone());
        assert_eq!(m.points(), pts.as_slice());
        assert_eq!(m.num_points(), 2);
        assert!(!m.is_empty());
    }

    #[test]
    fn empty_is_empty_and_matches_new_with_an_empty_vec() {
        assert!(MultiPoint3d::empty().is_empty());
        assert_eq!(MultiPoint3d::empty().num_points(), 0);
        assert_eq!(MultiPoint3d::empty(), MultiPoint3d::new(vec![]));
    }

    #[test]
    fn element_order_is_preserved() {
        let m = MultiPoint3d::new(vec![p(3.0, 0.0, 0.0), p(1.0, 0.0, 0.0)]);
        assert_eq!(m.points()[0], p(3.0, 0.0, 0.0));
        assert_eq!(m.points()[1], p(1.0, 0.0, 0.0));
    }

    #[test]
    fn duplicates_are_permitted() {
        let m = MultiPoint3d::new(vec![p(1.0, 1.0, 1.0), p(1.0, 1.0, 1.0)]);
        assert_eq!(m.num_points(), 2);
    }
}
