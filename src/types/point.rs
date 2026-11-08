// src/types/point.rs

//! Single-point geometry.

use crate::GeometryError;
use crate::types::Coord2d;

/// A single located point geometry.
///
/// A `Point` wraps exactly one `Coord2d` whose components are both finite.
///
/// # Invariants
///
/// - Both coordinates are finite (not NaN, not ±infinity). This is enforced
///   by [`Point::new`], and the inner field is private, so no `Point` with a
///   non-finite coordinate can exist. That includes values produced by
///   decoding the native binary format, which rebuilds through
///   [`Point::new`] and so goes through the same check.
///
/// # Why construction can fail
///
/// The exact predicates in `arithmetica` reject non-finite input. If a
/// non-finite coordinate reached a topology routine that maps that rejection
/// to a zero sign, the result would be a plausible but wrong answer
/// (for example "the point lies on the boundary"). Making the invalid state
/// unrepresentable removes that failure mode from every consumer at once,
/// instead of relying on each one to remember a check. This matches
/// `LineString::new`, `LinearRing::new`, and `Polygon::new`, which also
/// reject non-finite coordinates on construction.
///
/// The field is private, so a `Point` cannot be built around the check:
///
/// ```compile_fail
/// use geometres::{Coord2d, Point};
/// let _ = Point(Coord2d::new(f64::NAN, 0.0));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point(Coord2d);

impl Point {
    /// A point from a coordinate already known finite: one read from a
    /// validated geometry. Crate-private, so the type's invariant still
    /// holds for every value that exists.
    pub(crate) fn from_finite(coord: Coord2d) -> Self {
        Self(coord)
    }

    /// Construct a `Point` from a `Coord2d`.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError::NonFiniteCoordinate`] (with `coord_index` 0
    /// and `axis` 0 for x or 1 for y) if either component is NaN or
    /// infinite. When both are non-finite, x is reported.
    #[inline]
    pub fn new(coord: Coord2d) -> Result<Self, GeometryError> {
        if !coord.x.is_finite() {
            return Err(GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 0,
            });
        }
        if !coord.y.is_finite() {
            return Err(GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 1,
            });
        }
        Ok(Self(coord))
    }

    /// Return the x-coordinate of this point.
    #[inline]
    pub fn x(&self) -> f64 {
        self.0.x
    }

    /// Return the y-coordinate of this point.
    #[inline]
    pub fn y(&self) -> f64 {
        self.0.y
    }

    /// Return the underlying `Coord2d`.
    #[inline]
    pub fn coord(&self) -> Coord2d {
        self.0
    }
}

impl TryFrom<Coord2d> for Point {
    type Error = GeometryError;

    /// Equivalent to [`Point::new`].
    fn try_from(coord: Coord2d) -> Result<Self, Self::Error> {
        Self::new(coord)
    }
}

impl From<Point> for Coord2d {
    fn from(p: Point) -> Self {
        p.0
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

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    #[test]
    fn finite_coordinates_are_accepted_and_read_back_exactly() {
        let p = Point::new(c(1.5, -2.25)).unwrap();
        assert_eq!(p.x(), 1.5);
        assert_eq!(p.y(), -2.25);
        assert_eq!(p.coord(), c(1.5, -2.25));
        assert_eq!(Coord2d::from(p), c(1.5, -2.25));
    }

    #[test]
    fn extreme_but_finite_coordinates_are_accepted() {
        assert!(Point::new(c(f64::MAX, f64::MIN)).is_ok());
        assert!(Point::new(c(f64::MIN_POSITIVE / 2.0, -0.0)).is_ok());
    }

    #[test]
    fn every_non_finite_component_is_rejected_with_its_axis() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                Point::new(c(bad, 0.0)).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 0,
                    axis: 0
                }
            );
            assert_eq!(
                Point::new(c(0.0, bad)).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 0,
                    axis: 1
                }
            );
        }
    }

    #[test]
    fn x_is_reported_when_both_components_are_non_finite() {
        assert_eq!(
            Point::new(c(f64::NAN, f64::INFINITY)).unwrap_err(),
            GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 0
            }
        );
    }

    #[test]
    fn try_from_matches_new() {
        assert_eq!(Point::try_from(c(1.0, 2.0)), Point::new(c(1.0, 2.0)));
        assert!(Point::try_from(c(f64::NAN, 2.0)).is_err());
    }

    /// Decoding must go through `Point::new`.
    ///
    /// A coordinate arriving from outside cannot bypass validation. The
    /// native decoder is the way in from outside, so the property is
    /// asserted against it.
    mod decode_validation {
        use super::*;
        use crate::codec::{CodecError, decode, encode};

        /// Build a `Point` document by hand, so a coordinate that `Point::new`
        /// would reject can be put on the wire. Encoding a real `Point` cannot
        /// produce one, which is the whole point.
        fn point_bytes(x: f64, y: f64) -> Vec<u8> {
            let mut v = vec![crate::codec::FORMAT_VERSION, crate::codec::tag::POINT];
            v.extend_from_slice(&x.to_le_bytes());
            v.extend_from_slice(&y.to_le_bytes());
            v
        }

        #[test]
        fn finite_coordinates_decode() {
            let p = Point::new(c(1.0, 2.0)).unwrap();
            assert_eq!(decode::<Point>(&encode(&p)).unwrap(), p);
        }

        #[test]
        fn non_finite_coordinates_fail_to_decode() {
            for (x, y, axis) in [
                (f64::NAN, 0.0, 0),
                (0.0, f64::INFINITY, 1),
                (f64::NEG_INFINITY, f64::NAN, 0),
            ] {
                let err = decode::<Point>(&point_bytes(x, y)).unwrap_err();
                assert_eq!(
                    err,
                    CodecError::Geometry(GeometryError::NonFiniteCoordinate {
                        coord_index: 0,
                        axis
                    }),
                    "({x}, {y}) must be refused, not decoded"
                );
            }
        }

        #[test]
        fn a_refused_stream_yields_no_point_at_all() {
            // Not merely an error alongside a value: there is no Point.
            assert!(decode::<Point>(&point_bytes(f64::NAN, 0.0)).is_err());
        }
    }
}
