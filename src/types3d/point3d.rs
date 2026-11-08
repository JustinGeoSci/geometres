// src/types3d/point3d.rs

//! Single-point 3D geometry.

use crate::error::GeometryError;
use crate::types3d::Coord3d;

/// A single located 3D point geometry.
///
/// A `Point3d` wraps exactly one [`Coord3d`] whose three components are all
/// finite.
///
/// # Invariants
///
/// - All three coordinates are finite (not NaN, not ±infinity). This is
///   enforced by [`Point3d::new`], and the inner field is private, so no
///   `Point3d` with a non-finite coordinate can exist. That includes values
///   produced by decoding the native binary format, which rebuilds through
///   [`Point3d::new`] and so goes through the same check.
///
/// # Why construction can fail
///
/// The exact predicates in `arithmetica` reject non-finite input, so a
/// `Point3d` holding a NaN is a value nothing downstream can use. Making the
/// state unrepresentable is strictly better than making every consumer
/// remember to check for it. The 2D [`Point`](crate::types::Point) follows
/// the same reasoning.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3d(Coord3d);

impl Point3d {
    /// Construct a `Point3d` from a `Coord3d`.
    ///
    /// # Errors
    ///
    /// [`GeometryError::NonFiniteCoordinate`] if any component is NaN or
    /// infinite. `coord_index` is always 0, since a `Point3d` holds one
    /// coordinate, and `axis` is `0` for x, `1` for y, `2` for z. The first
    /// offending axis in that order is reported.
    #[inline]
    pub fn new(coord: Coord3d) -> Result<Self, GeometryError> {
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
        if !coord.z.is_finite() {
            return Err(GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 2,
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

    /// Return the z-coordinate of this point.
    #[inline]
    pub fn z(&self) -> f64 {
        self.0.z
    }

    /// Return the underlying `Coord3d`.
    #[inline]
    pub fn coord(&self) -> Coord3d {
        self.0
    }
}

impl TryFrom<Coord3d> for Point3d {
    type Error = GeometryError;

    /// Equivalent to [`Point3d::new`].
    #[inline]
    fn try_from(coord: Coord3d) -> Result<Self, Self::Error> {
        Self::new(coord)
    }
}

impl From<Point3d> for Coord3d {
    /// A `Point3d` always holds a valid coordinate, so this direction is
    /// infallible.
    #[inline]
    fn from(p: Point3d) -> Self {
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

    #[test]
    fn accessors_round_trip() {
        let p = Point3d::new(Coord3d::new(1.0, 2.0, 3.0)).unwrap();
        assert_eq!((p.x(), p.y(), p.z()), (1.0, 2.0, 3.0));
        assert_eq!(p.coord(), Coord3d::new(1.0, 2.0, 3.0));
        assert_eq!(Coord3d::from(p), Coord3d::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn extreme_but_finite_coordinates_are_accepted() {
        assert!(Point3d::new(Coord3d::new(f64::MAX, f64::MIN, 0.0)).is_ok());
        assert!(
            Point3d::new(Coord3d::new(
                f64::MIN_POSITIVE / 2.0,
                -0.0,
                f64::from_bits(1)
            ))
            .is_ok()
        );
    }

    /// Every non-finite component is rejected, and the reported axis is the
    /// offending one.
    ///
    /// Axis 2 is the case the 2D type cannot have. `error.rs` has a
    /// `Display` arm for it, so a non-finite z reports "non-finite z
    /// coordinate" rather than being mislabeled.
    #[test]
    fn every_non_finite_component_is_rejected_with_its_axis() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for (axis, coord) in [
                (0u8, Coord3d::new(bad, 0.0, 0.0)),
                (1, Coord3d::new(0.0, bad, 0.0)),
                (2, Coord3d::new(0.0, 0.0, bad)),
            ] {
                assert_eq!(
                    Point3d::new(coord).unwrap_err(),
                    GeometryError::NonFiniteCoordinate {
                        coord_index: 0,
                        axis
                    },
                    "{bad:?} on axis {axis} must be rejected"
                );
            }
        }
    }

    /// x is reported when several components are non-finite.
    #[test]
    fn the_first_offending_axis_is_reported() {
        assert_eq!(
            Point3d::new(Coord3d::new(f64::NAN, f64::INFINITY, f64::NAN)).unwrap_err(),
            GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 0
            }
        );
        assert_eq!(
            Point3d::new(Coord3d::new(0.0, f64::NAN, f64::NAN)).unwrap_err(),
            GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 1
            }
        );
    }

    #[test]
    fn try_from_matches_new() {
        let c = Coord3d::new(1.0, 2.0, 3.0);
        assert_eq!(Point3d::try_from(c), Point3d::new(c));
        assert!(Point3d::try_from(Coord3d::new(f64::NAN, 2.0, 3.0)).is_err());
    }
}
