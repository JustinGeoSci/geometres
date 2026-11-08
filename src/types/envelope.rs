// src/types/envelope.rs

//! Axis-aligned bounding rectangle ([`Envelope`]), the rectangle that grid
//! tessellation tiles and that spatial queries use as a window.
//!
//! # Deliberately minimal
//!
//! The type offers construction, accessors, two containment tests,
//! [`Envelope::from_points`] and [`Envelope::is_point`]. **It has no union, no
//! `from_geometry` and no envelope of a `Polygon`.** Those operations are easy
//! to add later from the outside, while an API shipped speculatively has to
//! be supported whether or not it turns out to be the right shape.
//!
//! The two envelope operations that **construct** coordinates, the center and
//! the expanded envelope, are not here. Every bound this type returns is an
//! input coordinate (the charter's rule for this crate).
//!
//! An envelope is **not a geometry**. It has no interior/boundary/exterior
//! semantics, takes no part in DE-9IM, and is not a degenerate `Polygon`. It is
//! a pair of corners, in the same role `Coord2d` plays for a point.

use crate::error::GeometryError;
use crate::types::Coord2d;

/// An axis-aligned rectangle given by its lower-left and upper-right corners.
///
/// # Invariants
///
/// - All four bounds are finite.
/// - `min_x <= max_x` and `min_y <= max_y`.
///
/// A **degenerate** envelope (zero width, zero height, or both) is valid. A
/// zero-width envelope is a legitimate query rectangle (a vertical line of
/// interest), and refusing it would make callers special-case a case the type
/// can represent perfectly well. Only non-finite or *inverted* bounds are
/// errors.
///
/// ```
/// use geometres::types::{Coord2d, Envelope};
///
/// let e = Envelope::new(Coord2d::new(0.0, 0.0), Coord2d::new(4.0, 2.0))?;
/// assert_eq!(e.width(), 4.0);
/// assert_eq!(e.height(), 2.0);
/// assert!(e.contains_coord(Coord2d::new(1.0, 1.0)));
///
/// // Degenerate is valid.
/// let line = Envelope::new(Coord2d::new(1.0, 0.0), Coord2d::new(1.0, 5.0))?;
/// assert_eq!(line.width(), 0.0);
///
/// // Inverted is not.
/// assert!(Envelope::new(Coord2d::new(4.0, 0.0), Coord2d::new(0.0, 2.0)).is_err());
/// # Ok::<(), geometres::GeometryError>(())
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Envelope {
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
}

impl Envelope {
    /// Construct an envelope from its lower-left and upper-right corners.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::NonFiniteCoordinate`] if any bound is NaN or
    ///   infinite. `coord_index` is 0 for `min` and 1 for `max`; `axis` is 0
    ///   for x and 1 for y.
    /// - [`GeometryError::InvertedEnvelope`] if `min_x > max_x` or
    ///   `min_y > max_y`.
    ///
    /// A zero-extent envelope, by contrast, is **accepted**:
    /// `min == max` is a valid point-sized rectangle, while `min > max` is a
    /// caller error that would otherwise make every containment test silently
    /// false.
    pub fn new(min: Coord2d, max: Coord2d) -> Result<Self, GeometryError> {
        for (index, c) in [(0usize, min), (1, max)] {
            if !c.x.is_finite() {
                return Err(GeometryError::NonFiniteCoordinate {
                    coord_index: index,
                    axis: 0,
                });
            }
            if !c.y.is_finite() {
                return Err(GeometryError::NonFiniteCoordinate {
                    coord_index: index,
                    axis: 1,
                });
            }
        }
        if min.x > max.x || min.y > max.y {
            return Err(GeometryError::InvertedEnvelope);
        }
        Ok(Self {
            min_x: min.x,
            min_y: min.y,
            max_x: max.x,
            max_y: max.y,
        })
    }

    /// The smallest envelope containing every point, or `Ok(None)` for an
    /// empty iterator.
    ///
    /// Each bound is the minimum or maximum of the input values on its
    /// axis, chosen by exact comparison, so every bound is an input
    /// coordinate and nothing is computed. On a tie (`-0.0` against `0.0`,
    /// which compare equal) the value met first is kept, on every platform,
    /// so the result's bit patterns are determined by the input order. A
    /// single point gives a point-sized envelope ([`Self::is_point`]).
    ///
    /// # Errors
    ///
    /// [`GeometryError::NonFiniteCoordinate`] on the first non-finite value
    /// met, with `coord_index` its position in the iterator (from 0) and
    /// `axis` 0 for x or 1 for y. Non-finiteness is checked point by point
    /// before the bounds are updated, so a NaN can never be absorbed by a
    /// comparison that is false for it.
    pub fn from_points<I>(points: I) -> Result<Option<Envelope>, GeometryError>
    where
        I: IntoIterator<Item = Coord2d>,
    {
        let mut bounds: Option<(f64, f64, f64, f64)> = None;
        for (coord_index, c) in points.into_iter().enumerate() {
            if !c.x.is_finite() {
                return Err(GeometryError::NonFiniteCoordinate {
                    coord_index,
                    axis: 0,
                });
            }
            if !c.y.is_finite() {
                return Err(GeometryError::NonFiniteCoordinate {
                    coord_index,
                    axis: 1,
                });
            }
            // Explicit comparisons rather than `f64::min`/`max`, so that on a
            // tie (`-0.0` against `0.0`) the value met first is kept on every
            // platform; `f64::min` may return either of two equal operands.
            bounds = Some(match bounds {
                None => (c.x, c.y, c.x, c.y),
                Some((min_x, min_y, max_x, max_y)) => (
                    if c.x < min_x { c.x } else { min_x },
                    if c.y < min_y { c.y } else { min_y },
                    if c.x > max_x { c.x } else { max_x },
                    if c.y > max_y { c.y } else { max_y },
                ),
            });
        }
        Ok(bounds.map(|(min_x, min_y, max_x, max_y)| Self {
            min_x,
            min_y,
            max_x,
            max_y,
        }))
    }

    /// `true` if the envelope has zero extent on both axes: `min == max`
    /// on x and on y, by exact comparison.
    ///
    /// A point-sized envelope is valid (see the type documentation); this
    /// reports it. A zero-width or zero-height envelope that is not also
    /// zero on the other axis is a line, not a point, and reports `false`.
    #[allow(clippy::float_cmp)] // stored bounds, compared exactly by design
    pub fn is_point(&self) -> bool {
        self.min_x == self.max_x && self.min_y == self.max_y
    }

    /// The lower-left corner.
    pub fn min(&self) -> Coord2d {
        Coord2d::new(self.min_x, self.min_y)
    }

    /// The upper-right corner.
    pub fn max(&self) -> Coord2d {
        Coord2d::new(self.max_x, self.max_y)
    }

    /// Extent along x.
    ///
    /// Exact when the subtraction is exact and correctly rounded otherwise;
    /// both bounds are finite, so the result is finite and non-negative.
    pub fn width(&self) -> f64 {
        self.max_x - self.min_x
    }

    /// Extent along y.
    pub fn height(&self) -> f64 {
        self.max_y - self.min_y
    }

    /// Is `c` inside or on the boundary of this envelope?
    ///
    /// **Closed on all four sides**: a coordinate exactly on an edge is
    /// contained. This matters for tessellation, where adjacent cells share an
    /// edge and a half-open rule would have to be chosen consistently. That
    /// decision belongs to the tessellator rather than to the rectangle.
    ///
    /// The comparisons are exact `f64` comparisons, not tolerance tests.
    pub fn contains_coord(&self, c: Coord2d) -> bool {
        c.x >= self.min_x && c.x <= self.max_x && c.y >= self.min_y && c.y <= self.max_y
    }

    /// Is `other` entirely within this envelope, boundary included?
    ///
    /// **Closed**, matching [`Self::contains_coord`]: an inner envelope that
    /// shares an edge with this one is contained. Used as an exact `O(1)`
    /// screen in front of the `O(n·m)` ring-containment test in
    /// `Polygon::new`: a ring whose bounding box is not inside another's
    /// cannot be inside it, and that rejection costs nothing.
    pub fn contains_envelope(&self, other: &Envelope) -> bool {
        other.min_x >= self.min_x
            && other.max_x <= self.max_x
            && other.min_y >= self.min_y
            && other.max_y <= self.max_y
    }

    /// Do the two envelopes share at least one point?
    ///
    /// **Closed**, matching [`Self::contains_coord`]: envelopes that touch
    /// only along an edge or at a corner do intersect.
    pub fn intersects(&self, other: &Envelope) -> bool {
        self.min_x <= other.max_x
            && other.min_x <= self.max_x
            && self.min_y <= other.max_y
            && other.min_y <= self.max_y
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

    fn env(a: f64, b: f64, x: f64, y: f64) -> Envelope {
        Envelope::new(c(a, b), c(x, y)).unwrap()
    }

    #[test]
    fn accessors_return_the_constructed_corners() {
        let e = env(-1.0, -2.0, 3.0, 4.0);
        assert_eq!(e.min(), c(-1.0, -2.0));
        assert_eq!(e.max(), c(3.0, 4.0));
        assert_eq!(e.width(), 4.0);
        assert_eq!(e.height(), 6.0);
    }

    #[test]
    fn a_zero_width_envelope_is_valid() {
        let e = env(1.0, 0.0, 1.0, 5.0);
        assert_eq!(e.width(), 0.0);
        assert_eq!(e.height(), 5.0);
        assert!(e.contains_coord(c(1.0, 2.5)));
    }

    #[test]
    fn a_point_sized_envelope_is_valid() {
        let e = env(2.0, 3.0, 2.0, 3.0);
        assert_eq!(e.width(), 0.0);
        assert_eq!(e.height(), 0.0);
        assert!(e.contains_coord(c(2.0, 3.0)));
        assert!(!e.contains_coord(c(2.0, 3.000_000_1)));
    }

    #[test]
    fn inverted_bounds_are_rejected_on_each_axis() {
        assert_eq!(
            Envelope::new(c(4.0, 0.0), c(0.0, 2.0)).unwrap_err(),
            GeometryError::InvertedEnvelope
        );
        assert_eq!(
            Envelope::new(c(0.0, 4.0), c(2.0, 0.0)).unwrap_err(),
            GeometryError::InvertedEnvelope
        );
    }

    #[test]
    fn non_finite_bounds_are_rejected_with_their_position() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                Envelope::new(c(bad, 0.0), c(1.0, 1.0)).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 0,
                    axis: 0
                }
            );
            assert_eq!(
                Envelope::new(c(0.0, bad), c(1.0, 1.0)).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 0,
                    axis: 1
                }
            );
            assert_eq!(
                Envelope::new(c(0.0, 0.0), c(bad, 1.0)).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 1,
                    axis: 0
                }
            );
            assert_eq!(
                Envelope::new(c(0.0, 0.0), c(1.0, bad)).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 1,
                    axis: 1
                }
            );
        }
    }

    /// Non-finiteness is checked before inversion.
    ///
    /// A NaN bound makes every comparison false, so `min > max` would not fire
    /// and an inverted-bounds error would be the wrong diagnosis.
    #[test]
    fn non_finite_is_reported_even_when_bounds_look_inverted() {
        assert_eq!(
            Envelope::new(c(f64::NAN, 5.0), c(0.0, 0.0)).unwrap_err(),
            GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 0
            }
        );
    }

    #[test]
    fn containment_is_closed_on_every_edge() {
        let e = env(0.0, 0.0, 2.0, 2.0);
        for p in [
            c(0.0, 0.0),
            c(2.0, 2.0),
            c(0.0, 2.0),
            c(2.0, 0.0),
            c(1.0, 0.0),
            c(0.0, 1.0),
        ] {
            assert!(e.contains_coord(p), "{p:?} is on the boundary");
        }
        assert!(!e.contains_coord(c(-0.000_000_1, 1.0)));
        assert!(!e.contains_coord(c(1.0, 2.000_000_1)));
    }

    #[test]
    fn intersection_is_closed_so_touching_counts() {
        let a = env(0.0, 0.0, 2.0, 2.0);
        assert!(a.intersects(&env(1.0, 1.0, 3.0, 3.0)), "overlapping");
        assert!(a.intersects(&env(2.0, 0.0, 4.0, 2.0)), "edge contact");
        assert!(a.intersects(&env(2.0, 2.0, 4.0, 4.0)), "corner contact");
        assert!(a.intersects(&env(0.5, 0.5, 1.5, 1.5)), "contained");
        assert!(!a.intersects(&env(2.5, 0.0, 4.0, 2.0)), "separated in x");
        assert!(!a.intersects(&env(0.0, 2.5, 2.0, 4.0)), "separated in y");
    }

    #[test]
    fn from_points_of_nothing_is_none() {
        assert_eq!(Envelope::from_points(std::iter::empty()), Ok(None));
    }

    #[test]
    fn from_points_of_one_point_is_that_point() {
        let e = Envelope::from_points([c(2.0, 3.0)]).unwrap().unwrap();
        assert_eq!(e, env(2.0, 3.0, 2.0, 3.0));
        assert!(e.is_point());
        assert_eq!(e.width(), 0.0);
    }

    #[test]
    fn from_points_takes_the_extremes_in_any_order() {
        let pts = [c(3.0, -1.0), c(-2.0, 7.0), c(0.5, 0.5), c(3.0, 7.0)];
        let e = Envelope::from_points(pts).unwrap().unwrap();
        assert_eq!(e, env(-2.0, -1.0, 3.0, 7.0));
        assert!(!e.is_point());
        // Every bound is one of the inputs, bit for bit.
        for v in [e.min().x, e.min().y, e.max().x, e.max().y] {
            assert!(
                pts.iter()
                    .any(|p| p.x.to_bits() == v.to_bits() || p.y.to_bits() == v.to_bits())
            );
        }
    }

    #[test]
    fn from_points_keeps_negative_zero_as_the_input_bit_pattern() {
        // `-0.0 == 0.0`, so the value met first is kept, in either order.
        let e = Envelope::from_points([c(-0.0, 1.0), c(0.0, 2.0)])
            .unwrap()
            .unwrap();
        assert_eq!(e.min().x.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(e.max().x.to_bits(), (-0.0_f64).to_bits());
        let e = Envelope::from_points([c(0.0, 1.0), c(-0.0, 2.0)])
            .unwrap()
            .unwrap();
        assert_eq!(e.min().x.to_bits(), 0.0_f64.to_bits());
        assert_eq!(e.max().x.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn from_points_reports_the_first_non_finite_value_by_index_and_axis() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                Envelope::from_points([c(0.0, 0.0), c(1.0, bad), c(bad, bad)]).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 1,
                    axis: 1
                }
            );
            assert_eq!(
                Envelope::from_points([c(bad, 0.0)]).unwrap_err(),
                GeometryError::NonFiniteCoordinate {
                    coord_index: 0,
                    axis: 0
                }
            );
        }
    }

    #[test]
    fn is_point_is_false_for_a_line_and_a_box() {
        assert!(!env(1.0, 0.0, 1.0, 5.0).is_point());
        assert!(!env(0.0, 2.0, 5.0, 2.0).is_point());
        assert!(!env(0.0, 0.0, 1.0, 1.0).is_point());
        assert!(env(-4.5, 9.0, -4.5, 9.0).is_point());
    }

    #[test]
    fn intersection_is_symmetric() {
        let pairs = [
            (env(0.0, 0.0, 2.0, 2.0), env(1.0, 1.0, 3.0, 3.0)),
            (env(0.0, 0.0, 2.0, 2.0), env(5.0, 5.0, 6.0, 6.0)),
            (env(0.0, 0.0, 2.0, 2.0), env(2.0, 2.0, 3.0, 3.0)),
        ];
        for (a, b) in pairs {
            assert_eq!(a.intersects(&b), b.intersects(&a));
        }
    }
}
