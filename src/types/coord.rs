// src/types/coord.rs

//! `Coord2d` scalar primitive and `CoordSeq` coordinate buffer.

use std::ops::Index;

/// A 2-dimensional coordinate pair.
///
/// `Coord2d` is a scalar primitive, **not a geometry**. It carries a pair of
/// `f64` values and makes no validity guarantees: callers may construct a
/// `Coord2d` with NaN, ±infinity, or any other non-finite value. Validation
/// of coordinate values is the responsibility of the geometry constructor
/// that accepts a sequence of `Coord2d` values (see `CoordSeq` and the
/// constructors on `LineString`, `LinearRing`, and `Polygon`).
///
/// Use `Coord2d` to pass coordinates *into* geometry constructors, and to
/// receive coordinate data *out of* geometry accessors. Do not interpret a
/// standalone `Coord2d` as having any geometric meaning.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coord2d {
    /// The x-axis (horizontal) component.
    pub x: f64,
    /// The y-axis (vertical) component.
    pub y: f64,
}

impl Coord2d {
    /// Construct a new `Coord2d` from x and y components.
    ///
    /// No validation is performed. Non-finite values are accepted here;
    /// they will be rejected when the `Coord2d` is used to construct a
    /// validated geometry type.
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

/// A read-only, ordered buffer of `Coord2d` values.
///
/// `CoordSeq` is a pure coordinate buffer with no closure semantics. It
/// stores `n` distinct vertices in the order provided and exposes only
/// read access after construction. The closure invariant (closing vertex
/// equals first vertex) is a property of `LinearRing`, not of `CoordSeq`
/// itself; `CoordSeq` carries no such constraint.
///
/// Mutation methods (`push`, `pop`, `clear`, and similar) are intentionally
/// absent from the public API. All coordinates are fixed at construction time.
#[derive(Clone, Debug, PartialEq)]
pub struct CoordSeq(Vec<Coord2d>);

impl CoordSeq {
    /// Construct a `CoordSeq` from an owned `Vec<Coord2d>`.
    ///
    /// The supplied coordinates are stored as-is. No validation is performed.
    pub fn new(coords: Vec<Coord2d>) -> Self {
        Self(coords)
    }

    /// Return the number of coordinates in the sequence.
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Return `true` if the sequence contains no coordinates.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Return a slice over the coordinate values.
    #[inline]
    pub fn as_slice(&self) -> &[Coord2d] {
        &self.0
    }

    /// Return an iterator over references to the coordinate values.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &Coord2d> {
        self.0.iter()
    }
}

impl Index<usize> for CoordSeq {
    type Output = Coord2d;

    #[inline]
    fn index(&self, i: usize) -> &Self::Output {
        &self.0[i]
    }
}

/// A bounding box `(min_x, min_y, max_x, max_y)`, each bound one of the
/// input values (selected by comparison, never computed).
pub(crate) type Bbox = (f64, f64, f64, f64);

/// The bounding box of a nonempty slice of finite coordinates, folded from
/// the first with `f64::min` and `f64::max`: the fold the boxes were always
/// computed with, so a kept box holds the same values a fresh pass would.
///
/// For an empty slice (unreachable from a constructed geometry) the result
/// is the inverted box `(+∞, +∞, −∞, −∞)`, which meets nothing.
pub(crate) fn bbox_of(cs: &[Coord2d]) -> Bbox {
    let Some(first) = cs.first() else {
        return (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
    };
    cs[1..]
        .iter()
        .fold((first.x, first.y, first.x, first.y), |e, c| {
            (e.0.min(c.x), e.1.min(c.y), e.2.max(c.x), e.3.max(c.y))
        })
}

/// The smallest box holding every box, by the same comparisons; `None`
/// for none.
///
/// Folding cached boxes gives the bounds a fold over all the coordinates
/// gives, as values: `min` and `max` are associative and commutative on
/// finite values, up to which of two equal operands is returned. The one
/// pair of equal values with different bits is `0.0` and `-0.0`, so a bound
/// may differ from the coordinate fold's only in the sign of a zero. Every
/// caller only compares the bounds, and the two zeros compare equal.
pub(crate) fn bbox_fold(boxes: impl IntoIterator<Item = Bbox>) -> Option<Bbox> {
    boxes
        .into_iter()
        .reduce(|a, b| (a.0.min(b.0), a.1.min(b.1), a.2.max(b.2), a.3.max(b.3)))
}

/// Validates a `CoordSeq` for non-finite values and consecutive duplicate vertices.
///
/// This is the shared validation pass used by all geometry constructors.
/// It checks non-finite coordinates first (full pass), then consecutive
/// duplicates. Callers that need a minimum-length check perform it separately
/// after this function returns `Ok`.
///
/// # Errors
/// - `NonFiniteCoordinate` if any `x` or `y` value is NaN or infinite.
/// - `ConsecutiveDuplicateVertices` if any two adjacent coordinates are identical.
pub(super) fn check_coords(coords: &CoordSeq) -> Result<(), crate::GeometryError> {
    for (i, c) in coords.iter().enumerate() {
        if !c.x.is_finite() {
            return Err(crate::GeometryError::NonFiniteCoordinate {
                coord_index: i,
                axis: 0,
            });
        }
        if !c.y.is_finite() {
            return Err(crate::GeometryError::NonFiniteCoordinate {
                coord_index: i,
                axis: 1,
            });
        }
    }
    for i in 0..coords.len().saturating_sub(1) {
        if coords[i] == coords[i + 1] {
            return Err(crate::GeometryError::ConsecutiveDuplicateVertices { index: i });
        }
    }
    Ok(())
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
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
    fn coord2d_new_round_trips_through_accessors() {
        let c = Coord2d::new(1.0, 2.0);
        assert_eq!(c.x, 1.0);
        assert_eq!(c.y, 2.0);
    }

    #[test]
    fn coord_seq_as_slice_returns_expected_coordinates() {
        let a = Coord2d::new(0.0, 0.0);
        let b = Coord2d::new(1.0, 0.0);
        let c = Coord2d::new(1.0, 1.0);
        let seq = CoordSeq::new(vec![a, b, c]);
        let sl = seq.as_slice();
        assert_eq!(sl.len(), 3);
        assert_eq!(sl[0], a);
        assert_eq!(sl[1], b);
        assert_eq!(sl[2], c);
    }

    #[test]
    fn coord_seq_index_matches_as_slice() {
        let a = Coord2d::new(3.0, 4.0);
        let seq = CoordSeq::new(vec![a]);
        assert_eq!(seq[0], a);
    }

    #[test]
    fn coord_seq_len_and_is_empty() {
        let empty = CoordSeq::new(vec![]);
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);

        let one = CoordSeq::new(vec![Coord2d::new(0.0, 0.0)]);
        assert!(!one.is_empty());
        assert_eq!(one.len(), 1);
    }

    #[test]
    fn coord_seq_iter_yields_all_coordinates() {
        let coords = vec![
            Coord2d::new(1.0, 2.0),
            Coord2d::new(3.0, 4.0),
            Coord2d::new(5.0, 6.0),
        ];
        let seq = CoordSeq::new(coords.clone());
        let collected: Vec<&Coord2d> = seq.iter().collect();
        assert_eq!(collected.len(), 3);
        for (i, c) in collected.into_iter().enumerate() {
            assert_eq!(*c, coords[i]);
        }
    }
}
