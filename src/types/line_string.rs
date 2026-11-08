// src/types/line_string.rs

//! Open coordinate-sequence geometry (`LineString`).

use std::sync::OnceLock;

use crate::GeometryError;
use crate::types::coord::{Bbox, Coord2d, CoordSeq, bbox_of, check_coords};

/// An open sequence of two or more `Coord2d` vertices.
///
/// A `LineString` represents a polyline: an ordered list of positions
/// connected by straight line segments. The first and last vertices are
/// distinct (an open curve). For a closed ring boundary, use `LinearRing`.
///
/// # Construction
///
/// `LineString::new` validates the coordinate sequence before accepting it.
/// The validation order is: non-finite coordinates, minimum length (≥ 2),
/// consecutive duplicate vertices.
///
/// # Invariants
///
/// - At least 2 vertices.
/// - No coordinate value is NaN or infinite.
/// - No two consecutive vertices are identical.
///
/// # The cached box
///
/// A line string keeps its bounding box once it is first asked for, so the
/// operations that screen a line by its box pay the pass over its vertices
/// once per line, not once per call. It is computed on first use, not at
/// construction, so building lines costs nothing more. The box is a
/// function of the coordinates; it takes no part in equality or in the
/// `Debug` output, which are those of the coordinates alone.
#[derive(Clone)]
pub struct LineString {
    coords: CoordSeq,
    bbox: OnceLock<Bbox>,
}

impl PartialEq for LineString {
    fn eq(&self, other: &Self) -> bool {
        self.coords == other.coords
    }
}

// The box is derived from the coordinates and left out on purpose, so the
// output is the coordinates' alone, as it was before the box was cached.
#[allow(clippy::missing_fields_in_debug)]
impl core::fmt::Debug for LineString {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("LineString")
            .field("coords", &self.coords)
            .finish()
    }
}

impl LineString {
    /// A line string from coordinates already known valid (finite, at least
    /// two, no two consecutive equal): a rearrangement of a validated line's
    /// coordinates, or a stretch of one, that the caller has argued keeps
    /// those properties. Crate-private.
    pub(crate) fn from_valid(coords: CoordSeq) -> Self {
        debug_assert!(coords.len() >= 2);
        Self {
            coords,
            bbox: OnceLock::new(),
        }
    }

    /// Construct a `LineString` from a coordinate sequence.
    ///
    /// Validation order: non-finite check (full sequence), minimum-length
    /// check (≥ 2), consecutive duplicate check.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::NonFiniteCoordinate`] if any coordinate value is
    ///   NaN or infinite. `coord_index` names the offending position; `axis`
    ///   is `0` for x, `1` for y.
    /// - [`GeometryError::InsufficientVertices`] if `coords.len() < 2`.
    /// - [`GeometryError::ConsecutiveDuplicateVertices`] if any two
    ///   consecutive coordinates are identical.
    pub fn new(coords: CoordSeq) -> Result<Self, GeometryError> {
        // 1. Non-finite check first: a NaN length check could mask a bad value.
        check_coords(&coords)?;
        // 2. Minimum length.
        if coords.len() < 2 {
            return Err(GeometryError::InsufficientVertices {
                required: 2,
                found: coords.len(),
            });
        }
        Ok(Self {
            coords,
            bbox: OnceLock::new(),
        })
    }

    /// The box `(min_x, min_y, max_x, max_y)` over every vertex, computed on
    /// first use and kept.
    pub(crate) fn bbox(&self) -> Bbox {
        *self.bbox.get_or_init(|| bbox_of(self.coords.as_slice()))
    }

    /// Return a reference to the underlying coordinate sequence.
    pub fn coords(&self) -> &CoordSeq {
        &self.coords
    }

    /// Return the number of vertices in this line string.
    pub fn num_points(&self) -> usize {
        self.coords.len()
    }

    /// Return the first (start) vertex.
    pub fn start(&self) -> Coord2d {
        self.coords[0]
    }

    /// Return the last (end) vertex.
    pub fn end(&self) -> Coord2d {
        self.coords[self.coords.len() - 1]
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

    fn seq(pts: &[(f64, f64)]) -> CoordSeq {
        CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())
    }

    #[test]
    fn rejects_fewer_than_two_coords() {
        let err = LineString::new(seq(&[(0.0, 0.0)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::InsufficientVertices {
                required: 2,
                found: 1
            }
        );
        let err2 = LineString::new(seq(&[])).unwrap_err();
        assert_eq!(
            err2,
            GeometryError::InsufficientVertices {
                required: 2,
                found: 0
            }
        );
    }

    #[test]
    fn rejects_consecutive_duplicate_vertices() {
        let err =
            LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0), (1.0, 1.0), (2.0, 2.0)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::ConsecutiveDuplicateVertices { index: 1 }
        );
    }

    #[test]
    fn rejects_non_finite_coordinate() {
        let err = LineString::new(seq(&[(f64::NAN, 0.0), (1.0, 1.0)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 0
            }
        );
        let err2 = LineString::new(seq(&[(0.0, f64::INFINITY), (1.0, 1.0)])).unwrap_err();
        assert_eq!(
            err2,
            GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 1
            }
        );
    }

    #[test]
    fn non_finite_checked_before_length() {
        // A single NaN coordinate must report NonFiniteCoordinate rather
        // than InsufficientVertices.
        let err = LineString::new(seq(&[(f64::NAN, 0.0)])).unwrap_err();
        assert!(matches!(err, GeometryError::NonFiniteCoordinate { .. }));
    }

    #[test]
    fn valid_line_string_accessors() {
        let ls = LineString::new(seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])).unwrap();
        assert_eq!(ls.num_points(), 3);
        assert_eq!(ls.start(), c(0.0, 0.0));
        assert_eq!(ls.end(), c(1.0, 1.0));
    }
}
