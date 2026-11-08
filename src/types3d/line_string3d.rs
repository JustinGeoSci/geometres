// src/types3d/line_string3d.rs

//! Open coordinate-sequence 3D geometry (`LineString3d`).

use crate::GeometryError;
use crate::types3d::coord3d::check_coords3d;
use crate::types3d::{Coord3d, CoordSeq3d};

/// An open sequence of two or more `Coord3d` vertices.
///
/// A separate type from `LineString`, not a `LineString` with optional Z.
///
/// # Invariants
///
/// - At least 2 vertices.
/// - No component of any coordinate is NaN or infinite.
/// - No two consecutive vertices are identical in all of x, y, and z.
#[derive(Clone, Debug, PartialEq)]
pub struct LineString3d {
    coords: CoordSeq3d,
}

impl LineString3d {
    /// Construct a `LineString3d` from a coordinate sequence.
    ///
    /// Validation order: non-finite check (full sequence), minimum-length
    /// check (≥ 2), consecutive duplicate check.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::NonFiniteCoordinate`] if any component is NaN or
    ///   infinite; `axis` is `0` for x, `1` for y, `2` for z.
    /// - [`GeometryError::InsufficientVertices`] if `coords.len() < 2`.
    /// - [`GeometryError::ConsecutiveDuplicateVertices`] if two consecutive
    ///   coordinates are identical.
    pub fn new(coords: CoordSeq3d) -> Result<Self, GeometryError> {
        check_coords3d(&coords)?;
        if coords.len() < 2 {
            return Err(GeometryError::InsufficientVertices {
                required: 2,
                found: coords.len(),
            });
        }
        Ok(Self { coords })
    }

    /// Return a reference to the underlying coordinate sequence.
    pub fn coords(&self) -> &CoordSeq3d {
        &self.coords
    }

    /// Return the number of vertices in this line string.
    pub fn num_points(&self) -> usize {
        self.coords.len()
    }

    /// Return the first (start) vertex.
    pub fn start(&self) -> Coord3d {
        self.coords[0]
    }

    /// Return the last (end) vertex.
    pub fn end(&self) -> Coord3d {
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

    fn seq(pts: &[(f64, f64, f64)]) -> CoordSeq3d {
        CoordSeq3d::new(pts.iter().map(|&(x, y, z)| Coord3d::new(x, y, z)).collect())
    }

    #[test]
    fn rejects_fewer_than_two_coords() {
        assert_eq!(
            LineString3d::new(seq(&[(0.0, 0.0, 0.0)])).unwrap_err(),
            GeometryError::InsufficientVertices {
                required: 2,
                found: 1
            }
        );
        assert_eq!(
            LineString3d::new(seq(&[])).unwrap_err(),
            GeometryError::InsufficientVertices {
                required: 2,
                found: 0
            }
        );
    }

    #[test]
    fn rejects_non_finite_z() {
        let err = LineString3d::new(seq(&[(0.0, 0.0, 0.0), (1.0, 1.0, f64::NAN)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::NonFiniteCoordinate {
                coord_index: 1,
                axis: 2
            }
        );
    }

    #[test]
    fn rejects_consecutive_duplicate() {
        let err = LineString3d::new(seq(&[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::ConsecutiveDuplicateVertices { index: 0 }
        );
    }

    #[test]
    fn accepts_vertical_segment_and_exposes_endpoints() {
        let ls = LineString3d::new(seq(&[(1.0, 1.0, 0.0), (1.0, 1.0, 5.0)])).unwrap();
        assert_eq!(ls.num_points(), 2);
        assert_eq!(ls.start(), Coord3d::new(1.0, 1.0, 0.0));
        assert_eq!(ls.end(), Coord3d::new(1.0, 1.0, 5.0));
        assert_eq!(ls.coords().len(), 2);
    }
}
