// src/types3d/coord3d.rs

//! `Coord3d` scalar primitive and `CoordSeq3d` coordinate buffer.

use std::ops::Index;

/// A 3-dimensional coordinate triple.
///
/// `Coord3d` is a scalar primitive, **not a geometry**, and is an independent
/// struct rather than an extension of `Coord2d`. Like `Coord2d` it makes no
/// validity guarantees: NaN and infinite components are accepted here and
/// rejected by the geometry constructors that consume a [`CoordSeq3d`].
///
/// All topology computation in `geometres` uses only the XY components.
/// The Z component is stored and round-tripped but is not used in any
/// predicate or topology operation in this version.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coord3d {
    /// The x-axis component.
    pub x: f64,
    /// The y-axis component.
    pub y: f64,
    /// The z-axis (vertical) component.
    pub z: f64,
}

impl Coord3d {
    /// Construct a new `Coord3d` from x, y, and z components.
    ///
    /// No validation is performed. Non-finite values are accepted here;
    /// they are rejected when the coordinate is used to construct a
    /// validated geometry type.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Projects to 2D by dropping Z. Used by `Polygon3d` for its certified
    /// winding and hole placement, which are decided on the XY projection.
    /// Kept `pub(crate)` so no 3D-to-2D conversion leaks into the public
    /// API.
    #[inline]
    pub(crate) fn to_2d(self) -> crate::types::Coord2d {
        crate::types::Coord2d::new(self.x, self.y)
    }
}

/// A read-only, ordered buffer of `Coord3d` values.
///
/// `CoordSeq3d` mirrors `CoordSeq`: a pure coordinate buffer with no closure
/// semantics and no mutation methods. All coordinates are fixed at
/// construction time.
#[derive(Clone, Debug, PartialEq)]
pub struct CoordSeq3d(Vec<Coord3d>);

impl CoordSeq3d {
    /// Construct a `CoordSeq3d` from an owned `Vec<Coord3d>`.
    ///
    /// The supplied coordinates are stored as-is. No validation is performed.
    pub fn new(coords: Vec<Coord3d>) -> Self {
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
    pub fn as_slice(&self) -> &[Coord3d] {
        &self.0
    }

    /// Return an iterator over references to the coordinate values.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &Coord3d> {
        self.0.iter()
    }
}

impl Index<usize> for CoordSeq3d {
    type Output = Coord3d;

    #[inline]
    fn index(&self, i: usize) -> &Self::Output {
        &self.0[i]
    }
}

/// Validates a `CoordSeq3d` for non-finite values and consecutive duplicates.
///
/// Shared validation pass for all 3D geometry constructors. Checks every
/// component of every coordinate for finiteness first (x, y, then z per
/// coordinate), then checks consecutive duplicates. Two 3D coordinates are
/// duplicates only if x, y, **and** z are all equal: points that differ only
/// in z are distinct 3D vertices even though they coincide in projection.
///
/// # Errors
/// - `NonFiniteCoordinate` (`axis` 0, 1, or 2) if any component is NaN or infinite.
/// - `ConsecutiveDuplicateVertices` if two adjacent coordinates are identical.
pub(super) fn check_coords3d(coords: &CoordSeq3d) -> Result<(), crate::GeometryError> {
    for (i, c) in coords.iter().enumerate() {
        for (axis, v) in [(0_u8, c.x), (1, c.y), (2, c.z)] {
            if !v.is_finite() {
                return Err(crate::GeometryError::NonFiniteCoordinate {
                    coord_index: i,
                    axis,
                });
            }
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
    fn coord3d_new_round_trips() {
        let c = Coord3d::new(1.0, 2.0, 3.0);
        assert_eq!((c.x, c.y, c.z), (1.0, 2.0, 3.0));
    }

    #[test]
    fn to_2d_drops_z() {
        let c = Coord3d::new(1.0, 2.0, 99.0);
        assert_eq!(c.to_2d(), crate::types::Coord2d::new(1.0, 2.0));
    }

    #[test]
    fn coord_seq3d_accessors() {
        let a = Coord3d::new(0.0, 0.0, 0.0);
        let b = Coord3d::new(1.0, 2.0, 3.0);
        let seq = CoordSeq3d::new(vec![a, b]);
        assert_eq!(seq.len(), 2);
        assert!(!seq.is_empty());
        assert_eq!(seq.as_slice(), &[a, b]);
        assert_eq!(seq[1], b);
        assert_eq!(seq.iter().count(), 2);
        assert!(CoordSeq3d::new(vec![]).is_empty());
    }

    #[test]
    fn check_coords3d_reports_each_axis() {
        for (axis, bad) in [
            (0_u8, Coord3d::new(f64::NAN, 0.0, 0.0)),
            (1, Coord3d::new(0.0, f64::INFINITY, 0.0)),
            (2, Coord3d::new(0.0, 0.0, f64::NEG_INFINITY)),
        ] {
            let seq = CoordSeq3d::new(vec![Coord3d::new(1.0, 1.0, 1.0), bad]);
            assert_eq!(
                check_coords3d(&seq),
                Err(crate::GeometryError::NonFiniteCoordinate {
                    coord_index: 1,
                    axis
                })
            );
        }
    }

    #[test]
    fn check_coords3d_duplicate_requires_all_three_components() {
        let a = Coord3d::new(1.0, 1.0, 1.0);
        let same_xy = Coord3d::new(1.0, 1.0, 2.0);
        assert!(check_coords3d(&CoordSeq3d::new(vec![a, same_xy])).is_ok());
        assert_eq!(
            check_coords3d(&CoordSeq3d::new(vec![a, a])),
            Err(crate::GeometryError::ConsecutiveDuplicateVertices { index: 0 })
        );
    }
}
