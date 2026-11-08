// src/types3d/linear_ring3d.rs

//! Closed 3D ring boundary geometry (`LinearRing3d`).
//!
//! Closure convention is identical to `LinearRing`: exactly `n` distinct
//! vertices are stored and the closing vertex is implied. An OGC-format
//! sequence whose first and last coordinates are identical in x, y, and z
//! has the duplicate closing vertex stripped silently.

use crate::GeometryError;
use crate::types3d::coord3d::check_coords3d;
use crate::types3d::{Coord3d, CoordSeq3d};

/// A closed ring defined by three or more distinct 3D vertices.
///
/// The ring is closed implicitly: an edge from the last vertex back to the
/// first is always present, and the closing vertex is not stored.
///
/// Winding order is not enforced by `LinearRing3d` itself; that is the
/// responsibility of `Polygon3d::new()`. [`LinearRing3d::newell_normal`]
/// exposes the quantity on which that check is based.
///
/// # Invariants
///
/// - At least 3 distinct vertices.
/// - No component of any coordinate is NaN or infinite.
/// - No two consecutive stored vertices are identical in all of x, y, z
///   (including the wrap-around pair after any closing-vertex strip).
#[derive(Clone, Debug, PartialEq)]
pub struct LinearRing3d {
    coords: CoordSeq3d,
}

impl LinearRing3d {
    /// Construct a `LinearRing3d` from a coordinate sequence.
    ///
    /// If `coords[0] == coords[n-1]` (OGC wire-format closing vertex), the
    /// last coordinate is stripped before validation. This is the only
    /// implicit removal.
    ///
    /// Validation order after stripping: non-finite check, minimum-length
    /// check (≥ 3 distinct vertices), wrap-around duplicate check.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::NonFiniteCoordinate`] if any component is NaN or
    ///   infinite; `axis` is `0` for x, `1` for y, `2` for z.
    /// - [`GeometryError::InsufficientVertices`] if fewer than 3 distinct
    ///   vertices remain after stripping.
    /// - [`GeometryError::ConsecutiveDuplicateVertices`] if any consecutive
    ///   pair, including the wrap-around pair, is identical.
    pub fn new(coords: CoordSeq3d) -> Result<Self, GeometryError> {
        let coords = if coords.len() >= 2 && coords[0] == coords[coords.len() - 1] {
            CoordSeq3d::new(coords.as_slice()[..coords.len() - 1].to_vec())
        } else {
            coords
        };
        Self::from_distinct(coords)
    }

    /// As [`LinearRing3d::new`], without stripping a closing vertex: the
    /// native decoder's constructor (see `LinearRing::from_distinct`).
    pub(crate) fn from_distinct(coords: CoordSeq3d) -> Result<Self, GeometryError> {
        check_coords3d(&coords)?;

        if coords.len() < 3 {
            return Err(GeometryError::InsufficientVertices {
                required: 3,
                found: coords.len(),
            });
        }

        // check_coords3d covered pairs 0..n-1; check the implied closing pair.
        let n = coords.len();
        if coords[n - 1] == coords[0] {
            return Err(GeometryError::ConsecutiveDuplicateVertices { index: n - 1 });
        }

        Ok(Self { coords })
    }

    /// Return a reference to the stored distinct vertices (closing vertex not included).
    pub fn coords(&self) -> &CoordSeq3d {
        &self.coords
    }

    /// Return the number of distinct stored vertices.
    pub fn num_vertices(&self) -> usize {
        self.coords.len()
    }

    /// Iterate over all `n` edges of the ring, including the implied closing edge.
    ///
    /// The iterator is `Clone`, as [`LinearRing::edges`](crate::types::LinearRing::edges)
    /// is, for the same reason: the certified winding sum walks the projected
    /// edges a second time when it falls back to exact arithmetic.
    pub fn edges(&self) -> impl Iterator<Item = (Coord3d, Coord3d)> + Clone + '_ {
        let n = self.coords.len();
        (0..n).map(move |i| (self.coords[i], self.coords[(i + 1) % n]))
    }

    /// Return the (unnormalized) Newell normal `[nx, ny, nz]` of the ring.
    ///
    /// For consecutive vertices `p_i`, `p_j` (with `j = i + 1` modulo `n`):
    ///
    /// ```text
    /// nx += (y_i - y_j) * (z_i + z_j)
    /// ny += (z_i - z_j) * (x_i + x_j)
    /// nz += (x_i - x_j) * (y_i + y_j)
    /// ```
    ///
    /// The vector's magnitude is twice the area of the ring's projection onto
    /// the plane it is normal to, and its direction follows the right-hand
    /// rule on the vertex order. The `nz` term telescopes algebraically to the
    /// 2D shoelace sum `Σ (x_i·y_j − x_j·y_i)` over the XY projection, so
    /// `nz > 0` exactly when the XY projection is counterclockwise viewed
    /// from +Z. Unlike a normal built from only the first three vertices, the
    /// Newell sum is well defined for rings with collinear leading vertices
    /// and averages out small out-of-plane noise.
    ///
    /// This computation uses native `f64` arithmetic and is approximate.
    /// Its `z` agrees in sign with `Polygon3d`'s certified winding only
    /// where that certification is not needed, when the projected area is
    /// large against the rounding of this sum, which grows with the
    /// coordinates' magnitude; use `Polygon3d::new` for winding. Nothing in
    /// this crate decides winding from this value.
    pub fn newell_normal(&self) -> [f64; 3] {
        let n = self.coords.len();
        let (mut nx, mut ny, mut nz) = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..n {
            let a = self.coords[i];
            let b = self.coords[(i + 1) % n];
            nx += (a.y - b.y) * (a.z + b.z);
            ny += (a.z - b.z) * (a.x + b.x);
            nz += (a.x - b.x) * (a.y + b.y);
        }
        [nx, ny, nz]
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

    fn square_ccw() -> CoordSeq3d {
        seq(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
        ])
    }

    #[test]
    fn accepts_ogc_format_strips_closing_vertex() {
        let ogc = seq(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0),
        ]);
        let ring = LinearRing3d::new(ogc).unwrap();
        assert_eq!(ring.num_vertices(), 4);
        assert_eq!(ring.coords()[0], Coord3d::new(0.0, 0.0, 0.0));
        assert_eq!(ring.coords()[3], Coord3d::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn closing_vertex_strip_requires_equal_z() {
        // First and last agree in x, y but differ in z: distinct 3D vertices,
        // so nothing is stripped and the wrap-around pair is not a duplicate.
        let pts = seq(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ]);
        let ring = LinearRing3d::new(pts).unwrap();
        assert_eq!(ring.num_vertices(), 4);
    }

    #[test]
    fn rejects_fewer_than_three_distinct_vertices() {
        let err = LinearRing3d::new(seq(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::InsufficientVertices {
                required: 3,
                found: 2
            }
        );
    }

    #[test]
    fn rejects_fewer_than_three_after_ogc_strip() {
        let err = LinearRing3d::new(seq(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)]))
            .unwrap_err();
        assert_eq!(
            err,
            GeometryError::InsufficientVertices {
                required: 3,
                found: 2
            }
        );
    }

    #[test]
    fn rejects_consecutive_duplicate_in_middle() {
        let err = LinearRing3d::new(seq(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
        ]))
        .unwrap_err();
        assert_eq!(
            err,
            GeometryError::ConsecutiveDuplicateVertices { index: 1 }
        );
    }

    #[test]
    fn rejects_non_finite_z() {
        let err = LinearRing3d::new(seq(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, f64::INFINITY),
            (1.0, 1.0, 0.0),
        ]))
        .unwrap_err();
        assert_eq!(
            err,
            GeometryError::NonFiniteCoordinate {
                coord_index: 1,
                axis: 2
            }
        );
    }

    #[test]
    fn edges_yields_n_edges_including_closing() {
        let ring = LinearRing3d::new(square_ccw()).unwrap();
        let edges: Vec<_> = ring.edges().collect();
        assert_eq!(edges.len(), 4);
        assert_eq!(
            edges[3],
            (Coord3d::new(0.0, 1.0, 0.0), Coord3d::new(0.0, 0.0, 0.0))
        );
    }

    #[test]
    fn newell_normal_of_ccw_square_is_plus_z() {
        let ring = LinearRing3d::new(square_ccw()).unwrap();
        // Twice the unit area, along +Z; exact in f64 for these values.
        assert_eq!(ring.newell_normal(), [0.0, 0.0, 2.0]);
    }

    #[test]
    fn newell_normal_of_vertical_ring_has_zero_z() {
        // Square in the plane y = 0 (normal along +/- Y).
        let ring = LinearRing3d::new(seq(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 0.0, 1.0),
            (0.0, 0.0, 1.0),
        ]))
        .unwrap();
        let n = ring.newell_normal();
        assert_eq!(n[2], 0.0);
        assert_eq!(n[1].abs(), 2.0);
    }
}
