// src/types/linear_ring.rs

//! Closed ring boundary geometry (`LinearRing`).
//!
//! # Closure convention
//!
//! Exactly `n` **distinct** vertices are stored. The closing vertex is
//! **implied**: it is the same as `coords[0]`. The closing edge from
//! `coords[n-1]` back to `coords[0]` is produced on demand by `edges()`.
//! This avoids storing a redundant vertex and makes `coords.len()` equal
//! to the number of unique positions.
//!
//! If the caller supplies an OGC-format sequence where `coords[0] == coords[n-1]`
//! (the common wire-format convention), the duplicate closing vertex is stripped
//! silently before construction. This is the only silent removal; consecutive
//! duplicates elsewhere in the sequence remain an error.

use std::sync::OnceLock;

use crate::GeometryError;
use crate::exact::kernel as exact;
use crate::types::coord::{Bbox, Coord2d, CoordSeq, bbox_of, check_coords};
use crate::types::envelope::Envelope;

/// A closed ring defined by three or more distinct vertices.
///
/// The ring is closed implicitly: an edge from the last vertex back to the
/// first vertex is always present. The closing vertex is not stored.
///
/// Winding order is not enforced by `LinearRing` itself; that is the
/// responsibility of `Polygon::new()`. `signed_area()` returns the
/// conventional shoelace value: positive for counterclockwise, negative for
/// clockwise.
///
/// # Invariants
///
/// - At least 3 distinct vertices.
/// - No coordinate value is NaN or infinite.
/// - No two consecutive stored vertices are identical.
///   (The pair `coords[n-1]` / `coords[0]` is checked after any OGC
///   closing vertex has been stripped.)
///
/// # The cached box
///
/// A ring keeps its bounding box once it is first asked for, so
/// [`LinearRing::envelope`] and every operation that screens a ring by its
/// box pay the pass over the vertices once per ring, not once per call.
/// It is computed on first use, not at construction, so an operation that
/// builds rings and never asks for their boxes pays nothing for it. The box
/// is a function of the coordinates; it takes no part in equality or in the
/// `Debug` output, which are those of the coordinates alone.
#[derive(Clone)]
pub struct LinearRing {
    coords: CoordSeq,
    bbox: OnceLock<Bbox>,
}

impl PartialEq for LinearRing {
    fn eq(&self, other: &Self) -> bool {
        self.coords == other.coords
    }
}

// The box is derived from the coordinates and left out on purpose, so the
// output is the coordinates' alone, as it was before the box was cached.
#[allow(clippy::missing_fields_in_debug)]
impl core::fmt::Debug for LinearRing {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("LinearRing")
            .field("coords", &self.coords)
            .finish()
    }
}

impl LinearRing {
    /// A ring from distinct vertices already known valid: a rotation of a
    /// validated ring's vertices, which keeps every property `new` checks
    /// (the same vertices, the same cyclic neighbors). Crate-private.
    pub(crate) fn from_valid(coords: CoordSeq) -> Self {
        debug_assert!(coords.len() >= 3);
        Self {
            coords,
            bbox: OnceLock::new(),
        }
    }

    /// Construct a `LinearRing` from a coordinate sequence.
    ///
    /// If `coords[0] == coords[n-1]` (OGC wire-format closing vertex), the
    /// last coordinate is stripped before validation. This is the only
    /// implicit removal; duplicates elsewhere in the sequence are errors.
    ///
    /// Validation order after stripping: non-finite check, minimum-length
    /// check (≥ 3 distinct vertices), consecutive duplicate check (including
    /// the wrap-around pair `coords[n-1]` / `coords[0]`).
    ///
    /// **Simplicity is not checked**: a constructed ring may touch or cross
    /// itself (the construction contract in [`types`](crate::types)).
    /// [`ops::is_simple`](crate::ops::is_simple) and
    /// [`ops::validity_issues`](crate::ops::validity_issues) decide it
    /// exactly, and [`ops::split_pinches`](crate::ops::split_pinches)
    /// repairs a pinch.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::NonFiniteCoordinate`] if any coordinate value is
    ///   NaN or infinite.
    /// - [`GeometryError::InsufficientVertices`] if fewer than 3 distinct
    ///   vertices remain after stripping.
    /// - [`GeometryError::ConsecutiveDuplicateVertices`] if any non-closing
    ///   consecutive pair is identical (or the wrap-around pair after the
    ///   closing vertex strip).
    pub fn new(coords: CoordSeq) -> Result<Self, GeometryError> {
        // Strip OGC closing vertex if present: coords[0] == coords[n-1].
        let coords = if coords.len() >= 2 && coords[0] == coords[coords.len() - 1] {
            CoordSeq::new(coords.as_slice()[..coords.len() - 1].to_vec())
        } else {
            coords
        };
        Self::from_distinct(coords)
    }

    /// As [`LinearRing::new`], without stripping a closing vertex: `coords`
    /// must already be the ring's distinct vertices, and a last vertex equal
    /// to the first is refused as the duplicate pair across the closing edge.
    ///
    /// The native decoder's constructor. The format writes a ring's distinct
    /// vertices with closure implied, so a stream carrying a closing vertex
    /// is not one the encoder writes; accepting it would give one ring two
    /// encodings, and re-encoding it would not reproduce the bytes read.
    pub(crate) fn from_distinct(coords: CoordSeq) -> Result<Self, GeometryError> {
        // Non-finite check first.
        check_coords(&coords)?;

        // Minimum length: 3 distinct vertices.
        if coords.len() < 3 {
            return Err(GeometryError::InsufficientVertices {
                required: 3,
                found: coords.len(),
            });
        }

        // Wrap-around duplicate check: the closing pair coords[n-1] / coords[0].
        // check_coords already handled consecutive pairs 0..n-1; check the
        // final implied pair here.
        let n = coords.len();
        if coords[n - 1] == coords[0] {
            return Err(GeometryError::ConsecutiveDuplicateVertices { index: n - 1 });
        }

        Ok(Self {
            coords,
            bbox: OnceLock::new(),
        })
    }

    /// Return a reference to the stored distinct vertices (closing vertex not included).
    pub fn coords(&self) -> &CoordSeq {
        &self.coords
    }

    /// Return the number of distinct stored vertices.
    ///
    /// The ring has exactly this many edges (including the implied closing edge).
    pub fn num_vertices(&self) -> usize {
        self.coords.len()
    }

    /// Iterate over all `n` edges of the ring, including the implied closing edge.
    ///
    /// Yields `(start, end)` coordinate pairs. The last yielded pair is
    /// `(coords[n-1], coords[0])`. The iterator is `Clone`, so a caller that
    /// must walk the edges twice (the certified area sum does, when it falls
    /// back to exact arithmetic) can do so without collecting them.
    pub fn edges(&self) -> impl Iterator<Item = (Coord2d, Coord2d)> + Clone + '_ {
        let n = self.coords.len();
        (0..n).map(move |i| (self.coords[i], self.coords[(i + 1) % n]))
    }

    /// Return the signed area of the ring (shoelace formula).
    ///
    /// The sign encodes orientation:
    ///
    /// - **Positive**: counterclockwise (CCW), the required orientation for
    ///   exterior polygon rings.
    /// - **Negative**: clockwise (CW), the required orientation for hole rings.
    /// - **Zero**: the ring has exactly zero area.
    ///
    /// The sign is **exact** for every finite result: the sum is computed in
    /// compensated arithmetic and redone in exact expansion arithmetic when
    /// its sign cannot be certified (see `ops::signed_area_ring`). The
    /// magnitude is accurate to about one rounding of the exact value, except
    /// where a product underflows below about `1e-292`. If an intermediate
    /// product overflows the result is non-finite.
    pub fn signed_area(&self) -> f64 {
        crate::ops::signed_area_ring(self)
    }

    /// The ring's axis-aligned bounding box.
    ///
    /// Cannot fail for a constructed ring: every coordinate is finite by
    /// `check_coords` and the minimum of a nonempty set is never greater than
    /// its maximum. The `Result` is propagated rather than unwrapped because
    /// the crate denies `unwrap` in production code, and because an `Err`
    /// here would mean a `LinearRing` invariant had been violated, which is
    /// worth surfacing rather than masking.
    ///
    /// # Errors
    ///
    /// Unreachable for a valid ring; see above.
    pub fn envelope(&self) -> Result<Envelope, GeometryError> {
        let (min_x, min_y, max_x, max_y) = self.bbox();
        Envelope::new(Coord2d::new(min_x, min_y), Coord2d::new(max_x, max_y))
    }

    /// The box `(min_x, min_y, max_x, max_y)`: the values
    /// [`LinearRing::envelope`] returns, computed on first use and kept.
    pub(crate) fn bbox(&self) -> Bbox {
        *self.bbox.get_or_init(|| bbox_of(self.coords.as_slice()))
    }

    /// Return `true` if the ring has no self-intersections.
    ///
    /// Uses an O(n²) edge-pair test. Adequate for validation and testing; a
    /// sweep-line check is introduced when performance matters.
    ///
    /// # Exactness
    ///
    /// Each pair is decided by `exact::kernel::segments_properly_cross`, which is
    /// four `arithmetica::orient2d` signs on the stored coordinates. A naive
    /// `f64` cross product would make the simplicity validator, and every
    /// validity check routed through it, inexact on the near-degenerate
    /// input where exactness matters most.
    ///
    /// Collinear overlapping edges are not reported as a self-intersection;
    /// consecutive ones are excluded by the no-consecutive-duplicate
    /// invariant and non-consecutive ones are a documented limit.
    pub fn is_simple(&self) -> bool {
        let edges: Vec<(Coord2d, Coord2d)> = self.edges().collect();
        let n = edges.len();
        for i in 0..n {
            for j in (i + 2)..n {
                // Skip adjacent edges (they share a vertex, not a crossing).
                if i == 0 && j == n - 1 {
                    continue;
                }
                if exact::segments_properly_cross(edges[i].0, edges[i].1, edges[j].0, edges[j].1) {
                    return false;
                }
            }
        }
        true
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

    /// A simple unit square, CCW.
    fn unit_square_ccw() -> CoordSeq {
        seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    }

    /// A simple unit square, CW.
    fn unit_square_cw() -> CoordSeq {
        seq(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])
    }

    #[test]
    fn accepts_ogc_format_strips_closing_vertex() {
        // OGC wire format: coords[0] == coords[n-1]
        let ogc = seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]);
        let ring = LinearRing::new(ogc).unwrap();
        // Stored vertices: 4, not 5.
        assert_eq!(ring.num_vertices(), 4);
        assert_eq!(ring.coords()[0], c(0.0, 0.0));
        assert_eq!(ring.coords()[3], c(0.0, 1.0));
    }

    #[test]
    fn rejects_fewer_than_three_distinct_vertices() {
        let err = LinearRing::new(seq(&[(0.0, 0.0), (1.0, 0.0)])).unwrap_err();
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
        // After stripping the closing vertex: only 2 remain.
        let err = LinearRing::new(seq(&[(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::InsufficientVertices {
                required: 3,
                found: 2
            }
        );
    }

    #[test]
    fn edges_yields_n_edges_including_closing() {
        let ring = LinearRing::new(unit_square_ccw()).unwrap();
        let edges: Vec<(Coord2d, Coord2d)> = ring.edges().collect();
        // 4 vertices → 4 edges.
        assert_eq!(edges.len(), 4);
        // Last edge is closing: (0,1) → (0,0).
        assert_eq!(edges[3], (c(0.0, 1.0), c(0.0, 0.0)));
        // First edge.
        assert_eq!(edges[0], (c(0.0, 0.0), c(1.0, 0.0)));
    }

    #[test]
    fn signed_area_positive_for_ccw_ring() {
        let ring = LinearRing::new(unit_square_ccw()).unwrap();
        assert!(
            ring.signed_area() > 0.0,
            "CCW ring must have positive signed area"
        );
    }

    #[test]
    fn signed_area_negative_for_cw_ring() {
        let ring = LinearRing::new(unit_square_cw()).unwrap();
        assert!(
            ring.signed_area() < 0.0,
            "CW ring must have negative signed area"
        );
    }

    #[test]
    fn signed_area_magnitude_correct() {
        let ring = LinearRing::new(unit_square_ccw()).unwrap();
        // Unit square: area = 1.0
        assert!((ring.signed_area() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_consecutive_duplicate_in_middle() {
        let err =
            LinearRing::new(seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 1.0)])).unwrap_err();
        assert_eq!(
            err,
            GeometryError::ConsecutiveDuplicateVertices { index: 1 }
        );
    }

    #[test]
    fn rejects_non_finite_coordinate() {
        let err =
            LinearRing::new(seq(&[(0.0, 0.0), (f64::INFINITY, 0.0), (1.0, 1.0)])).unwrap_err();
        assert!(matches!(
            err,
            GeometryError::NonFiniteCoordinate { axis: 0, .. }
        ));
    }

    /// The box of `cs` by a fresh pass, with explicit comparisons.
    fn fresh_box(cs: &[Coord2d]) -> [u64; 4] {
        let (mut x0, mut y0, mut x1, mut y1) = (cs[0].x, cs[0].y, cs[0].x, cs[0].y);
        for p in cs {
            x0 = x0.min(p.x);
            y0 = y0.min(p.y);
            x1 = x1.max(p.x);
            y1 = y1.max(p.y);
        }
        [x0.to_bits(), y0.to_bits(), x1.to_bits(), y1.to_bits()]
    }

    fn box_bits(b: Bbox) -> [u64; 4] {
        [b.0.to_bits(), b.1.to_bits(), b.2.to_bits(), b.3.to_bits()]
    }

    /// After every constructor, the box is the box of the coordinates: `new`,
    /// `from_valid` (a rotation), the decoder's `from_distinct`, and the
    /// line string's two constructors. The extreme values sit at every
    /// position, the last included.
    #[test]
    fn the_cached_box_is_the_box_of_the_coordinates() {
        use crate::types::{LineString, Polygon};
        let mut state = 0x0B0C_5EEDu64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        for case in 0..2000 {
            let n = 3 + case % 13;
            let pts: Vec<Coord2d> = (0..n)
                .map(|_| c((next() % 41) as f64 - 20.0, (next() % 41) as f64 - 20.0))
                .collect();
            let Ok(ring) = LinearRing::new(CoordSeq::new(pts.clone())) else {
                continue;
            };
            let cs = ring.coords().as_slice();
            assert_eq!(box_bits(ring.bbox()), fresh_box(cs), "case {case}");
            let e = ring.envelope().unwrap();
            assert_eq!(
                [e.min().x, e.min().y, e.max().x, e.max().y].map(f64::to_bits),
                fresh_box(cs)
            );
            let k = case % cs.len();
            let rotated: Vec<Coord2d> = cs[k..].iter().chain(&cs[..k]).copied().collect();
            let r2 = LinearRing::from_valid(CoordSeq::new(rotated.clone()));
            assert_eq!(box_bits(r2.bbox()), fresh_box(&rotated), "case {case}");
            if let Ok(p) = Polygon::new_simple(ring.clone()) {
                let back: Polygon = crate::codec::decode(&crate::codec::encode(&p)).unwrap();
                assert_eq!(
                    box_bits(back.exterior().bbox()),
                    fresh_box(cs),
                    "case {case}"
                );
            }
            let Ok(line) = LineString::new(CoordSeq::new(pts.clone())) else {
                continue;
            };
            assert_eq!(box_bits(line.bbox()), fresh_box(&pts), "case {case}");
            let back: LineString = crate::codec::decode(&crate::codec::encode(&line)).unwrap();
            assert_eq!(box_bits(back.bbox()), fresh_box(&pts), "case {case}");
            let l2 = LineString::from_valid(CoordSeq::new(pts.clone()));
            assert_eq!(box_bits(l2.bbox()), fresh_box(&pts), "case {case}");
        }
    }

    /// The box takes no part in equality or `Debug`, which stay those of
    /// the coordinates.
    #[test]
    fn equality_and_debug_are_the_coordinates_alone() {
        let ring = LinearRing::new(unit_square_ccw()).unwrap();
        assert_eq!(
            format!("{ring:?}"),
            format!("LinearRing {{ coords: {:?} }}", ring.coords())
        );
        let line = crate::types::LineString::new(unit_square_ccw()).unwrap();
        assert_eq!(
            format!("{line:?}"),
            format!("LineString {{ coords: {:?} }}", line.coords())
        );
        assert_eq!(ring, LinearRing::new(unit_square_ccw()).unwrap());
    }
}
