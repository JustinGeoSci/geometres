// src/types3d/polygon3d.rs

//! 3D polygon geometry (`Polygon3d`).
//!
//! # Winding order in 3D
//!
//! Winding is judged from the **Z component of each ring's Newell normal**,
//! i.e. counterclockwise or clockwise as seen looking down the −Z direction
//! from above (+Z). The exterior ring must be counterclockwise from above
//! (`nz > 0`) and each hole clockwise from above (`nz < 0`).
//!
//! # The sign is certified
//!
//! The Newell `z` is `Σ (xᵢ − xⱼ)(yᵢ + yⱼ)` over the edges, which telescopes
//! to the 2D shoelace sum `Σ (xᵢ·yⱼ − xⱼ·yᵢ)`: exactly twice the signed area
//! of the ring's XY projection. So it is computed here by the crate's one
//! certified shoelace routine, `ops::area::certified_twice_area_edges`, on
//! the projected edges. It is a compensated sum whose sign is accepted only
//! when an error bound proves it, and otherwise redone in exact expansion
//! arithmetic. The sign of `nz` is therefore the sign of the exact projected
//! area for every finite ring, however far from the origin it lies.
//!
//! The plain value [`LinearRing3d::newell_normal`]`[2]` is a native `f64`
//! sum whose rounding error grows with the coordinates' magnitude
//! (`ε·|Δx|·|y|` per edge, with no re-centering). At state-plane magnitudes
//! (about `1.4·10⁶` by `5.5·10⁵`) it can give a thin parcel-sized ring a
//! sign the arithmetic did not determine, or sum to an exact `0.0` and refuse
//! the ring as degenerate despite a real area. Both cases are pinned by tests
//! in this file.
//!
//! # Degeneracy
//!
//! A ring is rejected as [`GeometryError::DegenerateRing`] **if and only if
//! the exact area of its XY projection is zero**, the same criterion
//! `Polygon` uses on the 2D shoelace area: `nz`, whose sign is certified, is
//! compared with zero. This is a permanent contract with no tolerance; a
//! ring with any nonzero projected area, however thin, is not degenerate. A
//! ring in a vertical plane
//! has an exact projected area of zero and is refused, because its XY
//! projection has no area from which a +Z winding can be read.
//!
//! # Every ring's XY projection is a valid `LinearRing`
//!
//! Hole placement and every multipart rule are decided on the XY projection
//! with the 2D exact primitives, which take `LinearRing` values. So
//! `Polygon3d::new` projects each ring through `LinearRing::new` and refuses
//! a ring whose projection it rejects (two consecutive vertices equal in
//! `x` and `y`, a vertical edge) as `DegenerateRing`. Such a ring has a
//! zero-length projected edge, which no XY topology can work with.

use crate::GeometryError;
use crate::ops::certified_twice_area_edges;
use crate::types::{CoordSeq, LinearRing, Polygon, validate_hole_placement};
use crate::types3d::LinearRing3d;

/// A 3D polygon with an exterior ring and zero or more interior hole rings.
///
/// # Winding Order Convention (read this before using tilted geometry)
///
/// `Polygon3d::new()` enforces winding relative to the **+Z axis**:
///
/// - The **exterior** ring must be counterclockwise when viewed from above
///   (+Z), i.e. the Z component of its Newell normal must be positive.
/// - Each **hole** ring must be **clockwise when viewed from above**, i.e.
///   the Z component of its Newell normal must be negative.
///
/// The sign of that Z component is **certified**: it is the sign of the
/// exact area of the ring's XY projection (see the module documentation).
///
/// Tilted planes are a limitation. This convention is tied to the global +Z
/// axis, **not** to the polygon's own plane. If the geometry lies in a tilted
/// plane, "counterclockwise from above" may not match the orientation the
/// caller intends (for example, a ring that is counterclockwise as seen from
/// the polygon's own outward normal will be rejected if that normal points
/// into −Z). A ring lying in a vertical plane has `nz = 0` and is always
/// rejected with [`GeometryError::DegenerateRing`], because its XY projection
/// has no area from which a winding order can be read. Callers with tilted
/// geometry must order their vertices so the XY projection is counterclockwise
/// for exteriors and clockwise for holes.
///
/// # Planarity
///
/// Planarity of a ring is **not** verified. The Newell normal is defined for
/// non-planar rings (it is the area-weighted average direction), so a
/// non-planar ring is accepted if its normal passes the checks above.
///
/// # Hole placement
///
/// Every hole must lie inside the exterior and holes must be
/// interior-disjoint, decided exactly on the XY projection with the same
/// primitives `Polygon::new` uses. Touching at a vertex is allowed.
///
/// # Invariants
///
/// - The exterior and every hole are valid `LinearRing3d` values, and each
///   ring's XY projection is a valid `LinearRing`.
/// - The exterior's projected area is positive; every hole's is negative.
/// - No ring is degenerate: every ring's exact projected area is nonzero.
/// - Every hole lies inside the exterior, and holes are interior-disjoint,
///   on the XY projection.
#[derive(Clone, Debug, PartialEq)]
pub struct Polygon3d {
    exterior: LinearRing3d,
    holes: Vec<LinearRing3d>,
}

impl Polygon3d {
    /// Construct a `Polygon3d` from an exterior ring and zero or more hole rings.
    ///
    /// The exterior ring is checked first; holes are checked in order and the
    /// first offending hole's index is reported. For each ring, projection
    /// is checked before degeneracy, and degeneracy before winding order.
    /// Hole placement is checked last, after every ring has passed.
    ///
    /// **Winding is relative to +Z.** The exterior must be counterclockwise
    /// viewed from above (positive Z component of the Newell normal) and holes
    /// clockwise viewed from above (negative Z component). This is the
    /// convention of the global +Z axis rather than the polygon's own plane:
    /// for geometry in a tilted plane, "counterclockwise from above" may not
    /// match the caller's expectation. See the type-level documentation.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::DegenerateRing`] if any ring's exact projected area
    ///   is zero (collinear vertices, or a ring in a vertical plane), or if
    ///   any ring has two consecutive vertices equal in
    ///   `x` and `y`, so that its XY projection is not a valid `LinearRing`.
    /// - [`GeometryError::NonFiniteResult`] if a ring's projected area
    ///   overflows `f64`.
    /// - [`GeometryError::IncorrectWindingOrder`] if the exterior ring is
    ///   clockwise viewed from above.
    /// - [`GeometryError::IncorrectHoleWindingOrder`] if a hole is
    ///   counterclockwise viewed from above; `hole_index` names the first
    ///   offending hole.
    /// - [`GeometryError::HoleNotInsideExterior`] if a hole's XY projection
    ///   is not inside the exterior's; `hole_index` names the first.
    /// - [`GeometryError::OverlappingHoles`] if two holes' XY projections
    ///   have overlapping interiors; the first offending pair is named.
    pub fn new(exterior: LinearRing3d, holes: Vec<LinearRing3d>) -> Result<Self, GeometryError> {
        let exterior_xy = xy_projection(&exterior)?;
        let nz = checked_nz(&exterior_xy)?;
        if nz < 0.0 {
            return Err(GeometryError::IncorrectWindingOrder);
        }
        let mut holes_xy = Vec::with_capacity(holes.len());
        for (i, hole) in holes.iter().enumerate() {
            let hole_xy = xy_projection(hole)?;
            let nz = checked_nz(&hole_xy)?;
            if nz > 0.0 {
                return Err(GeometryError::IncorrectHoleWindingOrder { hole_index: i });
            }
            holes_xy.push(hole_xy);
        }
        validate_hole_placement(&exterior_xy, &holes_xy)?;
        Ok(Self { exterior, holes })
    }

    /// Construct a `Polygon3d` with no holes.
    ///
    /// Equivalent to `Polygon3d::new(exterior, vec![])`.
    ///
    /// # Errors
    ///
    /// See [`Polygon3d::new`].
    pub fn new_simple(exterior: LinearRing3d) -> Result<Self, GeometryError> {
        Self::new(exterior, vec![])
    }

    /// Return a reference to the exterior ring.
    pub fn exterior(&self) -> &LinearRing3d {
        &self.exterior
    }

    /// Return a slice over the hole rings.
    pub fn holes(&self) -> &[LinearRing3d] {
        &self.holes
    }

    /// Return the number of hole rings.
    pub fn num_holes(&self) -> usize {
        self.holes.len()
    }

    /// The polygon's XY projection as a validated 2D `Polygon`, for the
    /// multipart rule (`MultiPolygon3d::new`), which is decided in 2D.
    ///
    /// Every coordinate of the result is an input `x`, `y` pair. The
    /// projection passes `Polygon::new` for every valid `Polygon3d`: the
    /// winding sign is the same certified sum, the degeneracy test is the
    /// same (the exact projected area is nonzero), and hole placement was
    /// decided on this projection. The `Result` is
    /// kept rather than unwrapped so that, should that argument ever fail,
    /// the failure is an error and not a panic or an unchecked value.
    pub(crate) fn xy_projection(&self) -> Result<Polygon, GeometryError> {
        let exterior = xy_projection(&self.exterior)?;
        let holes = self
            .holes
            .iter()
            .map(xy_projection)
            .collect::<Result<Vec<_>, _>>()?;
        Polygon::new(exterior, holes)
    }
}

/// A ring's XY projection as a `LinearRing`, or `DegenerateRing` if
/// `LinearRing::new` refuses it.
///
/// Every coordinate is an input `x`, `y` pair: nothing is computed. The only
/// refusal `LinearRing::new` can raise here is a consecutive duplicate (the
/// 3D ring already has at least three finite, pairwise-distinct consecutive
/// vertices, so after dropping `z` only equality in `x` and `y` is new), and
/// such a ring has a zero-length projected edge: degenerate for every XY
/// purpose this type serves.
fn xy_projection(ring: &LinearRing3d) -> Result<LinearRing, GeometryError> {
    let coords = ring.coords().iter().map(|c| c.to_2d()).collect();
    LinearRing::new(CoordSeq::new(coords)).map_err(|_| GeometryError::DegenerateRing)
}

/// Twice the signed area of the ring's XY projection, the Z component of
/// its Newell normal, with a **certified sign**; or `DegenerateRing` if that
/// area is exactly zero.
///
/// `ring_xy` is the ring's projection; the sum runs over its edges, which
/// are the projected edges of the ring in order. The certified sum returns
/// exactly `0.0` if and only if the exact area is zero (see
/// [`crate::types::ring_is_degenerate`]), so `twice == 0.0` is the exact
/// test. A non-finite sum (overflow) is `NonFiniteResult`, as in
/// `Polygon::new`, rather than a value whose sign comparisons would all be
/// false.
fn checked_nz(ring_xy: &LinearRing) -> Result<f64, GeometryError> {
    let twice = certified_twice_area_edges(ring_xy.edges());
    if !twice.is_finite() {
        return Err(GeometryError::NonFiniteResult);
    }
    if twice == 0.0 {
        return Err(GeometryError::DegenerateRing);
    }
    Ok(twice)
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
    use crate::types3d::{Coord3d, CoordSeq3d};

    fn ring(pts: &[(f64, f64, f64)]) -> LinearRing3d {
        LinearRing3d::new(CoordSeq3d::new(
            pts.iter().map(|&(x, y, z)| Coord3d::new(x, y, z)).collect(),
        ))
        .unwrap()
    }

    fn ccw_square() -> LinearRing3d {
        ring(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
        ])
    }

    fn cw_square() -> LinearRing3d {
        ring(&[
            (0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 0.0),
            (1.0, 0.0, 0.0),
        ])
    }

    fn cw_inner() -> LinearRing3d {
        ring(&[
            (0.25, 0.25, 0.0),
            (0.25, 0.75, 0.0),
            (0.75, 0.75, 0.0),
            (0.75, 0.25, 0.0),
        ])
    }

    fn ccw_inner() -> LinearRing3d {
        ring(&[
            (0.25, 0.25, 0.0),
            (0.75, 0.25, 0.0),
            (0.75, 0.75, 0.0),
            (0.25, 0.75, 0.0),
        ])
    }

    #[test]
    fn accepts_ccw_exterior_with_cw_hole() {
        let p = Polygon3d::new(ccw_square(), vec![cw_inner()]).unwrap();
        assert_eq!(p.num_holes(), 1);
        assert_eq!(p.holes().len(), 1);
        assert_eq!(*p.exterior(), ccw_square());
    }

    #[test]
    fn rejects_cw_exterior_ring() {
        // A clockwise exterior is refused.
        let err = Polygon3d::new_simple(cw_square()).unwrap_err();
        assert_eq!(err, GeometryError::IncorrectWindingOrder);
    }

    #[test]
    fn rejects_ccw_hole_ring_with_index() {
        let err = Polygon3d::new(ccw_square(), vec![cw_inner(), ccw_inner()]).unwrap_err();
        assert_eq!(
            err,
            GeometryError::IncorrectHoleWindingOrder { hole_index: 1 }
        );
    }

    #[test]
    fn rejects_near_zero_normal_ring() {
        // Three collinear points in 3D: Newell normal is the zero vector.
        let r = ring(&[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2.0, 2.0, 2.0)]);
        assert_eq!(r.newell_normal(), [0.0, 0.0, 0.0]);
        assert_eq!(
            Polygon3d::new_simple(r).unwrap_err(),
            GeometryError::DegenerateRing
        );
    }

    #[test]
    fn rejects_degenerate_hole() {
        let r = ring(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]);
        assert_eq!(
            Polygon3d::new(ccw_square(), vec![r]).unwrap_err(),
            GeometryError::DegenerateRing
        );
    }

    #[test]
    fn vertical_ring_is_rejected_as_degenerate() {
        // Square in the plane y = 0: nz == 0 exactly.
        let r = ring(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 0.0, 1.0),
            (0.0, 0.0, 1.0),
        ]);
        assert_eq!(
            Polygon3d::new_simple(r).unwrap_err(),
            GeometryError::DegenerateRing
        );
    }

    #[test]
    fn z_values_do_not_change_winding_for_horizontal_offsets() {
        let r = ring(&[
            (0.0, 0.0, 7.0),
            (1.0, 0.0, 7.0),
            (1.0, 1.0, 7.0),
            (0.0, 1.0, 7.0),
        ]);
        assert!(Polygon3d::new_simple(r).is_ok());
    }

    #[test]
    fn tilted_plane_winding_follows_xy_projection() {
        // Plane z = x. XY projection is CCW for this order -> nz > 0 -> accepted.
        let ccw = ring(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
            (0.0, 1.0, 0.0),
        ]);
        assert!(Polygon3d::new_simple(ccw).is_ok());
        // Reversed order: CW from above -> rejected, documenting the limitation.
        let cw = ring(&[
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 1.0),
            (1.0, 0.0, 1.0),
            (0.0, 0.0, 0.0),
        ]);
        assert_eq!(
            Polygon3d::new_simple(cw).unwrap_err(),
            GeometryError::IncorrectWindingOrder
        );
    }

    // ===================================================================
    // Certified winding, and the defect it closes
    // ===================================================================

    /// Sum of 3D edge lengths of a ring, including the closing edge: the
    /// perimeter the old degeneracy threshold used.
    fn ring_perimeter(ring: &LinearRing3d) -> f64 {
        ring.edges()
            .map(|(a, b)| (b.x - a.x).hypot(b.y - a.y).hypot(b.z - a.z))
            .sum()
    }

    /// The former native `checked_nz`, reproduced here so the defect stays
    /// demonstrable: the native `f64` Newell `z` against a perimeter-based
    /// degeneracy threshold.
    fn native_checked_nz(r: &LinearRing3d) -> Result<f64, GeometryError> {
        let nz = r.newell_normal()[2];
        if (nz * 0.5).abs() < f64::EPSILON * ring_perimeter(r) * ring_perimeter(r) {
            return Err(GeometryError::DegenerateRing);
        }
        Ok(nz)
    }

    /// Two thin triangles at state-plane magnitudes (about `1.45·10⁶` by
    /// `5.5·10⁵`) whose exact projected areas are
    /// positive and 28 and 92 times the degeneracy threshold, with `z`
    /// varying. The native sum gives `nz = −2.98e-8` and `−5.96e-8`
    /// for them: the wrong sign, and far enough past the threshold to be
    /// asserted. Found by a search over near-collinear triangles at these
    /// magnitudes; the coordinates are written as shortest
    /// round-trip decimals and parse to the same `f64` values.
    fn state_plane_slivers() -> [LinearRing3d; 2] {
        [
            ring(&[
                (
                    1_449_864.908_798_434_8,
                    549_999.994_734_418_9,
                    718.970_082_803_438_8,
                ),
                (
                    1_449_729.817_614_475_3,
                    549_999.989_468_837_3,
                    644.439_972_030_012_4,
                ),
                (
                    1_449_999.999_962_733_4,
                    549_999.999_999_999_8,
                    878.229_677_691_913_4,
                ),
            ]),
            ring(&[
                (
                    1_450_364.025_170_730_7,
                    550_000.010_303_409_4,
                    652.734_002_687_942_4,
                ),
                (1_449_999.999_993_060_7, 550_000.0, 570.035_090_337_130_4),
                (
                    1_450_182.012_582_091,
                    550_000.005_151_704_7,
                    574.906_760_078_228_8,
                ),
            ]),
        ]
    }

    fn reversed(r: &LinearRing3d) -> LinearRing3d {
        let mut v: Vec<Coord3d> = r.coords().iter().copied().collect();
        v.reverse();
        LinearRing3d::new(CoordSeq3d::new(v)).unwrap()
    }

    #[test]
    fn state_plane_slivers_are_counterclockwise_and_their_reversals_clockwise() {
        for r in state_plane_slivers() {
            // Positive exact area: accepted as an exterior.
            assert!(Polygon3d::new_simple(r.clone()).is_ok());
            // Reversed: negative exact area, accepted as a hole of a square
            // around it and refused as an exterior.
            let rev = reversed(&r);
            assert_eq!(
                Polygon3d::new_simple(rev.clone()).unwrap_err(),
                GeometryError::IncorrectWindingOrder
            );
            let around = ring(&[
                (1_449_000.0, 549_000.0, 0.0),
                (1_451_000.0, 549_000.0, 0.0),
                (1_451_000.0, 551_000.0, 0.0),
                (1_449_000.0, 551_000.0, 0.0),
            ]);
            assert!(Polygon3d::new(around, vec![rev]).is_ok());
        }
    }

    /// The proof that the defect was real: on both slivers the native sum
    /// returns a negative `nz` past the threshold, so it judges a
    /// counterclockwise ring clockwise.
    #[test]
    fn the_native_sum_asserted_the_wrong_sign_on_the_slivers() {
        for r in state_plane_slivers() {
            let old = native_checked_nz(&r).unwrap();
            assert!(old < 0.0, "native nz = {old:e}");
            // The certified sum, on the same edges, is positive.
            let xy = xy_projection(&r).unwrap();
            assert!(checked_nz(&xy).unwrap() > 0.0);
        }
    }

    /// A thin triangle whose projected area is exactly `16 · 2⁻³³`, 190
    /// times the threshold, for which every term of the native sum rounds
    /// to a value that cancels exactly: it returns `nz = 0.0` and refuses
    /// the ring as degenerate. Both orders.
    #[test]
    fn the_native_sum_refused_a_ring_two_hundred_times_the_threshold() {
        let h = 2f64.powi(-33);
        let ccw = ring(&[
            (1_450_000.0, 550_000.0, 700.0),
            (1_450_016.0, 550_016.0, 650.0),
            (1_450_008.0, 550_008.0 + h, 720.0),
        ]);
        let cw = reversed(&ccw);
        for r in [&ccw, &cw] {
            assert_eq!(
                native_checked_nz(r).unwrap_err(),
                GeometryError::DegenerateRing
            );
            assert_eq!(r.newell_normal()[2], 0.0);
        }
        assert!(Polygon3d::new_simple(ccw).is_ok());
        assert_eq!(
            Polygon3d::new_simple(cw).unwrap_err(),
            GeometryError::IncorrectWindingOrder
        );
    }

    /// The 2D sliver fixture of `ops::area` lifted to 3D: a near-collinear
    /// triangle at `1e9` whose exact projected area is `±0.5`.
    #[test]
    fn a_sliver_at_1e9_has_its_exact_sign_in_either_order() {
        let o = 1.0e9_f64;
        let ccw = ring(&[
            (o, o, 1.0),
            (o + 4.0, o + 1.0, 2.0),
            (o + 3.0, o + 1.0, 3.0),
        ]);
        let cw = ring(&[
            (o, o, 1.0),
            (o + 3.0, o + 1.0, 3.0),
            (o + 4.0, o + 1.0, 2.0),
        ]);
        let xy = xy_projection(&ccw).unwrap();
        assert_eq!(checked_nz(&xy).unwrap(), 1.0);
        let xy = xy_projection(&cw).unwrap();
        assert_eq!(checked_nz(&xy).unwrap(), -1.0);
        assert!(Polygon3d::new_simple(ccw).is_ok());
        assert_eq!(
            Polygon3d::new_simple(cw).unwrap_err(),
            GeometryError::IncorrectWindingOrder
        );
    }

    #[test]
    fn an_exactly_collinear_ring_far_from_the_origin_is_degenerate() {
        let o = 1.0e9_f64;
        let r = ring(&[
            (o, o, 5.0),
            (o + 1.0, o + 1.0, 6.0),
            (o + 2.0, o + 2.0, 7.0),
        ]);
        assert_eq!(
            Polygon3d::new_simple(r).unwrap_err(),
            GeometryError::DegenerateRing
        );
    }

    #[test]
    fn an_overflowing_projected_area_is_non_finite_not_accepted() {
        let r = ring(&[
            (0.0, 0.0, 0.0),
            (1e200, 0.0, 0.0),
            (1e200, 1e200, 0.0),
            (0.0, 1e200, 0.0),
        ]);
        assert_eq!(
            Polygon3d::new_simple(r).unwrap_err(),
            GeometryError::NonFiniteResult
        );
    }

    // ===================================================================
    // Hole placement on the XY projection
    // ===================================================================

    #[test]
    fn a_hole_outside_the_exterior_is_refused() {
        // Clockwise square entirely to the right of the unit square.
        let outside = ring(&[
            (2.0, 0.0, 1.0),
            (2.0, 1.0, 1.0),
            (3.0, 1.0, 1.0),
            (3.0, 0.0, 1.0),
        ]);
        assert_eq!(
            Polygon3d::new(ccw_square(), vec![cw_inner(), outside]).unwrap_err(),
            GeometryError::HoleNotInsideExterior { hole_index: 1 }
        );
    }

    #[test]
    fn a_hole_crossing_the_exterior_is_refused() {
        let crossing = ring(&[
            (0.5, 0.25, 0.0),
            (0.5, 0.75, 0.0),
            (1.5, 0.75, 0.0),
            (1.5, 0.25, 0.0),
        ]);
        assert_eq!(
            Polygon3d::new(ccw_square(), vec![crossing]).unwrap_err(),
            GeometryError::HoleNotInsideExterior { hole_index: 0 }
        );
    }

    #[test]
    fn overlapping_holes_are_refused_with_the_pair() {
        let shifted = ring(&[
            (0.5, 0.5, 9.0),
            (0.5, 0.9, 9.0),
            (0.9, 0.9, 9.0),
            (0.9, 0.5, 9.0),
        ]);
        assert_eq!(
            Polygon3d::new(ccw_square(), vec![cw_inner(), shifted]).unwrap_err(),
            GeometryError::OverlappingHoles {
                first: 0,
                second: 1
            }
        );
    }

    #[test]
    fn two_separate_holes_at_different_heights_are_accepted() {
        let a = ring(&[
            (0.1, 0.1, 0.0),
            (0.1, 0.4, 0.0),
            (0.4, 0.4, 0.0),
            (0.4, 0.1, 0.0),
        ]);
        let b = ring(&[
            (0.6, 0.6, 50.0),
            (0.6, 0.9, 50.0),
            (0.9, 0.9, 50.0),
            (0.9, 0.6, 50.0),
        ]);
        let p = Polygon3d::new(ccw_square(), vec![a, b]).unwrap();
        assert_eq!(p.num_holes(), 2);
    }

    #[test]
    fn a_hole_touching_the_exterior_at_a_vertex_is_accepted() {
        let touching = ring(&[(0.0, 0.0, 0.0), (0.25, 0.5, 0.0), (0.5, 0.25, 0.0)]);
        assert!(Polygon3d::new(ccw_square(), vec![touching]).is_ok());
    }

    /// A vertical edge: two consecutive vertices equal in `x` and `y`. The
    /// ring is a valid `LinearRing3d` (the vertices differ in `z`), but its
    /// XY projection has a zero-length edge and is not a valid `LinearRing`.
    #[test]
    fn a_ring_with_a_vertical_edge_is_degenerate() {
        let r = ring(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 0.0, 5.0),
            (1.0, 1.0, 5.0),
            (0.0, 1.0, 5.0),
        ]);
        assert_eq!(
            Polygon3d::new_simple(r).unwrap_err(),
            GeometryError::DegenerateRing
        );
        // Also as a hole.
        let hole = ring(&[
            (0.25, 0.25, 0.0),
            (0.25, 0.75, 0.0),
            (0.25, 0.75, 3.0),
            (0.75, 0.75, 3.0),
            (0.75, 0.25, 3.0),
        ]);
        assert_eq!(
            Polygon3d::new(ccw_square(), vec![hole]).unwrap_err(),
            GeometryError::DegenerateRing
        );
    }

    #[test]
    fn the_xy_projection_is_the_polygon_without_z() {
        let p = Polygon3d::new(ccw_square(), vec![cw_inner()]).unwrap();
        let q = p.xy_projection().unwrap();
        assert_eq!(q.num_holes(), 1);
        let c: Vec<(f64, f64)> = q.exterior().coords().iter().map(|c| (c.x, c.y)).collect();
        assert_eq!(c, vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
    }

    #[test]
    fn new_simple_equivalent_to_new_with_empty_holes() {
        assert_eq!(
            Polygon3d::new_simple(ccw_square()).unwrap(),
            Polygon3d::new(ccw_square(), vec![]).unwrap()
        );
    }
}
