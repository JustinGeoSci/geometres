// src/types/polygon.rs

//! Polygon geometry (`Polygon`).
//!
//! A polygon consists of an exterior `LinearRing` bounding its area and zero or
//! more interior `LinearRing` holes. `Polygon::new()` enforces the OGC winding-order
//! convention: the exterior ring must be counterclockwise (positive signed area)
//! and each hole ring must be clockwise (negative signed area).
//!
//! # Degeneracy
//!
//! A ring is degenerate, and rejected as [`GeometryError::DegenerateRing`],
//! **if and only if its exact signed area is zero**: every vertex lies on
//! one line, or the ring's lobes cancel exactly (a symmetric bow tie). This
//! is a permanent contract, decided exactly; it involves no tolerance. A ring
//! with any nonzero area, however thin, is not degenerate. See
//! [`ring_is_degenerate`] for why the test is exact.
//!
//! A ring whose area computation overflows is rejected with
//! [`GeometryError::NonFiniteResult`] rather than judged on a non-finite area.

use crate::GeometryError;
use crate::exact::kernel as exact;
use crate::types::linear_ring::LinearRing;

/// A polygon with an exterior ring and zero or more interior hole rings.
///
/// # Winding Order Convention
///
/// `Polygon::new()` enforces the OGC/simple-features winding convention:
/// - The **exterior** ring must be **counterclockwise** (positive shoelace area).
/// - Each **hole** ring must be **clockwise** (negative shoelace area).
///
/// # Invariants
///
/// - The exterior ring is a valid `LinearRing` (≥ 3 distinct vertices, no
///   consecutive duplicates, all finite coordinates).
/// - The exterior ring has a positive signed area (counterclockwise).
/// - Every hole ring has a negative signed area (clockwise).
/// - No ring (exterior or hole) is degenerate: its exact signed area is nonzero.
/// - **Every hole lies inside the exterior**. Its boundary does not cross the
///   exterior's or share a stretch of it, and the region it bounds is inside
///   the exterior's. Touching the exterior at isolated points is allowed.
/// - **Holes are interior-disjoint**. No hole lies inside, crosses, or
///   coincides with another. Touching at isolated points is allowed.
///
/// The last two make the net signed area of every `Polygon` strictly
/// positive by construction, which `ops::area` relies on. Without them, a
/// polygon whose hole is twenty times its exterior would be constructible
/// and would report a positive area.
///
/// # What construction does not check
///
/// A `Polygon` that exists can still be invalid in the OGC sense, in two
/// ways, and [`ops::validity_issues`](crate::ops::validity_issues) reports
/// both:
///
/// - **A ring may touch or cross itself.** Ring simplicity is checked when a
///   polygon is prepared, not when it is built (a pinch can be split with
///   [`ops::split_pinches`](crate::ops::split_pinches)).
/// - **The interior may be disconnected.** Holes may touch the exterior and
///   each other at points, and a cycle of such touches (a hole touching the
///   exterior twice, or two holes touching each other and each touching the
///   exterior) cuts the interior into pieces. OGC calls that invalid; this
///   constructor builds it. The test that decides it, exactly, is in
///   `ops::validate`.
#[derive(Clone, Debug, PartialEq)]
pub struct Polygon {
    exterior: LinearRing,
    holes: Vec<LinearRing>,
}

impl Polygon {
    /// A polygon from rings already known to form a valid polygon: the rings
    /// of a validated polygon, each rotated, the holes permuted. Rotating a
    /// ring changes neither its point set nor its winding, and the hole
    /// conditions are pairwise and against the exterior, so they do not
    /// depend on order. Crate-private.
    pub(crate) fn from_valid(exterior: LinearRing, holes: Vec<LinearRing>) -> Self {
        Self { exterior, holes }
    }

    /// A polygon assembled *without* validation, for this crate's own unit
    /// tests only: compiled under `cfg(test)` and never in any other build, so
    /// no crate built on this one can obtain an unvalidated `Polygon`.
    #[cfg(test)]
    pub fn from_parts_unvalidated(exterior: LinearRing, holes: Vec<LinearRing>) -> Self {
        Self { exterior, holes }
    }

    /// Construct a `Polygon` from an exterior ring and zero or more hole rings.
    ///
    /// # Winding Order Enforcement
    ///
    /// Winding order is validated via the shoelace formula, whose sign is
    /// certified. Degenerate rings (exact area zero) are rejected before the
    /// winding check (see the module-level note on degeneracy). The exterior
    /// ring is checked first; holes are checked in order, and the first
    /// offending hole's index is reported.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::DegenerateRing`] if the exterior ring or any hole ring
    ///   has an exact signed area of zero.
    /// - [`GeometryError::IncorrectWindingOrder`] if the exterior ring is clockwise
    ///   (negative signed area).
    /// - [`GeometryError::IncorrectHoleWindingOrder`] if any hole ring is
    ///   counterclockwise (positive signed area); `hole_index` names the first
    ///   offending hole.
    /// - [`GeometryError::HoleNotInsideExterior`] if any hole is *not* inside
    ///   the exterior; `hole_index` names the first offending hole.
    /// - [`GeometryError::OverlappingHoles`] if any two holes have
    ///   overlapping interiors; the first offending pair is named.
    ///
    /// # Cost of the containment checks
    ///
    /// Each hole is screened against the exterior's bounding box first, which
    /// is `O(1)` and rejects the gross cases outright: a hole larger than
    /// its exterior never reaches the exact test. A hole whose box does fit is
    /// then tested exactly in `O(n·m)` edge pairs. Hole pairs are screened by
    /// box intersection the same way, so the `O(h²)` pairwise pass does
    /// exact work only for holes whose boxes overlap, which on parcel data is
    /// rare. The exact tests are `arithmetica::orient2d` on stored
    /// coordinates: no intersection point is constructed and no tolerance is
    /// applied.
    pub fn new(exterior: LinearRing, holes: Vec<LinearRing>) -> Result<Self, GeometryError> {
        validate_rings(&exterior, &holes)?;
        Ok(Self { exterior, holes })
    }

    /// As [`Polygon::new_simple`] for an exterior whose signed area the
    /// caller already holds, handing the exterior back when refused.
    /// Crate-private.
    ///
    /// `area` must be `exterior.signed_area()`. The checks `new_simple`
    /// makes are exactly those of `validate_exterior_ring`, and each is a
    /// test of that one value: finite, nonzero (a ring is degenerate exactly
    /// when its signed area is `0.0`, [`ring_is_degenerate`]), and not
    /// negative. Applying them to the value held gives the same verdict
    /// without computing the area a second time.
    pub(crate) fn from_exterior_with_area(
        exterior: LinearRing,
        area: f64,
    ) -> Result<Self, LinearRing> {
        debug_assert_eq!(area.to_bits(), exterior.signed_area().to_bits());
        // Finite, then not degenerate (`!= 0.0`), then not clockwise
        // (`>= 0.0`): together, finite and positive.
        if area.is_finite() && area > 0.0 {
            Ok(Self {
                exterior,
                holes: Vec::new(),
            })
        } else {
            Err(exterior)
        }
    }

    /// As [`Polygon::new`], the same checks in the same order, handing the
    /// exterior back when they refuse, so a caller that reports a refused
    /// exterior need not clone it beforehand. Crate-private.
    pub(crate) fn new_or_exterior(
        exterior: LinearRing,
        holes: Vec<LinearRing>,
    ) -> Result<Self, LinearRing> {
        match validate_rings(&exterior, &holes) {
            Ok(()) => Ok(Self { exterior, holes }),
            Err(_) => Err(exterior),
        }
    }

    /// Construct a `Polygon` with no holes.
    ///
    /// Equivalent to `Polygon::new(exterior, vec![])`.
    ///
    /// # Errors
    ///
    /// See [`Polygon::new`].
    pub fn new_simple(exterior: LinearRing) -> Result<Self, GeometryError> {
        Self::new(exterior, vec![])
    }

    /// Return a reference to the exterior ring.
    pub fn exterior(&self) -> &LinearRing {
        &self.exterior
    }

    /// Return a slice over the hole rings.
    pub fn holes(&self) -> &[LinearRing] {
        &self.holes
    }

    /// Return the number of hole rings.
    pub fn num_holes(&self) -> usize {
        self.holes.len()
    }
}

/// `true` if `ring` is degenerate: its **exact signed area is zero**. This is
/// the test `Polygon::new` applies to every ring. Exposed so that a caller
/// assembling rings into polygons can name the ring that fails, which
/// `GeometryError::DegenerateRing` does not.
///
/// A permanent contract: exactly the rings whose vertices are all collinear,
/// or whose lobes cancel exactly, are degenerate. No tolerance is involved,
/// so a ring with nonzero area is never refused, however thin.
///
/// # Why `signed_area(ring) == 0.0` is exact
///
/// [`LinearRing::signed_area`] returns exactly `0.0` if and only if the
/// exact area is zero:
///
/// - **A compensated result of `0.0` is never certified.** Certification
///   requires `r.abs() > bound` with `bound >= 0` (`sign_is_certified` in
///   `ops/area.rs`), so a zero always falls through to the exact expansion
///   path.
/// - **The expansion path returns `0.0` if and only if the exact sum is
///   zero** (`exact_twice_area_edges` in `ops/area.rs`). A non-overlapping
///   expansion that represents zero has every component zero, because the
///   components occupy disjoint magnitude ranges and cannot cancel. Its sum
///   is `0.0` or `-0.0`, both equal to `0.0`. An expansion that represents a
///   nonzero value has a nonzero largest component, and the other components
///   sum to less than one ulp of it, so summing smallest first cannot round
///   to zero.
///
/// The boundary tests (`tests/degeneracy_boundary.rs`) pin both halves.
///
/// A non-finite area (overflow) is not degenerate; `Polygon::new` refuses
/// such a ring with [`GeometryError::NonFiniteResult`] instead.
pub fn ring_is_degenerate(ring: &LinearRing) -> bool {
    ring.signed_area() == 0.0
}

/// Every check [`Polygon::new`] makes, in its order.
fn validate_rings(exterior: &LinearRing, holes: &[LinearRing]) -> Result<(), GeometryError> {
    validate_exterior_ring(exterior)?;
    for (i, hole) in holes.iter().enumerate() {
        validate_hole_ring(hole, i)?;
    }
    validate_hole_placement(exterior, holes)
}

/// Validate that an exterior ring is non-degenerate and counterclockwise.
///
/// Checks in order:
/// 1. the exact signed area is nonzero (degeneracy guard)
/// 2. `signed_area > 0` (counterclockwise)
fn validate_exterior_ring(ring: &LinearRing) -> Result<(), GeometryError> {
    let area = ring.signed_area();
    if !area.is_finite() {
        return Err(GeometryError::NonFiniteResult);
    }
    if ring_is_degenerate(ring) {
        return Err(GeometryError::DegenerateRing);
    }
    if area < 0.0 {
        return Err(GeometryError::IncorrectWindingOrder);
    }
    Ok(())
}

/// Validate that a hole ring is non-degenerate and clockwise.
///
/// Checks in order:
/// 1. the exact signed area is nonzero (degeneracy guard)
/// 2. `signed_area < 0` (clockwise)
fn validate_hole_ring(ring: &LinearRing, hole_index: usize) -> Result<(), GeometryError> {
    let area = ring.signed_area();
    if !area.is_finite() {
        return Err(GeometryError::NonFiniteResult);
    }
    if ring_is_degenerate(ring) {
        return Err(GeometryError::DegenerateRing);
    }
    if area > 0.0 {
        return Err(GeometryError::IncorrectHoleWindingOrder { hole_index });
    }
    Ok(())
}

/// Validate that every hole lies inside the exterior and that holes are
/// interior-disjoint.
///
/// Order matters for the error reported: every hole is first checked against
/// the exterior (so a hole outside the exterior is reported as such even if it
/// also overlaps another hole), then pairs are checked in index order.
///
/// # Why the envelope screen is sound
///
/// `Envelope::contains_envelope` is closed on all sides and compares stored
/// coordinates exactly, so a ring whose box is *not* inside the other's box
/// has at least one vertex strictly outside the other's box, hence strictly
/// outside the other ring. The screen can therefore only reject rings that
/// the exact test would also reject; it never accepts anything on its own.
/// Likewise two rings whose boxes do not intersect cannot have overlapping
/// interiors.
pub(crate) fn validate_hole_placement(
    exterior: &LinearRing,
    holes: &[LinearRing],
) -> Result<(), GeometryError> {
    if holes.is_empty() {
        return Ok(());
    }
    let ext_env = exterior.envelope()?;
    let mut hole_envs = Vec::with_capacity(holes.len());

    for (i, hole) in holes.iter().enumerate() {
        let env = hole.envelope()?;
        if !ext_env.contains_envelope(&env)
            || exact::place_ring(hole, exterior) != exact::RingPlacement::Inside
        {
            return Err(GeometryError::HoleNotInsideExterior { hole_index: i });
        }
        hole_envs.push(env);
    }

    for i in 0..holes.len() {
        for j in (i + 1)..holes.len() {
            if hole_envs[i].intersects(&hole_envs[j])
                && !exact::rings_interior_disjoint(&holes[i], &holes[j])
            {
                return Err(GeometryError::OverlappingHoles {
                    first: i,
                    second: j,
                });
            }
        }
    }
    Ok(())
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
    use crate::types::coord::{Coord2d, CoordSeq};

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    fn seq(pts: &[(f64, f64)]) -> CoordSeq {
        CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())
    }

    /// Unit square, CCW: a valid exterior ring.
    fn ccw_unit_square() -> LinearRing {
        LinearRing::new(seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])).unwrap()
    }

    /// Unit square, CW: a valid hole ring.
    fn cw_unit_square() -> LinearRing {
        LinearRing::new(seq(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])).unwrap()
    }

    /// Small CW square centered in the unit square: a valid hole.
    fn cw_inner_square() -> LinearRing {
        // CW: go clockwise so signed area < 0.
        LinearRing::new(seq(&[
            (0.25, 0.25),
            (0.25, 0.75),
            (0.75, 0.75),
            (0.75, 0.25),
        ]))
        .unwrap()
    }

    /// Small CCW square centered in the unit square: wrong winding for a hole.
    fn ccw_inner_square() -> LinearRing {
        // CCW: signed area > 0.
        LinearRing::new(seq(&[
            (0.25, 0.25),
            (0.75, 0.25),
            (0.75, 0.75),
            (0.25, 0.75),
        ]))
        .unwrap()
    }

    #[test]
    fn accepts_ccw_exterior_no_holes() {
        let poly = Polygon::new_simple(ccw_unit_square());
        assert!(poly.is_ok());
        let poly = poly.unwrap();
        assert_eq!(poly.num_holes(), 0);
    }

    #[test]
    fn accepts_ccw_exterior_with_cw_hole() {
        let poly = Polygon::new(ccw_unit_square(), vec![cw_inner_square()]);
        assert!(poly.is_ok(), "expected Ok, got: {:?}", poly.unwrap_err());
        let poly = poly.unwrap();
        assert_eq!(poly.num_holes(), 1);
    }

    #[test]
    fn rejects_cw_exterior_ring() {
        let err = Polygon::new_simple(cw_unit_square()).unwrap_err();
        assert_eq!(err, GeometryError::IncorrectWindingOrder);
    }

    #[test]
    fn rejects_ccw_hole_ring() {
        let err = Polygon::new(ccw_unit_square(), vec![ccw_inner_square()]).unwrap_err();
        assert_eq!(
            err,
            GeometryError::IncorrectHoleWindingOrder { hole_index: 0 }
        );
    }

    #[test]
    fn hole_index_identifies_first_bad_hole() {
        // First hole is valid (CW), second hole has wrong winding (CCW).
        let err = Polygon::new(
            ccw_unit_square(),
            vec![cw_inner_square(), ccw_inner_square()],
        )
        .unwrap_err();
        assert_eq!(
            err,
            GeometryError::IncorrectHoleWindingOrder { hole_index: 1 }
        );
    }

    #[test]
    fn rejects_degenerate_exterior_ring() {
        // Three collinear points: the exact area is zero.
        let ring = LinearRing::new(seq(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])).unwrap();
        let err = Polygon::new_simple(ring).unwrap_err();
        assert_eq!(err, GeometryError::DegenerateRing);
    }

    #[test]
    fn new_simple_equivalent_to_new_with_empty_holes() {
        let ring = ccw_unit_square();
        let p1 = Polygon::new_simple(ring.clone()).unwrap();
        let p2 = Polygon::new(ring, vec![]).unwrap();
        assert_eq!(p1, p2);
    }

    #[test]
    fn accessors_return_correct_values() {
        let exterior = ccw_unit_square();
        let hole = cw_inner_square();
        let poly = Polygon::new(exterior.clone(), vec![hole.clone()]).unwrap();
        assert_eq!(*poly.exterior(), exterior);
        assert_eq!(poly.holes(), &[hole]);
        assert_eq!(poly.num_holes(), 1);
    }

    /// A unit square translated to `(o, o)`; its area is exactly 1.
    fn square_at(o: f64, ccw: bool) -> LinearRing {
        let mut pts = [(o, o), (o + 1.0, o), (o + 1.0, o + 1.0), (o, o + 1.0)];
        if !ccw {
            pts.reverse();
        }
        LinearRing::new(CoordSeq::new(
            pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
        ))
        .unwrap()
    }

    #[test]
    fn valid_polygons_far_from_the_origin_are_accepted() {
        // The plain shoelace sum rejected these as degenerate near 1e9.
        for o in [1.0e6, 1.0e8, 1.0e9, 4.0e9, 1.0e12] {
            assert!(
                Polygon::new(square_at(o, true), vec![]).is_ok(),
                "offset {o}"
            );
        }
    }

    #[test]
    fn winding_is_judged_exactly_far_from_the_origin() {
        for o in [1.0e9, 1.0e12] {
            assert_eq!(
                Polygon::new(square_at(o, false), vec![]).unwrap_err(),
                GeometryError::IncorrectWindingOrder
            );
        }
    }

    #[test]
    fn a_thin_needle_far_from_the_origin_is_accepted() {
        // thickness/length = 1e-3 here: a valid ring with nonzero area; the
        // plain sum rejected needles like this at large offsets.
        for o in [0.0, 1.0e6, 1.0e7] {
            let pt = |x: f64, y: f64| Coord2d::new(o + x, o + y);
            let needle = LinearRing::new(CoordSeq::new(vec![
                pt(0.0, 0.0),
                pt(1000.0, 0.0),
                pt(1000.0, 1.0),
                pt(0.0, 1.0),
            ]))
            .unwrap();
            assert!(Polygon::new(needle, vec![]).is_ok(), "offset {o}");
        }
    }

    #[test]
    fn exactly_collinear_ring_is_degenerate_at_any_offset() {
        let o = 1.0e9_f64;
        let r = LinearRing::new(CoordSeq::new(vec![
            Coord2d::new(o, o),
            Coord2d::new(o + 1.0, o + 1.0),
            Coord2d::new(o + 2.0, o + 2.0),
        ]))
        .unwrap();
        assert_eq!(
            Polygon::new(r, vec![]).unwrap_err(),
            GeometryError::DegenerateRing
        );
    }

    #[test]
    fn an_overflowing_area_is_an_error_not_an_accepted_ring() {
        let r = LinearRing::new(CoordSeq::new(vec![
            Coord2d::new(0.0, 0.0),
            Coord2d::new(1e200, 0.0),
            Coord2d::new(1e200, 1e200),
            Coord2d::new(0.0, 1e200),
        ]))
        .unwrap();
        assert_eq!(
            Polygon::new(r, vec![]).unwrap_err(),
            GeometryError::NonFiniteResult
        );
    }
}
