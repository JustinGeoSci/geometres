// tests/degeneracy_boundary.rs

//! The constructed boundary set for the degeneracy criterion: a ring is
//! degenerate **if and only if its exact signed area is zero**.
//!
//! Every ring here sits at or near the boundary of a common tolerance rule,
//! `|signed_area| < f64::EPSILON × perimeter²`, against which the criterion
//! is compared. Each test states what that tolerance says as well, so the
//! cases where the two disagree are visible in the test's own text. The
//! fixtures use dyadic coordinates throughout, so every vertex is stored
//! exactly and every expected area below is the ring's exact area rather
//! than an approximation of it.

#![allow(clippy::unwrap_used, clippy::float_cmp)]

use geometres::GeometryError;
use geometres::ops::ring_perimeter;
use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon, ring_is_degenerate};

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(
        pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
    ))
    .unwrap()
}

/// Asserts one ring's verdict: `ring_is_degenerate` is `degenerate`, and
/// `Polygon::new` agrees for a counterclockwise ring, refusing it as
/// `DegenerateRing` exactly when the predicate holds.
#[track_caller]
fn check(r: &LinearRing, degenerate: bool, what: &str) {
    assert_eq!(ring_is_degenerate(r), degenerate, "{what}");
    if r.signed_area() >= 0.0 {
        let built = Polygon::new_simple(r.clone());
        if degenerate {
            assert_eq!(
                built.unwrap_err(),
                GeometryError::DegenerateRing,
                "{what}: Polygon::new must refuse"
            );
        } else if let Err(e) = built {
            // Not degenerate; a CCW ring may still be refused for another
            // reason, but never as DegenerateRing.
            assert_ne!(e, GeometryError::DegenerateRing, "{what}: wrong refusal");
        }
    }
}

// ===========================================================================
// Exactly zero area: degenerate
// ===========================================================================

/// Collinear rings, including repeated direction changes along the line,
/// large coordinates, an offset of 2^30, and a bow tie whose lobes cancel.
/// Every exact area is zero, so every one is degenerate (as it is under the
/// tolerance rule), and `signed_area` must return exactly `0.0` on each: the
/// zero half of the exactness argument in `ring_is_degenerate` (a
/// compensated `0.0` never certifies, and the exact expansion of a zero sum
/// is all zeros).
#[test]
fn exactly_collinear_rings_are_degenerate() {
    let o = 2f64.powi(30);
    let cases: [(&str, Vec<(f64, f64)>); 6] = [
        (
            "three on the x-axis",
            vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)],
        ),
        (
            "direction changes on y = x/4",
            vec![(0.0, 0.0), (8.0, 2.0), (4.0, 1.0), (12.0, 3.0), (2.0, 0.5)],
        ),
        (
            "direction changes on y = 3x, mixed magnitudes",
            vec![(0.0, 0.0), (1e8, 3e8), (3e7, 9e7), (7e7, 2.1e8), (2e7, 6e7)],
        ),
        (
            "collinear at offset 2^30",
            vec![(o, o), (o + 8.0, o + 2.0), (o + 4.0, o + 1.0)],
        ),
        (
            "bow tie: the lobes cancel exactly",
            vec![(0.0, 0.0), (10.0, 10.0), (10.0, 0.0), (0.0, 10.0)],
        ),
        (
            "seven vertices on y = x, repeatedly reversing direction",
            vec![
                (0.0, 0.0),
                (16.0, 16.0),
                (4.0, 4.0),
                (32.0, 32.0),
                (8.0, 8.0),
                (64.0, 64.0),
                (2.0, 2.0),
            ],
        ),
    ];
    for (what, pts) in cases {
        let r = ring(&pts);
        assert_eq!(r.signed_area(), 0.0, "{what}: the exact area is zero");
        check(&r, true, what);
    }
}

// ===========================================================================
// Needles near the origin: thickness 2^-44 to 2^-54 of their length
// ===========================================================================

/// Triangles `(0,0), (1024,0), (512, 1024·t)` for dyadic thickness `t`.
/// The exact area is `2^19·t` (the shoelace terms are exact products of
/// dyadic values), and the perimeter is `2048·(1 + δ)` with
/// `δ ≤ t²/2 ≤ 2^-87`, so far below half an ulp of `2048` that the
/// computed perimeter is exactly `2048.0` and the computed threshold is
/// exactly `ε·2048² = 2^-30`. The comparisons below hold with a
/// factor-of-two margin at least, and the `t = 2^-49` case lands exactly on
/// the computed threshold, where the strict `<` accepts: a tie at the
/// boundary goes to the ring, under the threshold as computed.
#[test]
fn needles_of_any_thickness_with_nonzero_area_are_valid() {
    // (t, what the tolerance rule says). Every one is valid, since the area
    // is nonzero.
    let cases: [(i32, bool, &str); 5] = [
        (-44, false, "t = 2^-44 ≈ 5.7e-14: 2^9 times the threshold"),
        (-47, false, "t = 2^-47 ≈ 7.1e-15: 4 times the threshold"),
        (
            -49,
            false,
            "t = 2^-49: area exactly the computed threshold; strict < accepts the tie",
        ),
        (-50, true, "t = 2^-50 ≈ 8.9e-16: half the threshold"),
        (
            -54,
            true,
            "t = 2^-54 ≈ 5.6e-17: refused by the threshold, real by the fact",
        ),
    ];
    for (k, old_tolerance_refused, what) in cases {
        let t = 2f64.powi(k);
        let r = ring(&[(0.0, 0.0), (1024.0, 0.0), (512.0, 1024.0 * t)]);
        assert_eq!(
            r.signed_area(),
            2f64.powi(19) * t,
            "{what}: the area is exact"
        );
        check(&r, false, what);
        // What the tolerance rule says, computed by its own formula: the
        // perimeter is exactly 2048.0 here (see above).
        let p = ring_perimeter(&r);
        assert_eq!(p, 2048.0, "{what}");
        assert_eq!(
            r.signed_area().abs() < f64::EPSILON * p * p,
            old_tolerance_refused,
            "{what}: the tolerance rule's verdict"
        );
    }
}

// ===========================================================================
// At offset 2^30: representability closes the gap for needles...
// ===========================================================================

/// The same needle shape at coordinate offset 2^30 cannot be stored: the
/// apex height `2^-34` is far below half an ulp of `2^30` (which is
/// `2^-23`), so `2^30 + 2^-34` rounds to `2^30` at construction and the
/// stored ring is exactly collinear, and degenerate. Representability, not
/// the degeneracy test, decides this ring.
#[test]
fn a_thin_needle_at_offset_two_to_30_collapses_to_collinear_when_stored() {
    let o = 2f64.powi(30);
    let apex_y = o + 1024.0 * 2f64.powi(-44);
    assert_eq!(apex_y, o, "the apex height is absorbed by the offset");
    let r = ring(&[(o, o), (o + 1024.0, o), (o + 512.0, apex_y)]);
    assert_eq!(r.signed_area(), 0.0);
    check(&r, true, "collapsed needle at offset 2^30");
}

// ===========================================================================
// ...but a determinant-1 sliver keeps the gap open at any magnitude
// ===========================================================================

/// The thinnest non-collinear ring that stored coordinates allow at offset
/// 2^30. Every coordinate in `[2^30, 2^31)` is a multiple of
/// `u = ulp(2^30) = 2^-22`, and the shoelace sum is translation-invariant
/// in exact arithmetic, so the exact doubled area of any stored ring there
/// is an integer multiple of `u² = 2^-44`. The multiple is the integer
/// determinant of the vertex offsets in units of `u`; consecutive pairs
/// with determinant ±1 give `|area| = u²/2 = 2^-45`, about `4e-19` of the
/// ring's length, a thinness no input file is likely to produce, yet real.
///
/// With edge offsets of ~`2^36·u = 2^14` the perimeter is ~`6.7e4`, so the
/// tolerance `ε·P²` is ~`1e-6`. The ring is some `3.6e9` times thinner
/// than the tolerance, which refuses it, and it is valid. The same integer
/// triangle placed at the origin is the near-origin variant.
#[test]
fn a_determinant_one_sliver_is_valid_at_any_magnitude() {
    let u = 2f64.powi(-22);
    let k = 2f64.powi(36);
    // Integer offsets (0,0), (k, k−1), (k−1, k−2): determinant
    // k(k−2) − (k−1)² = −1, so the doubled area is −u²·1 at offset 2^30.
    for (o, unit, what) in [
        (2f64.powi(30), u, "at offset 2^30, in units of ulp(2^30)"),
        (0.0, 1.0, "at the origin, in integer units"),
    ] {
        let r = ring(&[
            (o, o),
            (o + k * unit, o + (k - 1.0) * unit),
            (o + (k - 1.0) * unit, o + (k - 2.0) * unit),
        ]);
        let area = r.signed_area();
        assert_eq!(area, -(unit * unit) * 0.5, "{what}: area is −unit²/2");
        assert!(area != 0.0, "{what}: nonzero, so the fact says valid");
        check(&r, false, what);
    }
}

// ===========================================================================
// Just above the boundary: accepted by both
// ===========================================================================

/// The `t = 2^-48` needle has area `2^-29`, twice the tolerance rule's
/// `2^-30`, so it is valid under both; `Polygon::new` builds it.
#[test]
fn a_needle_above_the_old_tolerance_is_valid() {
    let r = ring(&[(0.0, 0.0), (1024.0, 0.0), (512.0, 1024.0 * 2f64.powi(-48))]);
    check(&r, false, "t = 2^-48");
    assert!(Polygon::new_simple(r).is_ok());
}

// ===========================================================================
// Non-finite area: no verdict change between the criteria
// ===========================================================================

/// An overflowing shoelace sum is a non-finite area. It is not degenerate
/// (`∞` and `NaN` are not equal to `0.0`); `Polygon::new` refuses it as
/// `NonFiniteResult` instead.
#[test]
fn an_overflowing_area_is_not_degenerate() {
    let r = ring(&[(0.0, 0.0), (1e300, 0.0), (1e300, 1e300), (0.0, 1e300)]);
    assert!(!r.signed_area().is_finite());
    assert!(!ring_is_degenerate(&r));
}

// ===========================================================================
// The 3D analogue, on the XY projection
// ===========================================================================

#[cfg(feature = "3d")]
mod three_d {
    use geometres::GeometryError;
    use geometres::types3d::{Coord3d, CoordSeq3d, LinearRing3d, Polygon3d};

    fn ring3(pts: &[(f64, f64, f64)]) -> LinearRing3d {
        LinearRing3d::new(CoordSeq3d::new(
            pts.iter().map(|&(x, y, z)| Coord3d::new(x, y, z)).collect(),
        ))
        .unwrap()
    }

    /// A vertical ring projects to zero area exactly: degenerate.
    #[test]
    fn a_vertical_ring_is_degenerate() {
        let r = ring3(&[
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

    /// The `t = 2^-50` needle with constant `z`: its XY projection is the
    /// 2D case above, with nonzero area, so it is a valid ring (the tolerance
    /// rule refuses it).
    #[test]
    fn a_thin_3d_needle_with_nonzero_projected_area_is_valid() {
        let h = 1024.0 * 2f64.powi(-50);
        let r = ring3(&[(0.0, 0.0, 7.0), (1024.0, 0.0, 7.0), (512.0, h, 7.0)]);
        assert!(Polygon3d::new_simple(r).is_ok());
    }
}
