// tests/polygon3d_winding_oracle.rs

//! `Polygon3d::new`'s winding decision agrees with
//! an exact oracle on 100,000+ random rings, including near-collinear rings
//! far from the origin.
//!
//! # The oracle
//!
//! Exact integer arithmetic, independent of the expansion arithmetic the
//! crate falls back to. Every coordinate is generated as `m · 2^s` with `m`
//! an integer of at most 51 bits and `s` one exponent per ring, so it is an
//! exact `f64` and the projected shoelace sum `Σ (xᵢ·yⱼ − xⱼ·yᵢ)` equals
//! `2^(2s) · Σ (mxᵢ·myⱼ − mxⱼ·myᵢ)`, whose integer part fits an `i128`
//! (products of 102 bits, at most 2·64 terms). The sign of that integer is
//! the sign of the exact projected area, and that is what the type's
//! certified sum must reproduce.
//!
//! # What is checked
//!
//! Through the public API only. `Polygon3d::new_simple` returns `Ok` for a
//! positive projected area above the degeneracy threshold,
//! `IncorrectWindingOrder` for a negative one, and `DegenerateRing` for a
//! zero area, an area below the threshold `ε·perimeter²`, or a projection
//! `LinearRing::new` refuses. The threshold is an `f64` policy on an
//! accurate magnitude, so near it the exact oracle cannot say which way the
//! policy falls; the test therefore checks:
//!
//! - whenever the result is `Ok` or `IncorrectWindingOrder`, the sign is the
//!   oracle's, unconditionally;
//! - whenever the oracle's area is zero, the result is `DegenerateRing`;
//! - whenever the result is `DegenerateRing` and the projection is a valid
//!   ring, the oracle's area is below twice the threshold;
//! - whenever the oracle's area is above twice the threshold and the
//!   projection is valid, the result is not `DegenerateRing`.

#![cfg(feature = "3d")]
#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation
)]

use geometres::GeometryError;
use geometres::types3d::{Coord3d, CoordSeq3d, LinearRing3d, Polygon3d};
use proptest::prelude::*;
use proptest::test_runner::{Config, TestRunner};

/// A ring given as integer mantissas at one binary exponent.
#[derive(Clone, Debug)]
struct IntRing {
    scale: i32,
    pts: Vec<(i64, i64, i64)>,
}

impl IntRing {
    fn ring3d(&self) -> LinearRing3d {
        let f = |m: i64| (m as f64) * 2f64.powi(self.scale);
        LinearRing3d::new(CoordSeq3d::new(
            self.pts
                .iter()
                .map(|&(x, y, z)| Coord3d::new(f(x), f(y), f(z)))
                .collect(),
        ))
        .unwrap()
    }

    /// `Σ (xᵢ·yⱼ − xⱼ·yᵢ)` over the mantissas, exactly.
    fn twice_area_mantissa(&self) -> i128 {
        let n = self.pts.len();
        (0..n)
            .map(|i| {
                let (ax, ay, _) = self.pts[i];
                let (bx, by, _) = self.pts[(i + 1) % n];
                i128::from(ax) * i128::from(by) - i128::from(bx) * i128::from(ay)
            })
            .sum()
    }

    /// The exact projected area as an `f64` for comparison against the
    /// threshold: `|mantissa| · 2^(2s) / 2`, rounded once (the comparison
    /// below allows a factor of two, so one rounding cannot matter).
    fn area_magnitude(&self) -> f64 {
        (self.twice_area_mantissa().unsigned_abs() as f64) * 2f64.powi(2 * self.scale) * 0.5
    }

    fn xy_projection_is_a_valid_ring(&self) -> bool {
        let n = self.pts.len();
        (0..n).all(|i| {
            let (ax, ay, _) = self.pts[i];
            let (bx, by, _) = self.pts[(i + 1) % n];
            (ax, ay) != (bx, by)
        })
    }

    /// The type's degeneracy threshold `ε·perimeter²`, with the 3D perimeter.
    fn threshold(&self) -> f64 {
        let r = self.ring3d();
        let p: f64 = r
            .edges()
            .map(|(a, b)| (b.x - a.x).hypot(b.y - a.y).hypot(b.z - a.z))
            .sum();
        f64::EPSILON * p * p
    }
}

/// Checks one ring against the oracle; `Err` is the failure description.
fn check(r: &IntRing) -> Result<(), TestCaseError> {
    let mantissa = r.twice_area_mantissa();
    let result = Polygon3d::new_simple(r.ring3d());
    match &result {
        Ok(_) => prop_assert!(
            mantissa > 0,
            "accepted with exact area {mantissa} <= 0: {r:?}"
        ),
        Err(GeometryError::IncorrectWindingOrder) => {
            prop_assert!(
                mantissa < 0,
                "clockwise with exact area {mantissa} >= 0: {r:?}"
            );
        }
        Err(GeometryError::DegenerateRing) => {
            if r.xy_projection_is_a_valid_ring() {
                prop_assert!(
                    mantissa == 0 || r.area_magnitude() < 2.0 * r.threshold(),
                    "degenerate with area {:e} vs threshold {:e}: {r:?}",
                    r.area_magnitude(),
                    r.threshold()
                );
            }
        }
        Err(other) => prop_assert!(false, "unexpected {other:?}: {r:?}"),
    }
    if mantissa == 0 {
        prop_assert_eq!(result.clone().err(), Some(GeometryError::DegenerateRing));
    }
    if r.xy_projection_is_a_valid_ring() && r.area_magnitude() > 2.0 * r.threshold() {
        prop_assert_ne!(result.err(), Some(GeometryError::DegenerateRing));
    }
    Ok(())
}

/// Rings of 3 to 12 vertices with mantissas across the whole 51-bit range
/// at a random exponent: arbitrary shapes, arbitrary magnitudes.
fn arbitrary_rings() -> impl Strategy<Value = IntRing> {
    let m = -(1i64 << 50)..(1i64 << 50);
    (
        -60i32..20,
        prop::collection::vec((m.clone(), m.clone(), m), 3..=12),
    )
        .prop_filter("consecutive 3D duplicates", |(_, pts)| {
            let n = pts.len();
            (0..n).all(|i| pts[i] != pts[(i + 1) % n])
        })
        .prop_map(|(scale, pts)| IntRing { scale, pts })
}

/// Near-collinear rings far from the origin: 3 to 6 vertices along a line
/// through a base point of magnitude up to 2^50, each displaced from the
/// line by at most a few units in the last place, with `z` varying freely.
/// Their exact projected areas are tiny against their coordinates, which is
/// where a native sum has no information about the sign.
fn near_collinear_far_rings() -> impl Strategy<Value = IntRing> {
    let big = -(1i64 << 50)..(1i64 << 50);
    let dir = -(1i64 << 20)..(1i64 << 20);
    let small = -4i64..=4;
    (
        -60i32..20,
        (big.clone(), big),
        (dir.clone(), dir),
        prop::collection::vec((-64i64..64, small, -(1i64 << 30)..(1i64 << 30)), 3..=6),
    )
        .prop_map(|(scale, (bx, by), (dx, dy), steps)| {
            // Point k is base + t_k · (dx, dy) + perpendicular wobble e_k.
            let pts = steps
                .iter()
                .map(|&(t, e, z)| {
                    // Perpendicular to (dx, dy) is (-dy, dx); scaled by e,
                    // which keeps the displacement at e units of the
                    // perpendicular's length, tiny against the base.
                    (
                        bx + t * dx - e * dy.signum(),
                        by + t * dy + e * dx.signum(),
                        z,
                    )
                })
                .collect::<Vec<_>>();
            IntRing { scale, pts }
        })
        .prop_filter("consecutive 3D duplicates", |r| {
            let n = r.pts.len();
            (0..n).all(|i| r.pts[i] != r.pts[(i + 1) % n])
        })
        .prop_filter("mantissa overflow", |r| {
            r.pts
                .iter()
                .all(|&(x, y, _)| x.abs() < (1 << 51) && y.abs() < (1 << 51))
        })
}

fn run(strategy: impl Strategy<Value = IntRing>, cases: u32) {
    let mut runner = TestRunner::new(Config {
        cases,
        ..Config::default()
    });
    runner.run(&strategy, |r| check(&r)).unwrap();
}

#[test]
fn winding_matches_the_exact_oracle_on_arbitrary_rings() {
    run(arbitrary_rings(), 60_000);
}

#[test]
fn winding_matches_the_exact_oracle_on_near_collinear_far_rings() {
    run(near_collinear_far_rings(), 60_000);
}

/// The two state-plane slivers an earlier winding test misjudged (see
/// `types3d::polygon3d`), written as mantissas at `2^-33`, checked against
/// the same oracle: a fixed point in the random suite. These mantissas run
/// to 54 bits, beyond the generators' 51, but each is an `f64` value's own
/// significand shifted, so it converts back exactly, and the products still
/// fit an `i128` with room to spare.
#[test]
fn the_state_plane_slivers_agree_with_the_oracle() {
    let from_f64 = |pts: &[(f64, f64, f64)], scale: i32| {
        let s = 2f64.powi(-scale);
        IntRing {
            scale,
            pts: pts
                .iter()
                .map(|&(x, y, z)| {
                    let m = |v: f64| {
                        let q = v * s;
                        assert_eq!(q, q.trunc(), "{v} is not a multiple of 2^{scale}");
                        q as i64
                    };
                    (m(x), m(y), m(z.trunc()))
                })
                .collect(),
        }
    };
    let a = from_f64(
        &[
            (1_449_864.908_798_434_8, 549_999.994_734_418_9, 718.0),
            (1_449_729.817_614_475_3, 549_999.989_468_837_3, 644.0),
            (1_449_999.999_962_733_4, 549_999.999_999_999_8, 878.0),
        ],
        -33,
    );
    let b = from_f64(
        &[
            (1_450_364.025_170_730_7, 550_000.010_303_409_4, 652.0),
            (1_449_999.999_993_060_7, 550_000.0, 570.0),
            (1_450_182.012_582_091, 550_000.005_151_704_7, 574.0),
        ],
        -33,
    );
    for r in [a, b] {
        assert!(r.twice_area_mantissa() > 0);
        check(&r).unwrap();
        let rev = IntRing {
            scale: r.scale,
            pts: r.pts.iter().rev().copied().collect(),
        };
        assert!(rev.twice_area_mantissa() < 0);
        check(&rev).unwrap();
    }
}
