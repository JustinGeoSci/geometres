// tests/measure.rs

//! The measurements against exact values.
//!
//! - `signed_area` and `area` on integer coordinates against the exact
//!   shoelace sum in `i128`. Where twice the area and every partial sum is
//!   an integer below 2⁵³, the certified sum is exact and the result is
//!   compared bit for bit. On slivers with coordinates near 2³⁰, where a
//!   plain `f64` shoelace sum can get the sign wrong, only the sign is
//!   compared, and it must be exact.
//! - `length` of axis-parallel lines, which is exact.
//! - `distance` is `0.0` exactly when the operands intersect and positive
//!   otherwise; `dwithin` at `0` is `intersects`, and `dwithin` is monotone
//!   in its distance.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss
)]

mod common;

use common::{Rng, geometries, jagged, line, ring};
use geometres::ops::{area, distance, dwithin, length, signed_area};
use geometres::predicates::intersects;
use geometres::types::{Geometry, LinearRing, Polygon};

/// Twice the signed area of a ring of integer coordinates, exactly.
fn twice_area(r: &LinearRing) -> i128 {
    r.edges()
        .map(|(a, b)| {
            let (ax, ay, bx, by) = (a.x as i128, a.y as i128, b.x as i128, b.y as i128);
            ax * by - bx * ay
        })
        .sum()
}

fn twice_area_polygon(p: &Polygon) -> i128 {
    twice_area(p.exterior()) + p.holes().iter().map(twice_area).sum::<i128>()
}

#[test]
fn area_is_exact_on_integer_polygons() {
    let mut rng = Rng::new(0x5EED_0201);
    for k in 0..2_000 {
        let n = 2 + (k % 40);
        let (x, y) = (rng.int(-4_000, 4_000), rng.int(-4_000, 4_000));
        let p = jagged(&mut rng, x, y, n, 3_000);
        let exact = twice_area_polygon(&p);
        assert!(exact > 0);
        assert_eq!(signed_area(&p), exact as f64 * 0.5, "{p:?}");
        assert_eq!(area(&p), exact as f64 * 0.5);
    }
}

#[test]
fn the_sign_of_a_sliver_area_is_exact() {
    // A triangle with vertices near 2³⁰ whose twice-area is a small
    // integer: the products in the shoelace sum are near 2⁶⁰, far past the
    // precision of an `f64`, so a plain sum loses the sign. The certified
    // sum must not.
    let mut rng = Rng::new(0x5EED_0202);
    let base = 1_073_741_824.0;
    let mut checked = 0;
    for _ in 0..2_000 {
        let (dx, dy) = (rng.int(1, 1 << 20), rng.int(1, 1 << 20));
        let k = rng.int(1, 1 << 10);
        // a, then b = a + (dx, dy), then c = a + k·(dx, dy) + (0, s): the
        // third vertex is `s` above the line through a and b.
        let a = (base + rng.int(0, 1 << 20), base + rng.int(0, 1 << 20));
        let b = (a.0 + dx, a.1 + dy);
        let s = rng.int(1, 3);
        let c3 = (a.0 + k * dx, a.1 + k * dy + s);
        let r = ring(&[a, b, c3]);
        let exact = twice_area(&r);
        assert!(exact > 0, "the construction puts c to the left of a → b");
        let Ok(p) = Polygon::new_simple(r) else {
            continue;
        };
        assert!(signed_area(&p) > 0.0, "{p:?}");
        assert!(area(&p) > 0.0);
        checked += 1;
    }
    assert!(checked > 1_900, "{checked}");
}

#[test]
fn length_of_axis_parallel_lines_is_exact() {
    let mut rng = Rng::new(0x5EED_0203);
    for _ in 0..1_000 {
        let mut pts = vec![(rng.int(-1_000, 1_000), rng.int(-1_000, 1_000))];
        let mut want = 0.0;
        for i in 0..=rng.below(20) {
            let (x, y) = pts[pts.len() - 1];
            let step = rng.int(1, 500) * if rng.one_in(2) { 1.0 } else { -1.0 };
            want += step.abs();
            pts.push(if i % 2 == 0 {
                (x + step, y)
            } else {
                (x, y + step)
            });
        }
        assert_eq!(length(&line(&pts)), want);
    }
}

#[test]
fn distance_is_zero_exactly_when_the_operands_intersect() {
    let gs = geometries(0x5EED_0204, 56);
    let (mut zero, mut positive) = (0, 0);
    for a in &gs {
        for b in &gs {
            let d = distance(a, b).unwrap();
            let meet = intersects(a, b).unwrap();
            assert_eq!(d == 0.0, meet, "distance {d}: {a:?}\n{b:?}");
            assert!(d >= 0.0 && d.is_finite());
            assert_eq!(d.to_bits(), distance(b, a).unwrap().to_bits(), "symmetry");
            if meet {
                zero += 1;
            } else {
                positive += 1;
            }
        }
    }
    assert!(
        zero > 200 && positive > 200,
        "zero {zero}, positive {positive}"
    );
}

#[test]
fn dwithin_is_intersects_at_zero_and_monotone() {
    let gs = geometries(0x5EED_0205, 42);
    for a in &gs {
        for b in &gs {
            assert_eq!(
                dwithin(a, b, 0.0).unwrap(),
                intersects(a, b).unwrap(),
                "{a:?}\n{b:?}"
            );
            let mut was = false;
            for d in [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 64.0] {
                let now = dwithin(a, b, d).unwrap();
                assert!(now || !was, "dwithin turned false at {d}: {a:?}\n{b:?}");
                was = now;
            }
            assert!(was, "everything is within 64 here: {a:?}\n{b:?}");
        }
    }
}

#[test]
fn dwithin_agrees_with_distance_away_from_the_threshold() {
    // `distance` carries a documented bound, `dwithin` is exact; where the
    // distance is far from `d` relative to that bound they must agree.
    let gs: Vec<Geometry> = geometries(0x5EED_0206, 42);
    for a in &gs {
        for b in &gs {
            let dist = distance(a, b).unwrap();
            for d in [0.5, 1.0, 1.5, 3.0] {
                if (dist - d).abs() > 1e-9 {
                    assert_eq!(dwithin(a, b, d).unwrap(), dist <= d, "{dist} vs {d}");
                }
            }
        }
    }
}
