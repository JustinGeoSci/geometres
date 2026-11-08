// src/ops/frechet_tests.rs

//! Tests for `ops::frechet`: the predicate equals an exact `i128` brute
//! force (itself checked against an exhaustive enumeration of couplings on
//! short lines), and the boundary cases answer exactly.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]

use super::{frechet_distance, frechet_within, hausdorff_distance};
use crate::GeometryError;
use crate::types::{Coord2d, CoordSeq, LineString};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
}

fn line(v: &[(i64, i64)]) -> LineString {
    LineString::new(CoordSeq::new(
        v.iter()
            .map(|&(x, y)| Coord2d::new(x as f64, y as f64))
            .collect(),
    ))
    .unwrap()
}

fn fline(v: &[(f64, f64)]) -> LineString {
    LineString::new(CoordSeq::new(
        v.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
    ))
    .unwrap()
}

/// A random line of `n` integer vertices (consecutive repeats removed) in
/// a window of side `s` at `(ox, oy)`.
fn random_line(rng: &mut Rng, n: usize, ox: i64, oy: i64, s: i64) -> Vec<(i64, i64)> {
    let mut v: Vec<(i64, i64)> = Vec::new();
    while v.len() < n {
        let p = (rng.int(ox, ox + s), rng.int(oy, oy + s));
        if v.last() != Some(&p) {
            v.push(p);
        }
    }
    v
}

fn sq(p: (i64, i64), q: (i64, i64)) -> i128 {
    let (dx, dy) = (i128::from(p.0 - q.0), i128::from(p.1 - q.1));
    dx * dx + dy * dy
}

/// The dynamic program in exact integers: is some coupling within `√dd`?
fn brute_within(p: &[(i64, i64)], q: &[(i64, i64)], dd: i128) -> bool {
    let (n, m) = (p.len(), q.len());
    let mut r = vec![vec![false; m]; n];
    for i in 0..n {
        for j in 0..m {
            let from = (i == 0 && j == 0)
                || (i > 0 && r[i - 1][j])
                || (j > 0 && r[i][j - 1])
                || (i > 0 && j > 0 && r[i - 1][j - 1]);
            r[i][j] = from && sq(p[i], q[j]) <= dd;
        }
    }
    r[n - 1][m - 1]
}

/// The exact squared discrete Fréchet distance, by enumerating every
/// coupling. That is the definition itself, practical for short lines.
fn enumerate(p: &[(i64, i64)], q: &[(i64, i64)]) -> i128 {
    fn walk(p: &[(i64, i64)], q: &[(i64, i64)], i: usize, j: usize, worst: i128, best: &mut i128) {
        let worst = worst.max(sq(p[i], q[j]));
        if worst >= *best {
            return;
        }
        if i + 1 == p.len() && j + 1 == q.len() {
            *best = worst;
            return;
        }
        if i + 1 < p.len() {
            walk(p, q, i + 1, j, worst, best);
        }
        if j + 1 < q.len() {
            walk(p, q, i, j + 1, worst, best);
        }
        if i + 1 < p.len() && j + 1 < q.len() {
            walk(p, q, i + 1, j + 1, worst, best);
        }
    }
    let mut best = i128::MAX;
    walk(p, q, 0, 0, 0, &mut best);
    best
}

/// `⌊√x⌋` for a non-negative integer.
fn isqrt(x: i128) -> i64 {
    let mut r = (x as f64).sqrt() as i128;
    while r * r > x {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= x {
        r += 1;
    }
    r as i64
}

// ===========================================================================
// The predicate is exact
// ===========================================================================

#[test]
fn the_predicate_equals_an_i128_brute_force() {
    let mut rng = Rng(0xF2EC_4E7A);
    let top = 1_i64 << 26;
    let (mut cases, mut trues, mut boundary, mut enumerated) = (0_u32, 0_u32, 0_u32, 0_u32);
    while cases < 10_000 {
        let s = [8, 1_000, top - 1][(cases % 3) as usize];
        let (ox, oy) = (rng.int(0, top - 1 - s), rng.int(0, top - 1 - s));
        let n = rng.int(2, 12) as usize;
        let m = rng.int(2, 12) as usize;
        let p = random_line(&mut rng, n, ox, oy, s);
        let q = random_line(&mut rng, m, ox, oy, s);
        // The exact value, by the program on squared values (bisection on
        // the candidate pair distances), then `d` aimed at its root.
        let mut cands: Vec<i128> = p
            .iter()
            .flat_map(|&a| q.iter().map(move |&b| sq(a, b)))
            .collect();
        cands.sort_unstable();
        let true_sq = *cands.iter().find(|&&c| brute_within(&p, &q, c)).unwrap();
        if n <= 6 && m <= 6 {
            assert_eq!(enumerate(&p, &q), true_sq, "the program is the definition");
            enumerated += 1;
        }
        let r = isqrt(true_sq);
        for k in [r - 1, r, r + 1, 0] {
            if k < 0 {
                continue;
            }
            let dd = i128::from(k) * i128::from(k);
            let want = true_sq <= dd;
            assert_eq!(brute_within(&p, &q, dd), want);
            let (a, b) = (line(&p), line(&q));
            assert_eq!(
                frechet_within(&a, &b, k as f64).unwrap(),
                want,
                "d = {k}: {p:?} {q:?}"
            );
            assert_eq!(frechet_within(&b, &a, k as f64).unwrap(), want, "symmetric");
            trues += u32::from(want);
        }
        // An exact boundary: the true value is the integer `r`.
        if r * r == true_sq as i64 && r > 0 {
            boundary += 1;
            let (a, b) = (line(&p), line(&q));
            assert!(frechet_within(&a, &b, r as f64).unwrap());
            assert!(!frechet_within(&a, &b, (r as f64).next_down()).unwrap());
        }
        cases += 1;
    }
    assert!(trues > 5_000, "true answers: {trues}");
    assert!(boundary > 200, "exact integer boundaries: {boundary}");
    assert!(enumerated > 1_000, "enumerated: {enumerated}");
}

/// Exactly: is `S ≤ d²`, for an integer `S` and a non-negative finite `d`?
fn exact_le(s: i128, d: f64) -> bool {
    if d == 0.0 {
        return s == 0;
    }
    let bits = d.to_bits();
    let exp = ((bits >> 52) & 0x7FF) as i64;
    let mant = i128::from(bits & ((1 << 52) - 1)) | if exp == 0 { 0 } else { 1 << 52 };
    let e = if exp == 0 { -1074 } else { exp - 1075 };
    // d = mant · 2^e, so d² = mant² · 2^(2e).
    let m2 = mant * mant;
    if e >= 0 {
        s <= m2 << (2 * e)
    } else {
        s << (-2 * e) <= m2
    }
}

#[test]
fn the_predicate_is_not_the_rounded_value_compared() {
    // Two-vertex lines with offsets near 2^27: the squared distance passes 2^53 and
    // rounds, so `frechet_distance(a, b) <= d` and the exact answer part
    // ways at some `d` next to the computed value.
    let mut rng = Rng(0x0F3A_77E1);
    let mut disagreements = 0_u32;
    for _ in 0..20_000 {
        let p = [
            (rng.int(0, 1 << 20), rng.int(0, 1 << 20)),
            (rng.int(0, 1 << 20), rng.int(0, 1 << 20)),
        ];
        let q = [
            (
                p[0].0 + rng.int(1 << 26, (1 << 27) - 1),
                p[0].1 + rng.int(1 << 26, (1 << 27) - 1),
            ),
            (p[1].0 + 1, p[1].1),
        ];
        if p[0] == p[1] || q[0] == q[1] {
            continue;
        }
        let (a, b) = (line(&p), line(&q));
        let value = frechet_distance(&a, &b).unwrap();
        let mut cands: Vec<i128> = p
            .iter()
            .flat_map(|&x| q.iter().map(move |&y| sq(x, y)))
            .collect();
        cands.sort_unstable();
        let true_sq = *cands.iter().find(|&&c| brute_within(&p, &q, c)).unwrap();
        for d in [value.next_down(), value, value.next_up()] {
            let want = exact_le(true_sq, d);
            assert_eq!(frechet_within(&a, &b, d).unwrap(), want);
            disagreements += u32::from((value <= d) != want);
        }
    }
    assert!(
        disagreements > 0,
        "the rounded comparison must differ somewhere"
    );
}

// ===========================================================================
// The boundary cases
// ===========================================================================

#[test]
fn the_best_coupling_need_not_be_the_diagonal() {
    // The diagonal coupling's worst pair is 6 apart; the best coupling
    // holds p's first four vertices on q's first, then walks q: worst 3.
    let p = line(&[(0, 0), (1, 0), (2, 0), (3, 0), (10, 0)]);
    let q = line(&[(0, 0), (7, 0), (8, 0), (9, 0), (10, 0)]);
    assert!(frechet_within(&p, &q, 3.0).unwrap());
    assert!(!frechet_within(&p, &q, 3.0_f64.next_down()).unwrap());
    assert_eq!(frechet_distance(&p, &q).unwrap(), 3.0);
    // Lines of different lengths have no all-diagonal coupling at all:
    // hold r's first vertex through p's fourth, then step to both ends.
    let r = line(&[(0, 0), (10, 0)]);
    assert!(frechet_within(&p, &r, 3.0).unwrap());
    assert!(!frechet_within(&p, &r, 3.0_f64.next_down()).unwrap());
    assert_eq!(frechet_distance(&p, &r).unwrap(), 3.0);
}

#[test]
fn a_reversed_line_is_far_and_hausdorff_zero() {
    let p = line(&[(0, 0), (3, 0), (3, 4)]);
    let q = line(&[(3, 4), (3, 0), (0, 0)]);
    assert_eq!(frechet_distance(&p, &q).unwrap(), 5.0);
    assert!(frechet_within(&p, &q, 5.0).unwrap());
    assert!(!frechet_within(&p, &q, 5.0_f64.next_down()).unwrap());
    assert_eq!(hausdorff_distance(&p, &q).unwrap(), 0.0);
}

#[test]
fn the_value_is_zero_exactly_for_equal_lines() {
    let p = fline(&[(1e9, 1e9), (1e9 + 0.5, 1e9), (1e9 + 0.5, 1e9 + 3.25)]);
    assert_eq!(frechet_distance(&p, &p.clone()).unwrap(), 0.0);
    assert!(frechet_within(&p, &p, 0.0).unwrap());
    // One vertex moved by one ulp: not zero, and positive.
    let q = fline(&[
        (1e9, 1e9),
        ((1e9_f64 + 0.5).next_up(), 1e9),
        (1e9 + 0.5, 1e9 + 3.25),
    ]);
    assert!(frechet_distance(&p, &q).unwrap() > 0.0);
    assert!(!frechet_within(&p, &q, 0.0).unwrap());
    // A collinear vertex inserted: Hausdorff 0, Fréchet not.
    let r = fline(&[(0.0, 0.0), (10.0, 0.0)]);
    let s = fline(&[(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]);
    assert_eq!(hausdorff_distance(&r, &s).unwrap(), 0.0);
    assert_eq!(frechet_distance(&r, &s).unwrap(), 5.0);
}

#[test]
fn invalid_tolerances_are_refused() {
    let p = line(&[(0, 0), (1, 0)]);
    for d in [-1.0, (-0.0_f64).next_down(), f64::NAN, f64::INFINITY] {
        assert!(matches!(
            frechet_within(&p, &p, d),
            Err(GeometryError::InvalidTolerance { .. })
        ));
    }
}
