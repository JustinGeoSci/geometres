// fuzz/fuzz_targets/dwithin.rs

//! Fuzz target: `dwithin` against an exact integer brute force.
//!
//! The input is decoded into two small geometries (points, multipoints,
//! line strings or multilinestrings) with integer coordinates on a grid in
//! `[0, 48]`, where touches, shared vertices, overlaps and exact
//! ties at integer distances are common, and an integer `d` in `[0, 64]`.
//! Puntal and lineal types keep the oracle purely in integers: the oracle
//! decides intersection itself (point-on-segment and proper segment
//! crossing by `i128` signs), then takes the minimum element distance by
//! exact squared comparisons. The oracle never calls back into the crate.
//!
//! The contract checked, for every constructed pair:
//!
//! - `dwithin(a, b, d)` equals the exact brute force, for `d` and for
//!   `d.next_down()` (one ulp below an exact integer tie must flip exactly
//!   when the tie is exact);
//! - `distance(a, b)` is `0.0` exactly when the oracle sees an
//!   intersection, and otherwise lies within `64·u·M` of the oracle's
//!   squared distance (checked by exact `i128` comparison of the bracketing
//!   squares).

#![no_main]

use geometres::ops::{distance, dwithin};
use geometres::types::{
    Coord2d, CoordSeq, Geometry, LineString, MultiLineString, MultiPoint, Point,
};
use libfuzzer_sys::fuzz_target;

/// A cursor over the fuzz bytes; runs out gracefully (zeros).
struct Bytes<'a> {
    data: &'a [u8],
    at: usize,
}

impl Bytes<'_> {
    fn u8(&mut self) -> u8 {
        let b = self.data.get(self.at).copied().unwrap_or(0);
        self.at += 1;
        b
    }

    fn coord(&mut self) -> (i64, i64) {
        (i64::from(self.u8() % 49), i64::from(self.u8() % 49))
    }

    fn geometry(&mut self) -> Option<(Geometry, Vec<Elem>)> {
        let kind = self.u8() % 4;
        let mut elems = Vec::new();
        let g = match kind {
            0 => {
                let (x, y) = self.coord();
                elems.push(Elem::P(x, y));
                Geometry::Point(Point::new(Coord2d::new(x as f64, y as f64)).ok()?)
            }
            1 => {
                let n = 1 + usize::from(self.u8() % 4);
                let mut pts = Vec::new();
                for _ in 0..n {
                    let (x, y) = self.coord();
                    elems.push(Elem::P(x, y));
                    pts.push(Point::new(Coord2d::new(x as f64, y as f64)).ok()?);
                }
                Geometry::MultiPoint(MultiPoint::new(pts))
            }
            2 => {
                let l = self.line(&mut elems)?;
                Geometry::LineString(l)
            }
            _ => {
                let a = self.line(&mut elems)?;
                let b = self.line(&mut elems)?;
                Geometry::MultiLineString(MultiLineString::new(vec![a, b]))
            }
        };
        Some((g, elems))
    }

    fn line(&mut self, elems: &mut Vec<Elem>) -> Option<LineString> {
        let n = 2 + usize::from(self.u8() % 4);
        let pts: Vec<(i64, i64)> = (0..n).map(|_| self.coord()).collect();
        for w in pts.windows(2) {
            elems.push(Elem::S(w[0].0, w[0].1, w[1].0, w[1].1));
        }
        let coords: Vec<Coord2d> = pts
            .iter()
            .map(|&(x, y)| Coord2d::new(x as f64, y as f64))
            .collect();
        LineString::new(CoordSeq::new(coords)).ok()
    }
}

/// An oracle element with integer coordinates.
#[derive(Clone, Copy)]
enum Elem {
    P(i64, i64),
    S(i64, i64, i64, i64),
}

fn sign(v: i128) -> i32 {
    match v {
        0 => 0,
        v if v > 0 => 1,
        _ => -1,
    }
}

/// Whether the point lies on the closed segment, decided exactly.
fn on_segment(px: i64, py: i64, ax: i64, ay: i64, bx: i64, by: i64) -> bool {
    let cross = i128::from(bx - ax) * i128::from(py - ay)
        - i128::from(by - ay) * i128::from(px - ax);
    if cross != 0 {
        return false;
    }
    let dot = i128::from(px - ax) * i128::from(bx - ax)
        + i128::from(py - ay) * i128::from(by - ay);
    let len2 = i128::from(bx - ax) * i128::from(bx - ax)
        + i128::from(by - ay) * i128::from(by - ay);
    if len2 == 0 {
        return px == ax && py == ay;
    }
    (0..=len2).contains(&dot)
}

/// Whether two closed segments meet, decided exactly by orientation signs
/// and the collinear-endpoint cases.
fn segments_meet(a: Elem, b: Elem) -> bool {
    let (Elem::S(a1x, a1y, a2x, a2y), Elem::S(b1x, b1y, b2x, b2y)) = (a, b) else {
        unreachable!("called on segments only");
    };
    let orient = |ox: i64, oy: i64, px: i64, py: i64, qx: i64, qy: i64| {
        sign(
            i128::from(px - ox) * i128::from(qy - oy)
                - i128::from(py - oy) * i128::from(qx - ox),
        )
    };
    let d1 = orient(b1x, b1y, b2x, b2y, a1x, a1y);
    let d2 = orient(b1x, b1y, b2x, b2y, a2x, a2y);
    let d3 = orient(a1x, a1y, a2x, a2y, b1x, b1y);
    let d4 = orient(a1x, a1y, a2x, a2y, b2x, b2y);
    if d1 * d2 < 0 && d3 * d4 < 0 {
        return true;
    }
    on_segment(a1x, a1y, b1x, b1y, b2x, b2y)
        || on_segment(a2x, a2y, b1x, b1y, b2x, b2y)
        || on_segment(b1x, b1y, a1x, a1y, a2x, a2y)
        || on_segment(b2x, b2y, a1x, a1y, a2x, a2y)
}

/// Whether two elements meet, decided exactly.
fn elems_meet(a: Elem, b: Elem) -> bool {
    match (a, b) {
        (Elem::P(px, py), Elem::P(qx, qy)) => px == qx && py == qy,
        (Elem::P(px, py), Elem::S(ax, ay, bx, by))
        | (Elem::S(ax, ay, bx, by), Elem::P(px, py)) => on_segment(px, py, ax, ay, bx, by),
        (Elem::S(..), Elem::S(..)) => segments_meet(a, b),
    }
}

/// Exact comparison of two rationals: `p/q ≤ r/s` (or `<` when `strict`),
/// with every product inside `i128` for this grid and these scales.
fn le(p: i128, q: i128, r: i128, s: i128, strict: bool) -> bool {
    if strict { p * s < r * q } else { p * s <= r * q }
}

/// Exact: `dist(p, seg)² ≤ num/den` (or `<` when `strict`): the clamped
/// projection case split, with the interior case as the rational `cross²/len2`.
#[allow(clippy::too_many_arguments)]
fn ps_sq_le(
    px: i64,
    py: i64,
    ax: i64,
    ay: i64,
    bx: i64,
    by: i64,
    num: i128,
    den: i128,
    strict: bool,
) -> bool {
    let (ux, uy) = (i128::from(bx - ax), i128::from(by - ay));
    let (wx, wy) = (i128::from(px - ax), i128::from(py - ay));
    let dot = wx * ux + wy * uy;
    let len2 = ux * ux + uy * uy;
    if dot <= 0 || len2 == 0 {
        return le(wx * wx + wy * wy, 1, num, den, strict);
    }
    if dot >= len2 {
        let (vx, vy) = (i128::from(px - bx), i128::from(py - by));
        return le(vx * vx + vy * vy, 1, num, den, strict);
    }
    let cross = ux * wy - uy * wx;
    le(cross * cross, len2, num, den, strict)
}

/// Exact pair test: distance² ≤ num/den (or `<`), for elements of
/// geometries the oracle has already found non-intersecting (so the
/// endpoint reduction is complete for segment pairs).
fn pair_sq_le(a: Elem, b: Elem, num: i128, den: i128, strict: bool) -> bool {
    match (a, b) {
        (Elem::P(px, py), Elem::P(qx, qy)) => {
            let (dx, dy) = (i128::from(px - qx), i128::from(py - qy));
            le(dx * dx + dy * dy, 1, num, den, strict)
        }
        (Elem::P(px, py), Elem::S(ax, ay, bx, by))
        | (Elem::S(ax, ay, bx, by), Elem::P(px, py)) => {
            ps_sq_le(px, py, ax, ay, bx, by, num, den, strict)
        }
        (Elem::S(a1x, a1y, a2x, a2y), Elem::S(b1x, b1y, b2x, b2y)) => {
            ps_sq_le(a1x, a1y, b1x, b1y, b2x, b2y, num, den, strict)
                || ps_sq_le(a2x, a2y, b1x, b1y, b2x, b2y, num, den, strict)
                || ps_sq_le(b1x, b1y, a1x, a1y, a2x, a2y, num, den, strict)
                || ps_sq_le(b2x, b2y, a1x, a1y, a2x, a2y, num, den, strict)
        }
    }
}

fuzz_target!(|data: &[u8]| {
    let mut bytes = Bytes { data, at: 0 };
    let d = i64::from(bytes.u8() % 65);
    let Some((ga, ea)) = bytes.geometry() else { return };
    let Some((gb, eb)) = bytes.geometry() else { return };

    let meets = ea.iter().any(|&a| eb.iter().any(|&b| elems_meet(a, b)));
    let within = |num: i128, den: i128, strict: bool| {
        meets
            || ea
                .iter()
                .any(|&a| eb.iter().any(|&b| pair_sq_le(a, b, num, den, strict)))
    };

    // dwithin at the integer d, exact both ways.
    let got = dwithin(&ga, &gb, d as f64).expect("grid inputs are valid");
    let d2 = i128::from(d) * i128::from(d);
    assert_eq!(got, within(d2, 1, false), "dwithin at d = {d}");

    // One ulp below d: always monotone, and an exact tie flips. When the
    // minimum distance is exactly d (within d² but not strictly within),
    // one ulp below must answer false.
    if d > 0 {
        let got_below =
            dwithin(&ga, &gb, (d as f64).next_down()).expect("grid inputs are valid");
        assert!(!got_below || got, "dwithin is monotone in d");
        if got && !within(d2, 1, true) {
            assert!(!got_below, "an exact tie at d = {d} must flip one ulp below");
        }
    }

    // distance: zero exactly on contact; otherwise within the bound, checked
    // by exact bracketing squares: within (dist + δ) and not strictly inside
    // (dist − δ), both as exact rationals. The f64 brackets are scaled by
    // 2⁴⁶ and truncated outward, which widens each bracket by under
    // 2⁻⁴⁵ ≈ 0.06·δ. That widening is lenient only and never raises a false alarm.
    let dist = distance(&ga, &gb).expect("non-empty by construction");
    assert_eq!(dist == 0.0, meets, "distance zero iff the oracle meets");
    if !meets {
        let m = 68.0; // grid magnitude bound (coords ≤ 48; M ≤ 48√2 < 68)
        let delta = 64.0 * (f64::EPSILON / 2.0) * m;
        let scale = 1_i128 << 46;
        #[allow(clippy::cast_possible_truncation)]
        let up = ((dist + delta) * scale as f64).ceil() as i128 + 1;
        assert!(
            within(up * up, scale * scale, false),
            "distance {dist} more than δ above the true distance"
        );
        #[allow(clippy::cast_possible_truncation)]
        let down = ((dist - delta) * scale as f64).floor() as i128 - 1;
        if down > 0 {
            assert!(
                !within(down * down, scale * scale, true),
                "distance {dist} more than δ below the true distance"
            );
        }
    }
});
