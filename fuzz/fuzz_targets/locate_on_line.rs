// fuzz/fuzz_targets/locate_on_line.rs

//! Fuzz target: `PreparedPolygonal::locate_on_line` against `locate`.
//!
//! The input decodes into a mask of up to three parts in separate cells of
//! an integer grid. Rings are star-shaped with up to 90 vertices, so both
//! small and indexed masks occur, or raw grid points, which usually cross
//! themselves and make delegating masks. Each part may have a small square
//! hole. Then a row `y` and a list of `x` values on a half-unit grid, in the
//! input's order, repeats allowed, with signed zeros. The
//! contract checked: for every `i`, `locate_on_line(y, xs)[i]` equals
//! `locate((xs[i], y))`, and neither panics.

#![no_main]

use geometres::prepared::PreparedPolygonal;
use geometres::types::{Coord2d, CoordSeq, LinearRing, MultiPolygon, Polygon};
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

    fn grid(&mut self) -> f64 {
        f64::from(self.u8() % 25)
    }

    /// A value on the half-unit grid over the masks' span, `-0.0` included.
    fn half(&mut self) -> f64 {
        let v = f64::from(self.u8() % 200) * 0.5 - 14.0;
        if v == 0.0 && self.u8() % 2 == 0 { -0.0 } else { v }
    }
}

/// A ring around `(cx, 12)`: star-shaped (vertices in angle order, so
/// usually simple, up to 90 of them, which makes indexed masks), or raw grid
/// points (usually self-crossing, which makes delegating masks).
fn gen_ring(b: &mut Bytes<'_>, cx: f64) -> Option<LinearRing> {
    let mut pts: Vec<Coord2d> = Vec::new();
    if b.u8() % 4 == 0 {
        let n = 3 + usize::from(b.u8() % 10);
        for _ in 0..n {
            pts.push(Coord2d::new(cx - 12.0 + b.grid(), b.grid()));
        }
    } else {
        let n = 3 + usize::from(b.u8() % 88);
        for k in 0..n {
            let t = k as f64 / n as f64 * core::f64::consts::TAU;
            let r = 2.0 + f64::from(b.u8() % 10);
            let p = Coord2d::new((cx + r * t.cos()).round(), (12.0 + r * t.sin()).round());
            if pts.last() != Some(&p) {
                pts.push(p);
            }
        }
        while pts.len() > 1 && pts.first() == pts.last() {
            pts.pop();
        }
    }
    LinearRing::new(CoordSeq::new(pts)).ok()
}

/// A part in its own cell (`cx` apart from the others), with an optional
/// small square hole at its center.
fn gen_part(b: &mut Bytes<'_>, cx: f64) -> Option<Polygon> {
    let ext = gen_ring(b, cx)?;
    if b.u8() % 2 == 0 {
        let s = f64::from(b.u8() % 2) + 0.5;
        let hole = LinearRing::new(CoordSeq::new(vec![
            Coord2d::new(cx - s, 12.0 - s),
            Coord2d::new(cx - s, 12.0 + s),
            Coord2d::new(cx + s, 12.0 + s),
            Coord2d::new(cx + s, 12.0 - s),
        ]))
        .ok()?;
        if let Ok(p) = Polygon::new(ext.clone(), vec![hole]) {
            return Some(p);
        }
    }
    Polygon::new_simple(ext).ok()
}

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let parts: Vec<Polygon> = (0..1 + b.u8() % 3)
        .filter_map(|k| gen_part(&mut b, 30.0 * f64::from(k)))
        .collect();
    let Ok(mask) = MultiPolygon::new(parts) else {
        return;
    };
    let Ok(prepared) = PreparedPolygonal::new(&mask) else {
        return;
    };
    // The row: a half-unit grid over the masks' height.
    let y = f64::from(b.u8() % 56) * 0.5 - 2.0;
    let n = usize::from(b.u8() % 64);
    let xs: Vec<f64> = (0..n).map(|_| b.half()).collect();
    let got = prepared
        .locate_on_line(y, &xs)
        .expect("finite input is never refused");
    assert_eq!(got.len(), xs.len());
    for (i, &x) in xs.iter().enumerate() {
        assert_eq!(
            got[i],
            prepared.locate(Coord2d::new(x, y)),
            "x {x}, y {y}, mask {mask:?}"
        );
    }
});
