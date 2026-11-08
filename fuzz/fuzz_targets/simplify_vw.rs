// fuzz/fuzz_targets/simplify_vw.rs

//! Fuzz target: Visvalingam–Whyatt simplification is vouched, and removes
//! only input vertices.
//!
//! The input is decoded into a small geometry with integer coordinates and
//! a threshold, an area up to about 2^24. The geometry is a line (repeats
//! and axis-aligned runs common), a star-shaped polygon (spiky, optionally
//! with a star-shaped hole), or two such polygons as a multipolygon.
//!
//! Checked on every input:
//!
//! - a returned result has no validity issue, by `validity_issues`;
//! - a refusal is one the documentation lists: the validity refusal or a
//!   constructor's error;
//! - every line and ring of a returned result is an ordered subsequence of
//!   the corresponding input line or ring (cyclically, for a ring), with a
//!   line's endpoints kept and every ring keeping at least three vertices.

#![no_main]

use geometres::GeometryError;
use geometres::ops::{simplify_vw, validity_issues};
use geometres::types::{
    Coord2d, CoordSeq, LineString, LinearRing, MultiPolygon, Polygon, UnvouchedReason,
};
use libfuzzer_sys::fuzz_target;

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

    fn u16(&mut self) -> u16 {
        u16::from_le_bytes([self.u8(), self.u8()])
    }
}

fn c(x: i64, y: i64) -> Coord2d {
    Coord2d::new(x as f64, y as f64)
}

/// A star-shaped ring around `(cx, cy)`: `n` vertices at increasing angles
/// (integer-rounded positions on a circle of 1,000 radius steps), radii in
/// `[lo, hi]`. Counterclockwise unless `cw`. `None` if rounding made it
/// unusable.
fn star(b: &mut Bytes<'_>, cx: i64, cy: i64, lo: i64, hi: i64, cw: bool) -> Option<LinearRing> {
    let n = 3 + usize::from(b.u8() % 14);
    let mut v: Vec<Coord2d> = Vec::with_capacity(n);
    for k in 0..n {
        let angle = (k as f64 + f64::from(b.u8()) / 300.0) * std::f64::consts::TAU / n as f64;
        let r = lo + i64::from(b.u16()) % (hi - lo + 1);
        let p = c(
            cx + (angle.cos() * r as f64).round() as i64,
            cy + (angle.sin() * r as f64).round() as i64,
        );
        if v.last() != Some(&p) {
            v.push(p);
        }
    }
    while v.len() > 1 && v[0] == v[v.len() - 1] {
        v.pop();
    }
    if cw {
        v.reverse();
    }
    LinearRing::new(CoordSeq::new(v)).ok()
}

fn polygon(b: &mut Bytes<'_>, cx: i64, cy: i64) -> Option<Polygon> {
    let ext = star(b, cx, cy, 400, 1_000, false)?;
    let holes = if b.u8().is_multiple_of(2) {
        vec![star(b, cx, cy, 50, 300, true)?]
    } else {
        Vec::new()
    };
    Polygon::new(ext, holes).ok()
}

/// `out` is an ordered subsequence of `input`; cyclically when `ring`.
fn subsequence(out: &[Coord2d], input: &[Coord2d], ring: bool) -> bool {
    let starts: Vec<usize> = if ring {
        (0..input.len()).filter(|&s| input[s] == out[0]).collect()
    } else {
        vec![0]
    };
    starts.into_iter().any(|s| {
        let mut k = 0;
        for j in 0..input.len() {
            if k < out.len() && input[(s + j) % input.len()] == out[k] {
                k += 1;
            }
        }
        k == out.len()
    })
}

fn check_polygon(out: &Polygon, input: &Polygon) {
    let rings = |p: &Polygon| -> Vec<Vec<Coord2d>> {
        std::iter::once(p.exterior())
            .chain(p.holes())
            .map(|r| r.coords().as_slice().to_vec())
            .collect()
    };
    let (o, i) = (rings(out), rings(input));
    assert_eq!(o.len(), i.len(), "rings are kept, not dropped");
    for (a, b) in o.iter().zip(&i) {
        assert!(a.len() >= 3);
        assert!(
            subsequence(a, b, true),
            "a ring gained or reordered a vertex"
        );
    }
}

fn refusal_is_documented(e: &GeometryError) {
    match e {
        GeometryError::Unvouched {
            reason: UnvouchedReason::SimplifiedInvalid { .. },
            ..
        }
        | GeometryError::DegenerateRing
        | GeometryError::IncorrectWindingOrder
        | GeometryError::IncorrectHoleWindingOrder { .. }
        | GeometryError::HoleNotInsideExterior { .. }
        | GeometryError::OverlappingHoles { .. }
        | GeometryError::OverlappingParts { .. } => {}
        other => panic!("undocumented refusal: {other:?}"),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let threshold = f64::from(b.u16()) * f64::from(b.u8());
    match b.u8() % 3 {
        0 => {
            let n = 2 + usize::from(b.u8() % 30);
            let mut v: Vec<Coord2d> = Vec::new();
            for _ in 0..n {
                let p = match (b.u8() % 4, v.last()) {
                    (0, Some(q)) => Coord2d::new(q.x, f64::from(b.u16())),
                    (1, Some(q)) => Coord2d::new(f64::from(b.u16()), q.y),
                    _ => Coord2d::new(f64::from(b.u16()), f64::from(b.u16())),
                };
                if v.last() != Some(&p) {
                    v.push(p);
                }
            }
            let Ok(l) = LineString::new(CoordSeq::new(v)) else {
                return;
            };
            match simplify_vw(&l, threshold) {
                Ok(out) => {
                    assert!(validity_issues(&out).is_empty());
                    let (o, i) = (out.coords().as_slice(), l.coords().as_slice());
                    assert_eq!(o[0], i[0]);
                    assert_eq!(o[o.len() - 1], i[i.len() - 1]);
                    assert!(subsequence(o, i, false));
                }
                Err(e) => refusal_is_documented(&e),
            }
        }
        1 => {
            let Some(p) = polygon(&mut b, 5_000, 5_000) else {
                return;
            };
            match simplify_vw(&p, threshold) {
                Ok(out) => {
                    assert!(validity_issues(&out).is_empty());
                    check_polygon(&out, &p);
                }
                Err(e) => refusal_is_documented(&e),
            }
        }
        _ => {
            let gap = 2_000 + i64::from(b.u8()) * 4;
            let (Some(p), Some(q)) = (
                polygon(&mut b, 5_000, 5_000),
                polygon(&mut b, 5_000 + gap, 5_000),
            ) else {
                return;
            };
            let Ok(mp) = MultiPolygon::new(vec![p, q]) else {
                return;
            };
            match simplify_vw(&mp, threshold) {
                Ok(out) => {
                    assert!(validity_issues(&out).is_empty());
                    for (o, i) in out.polygons().iter().zip(mp.polygons()) {
                        check_polygon(o, i);
                    }
                }
                Err(e) => refusal_is_documented(&e),
            }
        }
    }
});
