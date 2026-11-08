// fuzz/fuzz_targets/hausdorff_frechet.rs

//! Fuzz target: the discrete Hausdorff and Fréchet predicates against an
//! exact integer brute force.
//!
//! The input is decoded into two small geometries with integer coordinates
//! below 2²⁰ (points, lines and rectangles, with repeats and axis-aligned
//! runs common) and an integer distance `d`, so `d²` is an exact integer.
//!
//! Checked on every input:
//!
//! - `hausdorff_within(a, b, d)` equals an `i128` brute force over every
//!   vertex and every linework element (a polygon's rings, a line's
//!   segments, a point), and is symmetric;
//! - `hausdorff_distance(a, b)` is `0.0` exactly when the brute force says
//!   every vertex lies on the other's linework;
//! - for two lines, `frechet_within(a, b, d)` equals the coupling program
//!   in `i128`, and `frechet_distance(a, b)` is `0.0` exactly when the
//!   lines are equal.

#![no_main]

use geometres::ops::{frechet_distance, frechet_within, hausdorff_distance, hausdorff_within};
use geometres::types::{
    Coord2d, CoordSeq, Geometry, LineString, LinearRing, MultiPoint, Point, Polygon,
};
use libfuzzer_sys::fuzz_target;

type I2 = (i128, i128);

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

    fn u32(&mut self) -> u32 {
        u32::from_le_bytes([self.u8(), self.u8(), self.u8(), self.u8()])
    }

    fn coord(&mut self, prev: Option<(i64, i64)>) -> (i64, i64) {
        let fresh = |b: &mut Self| i64::from(b.u32() % (1 << 20));
        match (self.u8() % 4, prev) {
            (1, Some((x, _))) => (x, fresh(self)),
            (2, Some((_, y))) => (fresh(self), y),
            (3, Some((x, y))) => (x + i64::from(self.u8() % 8), y + i64::from(self.u8() % 8)),
            _ => (fresh(self), fresh(self)),
        }
    }
}

/// A geometry and, when it is a line, its vertex sequence.
struct Operand {
    geom: Geometry,
    line: Option<Vec<(i64, i64)>>,
    vertices: Vec<I2>,
    linework: Vec<(I2, I2)>,
}

fn c(p: (i64, i64)) -> Coord2d {
    Coord2d::new(p.0 as f64, p.1 as f64)
}

fn i2(p: (i64, i64)) -> I2 {
    (i128::from(p.0), i128::from(p.1))
}

fn decode(b: &mut Bytes<'_>) -> Option<Operand> {
    match b.u8() % 3 {
        0 => {
            let n = 1 + usize::from(b.u8() % 6);
            let mut pts = Vec::new();
            for _ in 0..n {
                let p = b.coord(pts.last().copied());
                pts.push(p);
            }
            Some(Operand {
                geom: Geometry::MultiPoint(MultiPoint::new(
                    pts.iter().map(|&p| Point::new(c(p)).unwrap()).collect(),
                )),
                line: None,
                vertices: pts.iter().map(|&p| i2(p)).collect(),
                linework: pts.iter().map(|&p| (i2(p), i2(p))).collect(),
            })
        }
        1 => {
            let n = 2 + usize::from(b.u8() % 10);
            let mut v: Vec<(i64, i64)> = Vec::new();
            for _ in 0..n {
                let p = b.coord(v.last().copied());
                if v.last() != Some(&p) {
                    v.push(p);
                }
            }
            if v.len() < 2 {
                return None;
            }
            let geom = Geometry::LineString(
                LineString::new(CoordSeq::new(v.iter().map(|&p| c(p)).collect())).unwrap(),
            );
            Some(Operand {
                geom,
                vertices: v.iter().map(|&p| i2(p)).collect(),
                linework: v.windows(2).map(|w| (i2(w[0]), i2(w[1]))).collect(),
                line: Some(v),
            })
        }
        _ => {
            let a = b.coord(None);
            let (w, h) = (1 + i64::from(b.u8()), 1 + i64::from(b.u8()));
            let r = [a, (a.0 + w, a.1), (a.0 + w, a.1 + h), (a.0, a.1 + h)];
            let ring = LinearRing::new(CoordSeq::new(r.iter().map(|&p| c(p)).collect())).unwrap();
            Some(Operand {
                geom: Geometry::Polygon(Polygon::new(ring, vec![]).unwrap()),
                line: None,
                vertices: r.iter().map(|&p| i2(p)).collect(),
                linework: (0..4).map(|k| (i2(r[k]), i2(r[(k + 1) % 4]))).collect(),
            })
        }
    }
}

/// Exactly: is `v` within `√dd` of the closed segment `(a, b)` (a point
/// when `a == b`)?
fn within(v: I2, (a, b): (I2, I2), dd: i128) -> bool {
    let (ux, uy) = (b.0 - a.0, b.1 - a.1);
    let (wx, wy) = (v.0 - a.0, v.1 - a.1);
    let dot = wx * ux + wy * uy;
    let len2 = ux * ux + uy * uy;
    if len2 == 0 || dot <= 0 {
        return wx * wx + wy * wy <= dd;
    }
    if dot >= len2 {
        let (zx, zy) = (v.0 - b.0, v.1 - b.1);
        return zx * zx + zy * zy <= dd;
    }
    let cross = ux * wy - uy * wx;
    cross * cross <= dd * len2
}

fn hausdorff_brute(a: &Operand, b: &Operand, dd: i128) -> bool {
    let dir = |x: &Operand, y: &Operand| {
        x.vertices
            .iter()
            .all(|&v| y.linework.iter().any(|&e| within(v, e, dd)))
    };
    dir(a, b) && dir(b, a)
}

fn frechet_brute(p: &[(i64, i64)], q: &[(i64, i64)], dd: i128) -> bool {
    let sq = |a: (i64, i64), b: (i64, i64)| {
        let (dx, dy) = (i128::from(a.0 - b.0), i128::from(a.1 - b.1));
        dx * dx + dy * dy
    };
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

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let k = i64::from(b.u32() % (1 << 22));
    let (Some(x), Some(y)) = (decode(&mut b), decode(&mut b)) else {
        return;
    };
    let d = k as f64;
    let dd = i128::from(k) * i128::from(k);

    let want = hausdorff_brute(&x, &y, dd);
    assert_eq!(
        hausdorff_within(&x.geom, &y.geom, d).unwrap(),
        want,
        "hausdorff_within"
    );
    assert_eq!(
        hausdorff_within(&y.geom, &x.geom, d).unwrap(),
        want,
        "symmetric"
    );
    let zero = hausdorff_brute(&x, &y, 0);
    assert_eq!(
        hausdorff_distance(&x.geom, &y.geom).unwrap() == 0.0,
        zero,
        "hausdorff zero"
    );

    if let (Some(p), Some(q), Geometry::LineString(la), Geometry::LineString(lb)) =
        (&x.line, &y.line, &x.geom, &y.geom)
    {
        let want = frechet_brute(p, q, dd);
        assert_eq!(frechet_within(la, lb, d).unwrap(), want, "frechet_within");
        assert_eq!(
            frechet_within(lb, la, d).unwrap(),
            want,
            "frechet symmetric"
        );
        assert_eq!(
            frechet_distance(la, lb).unwrap() == 0.0,
            p == q,
            "frechet zero"
        );
    }
});
