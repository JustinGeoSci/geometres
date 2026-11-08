// fuzz/fuzz_targets/delaunay.rs

//! Fuzz target: the Delaunay triangulation against an exact integer brute
//! force.
//!
//! The input is decoded into up to 64 points with integer coordinates below
//! 2²⁰, so the incircle and orientation determinants fit `i128` and the
//! oracle is exact with no arithmetic beyond integers. Duplicates and
//! collinear runs are common by construction (one byte in four reuses or
//! aligns an earlier point).
//!
//! Checked on every input:
//!
//! - **empty circumcircles:** no point lies strictly inside any triangle's
//!   circumcircle (`i128` determinant);
//! - **orientation:** every triangle strictly counterclockwise (`i128`
//!   cross), its first vertex the lexicographically smallest;
//! - **structure:** every distinct point a vertex when triangles exist;
//!   Euler's `edges = points + triangles − 1`; for a degenerate set, the
//!   chain between lexicographically consecutive points;
//! - **canonicality:** the reversed input gives the identical output.

#![no_main]

use geometres::ops::delaunay;
use geometres::types::{Coord2d, MultiPoint, Point};
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

    fn u32(&mut self) -> u32 {
        u32::from_le_bytes([self.u8(), self.u8(), self.u8(), self.u8()])
    }
}

fn decode(data: &[u8]) -> Vec<(i64, i64)> {
    let mut b = Bytes { data, at: 0 };
    let n = 3 + usize::from(b.u8()) % 62;
    let mut pts: Vec<(i64, i64)> = Vec::with_capacity(n);
    for _ in 0..n {
        let kind = b.u8() % 4;
        let p = match (kind, pts.last().copied()) {
            // Reuse an earlier point: duplicates.
            (0, Some(_)) => pts[usize::from(b.u8()) % pts.len()],
            // Align with the previous point: collinear runs and shared axes.
            (1, Some((x, _))) => (x, i64::from(b.u32() % (1 << 20))),
            (2, Some((_, y))) => (i64::from(b.u32() % (1 << 20)), y),
            _ => (
                i64::from(b.u32() % (1 << 20)),
                i64::from(b.u32() % (1 << 20)),
            ),
        };
        pts.push(p);
    }
    pts
}

fn cross(o: (i64, i64), a: (i64, i64), b: (i64, i64)) -> i128 {
    let (ox, oy) = o;
    i128::from(a.0 - ox) * i128::from(b.1 - oy) - i128::from(a.1 - oy) * i128::from(b.0 - ox)
}

/// The incircle determinant for counterclockwise `a`, `b`, `c`: positive if
/// `d` is strictly inside. Coordinates below 2²⁰, so each term is below
/// 2^(42+21+21+2) and the sum fits `i128` with room to spare.
fn in_circle(a: (i64, i64), b: (i64, i64), c: (i64, i64), d: (i64, i64)) -> i128 {
    let row = |p: (i64, i64)| {
        let (dx, dy) = (i128::from(p.0 - d.0), i128::from(p.1 - d.1));
        (dx, dy, dx * dx + dy * dy)
    };
    let (ax, ay, aw) = row(a);
    let (bx, by, bw) = row(b);
    let (cx, cy, cw) = row(c);
    ax * (by * cw - cy * bw) - ay * (bx * cw - cx * bw) + aw * (bx * cy - cx * by)
}

fn as_multipoint(pts: &[(i64, i64)]) -> MultiPoint {
    MultiPoint::new(
        pts.iter()
            .map(|&(x, y)| Point::new(Coord2d::new(x as f64, y as f64)).unwrap())
            .collect(),
    )
}

#[expect(clippy::cast_possible_truncation, reason = "coordinates are small integers")]
fn as_int(c: Coord2d) -> (i64, i64) {
    (c.x as i64, c.y as i64)
}

fuzz_target!(|data: &[u8]| {
    let pts = decode(data);
    let t = delaunay(&as_multipoint(&pts)).unwrap();

    let distinct: Vec<(i64, i64)> = t.points().iter().map(|&c| as_int(c)).collect();
    let n = distinct.len();
    let tris: Vec<[(i64, i64); 3]> = t
        .triangles()
        .iter()
        .map(|tr| tr.map(|v| distinct[v as usize]))
        .collect();

    let mut seen = vec![false; n];
    for tr in t.triangles() {
        let [a, b, c] = tr.map(|v| distinct[v as usize]);
        assert!(cross(a, b, c) > 0, "triangle not counterclockwise");
        let va = t.points()[tr[0] as usize];
        for &v in &tr[1..] {
            let q = t.points()[v as usize];
            assert!(
                (va.x, va.y) < (q.x, q.y),
                "triangle does not start at its smallest vertex"
            );
        }
        for &p in &distinct {
            if p != a && p != b && p != c {
                assert!(in_circle(a, b, c, p) <= 0, "point inside a circumcircle");
            }
        }
        for &v in tr {
            seen[v as usize] = true;
        }
    }
    let edges = t.edges().line_strings().len();
    if tris.is_empty() {
        // Degenerate: the chain between lexicographically consecutive points.
        assert_eq!(edges, n.saturating_sub(1));
    } else {
        assert!(seen.into_iter().all(|s| s), "a point is not a vertex");
        assert_eq!(edges, n + tris.len() - 1, "Euler count");
    }

    // Canonicality: the reversed input gives the identical output.
    let rev: Vec<(i64, i64)> = pts.iter().rev().copied().collect();
    let t2 = delaunay(&as_multipoint(&rev)).unwrap();
    let tris2: Vec<[(i64, i64); 3]> = t2
        .triangles()
        .iter()
        .map(|tr| tr.map(|v| as_int(t2.points()[v as usize])))
        .collect();
    assert_eq!(tris, tris2, "a permutation changed the triangulation");
});
