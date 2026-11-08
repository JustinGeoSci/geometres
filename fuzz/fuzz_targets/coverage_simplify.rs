// fuzz/fuzz_targets/coverage_simplify.rs

//! Fuzz target: coverage simplification keeps a coverage a coverage.
//!
//! The input is decoded into a small grid coverage (up to 4 × 4 cells of
//! side 64) with integer vertices: interior grid corners jittered, every
//! grid edge a chain of up to 5 jittered points shared by its two cells,
//! some cells dropped. Inputs that are not valid coverages are skipped. A
//! threshold (an area up to about 2^24) and the mode (all edges, or inner)
//! are decoded too.
//!
//! Checked on every valid input:
//!
//! - the call succeeds: a valid coverage is never refused;
//! - the result passes `coverage_validate`, with one polygon per input;
//! - every output ring is an ordered (cyclic) subsequence of its input ring;
//! - shared stretches stay identical: for each pair of polygons, the
//!   vertices of their shared input segments that survive are the same on
//!   both sides;
//! - in inner mode, every vertex of an unshared input segment survives.

#![no_main]

use std::collections::{BTreeMap, BTreeSet};

use geometres::ops::{coverage_simplify, coverage_simplify_inner, coverage_validate};
use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon};
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

    fn jit(&mut self, r: i64) -> i64 {
        i64::from(self.u8()) % (2 * r + 1) - r
    }
}

type P = (i64, i64);
type K = (u64, u64);

fn key(c: Coord2d) -> K {
    ((c.x + 0.0).to_bits(), (c.y + 0.0).to_bits())
}

fn chain(b: &mut Bytes<'_>, a: P, z: P) -> Vec<P> {
    let k = i64::from(b.u8() % 6);
    let mut v = vec![a];
    for q in 1..=k {
        let (t, off) = (q as f64 / (k + 1) as f64, b.jit(5) as f64);
        let (dx, dy) = ((z.0 - a.0) as f64, (z.1 - a.1) as f64);
        let len = dx.hypot(dy).max(1.0);
        v.push((
            (a.0 as f64 + t * dx - off * dy / len).round() as i64,
            (a.1 as f64 + t * dy + off * dx / len).round() as i64,
        ));
    }
    v.push(z);
    v
}

fn decode(b: &mut Bytes<'_>) -> Vec<Polygon> {
    let (n, m) = (1 + usize::from(b.u8() % 4), 1 + usize::from(b.u8() % 4));
    let s = 64_i64;
    let mut vx = vec![vec![(0, 0); m + 1]; n + 1];
    for (i, row) in vx.iter_mut().enumerate() {
        for (q, v) in row.iter_mut().enumerate() {
            let interior = i > 0 && i < n && q > 0 && q < m;
            let (dx, dy) = if interior { (b.jit(7), b.jit(7)) } else { (0, 0) };
            *v = (i as i64 * s + dx, q as i64 * s + dy);
        }
    }
    let h: Vec<Vec<Vec<P>>> = (0..n)
        .map(|i| (0..=m).map(|q| chain(b, vx[i][q], vx[i + 1][q])).collect())
        .collect();
    let v: Vec<Vec<Vec<P>>> = (0..=n)
        .map(|i| (0..m).map(|q| chain(b, vx[i][q], vx[i][q + 1])).collect())
        .collect();
    let mut out = Vec::new();
    for i in 0..n {
        for q in 0..m {
            if b.u8().is_multiple_of(7) {
                continue;
            }
            let mut r: Vec<P> = Vec::new();
            r.extend(&h[i][q][..h[i][q].len() - 1]);
            r.extend(&v[i + 1][q][..v[i + 1][q].len() - 1]);
            let top: Vec<P> = h[i][q + 1].iter().rev().copied().collect();
            r.extend(&top[..top.len() - 1]);
            let left: Vec<P> = v[i][q].iter().rev().copied().collect();
            r.extend(&left[..left.len() - 1]);
            r.dedup();
            while r.len() > 1 && r[0] == r[r.len() - 1] {
                r.pop();
            }
            let pts = r.iter().map(|&(x, y)| Coord2d::new(x as f64, y as f64)).collect();
            if let Ok(ring) = LinearRing::new(CoordSeq::new(pts))
                && let Ok(p) = Polygon::new(ring, vec![])
            {
                out.push(p);
            }
        }
    }
    out
}

fn rings(p: &Polygon) -> Vec<Vec<Coord2d>> {
    std::iter::once(p.exterior())
        .chain(p.holes())
        .map(|r| r.coords().as_slice().to_vec())
        .collect()
}

fn seg(a: Coord2d, b: Coord2d) -> (K, K) {
    let (x, y) = (key(a), key(b));
    if x <= y { (x, y) } else { (y, x) }
}

fn cyclic_subsequence(out: &[Coord2d], inp: &[Coord2d]) -> bool {
    (0..inp.len()).filter(|&s| inp[s] == out[0]).any(|s| {
        let mut q = 0;
        for t in 0..inp.len() {
            if q < out.len() && inp[(s + t) % inp.len()] == out[q] {
                q += 1;
            }
        }
        q == out.len()
    })
}

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let threshold = f64::from(b.u8()) * f64::from(u16::from_le_bytes([b.u8(), b.u8()]));
    let inner = b.u8() % 2 == 1;
    let polys = decode(&mut b);
    if polys.is_empty() || !coverage_validate(&polys).is_empty() {
        return;
    }
    let out = if inner {
        coverage_simplify_inner(&polys, threshold)
    } else {
        coverage_simplify(&polys, threshold)
    }
    .expect("a valid coverage is never refused");
    assert_eq!(out.len(), polys.len());
    assert!(coverage_validate(&out).is_empty(), "the result is a coverage");

    let mut users: BTreeMap<(K, K), BTreeSet<usize>> = BTreeMap::new();
    for (pi, p) in polys.iter().enumerate() {
        for r in rings(p) {
            for t in 0..r.len() {
                users.entry(seg(r[t], r[(t + 1) % r.len()])).or_default().insert(pi);
            }
        }
    }
    let kept: Vec<BTreeSet<K>> = out
        .iter()
        .map(|p| rings(p).iter().flatten().map(|&c| key(c)).collect())
        .collect();
    for (pi, (p, q)) in polys.iter().zip(&out).enumerate() {
        for (a, o) in rings(p).iter().zip(&rings(q)) {
            assert!(o.len() >= 3);
            assert!(cyclic_subsequence(o, a), "output vertices are input vertices, in order");
            if inner {
                for t in 0..a.len() {
                    let s = seg(a[t], a[(t + 1) % a.len()]);
                    if users[&s].len() == 1 {
                        assert!(kept[pi].contains(&s.0) && kept[pi].contains(&s.1));
                    }
                }
            }
        }
    }
    let mut shared: BTreeMap<(usize, usize), BTreeSet<K>> = BTreeMap::new();
    for (&(x, y), u) in &users {
        if u.len() == 2 {
            let v: Vec<usize> = u.iter().copied().collect();
            let e = shared.entry((v[0], v[1])).or_default();
            e.insert(x);
            e.insert(y);
        }
    }
    for (&(i, j), verts) in &shared {
        let si: Vec<&K> = verts.iter().filter(|v| kept[i].contains(v)).collect();
        let sj: Vec<&K> = verts.iter().filter(|v| kept[j].contains(v)).collect();
        assert_eq!(si, sj, "shared vertices differ between {i} and {j}");
    }
});
