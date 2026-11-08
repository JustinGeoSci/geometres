// src/ops/coverage_simplify_tests.rs

//! Tests for `ops::coverage_simplify`: results are coverages, shared
//! edges match, nodes stay, inner mode keeps the outline, the veto catches
//! a jump over an island, and the closing vouch refuses a broken rebuild.
//!
//! The property checks here are independent of the module's edge graph:
//! nodes and shared stretches are recomputed from the input's segments
//! directly, and outputs are compared vertex by vertex with the inputs.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

use std::collections::{BTreeMap, BTreeSet};

use super::coverage_simplify::{EdgeGraph, vouch};
use super::{CoverageIssue, coverage_simplify, coverage_simplify_inner, coverage_validate};
use crate::GeometryError;
use crate::types::{Coord2d, CoordSeq, LinearRing, Polygon, UnvouchedReason};

// ===========================================================================
// Fixtures
// ===========================================================================

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

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn ring(v: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(v.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn poly(ext: &[(f64, f64)]) -> Polygon {
    Polygon::new(ring(ext), vec![]).unwrap()
}

type P = (i64, i64);

/// A jittered chain from `a` to `b`: `k` interior points along the segment,
/// each pushed sideways by up to `p`, rounded to integers.
fn chain(rng: &mut Rng, a: P, b: P, k: usize, p: i64) -> Vec<P> {
    let mut v = vec![a];
    let (dx, dy) = ((b.0 - a.0) as f64, (b.1 - a.1) as f64);
    let len = dx.hypot(dy).max(1.0);
    for q in 1..=k {
        let t = q as f64 / (k + 1) as f64;
        let off = rng.int(-p, p) as f64;
        let x = a.0 as f64 + t * dx - off * dy / len;
        let y = a.1 as f64 + t * dy + off * dx / len;
        v.push((x.round() as i64, y.round() as i64));
    }
    v.push(b);
    v
}

/// An `n × m` grid coverage of cells of side `s`: interior grid vertices
/// jittered by up to `s/8`, every grid edge a jittered chain of `k`
/// interior points shared by its two cells, some cells dropped (gaps), and
/// one cell with a hole filled by an island. Cells whose ring came out
/// unusable are skipped; the caller keeps only valid coverages.
fn grid(rng: &mut Rng, n: usize, m: usize, s: i64, k: usize) -> Vec<Polygon> {
    let j = s / 8;
    let mut vx = vec![vec![(0, 0); m + 1]; n + 1];
    for (i, row) in vx.iter_mut().enumerate() {
        for (q, v) in row.iter_mut().enumerate() {
            let interior = i > 0 && i < n && q > 0 && q < m;
            let (dx, dy) = if interior {
                (rng.int(-j, j), rng.int(-j, j))
            } else {
                (0, 0)
            };
            *v = (i as i64 * s + dx, q as i64 * s + dy);
        }
    }
    let p = s / 10;
    let kk = |rng: &mut Rng| rng.int(0, k as i64) as usize;
    let mut h = vec![vec![Vec::new(); m + 1]; n];
    for (i, col) in h.iter_mut().enumerate() {
        for (q, e) in col.iter_mut().enumerate() {
            let kq = kk(rng);
            *e = chain(rng, vx[i][q], vx[i + 1][q], kq, p);
        }
    }
    let mut v = vec![vec![Vec::new(); m]; n + 1];
    for (i, col) in v.iter_mut().enumerate() {
        for (q, e) in col.iter_mut().enumerate() {
            let kq = kk(rng);
            *e = chain(rng, vx[i][q], vx[i][q + 1], kq, p);
        }
    }
    let island_cell = (
        rng.int(0, n as i64 - 1) as usize,
        rng.int(0, m as i64 - 1) as usize,
    );
    let mut out = Vec::new();
    for i in 0..n {
        for q in 0..m {
            if (i, q) != island_cell && rng.int(0, 9) == 0 {
                continue; // a gap
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
            let pts: Vec<(f64, f64)> = r.iter().map(|&(x, y)| (x as f64, y as f64)).collect();
            let Ok(ext) =
                LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect()))
            else {
                continue;
            };
            if (i, q) == island_cell {
                // A hole around the cell's center, filled by an island.
                let (cx, cy) = (i as i64 * s + s / 2, q as i64 * s + s / 2);
                let rr = s / 6;
                let mut hole: Vec<P> = vec![
                    (cx, cy - rr),
                    (cx + rr / 3, cy - rr / 3 + rng.int(-1, 1)),
                    (cx + rr, cy),
                    (cx, cy + rr),
                    (cx - rr, cy + rng.int(-2, 2)),
                ];
                hole.dedup();
                let island: Vec<(f64, f64)> =
                    hole.iter().map(|&(x, y)| (x as f64, y as f64)).collect();
                let hole_cw: Vec<(f64, f64)> = island.iter().rev().copied().collect();
                if let Ok(p) = Polygon::new(ext.clone(), vec![ring(&hole_cw)]) {
                    out.push(p);
                    out.push(poly(&island));
                    continue;
                }
            }
            if let Ok(p) = Polygon::new(ext, vec![]) {
                out.push(p);
            }
        }
    }
    out
}

/// Valid grid coverages, generated until `count` are found.
fn valid_grids(seed: u64, count: usize) -> Vec<Vec<Polygon>> {
    let mut rng = Rng(seed);
    let mut out = Vec::new();
    let mut tries = 0;
    while out.len() < count {
        tries += 1;
        assert!(
            tries < count * 4,
            "the generator should mostly give coverages"
        );
        let (n, m) = (rng.int(1, 6) as usize, rng.int(1, 6) as usize);
        let g = grid(&mut rng, n, m, 100, 6);
        if !g.is_empty() && coverage_validate(&g).is_empty() {
            out.push(g);
        }
    }
    out
}

// ===========================================================================
// Independent property checks
// ===========================================================================

type K = (u64, u64);

fn k(p: Coord2d) -> K {
    ((p.x + 0.0).to_bits(), (p.y + 0.0).to_bits())
}

fn rings(p: &Polygon) -> Vec<Vec<Coord2d>> {
    std::iter::once(p.exterior())
        .chain(p.holes())
        .map(|r| r.coords().as_slice().to_vec())
        .collect()
}

fn seg(a: Coord2d, b: Coord2d) -> (K, K) {
    let (x, y) = (k(a), k(b));
    if x <= y { (x, y) } else { (y, x) }
}

/// Which polygons use each undirected input segment.
fn segment_users(polys: &[Polygon]) -> BTreeMap<(K, K), BTreeSet<usize>> {
    let mut m: BTreeMap<(K, K), BTreeSet<usize>> = BTreeMap::new();
    for (pi, p) in polys.iter().enumerate() {
        for r in rings(p) {
            for t in 0..r.len() {
                m.entry(seg(r[t], r[(t + 1) % r.len()]))
                    .or_default()
                    .insert(pi);
            }
        }
    }
    m
}

/// The nodes, from the segments alone: coordinates where three or more
/// distinct segments meet, or where a ring's incoming and outgoing segments
/// have different users.
fn nodes(polys: &[Polygon]) -> BTreeSet<K> {
    let users = segment_users(polys);
    let mut degree: BTreeMap<K, usize> = BTreeMap::new();
    for &(a, b) in users.keys() {
        *degree.entry(a).or_default() += 1;
        *degree.entry(b).or_default() += 1;
    }
    let mut out: BTreeSet<K> = degree
        .iter()
        .filter(|&(_, &d)| d >= 3)
        .map(|(&p, _)| p)
        .collect();
    for p in polys {
        for r in rings(p) {
            let n = r.len();
            for t in 0..n {
                let (a, v, b) = (r[(t + n - 1) % n], r[t], r[(t + 1) % n]);
                if users[&seg(a, v)] != users[&seg(v, b)] {
                    out.insert(k(v));
                }
            }
        }
    }
    out
}

/// `out` is an ordered cyclic subsequence of `inp`.
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

/// Every coverage check on one result.
fn check_result(input: &[Polygon], out: &[Polygon]) {
    assert_eq!(out.len(), input.len(), "one polygon per input polygon");
    assert!(
        coverage_validate(out).is_empty(),
        "the result is a coverage"
    );
    let node_set = nodes(input);
    for (p, q) in input.iter().zip(out) {
        let (ri, ro) = (rings(p), rings(q));
        assert_eq!(ri.len(), ro.len(), "no ring is dropped");
        for (a, b) in ri.iter().zip(&ro) {
            assert!(b.len() >= 3);
            assert!(
                cyclic_subsequence(b, a),
                "output vertices are input vertices, in order"
            );
            let kept: BTreeSet<K> = b.iter().map(|&x| k(x)).collect();
            for &v in a {
                if node_set.contains(&k(v)) {
                    assert!(kept.contains(&k(v)), "a node was removed");
                }
            }
        }
    }
    // Shared stretches: for each pair of polygons, the vertices of their
    // shared input segments that survive are the same on both sides.
    let users = segment_users(input);
    let mut shared: BTreeMap<(usize, usize), BTreeSet<K>> = BTreeMap::new();
    for (&(a, b), u) in &users {
        if u.len() == 2 {
            let v: Vec<usize> = u.iter().copied().collect();
            let e = shared.entry((v[0], v[1])).or_default();
            e.insert(a);
            e.insert(b);
        }
    }
    let kept_of =
        |i: usize| -> BTreeSet<K> { rings(&out[i]).iter().flatten().map(|&x| k(x)).collect() };
    for (&(i, j), verts) in &shared {
        let (ki, kj) = (kept_of(i), kept_of(j));
        let si: BTreeSet<&K> = verts.iter().filter(|v| ki.contains(v)).collect();
        let sj: BTreeSet<&K> = verts.iter().filter(|v| kj.contains(v)).collect();
        assert_eq!(si, sj, "polygons {i} and {j} keep the same shared vertices");
    }
}

fn vertex_count(polys: &[Polygon]) -> usize {
    polys
        .iter()
        .map(|p| rings(p).iter().map(Vec::len).sum::<usize>())
        .sum()
}

// ===========================================================================
// The result is a coverage, and shared edges match
// ===========================================================================

#[test]
fn results_are_coverages_with_matching_shared_edges() {
    let grids = valid_grids(0xC0F3_5150, 160);
    let mut removed = 0_usize;
    for g in &grids {
        let before = vertex_count(g);
        for t in [0.0, 1.0, 10.0, 60.0, 300.0, 2_000.0, 1e9] {
            let out = coverage_simplify(g, t).unwrap();
            check_result(g, &out);
            removed += before - vertex_count(&out);
        }
    }
    assert!(
        removed > 20_000,
        "the thresholds actually simplify: {removed}"
    );
}

#[test]
fn threshold_zero_removes_only_collinear_vertices() {
    // Integer coordinates: the f64 area is exact, so a removal at 0 is an
    // exactly collinear vertex.
    for g in valid_grids(0x2E40_0001, 80) {
        let out = coverage_simplify(&g, 0.0).unwrap();
        check_result(&g, &out);
        for (p, q) in g.iter().zip(&out) {
            for (a, b) in rings(p).iter().zip(&rings(q)) {
                let kept: BTreeSet<K> = b.iter().map(|&x| k(x)).collect();
                let n = a.len();
                for t in 0..n {
                    if kept.contains(&k(a[t])) {
                        continue;
                    }
                    // The nearest kept neighbors on either side.
                    let mut lo = (t + n - 1) % n;
                    while !kept.contains(&k(a[lo])) {
                        lo = (lo + n - 1) % n;
                    }
                    let mut hi = (t + 1) % n;
                    while !kept.contains(&k(a[hi])) {
                        hi = (hi + 1) % n;
                    }
                    let (u, v, w) = (a[lo], a[t], a[hi]);
                    let cross = (v.x - u.x) * (w.y - u.y) - (v.y - u.y) * (w.x - u.x);
                    assert_eq!(cross, 0.0, "a non-collinear vertex was removed at 0");
                }
            }
        }
    }
    // Collinear midpoints on a shared edge go; the corners stay.
    let a = poly(&[
        (0.0, 0.0),
        (5.0, 0.0),
        (10.0, 0.0),
        (10.0, 5.0),
        (10.0, 10.0),
        (0.0, 10.0),
    ]);
    let b = poly(&[
        (10.0, 0.0),
        (20.0, 0.0),
        (20.0, 10.0),
        (10.0, 10.0),
        (10.0, 5.0),
    ]);
    let out = coverage_simplify(&[a, b], 0.0).unwrap();
    assert_eq!(out[0].exterior().num_vertices(), 4);
    assert_eq!(out[1].exterior().num_vertices(), 4);
}

#[test]
fn a_shared_edge_is_simplified_once_for_both_polygons() {
    let left = poly(&[(0.0, 0.0), (1.0, 0.0), (1.01, 0.5), (1.0, 1.0), (0.0, 1.0)]);
    let right = poly(&[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.01, 0.5)]);
    let out = coverage_simplify(&[left.clone(), right.clone()], 0.01).unwrap();
    assert_eq!(out[0].exterior().num_vertices(), 4);
    assert_eq!(out[1].exterior().num_vertices(), 4);
    check_result(&[left, right], &out);
}

#[test]
fn the_edge_graph_cuts_at_nodes() {
    // A 2 × 2 grid of unit cells: nodes are the center and the four edge
    // midpoints of the outline; the corners are not nodes. Eight edges:
    // four shared spokes, four outline pieces through the corners.
    let cells = [
        poly(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
        poly(&[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)]),
        poly(&[(0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]),
        poly(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]),
    ];
    let g = EdgeGraph::build(&cells).unwrap();
    assert_eq!(g.edges.len(), 8);
    assert_eq!(g.edges.iter().filter(|e| e.users == 2).count(), 4);
    assert_eq!(g.edges.iter().filter(|e| e.users == 1).count(), 4);
    let expected: BTreeSet<K> = [(1.0, 1.0), (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0)]
        .iter()
        .map(|&(x, y)| k(c(x, y)))
        .collect();
    for e in &g.edges {
        assert!(!e.closed());
        assert!(expected.contains(&k(e.verts[0])));
        assert!(expected.contains(&k(e.verts[e.verts.len() - 1])));
    }
    assert_eq!(nodes(&cells), expected, "the independent node set agrees");
    // A lone square has no node: one closed edge.
    let g = EdgeGraph::build(&cells[..1]).unwrap();
    assert_eq!(g.edges.len(), 1);
    assert!(g.edges[0].closed());
}

// ===========================================================================
// Inner mode
// ===========================================================================

#[test]
fn inner_mode_leaves_the_outline_unchanged() {
    for g in valid_grids(0x1AAE_0002, 100) {
        let users = segment_users(&g);
        for t in [10.0, 300.0, 1e9] {
            let out = coverage_simplify_inner(&g, t).unwrap();
            check_result(&g, &out);
            for (p, q) in g.iter().zip(&out) {
                for (a, b) in rings(p).iter().zip(&rings(q)) {
                    let kept: BTreeSet<K> = b.iter().map(|&x| k(x)).collect();
                    let n = a.len();
                    for t in 0..n {
                        let s = seg(a[t], a[(t + 1) % n]);
                        if users[&s].len() == 1 {
                            assert!(
                                kept.contains(&s.0) && kept.contains(&s.1),
                                "an outline vertex moved"
                            );
                        }
                    }
                }
            }
        }
    }
}

// ===========================================================================
// The veto, and refusals
// ===========================================================================

#[test]
fn a_removal_that_would_cover_a_neighbor_is_vetoed() {
    // A's notch vertex (20, 25) has effective area 75, the only one below
    // the thresholds used (the next smallest is 112.5), and B sits inside
    // the notch its removal would close: check (b) must reject it, since B
    // touches neither the chord nor the run.
    let a = poly(&[
        (0.0, 0.0),
        (40.0, 0.0),
        (40.0, 40.0),
        (25.0, 40.0),
        (20.0, 25.0),
        (15.0, 40.0),
        (0.0, 40.0),
    ]);
    let b = poly(&[(19.0, 35.0), (21.0, 35.0), (21.0, 37.0), (19.0, 37.0)]);
    let input = [a, b];
    for t in [80.0, 100.0, 112.0] {
        let out = coverage_simplify(&input, t).unwrap();
        check_result(&input, &out);
        assert!(
            out[0]
                .exterior()
                .coords()
                .iter()
                .any(|&p| p == c(20.0, 25.0)),
            "the notch vertex stays at threshold {t}"
        );
    }
    // Without B the same removal is allowed, and the top vertices it
    // leaves collinear go too: the rectangle.
    let out = coverage_simplify(&input[..1], 100.0).unwrap();
    assert_eq!(out[0].exterior().num_vertices(), 4);
}

#[test]
fn invalid_input_and_thresholds_are_refused() {
    let a = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
    let b = poly(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]);
    assert!(matches!(
        coverage_simplify(&[a.clone(), b], 1.0),
        Err(GeometryError::InvalidCoverage {
            first: CoverageIssue::Overlap { .. }
        })
    ));
    for t in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            coverage_simplify(std::slice::from_ref(&a), t),
            Err(GeometryError::InvalidTolerance { .. })
        ));
    }
    assert!(coverage_simplify(&[], 1.0).unwrap().is_empty());
}

#[test]
fn the_vouch_refuses_a_broken_rebuild() {
    // Injected: two overlapping squares, as a broken rebuild would hand to
    // the closing check. Refused, never returned.
    let a = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
    let b = poly(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]);
    assert!(matches!(
        vouch(vec![a.clone(), b]),
        Err(GeometryError::Unvouched {
            reason: UnvouchedReason::NotACoverage { .. },
            ..
        })
    ));
    assert_eq!(vouch(vec![a.clone()]).unwrap(), vec![a]);
}

#[test]
fn the_result_is_deterministic_and_order_independent_per_polygon() {
    // Reordering the input polygons gives each polygon the same result.
    for g in valid_grids(0x0DE7_0003, 40) {
        let out = coverage_simplify(&g, 300.0).unwrap();
        let rev: Vec<Polygon> = g.iter().rev().cloned().collect();
        let out_rev = coverage_simplify(&rev, 300.0).unwrap();
        let back: Vec<Polygon> = out_rev.into_iter().rev().collect();
        assert_eq!(out, back);
    }
}
