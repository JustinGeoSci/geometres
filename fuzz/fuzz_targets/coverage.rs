// fuzz/fuzz_targets/coverage.rs

//! Fuzz target: coverage validation and coverage union.
//!
//! The input builds a small coverage: a 3 × 3 grid of cells of side 8 whose
//! interior vertices are jittered by the input (shared by their cells, so
//! the grid is valid). The grid is then perturbed by moving one vertex of
//! one cell, which opens gaps, overlaps, T-junctions and vertices on edges.
//! Some cells are removed by a bitmask, and the gap width lies in `[0, 6]`.
//!
//! The contract checked, for every constructed input:
//!
//! - `coverage_validate_with_gap` equals a brute force over every pair of
//!   polygons and every pair of edges, written here with the public exact
//!   API only: `first_self_contact`, `relate` (overlaps, and collinear
//!   overlaps as a line-dimension interior intersection of two segments),
//!   `exact::on_segment`, and for gaps `intersects` and `dwithin` on the
//!   two edges as lines;
//! - whenever the coverage is valid, `coverage_union` succeeds, its result
//!   is valid, every output vertex is an input vertex, and its area equals
//!   the sum of the input areas.

#![no_main]

use std::collections::BTreeSet;

use geometres::exact::{first_self_contact, on_segment};
use geometres::ops::{
    CoverageIssue, area, coverage_union, coverage_validate_with_gap, dwithin, is_valid,
};
use geometres::predicates::{intersects, relate};
use geometres::topology::label::Dimension;
use geometres::types::{Coord2d, CoordSeq, LineString, LinearRing, Polygon};
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

    /// An offset in `[-k, k]`.
    fn offset(&mut self, k: i64) -> i64 {
        i64::from(self.u8()) % (2 * k + 1) - k
    }
}

fn cell(v: [(i64, i64); 4]) -> Option<Polygon> {
    let ring = LinearRing::new(CoordSeq::new(
        v.iter().map(|&(x, y)| Coord2d::new(x as f64, y as f64)).collect(),
    ))
    .ok()?;
    Polygon::new_simple(ring).ok()
}

fn build(b: &mut Bytes) -> Option<(Vec<Polygon>, f64)> {
    const N: usize = 3;
    const S: i64 = 8;
    let mut vx = [[(0_i64, 0_i64); N + 1]; N + 1];
    for (i, row) in vx.iter_mut().enumerate() {
        for (k, v) in row.iter_mut().enumerate() {
            let interior = i > 0 && i < N && k > 0 && k < N;
            let (dx, dy) = if interior { (b.offset(3), b.offset(3)) } else { (0, 0) };
            *v = (i as i64 * S + dx, k as i64 * S + dy);
        }
    }
    let mut cells: Vec<[(i64, i64); 4]> = Vec::new();
    for i in 0..N {
        for k in 0..N {
            cells.push([vx[i][k], vx[i + 1][k], vx[i + 1][k + 1], vx[i][k + 1]]);
        }
    }
    // Perturb one vertex of one cell.
    let which = usize::from(b.u8()) % cells.len();
    let corner = usize::from(b.u8()) % 4;
    let (dx, dy) = (b.offset(5), b.offset(5));
    cells[which][corner].0 += dx;
    cells[which][corner].1 += dy;
    // Remove cells by bitmask (bits set = removed), keeping at least one.
    let mask = u16::from(b.u8()) | (u16::from(b.u8()) << 8);
    let width = f64::from(b.u8() % 25) * 0.25;
    let polys: Vec<Polygon> = cells
        .iter()
        .enumerate()
        .filter(|(i, _)| mask & (1 << i) == 0)
        .filter_map(|(_, &v)| cell(v))
        .collect();
    (!polys.is_empty()).then_some((polys, width))
}

/// One ring edge, with where it came from.
#[derive(Clone, Copy)]
struct E {
    polygon: usize,
    ring: usize,
    edge: usize,
    a: Coord2d,
    b: Coord2d,
}

fn seg(e: E) -> LineString {
    LineString::new(CoordSeq::new(vec![e.a, e.b])).expect("ring edges are distinct points")
}

fn brute(polys: &[Polygon], width: f64) -> Vec<CoverageIssue> {
    let mut found = BTreeSet::new();
    let mut edges = Vec::new();
    for (polygon, p) in polys.iter().enumerate() {
        for (ring, r) in core::iter::once(p.exterior()).chain(p.holes()).enumerate() {
            if first_self_contact(r).is_some() {
                found.insert(CoverageIssue::NonSimpleRing { polygon, ring });
            }
            for (edge, (a, b)) in r.edges().enumerate() {
                edges.push(E { polygon, ring, edge, a, b });
            }
        }
    }
    for a in 0..polys.len() {
        for b in a + 1..polys.len() {
            let ii = relate(&polys[a], &polys[b]).expect("polygons relate").matrix().cells()[0];
            if ii != Dimension::Empty {
                found.insert(CoverageIssue::Overlap { a, b });
            }
        }
    }
    let un = |e: E, other: usize| CoverageIssue::UnmatchedEdge {
        polygon: e.polygon,
        ring: e.ring,
        edge: e.edge,
        other,
    };
    let inside = |p: Coord2d, e: E| p != e.a && p != e.b && on_segment(p, e.a, e.b).unwrap_or(false);
    for (i, &s) in edges.iter().enumerate() {
        for &t in &edges[i + 1..] {
            let same = (s.a == t.a && s.b == t.b) || (s.a == t.b && s.b == t.a);
            if s.polygon == t.polygon || same {
                continue;
            }
            let (ls, lt) = (seg(s), seg(t));
            let collinear = relate(&ls, &lt).expect("lines relate").matrix().cells()[0]
                == Dimension::Line;
            if collinear {
                found.insert(un(s, t.polygon));
                found.insert(un(t, s.polygon));
            }
            if inside(t.a, s) || inside(t.b, s) {
                found.insert(un(s, t.polygon));
            }
            if inside(s.a, t) || inside(s.b, t) {
                found.insert(un(t, s.polygon));
            }
            if width > 0.0
                && !intersects(&ls, &lt).expect("lines intersect")
                && dwithin(&ls, &lt, width).expect("valid width")
            {
                found.insert(CoverageIssue::NarrowGap {
                    polygon: s.polygon,
                    ring: s.ring,
                    edge: s.edge,
                    other: t.polygon,
                    other_ring: t.ring,
                    other_edge: t.edge,
                });
            }
        }
    }
    found.into_iter().collect()
}

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let Some((polys, width)) = build(&mut b) else { return };

    let got = coverage_validate_with_gap(&polys, width).expect("width is valid");
    assert_eq!(got, brute(&polys, width), "validation equals brute force");

    let valid = !got.iter().any(|i| !matches!(i, CoverageIssue::NarrowGap { .. }));
    if !valid {
        assert!(coverage_union(&polys).is_err(), "an invalid coverage is refused");
        return;
    }
    let u = coverage_union(&polys).expect("a valid coverage unions");
    assert!(is_valid(&u), "the union is valid");
    let inputs: BTreeSet<(u64, u64)> = polys
        .iter()
        .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes()))
        .flat_map(|r| r.coords().as_slice().to_vec())
        .map(|c| ((c.x + 0.0).to_bits(), (c.y + 0.0).to_bits()))
        .collect();
    for p in u.polygons() {
        for r in core::iter::once(p.exterior()).chain(p.holes()) {
            for c in r.coords().as_slice() {
                assert!(inputs.contains(&((c.x + 0.0).to_bits(), (c.y + 0.0).to_bits())));
            }
        }
    }
    let a_in: f64 = polys.iter().map(area).sum();
    let a_out: f64 = u.polygons().iter().map(area).sum();
    // Integer coordinates below 2^5: every area is an exact half-integer.
    assert_eq!(a_in, a_out, "areas match");
});
