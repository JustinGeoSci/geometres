// src/index/tests.rs

//! Tests that the index API is exact and that `EdgeRef` numbering matches
//! the order of the geometry types. The tests of the seam itself are in
//! `seam.rs`.
//!
//! The candidates of both indexes are compared with brute force over the
//! same envelopes and edges, on generated layouts that have what a real
//! polygon layer has: shared edges, touching corners, holes, coordinates
//! far from the origin, and repeated entries.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss
)]

use super::*;
use crate::types::{CoordSeq, LinearRing};
use std::collections::BTreeSet;

// ===========================================================================
// Builders and readers
// ===========================================================================

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn square(x: f64, y: f64, s: f64) -> Polygon {
    Polygon::new_simple(ring(&[(x, y), (x + s, y), (x + s, y + s), (x, y + s)])).unwrap()
}

fn env(a: f64, b: f64, x: f64, y: f64) -> Envelope {
    Envelope::new(c(a, b), c(x, y)).unwrap()
}

// ===========================================================================
// The brute-force oracles
// ===========================================================================

fn brute_parts(parts: &[Polygon], window: &Envelope) -> BTreeSet<usize> {
    parts
        .iter()
        .enumerate()
        .filter(|(_, p)| p.exterior().envelope().unwrap().intersects(window))
        .map(|(i, _)| i)
        .collect()
}

fn edge_envelope(a: Coord2d, b: Coord2d) -> Envelope {
    Envelope::from_points([a, b]).unwrap().unwrap()
}

fn brute_edges(parts: &[Polygon], window: &Envelope) -> BTreeSet<EdgeRef> {
    let mut out = BTreeSet::new();
    for (part, p) in parts.iter().enumerate() {
        let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
        for (ring, r) in rings.enumerate() {
            for (edge, (a, b)) in r.edges().enumerate() {
                if edge_envelope(a, b).intersects(window) {
                    out.insert(EdgeRef { part, ring, edge });
                }
            }
        }
    }
    out
}

fn part_candidates(index: &PartIndex, window: &Envelope) -> Vec<usize> {
    let mut v = Vec::new();
    index.candidates(window, |i| v.push(i));
    v
}

fn edge_candidates(index: &EdgeIndex, window: &Envelope) -> Vec<EdgeRef> {
    let mut v = Vec::new();
    index.candidates(window, |e| v.push(e));
    v
}

/// Both indexes against brute force for every window; each candidate is
/// reported once.
fn check_windows(parts: &[Polygon], windows: &[Envelope]) {
    let pi = PartIndex::new(parts).unwrap();
    let ei = EdgeIndex::new(parts).unwrap();
    assert_eq!(pi.len(), parts.len());
    let edges: usize = parts
        .iter()
        .map(|p| {
            p.exterior().num_vertices() + p.holes().iter().map(|h| h.num_vertices()).sum::<usize>()
        })
        .sum();
    assert_eq!(ei.len(), edges);
    for w in windows {
        let got = part_candidates(&pi, w);
        let set: BTreeSet<usize> = got.iter().copied().collect();
        assert_eq!(set.len(), got.len(), "a part was reported twice for {w:?}");
        assert_eq!(set, brute_parts(parts, w), "parts for {w:?}");

        let got = edge_candidates(&ei, w);
        let set: BTreeSet<EdgeRef> = got.iter().copied().collect();
        assert_eq!(set.len(), got.len(), "an edge was reported twice for {w:?}");
        assert_eq!(set, brute_edges(parts, w), "edges for {w:?}");
    }
}

// ===========================================================================
// Generated layouts
// ===========================================================================

/// A small deterministic generator (xorshift64*), so the layout is the same
/// on every run.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// A value in `0..n`.
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// A layout with the cases a real polygon layer has: lots that share whole
/// edges and single corners, lots with one or two holes, irregular
/// vertices on dyadic coordinates far from the origin, and repeated
/// entries. Lot `(i, j)` occupies the cell `[i, i + 1] × [j, j + 1]` scaled
/// by 64 and offset by `2²⁰`; its corners are the cell's, so neighbors share
/// edges, and every fifth lot is shrunk so it meets its neighbors at most at
/// a corner or not at all.
fn generated_layout() -> Vec<Polygon> {
    const N: u32 = 24;
    const S: f64 = 64.0;
    const O: f64 = 1_048_576.0;
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut out = Vec::new();
    for j in 0..N {
        for i in 0..N {
            let (x0, y0) = (O + f64::from(i) * S, O + f64::from(j) * S);
            let shrink = if (i + j) % 5 == 0 { 0.25 * S } else { 0.0 };
            let (a, b) = (x0 + shrink, y0 + shrink);
            let (cx, cy) = (x0 + S - shrink, y0 + S - shrink);
            // Extra vertices along the bottom and top edges, at dyadic
            // offsets, keep the ring simple and counterclockwise.
            let k = rng.below(4) as usize;
            let mut pts = vec![c(a, b)];
            for t in 1..=k {
                pts.push(c(a + (cx - a) * (t as f64) / (k as f64 + 1.0), b));
            }
            pts.push(c(cx, b));
            pts.push(c(cx, cy));
            pts.push(c(
                a + 0.5 * (cx - a),
                cy + if i % 2 == 0 { 0.0 } else { -0.125 },
            ));
            pts.push(c(a, cy));
            let ext = LinearRing::new(CoordSeq::new(pts)).unwrap();
            let holes = match rng.below(3) {
                0 => vec![],
                1 => vec![ring(&[
                    (x0 + 24.0, y0 + 24.0),
                    (x0 + 24.0, y0 + 40.0),
                    (x0 + 40.0, y0 + 40.0),
                    (x0 + 40.0, y0 + 24.0),
                ])],
                _ => vec![
                    ring(&[
                        (x0 + 20.0, y0 + 20.0),
                        (x0 + 20.0, y0 + 28.0),
                        (x0 + 28.0, y0 + 28.0),
                        (x0 + 28.0, y0 + 20.0),
                    ]),
                    ring(&[
                        (x0 + 36.0, y0 + 36.0),
                        (x0 + 36.0, y0 + 44.0),
                        (x0 + 44.0, y0 + 44.0),
                    ]),
                ],
            };
            out.push(Polygon::new(ext, holes).unwrap());
        }
    }
    // Repeated entries, as a layer with duplicated records has.
    let repeats: Vec<Polygon> = out.iter().step_by(37).cloned().collect();
    out.extend(repeats);
    out
}

/// Every lot's envelope as a window over the whole layout, then windows
/// that touch lots on a single coordinate.
#[test]
fn generated_layout_candidates_equal_brute_force() {
    let parts = generated_layout();
    assert_eq!(parts.len(), 576 + 16);
    let windows: Vec<Envelope> = parts
        .iter()
        .map(|p| p.exterior().envelope().unwrap())
        .collect();
    check_windows(&parts, &windows);
    let touching: Vec<Envelope> = windows
        .iter()
        .take(200)
        .map(|w| env(w.max().x, w.max().y, w.max().x + 10.0, w.max().y + 10.0))
        .collect();
    check_windows(&parts, &touching);
}

// ===========================================================================
// Structure and numbering
// ===========================================================================

#[test]
fn empty_indexes_have_no_candidates_and_no_bounds() {
    let pi = PartIndex::new(&[]).unwrap();
    let ei = EdgeIndex::new(&[]).unwrap();
    assert!(pi.is_empty() && ei.is_empty());
    assert_eq!((pi.len(), ei.len()), (0, 0));
    assert_eq!((pi.bounds(), ei.bounds()), (None, None));
    let w = env(-1e300, -1e300, 1e300, 1e300);
    assert!(part_candidates(&pi, &w).is_empty());
    assert!(edge_candidates(&ei, &w).is_empty());
}

#[test]
fn edge_refs_follow_exterior_then_holes_in_edges_order() {
    let outer = ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
    let hole_a = ring(&[(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)]);
    let hole_b = ring(&[(5.0, 5.0), (5.0, 6.0), (6.0, 6.0)]);
    let p = Polygon::new(outer, vec![hole_a, hole_b]).unwrap();
    let parts = [square(100.0, 100.0, 1.0), p];
    let ei = EdgeIndex::new(&parts).unwrap();
    assert_eq!(ei.len(), 4 + 4 + 4 + 3);

    // The closing edge of hole B, `(6,6) → (5,5)`, is edge 2 of ring 2.
    let w = env(5.5, 5.5, 5.5, 5.5);
    let got = edge_candidates(&ei, &w);
    assert_eq!(
        got,
        [EdgeRef {
            part: 1,
            ring: 2,
            edge: 2
        }]
    );

    // A window on the exterior's bottom edge, `(0,0) → (10,0)`, edge 0 of
    // ring 0, and nothing else.
    let w = env(4.0, 0.0, 6.0, 0.0);
    assert_eq!(
        edge_candidates(&ei, &w),
        [EdgeRef {
            part: 1,
            ring: 0,
            edge: 0
        }]
    );

    // Hole A's first edge is `(1,1) → (1,3)`: edge 0 of ring 1.
    let w = env(1.0, 2.0, 1.0, 2.0);
    assert_eq!(
        edge_candidates(&ei, &w),
        [EdgeRef {
            part: 1,
            ring: 1,
            edge: 0
        }]
    );

    // The numbering agrees with the types: rebuild each EdgeRef's endpoints
    // from `exterior`/`holes` and `edges`, and check its box meets the
    // window that found it.
    for w in [env(0.0, 0.0, 10.0, 10.0), env(2.0, 2.0, 2.0, 2.0)] {
        for e in edge_candidates(&ei, &w) {
            let part = &parts[e.part];
            let r = if e.ring == 0 {
                part.exterior()
            } else {
                &part.holes()[e.ring - 1]
            };
            let (a, b) = r.edges().nth(e.edge).unwrap();
            assert!(edge_envelope(a, b).intersects(&w), "{e:?}");
        }
    }
    assert_eq!(ei.bounds(), Some(env(0.0, 0.0, 101.0, 101.0)));
}

#[test]
fn part_index_reports_touching_parts_and_bounds() {
    let parts = [
        square(0.0, 0.0, 1.0),
        square(5.0, 5.0, 1.0),
        square(1.0, 0.0, 1.0),
    ];
    let pi = PartIndex::new(&parts).unwrap();
    assert_eq!(pi.len(), 3);
    let mut got = part_candidates(&pi, &env(1.0, 0.5, 1.5, 0.5));
    got.sort_unstable();
    assert_eq!(got, [0, 2]);
    assert_eq!(part_candidates(&pi, &env(6.0, 6.0, 7.0, 7.0)), [1]);
    assert!(part_candidates(&pi, &env(6.1, 6.0, 7.0, 7.0)).is_empty());
    assert_eq!(pi.bounds(), Some(env(0.0, 0.0, 6.0, 6.0)));
}

#[test]
fn a_part_with_holes_is_boxed_by_its_exterior() {
    let outer = ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
    let hole = ring(&[(4.0, 4.0), (4.0, 6.0), (6.0, 6.0), (6.0, 4.0)]);
    let parts = [Polygon::new(outer, vec![hole]).unwrap()];
    let pi = PartIndex::new(&parts).unwrap();
    // A window inside the hole still meets the part's box: the index
    // filters, the caller decides.
    assert_eq!(part_candidates(&pi, &env(5.0, 5.0, 5.0, 5.0)), [0]);
}

#[test]
fn too_many_entries_is_an_error_not_a_truncation() {
    assert_eq!(check_count(0), Ok(0));
    assert_eq!(check_count(u32::MAX as usize), Ok(u32::MAX));
    let over = u32::MAX as usize + 1;
    assert_eq!(
        check_count(over),
        Err(GeometryError::TooManyEntries { count: over })
    );
    assert_eq!(
        format!("{}", GeometryError::TooManyEntries { count: over }),
        "spatial index cannot hold 4294967296 entries (the limit is 4294967295)"
    );
}

#[test]
fn generated_squares_match_brute_force() {
    // A 30 × 30 grid of touching squares, queried with points, lines and
    // boxes on and between the lattice lines.
    let parts: Vec<Polygon> = (0..900)
        .map(|i| square(f64::from(i % 30), f64::from(i / 30), 1.0))
        .collect();
    let mut windows = Vec::new();
    for i in 0..30 {
        let x = f64::from(i);
        windows.push(env(x, x, x, x));
        windows.push(env(x + 0.5, 0.0, x + 0.5, 30.0));
        windows.push(env(x, 3.0, x + 2.0, 3.0));
        windows.push(env(x - 0.25, x - 0.25, x + 0.75, x + 0.75));
    }
    check_windows(&parts, &windows);
}
