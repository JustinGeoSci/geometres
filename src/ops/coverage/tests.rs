// src/ops/coverage/tests.rs

//! Tests of the coverage operations: the angular order against an
//! independent exact oracle, coverage validation against brute force, the
//! gap width at its exact boundary, the properties of coverage union, and
//! polygonize on noded figures.
#![allow(clippy::unwrap_used)]
// test-lint-posture: see the workspace lint comment in Cargo.toml
#![allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_lossless,
    clippy::too_many_lines
)]

use core::cmp::Ordering;
use std::collections::BTreeSet;

use super::graph::{coord_exact_cmp, edge_key};
use super::validate::{Edge, collinear_overlap, edges_of, gap_window, touch};
use super::*;
use crate::accum::Accum;
use crate::ops::area;
use crate::ops::distance::{Element, pair_within};
use crate::ops::is_valid;
use crate::predicates::relate;
use crate::topology::label::Dimension;
use crate::types::{Coord2d, CoordSeq, Geometry, LineString, LinearRing, MultiLineString, Polygon};

// ===========================================================================
// Fixtures
// ===========================================================================

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn poly(ext: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Polygon {
    Polygon::new(ring(ext), holes.iter().map(|h| ring(h)).collect()).unwrap()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    poly(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)], &[])
}

fn line(pts: &[(f64, f64)]) -> LineString {
    LineString::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
}

/// An `n × m` grid of cells of side `s`, every interior grid vertex jittered
/// by an integer in `[-j, j]`. Each such vertex is shared by all four of its
/// cells, so the result is a valid coverage while `j < s/2`. Cells listed row by row.
fn jittered_grid(rng: &mut Rng, n: usize, m: usize, s: i64, j: i64) -> Vec<Polygon> {
    let mut vx = vec![vec![(0.0, 0.0); m + 1]; n + 1];
    for (i, row) in vx.iter_mut().enumerate() {
        for (k, v) in row.iter_mut().enumerate() {
            let interior = i > 0 && i < n && k > 0 && k < m;
            let (dx, dy) = if interior {
                (rng.int(-j, j), rng.int(-j, j))
            } else {
                (0, 0)
            };
            *v = ((i as i64 * s + dx) as f64, (k as i64 * s + dy) as f64);
        }
    }
    let mut out = Vec::new();
    for i in 0..n {
        for k in 0..m {
            out.push(poly(
                &[vx[i][k], vx[i + 1][k], vx[i + 1][k + 1], vx[i][k + 1]],
                &[],
            ));
        }
    }
    out
}

// ===========================================================================
// The angular order
// ===========================================================================

/// The independent oracle: the "diamond angle" of an integer direction, an
/// exact rational strictly increasing with the true angle on `[0, 2π)`,
/// returned as `(quadrant, numerator, denominator)`.
fn diamond(dx: i64, dy: i64) -> (i64, i128, i128) {
    let (dx, dy) = (i128::from(dx), i128::from(dy));
    if dx > 0 && dy >= 0 {
        (0, dy, dx + dy)
    } else if dx <= 0 && dy > 0 {
        (1, -dx, -dx + dy)
    } else if dx < 0 && dy <= 0 {
        (2, -dy, -dx - dy)
    } else {
        (3, dx, dx - dy)
    }
}

fn diamond_cmp(a: (i64, i64), b: (i64, i64)) -> Ordering {
    let (qa, na, da) = diamond(a.0, a.1);
    let (qb, nb, db) = diamond(b.0, b.1);
    qa.cmp(&qb).then((na * db).cmp(&(nb * da)))
}

#[test]
fn the_angular_order_matches_the_exact_diamond_angle_on_random_fans() {
    let mut rng = Rng(0xA4C1_E5F0_0D12_3457);
    let mut checked = 0_u32;
    for _ in 0..400 {
        let (vx, vy) = (rng.int(-1_000, 1_000), rng.int(-1_000, 1_000));
        // Directions include axis-aligned, collinear, opposite and repeated
        // ones: small integer components.
        let dirs: Vec<(i64, i64)> = (0..12)
            .map(|_| {
                loop {
                    let d = (
                        rng.int(-3, 3) * rng.int(1, 4),
                        rng.int(-3, 3) * rng.int(1, 4),
                    );
                    if d != (0, 0) {
                        break d;
                    }
                }
            })
            .collect();
        let v = c(vx as f64, vy as f64);
        for &a in &dirs {
            for &b in &dirs {
                let pa = c((vx + a.0) as f64, (vy + a.1) as f64);
                let pb = c((vx + b.0) as f64, (vy + b.1) as f64);
                assert_eq!(
                    angular::ccw_cmp(v, pa, pb),
                    diamond_cmp(a, b),
                    "directions {a:?} and {b:?}"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 50_000);
}

#[test]
fn opposite_axis_and_signed_zero_directions_order_exactly() {
    let v = c(0.0, 0.0);
    let east = c(1.0, 0.0);
    let north = c(0.0, 1.0);
    let west = c(-1.0, 0.0);
    let south = c(0.0, -1.0);
    // The cyclic order from +x.
    let order = [
        east,
        c(1.0, 1.0),
        north,
        c(-1.0, 1.0),
        west,
        c(-1.0, -1.0),
        south,
        c(1.0, -1.0),
    ];
    for (i, &a) in order.iter().enumerate() {
        for (k, &b) in order.iter().enumerate() {
            assert_eq!(angular::ccw_cmp(v, a, b), i.cmp(&k), "{a:?} vs {b:?}");
        }
    }
    // Opposite directions are never Equal (the half classification).
    assert_eq!(angular::ccw_cmp(v, east, west), Ordering::Less);
    assert_eq!(angular::ccw_cmp(v, north, south), Ordering::Less);
    // Signed zeros: a direction with a -0.0 component is the same direction.
    let vz = c(-0.0, 0.0);
    assert_eq!(
        angular::ccw_cmp(vz, c(0.0, 1.0), c(-0.0, 2.0)),
        Ordering::Equal
    );
    assert_eq!(
        angular::ccw_cmp(vz, c(1.0, -0.0), c(-1.0, 0.0)),
        Ordering::Less
    );
    assert_eq!(
        angular::ccw_cmp(c(0.0, -0.0), c(-1.0, 0.0), c(1.0, -0.0)),
        Ordering::Greater
    );
}

/// Writes random non-integer fans, sorted by the angular order, to
/// `angular_cases.csv` in the system temporary directory, for an independent
/// check of every adjacent pair with exact rational diamond angles.
#[test]
#[ignore = "writes the angular-order case file angular_cases.csv to the temporary directory"]
fn write_angular_cases_csv() {
    use std::io::Write as _;
    let path = std::env::temp_dir().join("angular_cases.csv");
    let mut out = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
    let mut rng = Rng(0x0A46_0BAC_E5F0_0001);
    let f = |rng: &mut Rng, scale: f64| {
        ((rng.next() >> 11) as f64 / (1_u64 << 53) as f64 - 0.5) * scale
    };
    for _ in 0..5_000 {
        let scale = 2.0_f64.powi(rng.int(-40, 40) as i32);
        let v = c(f(&mut rng, scale * 8.0), f(&mut rng, scale * 8.0));
        let mut pts: Vec<Coord2d> = (0..10)
            .map(|_| {
                let mut p = c(v.x + f(&mut rng, scale), v.y + f(&mut rng, scale));
                // Exact axis and opposite directions, half the time.
                if rng.next().is_multiple_of(4) {
                    p = c(p.x, v.y);
                }
                if rng.next().is_multiple_of(4) {
                    p = c(v.x, p.y);
                }
                p
            })
            .filter(|p| *p != v)
            .collect();
        pts.sort_by(|&a, &b| angular::ccw_cmp(v, a, b));
        let fields: Vec<String> = core::iter::once(v)
            .chain(pts)
            .map(|p| format!("{:016x},{:016x}", p.x.to_bits(), p.y.to_bits()))
            .collect();
        writeln!(out, "{}", fields.join(";")).unwrap();
    }
    out.flush().unwrap();
}

// ===========================================================================
// Coverage validation: the brute-force oracle
// ===========================================================================

/// Every issue, by brute force over all polygon pairs and all edge pairs,
/// with no index. Uses the same exact per-pair rules.
fn brute_issues(polys: &[Polygon], gap: Option<f64>) -> Vec<CoverageIssue> {
    let mut found = BTreeSet::new();
    for (polygon, p) in polys.iter().enumerate() {
        for (ring, r) in core::iter::once(p.exterior()).chain(p.holes()).enumerate() {
            if crate::exact::first_self_contact(r).is_some() {
                found.insert(CoverageIssue::NonSimpleRing { polygon, ring });
            }
        }
    }
    for a in 0..polys.len() {
        for b in a + 1..polys.len() {
            if relate(&polys[a], &polys[b]).unwrap().matrix().cells()[0] != Dimension::Empty {
                found.insert(CoverageIssue::Overlap { a, b });
            }
        }
    }
    let edges = edges_of(polys);
    let inside =
        |p: Coord2d, e: Edge| p != e.a && p != e.b && crate::exact::kernel::on_segment(p, e.a, e.b);
    let un = |e: Edge, other: usize| CoverageIssue::UnmatchedEdge {
        polygon: e.polygon,
        ring: e.ring,
        edge: e.edge,
        other,
    };
    for (i, &s) in edges.iter().enumerate() {
        for &t in &edges[i + 1..] {
            if s.polygon == t.polygon || edge_key(s.a, s.b) == edge_key(t.a, t.b) {
                continue;
            }
            if collinear_overlap(s, t) {
                found.insert(un(s, t.polygon));
                found.insert(un(t, s.polygon));
            }
            if inside(t.a, s) || inside(t.b, s) {
                found.insert(un(s, t.polygon));
            }
            if inside(s.a, t) || inside(s.b, t) {
                found.insert(un(t, s.polygon));
            }
            if let Some(w) = gap.filter(|&w| w > 0.0)
                && !touch(s, t)
                && pair_within(Element::Segment(s.a, s.b), Element::Segment(t.a, t.b), w)
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

fn un(polygon: usize, ring: usize, edge: usize, other: usize) -> CoverageIssue {
    CoverageIssue::UnmatchedEdge {
        polygon,
        ring,
        edge,
        other,
    }
}

// ===========================================================================
// Coverage validation: hand cases
// ===========================================================================

#[test]
fn a_t_junction_is_four_unmatched_edges() {
    let a = rect(0.0, 0.0, 2.0, 1.0); // bottom edge 0: (0,0)→(2,0)
    let b = rect(0.0, -1.0, 1.0, 0.0); // top edge 2: (1,0)→(0,0)
    let cc = rect(1.0, -1.0, 2.0, 0.0); // top edge 2: (2,0)→(1,0)
    let got = coverage_validate(&[a, b, cc]);
    assert_eq!(
        got,
        vec![
            un(0, 0, 0, 1),
            un(0, 0, 0, 2),
            un(1, 0, 2, 0),
            un(2, 0, 2, 0)
        ]
    );
}

#[test]
fn an_overlap_is_reported_once_and_nothing_else() {
    let got = coverage_validate(&[rect(0.0, 0.0, 2.0, 2.0), rect(1.0, 1.0, 3.0, 3.0)]);
    assert_eq!(got, vec![CoverageIssue::Overlap { a: 0, b: 1 }]);
}

#[test]
fn a_sliver_overlap_is_an_overlap_and_unmatched_edges() {
    let a = rect(0.0, 0.0, 2.0, 1.0);
    let b = poly(&[(0.0, 0.875), (2.0, 1.0), (2.0, 2.0), (0.0, 2.0)], &[]);
    let polys = [a, b];
    let got = coverage_validate(&polys);
    assert!(got.contains(&CoverageIssue::Overlap { a: 0, b: 1 }));
    assert!(
        got.contains(&un(0, 0, 3, 1)),
        "A's left edge holds B's vertex"
    );
    assert!(got.contains(&un(1, 0, 3, 0)), "B's left edge overlaps A's");
    assert_eq!(got, brute_issues(&polys, None));
}

#[test]
fn a_valid_grid_and_a_filled_hole_are_valid() {
    let mut rng = Rng(7);
    assert!(coverage_validate(&jittered_grid(&mut rng, 4, 3, 10, 3)).is_empty());
    let holed = poly(
        &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        &[&[(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)]],
    );
    let filler = rect(1.0, 1.0, 3.0, 3.0);
    assert!(coverage_validate(&[holed, filler]).is_empty());
}

#[test]
fn a_vertex_on_another_polygons_edge_is_unmatched() {
    // B's top-left corner touches A's bottom edge at its midpoint only.
    let a = rect(0.0, 0.0, 2.0, 1.0);
    let b = rect(1.0, -1.0, 3.0, 0.0);
    // B shares A's bottom edge partly, collinearly: also an overlap stretch.
    let got = coverage_validate(&[a.clone(), b]);
    assert!(got.contains(&un(0, 0, 0, 1)));
    // A pure point touch: a triangle tip on A's bottom edge.
    let tip = poly(&[(1.0, 0.0), (0.5, -1.0), (1.5, -1.0)], &[]);
    let got = coverage_validate(&[a, tip]);
    assert_eq!(got, vec![un(0, 0, 0, 1)]);
}

#[test]
fn a_non_simple_ring_is_reported() {
    let bowtie =
        Polygon::new_simple(ring(&[(0.0, 0.0), (5.0, 4.0), (0.0, 3.0), (4.0, 0.0)])).unwrap();
    let got = coverage_validate(&[bowtie, rect(10.0, 10.0, 11.0, 11.0)]);
    assert_eq!(
        got,
        vec![CoverageIssue::NonSimpleRing {
            polygon: 0,
            ring: 0
        }]
    );
    assert!(matches!(
        coverage_union(&[Polygon::new_simple(ring(&[
            (0.0, 0.0),
            (5.0, 4.0),
            (0.0, 3.0),
            (4.0, 0.0)
        ]))
        .unwrap()]),
        Err(GeometryError::InvalidCoverage {
            first: CoverageIssue::NonSimpleRing {
                polygon: 0,
                ring: 0
            }
        })
    ));
}

#[test]
fn perturbed_grids_agree_with_brute_force() {
    let mut rng = Rng(0xC0FE_C0FE_0000_0001);
    let mut with_issues = 0_u32;
    for _ in 0..120 {
        let mut polys = jittered_grid(&mut rng, 3, 3, 10, 3);
        // Move one vertex of one cell only: a gap, an overlap, a T, or a
        // vertex on a neighbor's edge.
        let k = (rng.next() % polys.len() as u64) as usize;
        let mut v = polys[k].exterior().coords().as_slice().to_vec();
        let at = (rng.next() % 4) as usize;
        v[at] = c(
            v[at].x + rng.int(-4, 4) as f64,
            v[at].y + rng.int(-4, 4) as f64,
        );
        if let Ok(r) = LinearRing::new(CoordSeq::new(v))
            && let Ok(p) = Polygon::new_simple(r)
        {
            polys[k] = p;
        }
        let got = coverage_validate(&polys);
        assert_eq!(got, brute_issues(&polys, None));
        let gw = coverage_validate_with_gap(&polys, 2.0).unwrap();
        assert_eq!(gw, brute_issues(&polys, Some(2.0)));
        with_issues += u32::from(!got.is_empty());
    }
    assert!(
        with_issues > 40,
        "the perturbations exercise the issues: {with_issues}"
    );
}

// ===========================================================================
// The gap width
// ===========================================================================

fn gaps(issues: &[CoverageIssue]) -> usize {
    issues
        .iter()
        .filter(|i| matches!(i, CoverageIssue::NarrowGap { .. }))
        .count()
}

#[test]
fn a_gap_exactly_the_width_is_reported_and_one_ulp_wider_is_not() {
    let a = rect(0.0, 0.0, 10.0, 10.0);
    let b = rect(13.0, 0.0, 20.0, 10.0);
    let got = coverage_validate_with_gap(&[a.clone(), b], 3.0).unwrap();
    // Exactly 3 apart: A's right edge (1) and B's left edge (3), and the
    // edges at the gap's two mouths (A's bottom (0) and top (2) against
    // B's bottom (0) and top (2), and against each other's sides), whose
    // nearest points are the gap's corners, also exactly 3 apart.
    let ng = |edge: usize, other_edge: usize| CoverageIssue::NarrowGap {
        polygon: 0,
        ring: 0,
        edge,
        other: 1,
        other_ring: 0,
        other_edge,
    };
    assert_eq!(
        got,
        vec![
            ng(0, 0),
            ng(0, 3),
            ng(1, 0),
            ng(1, 2),
            ng(1, 3),
            ng(2, 2),
            ng(2, 3)
        ]
    );
    let wider = rect(13.0_f64.next_up(), 0.0, 20.0, 10.0);
    assert!(
        coverage_validate_with_gap(&[a.clone(), wider], 3.0)
            .unwrap()
            .is_empty()
    );
    // And the width one ulp below the gap does not report it either.
    let b = rect(13.0, 0.0, 20.0, 10.0);
    assert!(
        coverage_validate_with_gap(&[a, b], 3.0_f64.next_down())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn touching_and_matched_edges_are_not_gaps() {
    let mut rng = Rng(11);
    let grid = jittered_grid(&mut rng, 3, 3, 10, 3);
    // Every pair of edges of a valid grid either matches, touches, or is
    // farther apart than the width; with a width below the cell size,
    // nothing is reported.
    assert_eq!(gaps(&coverage_validate_with_gap(&grid, 1.0).unwrap()), 0);
    // A corner touch is not a gap: the four pairs of edges meeting at the
    // corner touch, and every other pair is at least 1 apart.
    let corner = [rect(0.0, 0.0, 1.0, 1.0), rect(1.0, 1.0, 2.0, 2.0)];
    assert_eq!(gaps(&coverage_validate_with_gap(&corner, 0.5).unwrap()), 0);
    // At width 1 the non-touching pairs 1 apart are reported, and still none
    // of the touching pairs (A's right and top edges against B's bottom and
    // left edges).
    let at1 = coverage_validate_with_gap(&corner, 1.0).unwrap();
    assert!(gaps(&at1) > 0);
    for i in &at1 {
        if let CoverageIssue::NarrowGap {
            edge, other_edge, ..
        } = *i
        {
            assert!(
                !([1, 2].contains(&edge) && [0, 3].contains(&other_edge)),
                "a touching pair was reported: {edge}/{other_edge}"
            );
        }
    }
    // width = 0 reports nothing; a bad width is refused.
    assert!(coverage_validate_with_gap(&grid, 0.0).unwrap().is_empty());
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            coverage_validate_with_gap(&grid, bad),
            Err(GeometryError::InvalidTolerance { .. })
        ));
    }
}

#[test]
fn vertices_moved_inward_by_a_known_amount_open_a_gap_of_that_width() {
    // Two cells sharing the edge x = 10; move both of the left cell's
    // vertices on it inward by exactly `k`: a uniform gap of width `k`.
    // (Moving one vertex makes a wedge narrowing to zero where the cells
    // still touch, narrower than `k` near the far corner: the neighbor's
    // corner (10, 10) is 10/√101 ≈ 0.995 from the slanted edge when k = 1.)
    for k in [1.0, 2.0, 4.0] {
        let left = rect(0.0, 0.0, 10.0 - k, 10.0);
        let right = rect(10.0, 0.0, 20.0, 10.0);
        let polys = [left, right];
        let at = coverage_validate_with_gap(&polys, k).unwrap();
        let below = coverage_validate_with_gap(&polys, k.next_down()).unwrap();
        assert!(gaps(&at) > 0, "reported at width {k}");
        assert_eq!(gaps(&below), 0, "not reported just below {k}");
        assert_eq!(at, brute_issues(&polys, Some(k)));
    }
    // The wedge: reported below 1, exactly as brute force says.
    let wedge = [
        poly(&[(0.0, 0.0), (10.0, 0.0), (9.0, 10.0), (0.0, 10.0)], &[]),
        rect(10.0, 0.0, 20.0, 10.0),
    ];
    let w = 0.996;
    let got = coverage_validate_with_gap(&wedge, w).unwrap();
    assert!(
        gaps(&got) > 0,
        "the wedge is 10/√101 ≈ 0.99504 wide at (10, 10)"
    );
    assert_eq!(got, brute_issues(&wedge, Some(w)));
    // And at 0.995, just below the irrational width, nothing is reported.
    assert_eq!(gaps(&coverage_validate_with_gap(&wedge, 0.995).unwrap()), 0);
}

/// V5's catch: the gap candidate window contains the exact grown box. Every
/// bound is checked exactly by TwoSum: `fl(v ± w) = s` and the exact error
/// `e`, so the real `v − w` is `s + e`, and `lo ≤ v − w ⟺ (s − lo) + e ≥ 0`
/// with `s − lo` exact (adjacent floats).
#[test]
fn gap_windows_contain_the_exact_grown_box() {
    fn two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        let bb = s - a;
        (s, (a - (s - bb)) + (b - bb))
    }
    let mut rng = Rng(0x6A9_0001_D00D_F00D);
    let mut inward = 0_u32;
    for _ in 0..20_000 {
        let frac = |rng: &mut Rng| 1.0 + (rng.next() >> 12) as f64 / (1_u64 << 52) as f64;
        let e1 = rng.int(-50, 50) as i32;
        let a = c(
            frac(&mut rng) * 2.0_f64.powi(e1),
            -frac(&mut rng) * 2.0_f64.powi(e1),
        );
        let b = c(a.x + frac(&mut rng), a.y + frac(&mut rng));
        let w = frac(&mut rng) * 2.0_f64.powi(rng.int(-50, 50) as i32);
        let win = gap_window(a, b, w);
        let (s, e) = two_sum(a.x.min(b.x), -w);
        assert!(
            s - win.min().x >= -e,
            "x lower bound above the exact grown box"
        );
        inward += u32::from(e < 0.0);
        let (s, e) = two_sum(a.y.min(b.y), -w);
        assert!(
            s - win.min().y >= -e,
            "y lower bound above the exact grown box"
        );
        let (s, e) = two_sum(a.x.max(b.x), w);
        assert!(
            win.max().x - s >= e,
            "x upper bound below the exact grown box"
        );
        let (s, e) = two_sum(a.y.max(b.y), w);
        assert!(
            win.max().y - s >= e,
            "y upper bound below the exact grown box"
        );
    }
    assert!(inward > 1_000, "the inward-rounding case occurs: {inward}");
}

// ===========================================================================
// Coverage union
// ===========================================================================

fn input_vertices(polys: &[Polygon]) -> BTreeSet<(u64, u64)> {
    edges_of(polys)
        .iter()
        .map(|e| super::graph::key(e.a))
        .collect()
}

fn check_union(polys: &[Polygon]) -> crate::types::MultiPolygon {
    let u = coverage_union(polys).unwrap();
    assert!(is_valid(&u), "the union is valid: {u:?}");
    let inputs = input_vertices(polys);
    for p in u.polygons() {
        for r in core::iter::once(p.exterior()).chain(p.holes()) {
            for v in r.coords().as_slice() {
                assert!(
                    inputs.contains(&super::graph::key(*v)),
                    "{v:?} is an input vertex"
                );
            }
        }
    }
    let mut sum_in = Accum::new();
    for p in polys {
        sum_in.add(area(p));
    }
    let mut sum_out = Accum::new();
    for p in u.polygons() {
        sum_out.add(area(p));
    }
    let (i, o) = (sum_in.value(), sum_out.value());
    assert!(
        (i - o).abs() <= 1e-12 * i.abs().max(1.0),
        "areas {i} vs {o}"
    );
    u
}

#[test]
fn a_grid_unions_to_one_square_keeping_collinear_vertices() {
    let polys: Vec<Polygon> = (0..3)
        .flat_map(|i| (0..3).map(move |k| rect(i as f64, k as f64, i as f64 + 1.0, k as f64 + 1.0)))
        .collect();
    let u = check_union(&polys);
    assert_eq!(u.polygons().len(), 1);
    assert_eq!(u.polygons()[0].holes().len(), 0);
    assert_eq!(
        u.polygons()[0].exterior().num_vertices(),
        12,
        "collinear vertices kept"
    );
}

#[test]
fn corner_touching_parts_stay_separate() {
    let u = check_union(&[rect(0.0, 0.0, 1.0, 1.0), rect(1.0, 1.0, 2.0, 2.0)]);
    assert_eq!(u.polygons().len(), 2);
}

#[test]
fn a_hole_touching_the_exterior_at_a_vertex_is_split_into_two_rings() {
    // A 4 × 4 block of unit cells without cell (0,0), a notch whose inner
    // corner is (1,1), and without cell (1,1), a hole touching the notch
    // corner at (1,1).
    let polys: Vec<Polygon> = (0..4)
        .flat_map(|i| (0..4).map(move |k| (i, k)))
        .filter(|&(i, k)| (i, k) != (0, 0) && (i, k) != (1, 1))
        .map(|(i, k)| rect(i as f64, k as f64, i as f64 + 1.0, k as f64 + 1.0))
        .collect();
    let u = check_union(&polys);
    assert_eq!(u.polygons().len(), 1);
    assert_eq!(u.polygons()[0].holes().len(), 1);
}

#[test]
fn an_island_in_a_hole_is_its_own_part() {
    // A 5 × 5 block without the ring of eight cells around the center: an
    // outer part with a 3 × 3 hole, and the center cell as an island.
    let polys: Vec<Polygon> = (0..5)
        .flat_map(|i| (0..5).map(move |k| (i, k)))
        .filter(|&(i, k)| {
            let ring8 = (1..=3).contains(&i) && (1..=3).contains(&k) && (i, k) != (2, 2);
            !ring8
        })
        .map(|(i, k)| rect(i as f64, k as f64, i as f64 + 1.0, k as f64 + 1.0))
        .collect();
    let u = check_union(&polys);
    assert_eq!(u.polygons().len(), 2);
    let holes: Vec<usize> = u.polygons().iter().map(|p| p.holes().len()).collect();
    assert_eq!(holes.iter().sum::<usize>(), 1);
}

#[test]
fn a_filled_hole_dissolves_away() {
    let holed = poly(
        &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
        &[&[(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)]],
    );
    let u = check_union(&[holed, rect(1.0, 1.0, 3.0, 3.0)]);
    assert_eq!(u.polygons().len(), 1);
    assert!(u.polygons()[0].holes().is_empty());
}

#[test]
fn random_jittered_grids_with_cells_removed_union_validly() {
    let mut rng = Rng(0x00C0_7E2A_6E00_0042);
    for _ in 0..80 {
        let grid = jittered_grid(&mut rng, 5, 4, 10, 4);
        let kept: Vec<Polygon> = grid
            .into_iter()
            .filter(|_| !rng.next().is_multiple_of(4))
            .collect();
        if kept.is_empty() {
            continue;
        }
        check_union(&kept);
    }
}

#[test]
fn an_invalid_coverage_is_refused_and_empty_input_is_empty() {
    let r = coverage_union(&[rect(0.0, 0.0, 2.0, 2.0), rect(1.0, 1.0, 3.0, 3.0)]);
    assert!(matches!(
        r,
        Err(GeometryError::InvalidCoverage {
            first: CoverageIssue::Overlap { a: 0, b: 1 }
        })
    ));
    assert!(coverage_union(&[]).unwrap().polygons().is_empty());
}

// ===========================================================================
// Polygonize
// ===========================================================================

fn rings_as_lines(polys: &[Polygon]) -> MultiLineString {
    let mut out = Vec::new();
    for p in polys {
        for r in core::iter::once(p.exterior()).chain(p.holes()) {
            let mut v = r.coords().as_slice().to_vec();
            v.push(v[0]);
            out.push(LineString::new(CoordSeq::new(v)).unwrap());
        }
    }
    MultiLineString::new(out)
}

fn normalized(polys: &[Polygon]) -> Geometry {
    crate::ops::normalize(&Geometry::GeometryCollection(
        crate::types::GeometryCollection::new(
            polys.iter().cloned().map(Geometry::Polygon).collect(),
        )
        .unwrap(),
    ))
}

#[test]
fn a_noded_grid_polygonizes_to_its_cells() {
    let mut rng = Rng(3);
    let grid = jittered_grid(&mut rng, 4, 3, 10, 3);
    // Every shared edge appears twice in the input (once per cell):
    // identical segments are one edge.
    let out = polygonize_noded(&rings_as_lines(&grid)).unwrap();
    assert_eq!(normalized(&out.polygons), normalized(&grid));
    assert!(out.dangles.line_strings().is_empty());
    assert!(out.cut_edges.line_strings().is_empty());
    assert!(out.invalid_rings.line_strings().is_empty());
}

#[test]
fn dangles_cut_edges_and_nested_islands_are_classified() {
    let lines = MultiLineString::new(vec![
        // The outer square, with a vertex at (4, 2) for the bridge.
        line(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 2.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]),
        // An island, and an island inside it.
        line(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]),
        line(&[(1.5, 1.5), (2.5, 1.5), (2.5, 2.5), (1.5, 2.5), (1.5, 1.5)]),
        // A dangle off the outer square's corner.
        line(&[(4.0, 0.0), (6.0, 0.0)]),
        // A cut edge bridging to a triangle.
        line(&[(4.0, 2.0), (8.0, 2.0)]),
        line(&[(8.0, 2.0), (10.0, 1.0), (10.0, 3.0), (8.0, 2.0)]),
    ]);
    let out = polygonize_noded(&lines).unwrap();
    let expected = [
        poly(
            &[(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (4.0, 4.0), (0.0, 4.0)],
            &[&[(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)]],
        ),
        poly(
            &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
            &[&[(1.5, 1.5), (1.5, 2.5), (2.5, 2.5), (2.5, 1.5)]],
        ),
        rect(1.5, 1.5, 2.5, 2.5),
        poly(&[(8.0, 2.0), (10.0, 1.0), (10.0, 3.0)], &[]),
    ];
    assert_eq!(normalized(&out.polygons), normalized(&expected));
    assert_eq!(out.dangles.line_strings().len(), 1);
    assert_eq!(out.cut_edges.line_strings().len(), 1);
    assert_eq!(out.cut_edges.line_strings()[0].coords().as_slice().len(), 2);
    assert!(
        out.invalid_rings.line_strings().is_empty(),
        "noded input has no invalid ring"
    );
}

#[test]
fn a_face_with_a_pinched_hole_is_one_polygon_with_a_touching_hole() {
    // A square with a diamond inside touching it at (0, 2).
    let lines = MultiLineString::new(vec![
        line(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]),
        line(&[(0.0, 2.0), (2.0, 1.0), (3.0, 2.0), (2.0, 3.0), (0.0, 2.0)]),
    ]);
    let out = polygonize_noded(&lines).unwrap();
    assert_eq!(out.polygons.len(), 2);
    let holes: usize = out.polygons.iter().map(|p| p.holes().len()).sum();
    assert_eq!(holes, 1);
    for p in &out.polygons {
        assert!(is_valid(p));
    }
}

#[test]
fn a_tree_is_all_dangles() {
    let out = polygonize_noded(&MultiLineString::new(vec![
        line(&[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0)]),
        line(&[(1.0, 0.0), (1.0, -1.0)]),
    ]))
    .unwrap();
    assert!(out.polygons.is_empty());
    let n: usize = out
        .dangles
        .line_strings()
        .iter()
        .map(|l| l.coords().as_slice().len() - 1)
        .sum();
    assert_eq!(n, 3, "every segment is a dangle");
}

/// The first offending pair, by brute force over the sorted distinct
/// segments.
fn brute_first_unnoded(segs: &[(Coord2d, Coord2d)]) -> Option<(usize, usize)> {
    for i in 0..segs.len() {
        for j in i + 1..segs.len() {
            if !polygonize::noded_pair(segs[i], segs[j]) {
                return Some((i, j));
            }
        }
    }
    None
}

#[test]
fn unnoded_input_is_refused_naming_a_brute_force_pair() {
    let cases = [
        // A proper crossing.
        MultiLineString::new(vec![
            line(&[(0.0, 0.0), (2.0, 2.0)]),
            line(&[(0.0, 2.0), (2.0, 0.0)]),
        ]),
        // An endpoint in another segment's interior (a T).
        MultiLineString::new(vec![
            line(&[(0.0, 0.0), (2.0, 0.0)]),
            line(&[(1.0, 0.0), (1.0, 1.0)]),
        ]),
        // A collinear overlap.
        MultiLineString::new(vec![
            line(&[(0.0, 0.0), (2.0, 0.0)]),
            line(&[(1.0, 0.0), (3.0, 0.0)]),
        ]),
        // A collinear overlap from a shared endpoint, either way round.
        MultiLineString::new(vec![
            line(&[(0.0, 0.0), (2.0, 1.0)]),
            line(&[(0.0, 0.0), (4.0, 2.0)]),
        ]),
        MultiLineString::new(vec![
            line(&[(4.0, 2.0), (2.0, 1.0)]),
            line(&[(0.0, 0.0), (4.0, 2.0)]),
        ]),
    ];
    for lines in &cases {
        let Err(GeometryError::NotNoded { first_pair }) = polygonize_noded(lines) else {
            panic!("refused");
        };
        let mut segs: Vec<(Coord2d, Coord2d)> = lines
            .line_strings()
            .iter()
            .flat_map(|l| {
                l.coords()
                    .as_slice()
                    .windows(2)
                    .map(|w| {
                        if coord_exact_cmp(w[0], w[1]).is_le() {
                            (w[0], w[1])
                        } else {
                            (w[1], w[0])
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        segs.sort_by(|s, t| coord_exact_cmp(s.0, t.0).then(coord_exact_cmp(s.1, t.1)));
        let (i, j) = brute_first_unnoded(&segs).unwrap();
        assert_eq!(first_pair, (segs[i], segs[j]));
    }
    // Shared endpoints only: noded.
    assert!(
        polygonize_noded(&MultiLineString::new(vec![
            line(&[(0.0, 0.0), (1.0, 0.0)]),
            line(&[(1.0, 0.0), (1.0, 1.0)]),
        ]))
        .is_ok()
    );
}

// ===========================================================================
// Polygonize internals: sorted vectors against ordered-map references
// ===========================================================================

/// The distinct segments as an ordered set computes them: first seen kept,
/// oriented, then a stable sort by `coord_exact_cmp`.
fn reference_distinct_segments(segs: &[(Coord2d, Coord2d)]) -> Vec<(Coord2d, Coord2d)> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for &(a, b) in segs {
        if seen.insert(edge_key(a, b)) {
            out.push(if coord_exact_cmp(a, b).is_le() {
                (a, b)
            } else {
                (b, a)
            });
        }
    }
    out.sort_by(|s, t| coord_exact_cmp(s.0, t.0).then(coord_exact_cmp(s.1, t.1)));
    out
}

/// The vertex numbering as a first-seen map computes it.
fn reference_number_vertices(
    segments: &[(Coord2d, Coord2d)],
) -> (Vec<Coord2d>, Vec<(usize, usize)>) {
    let mut ids = std::collections::BTreeMap::new();
    let mut coords: Vec<Coord2d> = Vec::new();
    let mut edges = Vec::new();
    for &(a, b) in segments {
        let mut id = |c: Coord2d| {
            *ids.entry(super::graph::key(c)).or_insert_with(|| {
                coords.push(c);
                coords.len() - 1
            })
        };
        edges.push((id(a), id(b)));
    }
    (coords, edges)
}

/// Random segments on a small lattice with negative values and signed
/// zeros, repeated and reversed, so that duplicates, shared endpoints,
/// touching boxes and every kind of offending pair occur.
fn lattice_segments(rng: &mut Rng, n: usize) -> Vec<(Coord2d, Coord2d)> {
    let v = |rng: &mut Rng| -> f64 {
        let x = rng.int(-3, 3) as f64;
        if x == 0.0 && rng.next().is_multiple_of(2) {
            -0.0
        } else {
            x
        }
    };
    let mut out = Vec::new();
    while out.len() < n {
        let (a, b) = (c(v(rng), v(rng)), c(v(rng), v(rng)));
        if a == b {
            continue;
        }
        out.push((a, b));
        if rng.next().is_multiple_of(4) {
            out.push((b, a));
        }
    }
    out
}

fn bits(s: &[(Coord2d, Coord2d)]) -> Vec<[u64; 4]> {
    s.iter()
        .map(|(a, b)| [a.x.to_bits(), a.y.to_bits(), b.x.to_bits(), b.y.to_bits()])
        .collect()
}

fn as_borrowed_lines(segs: &[(Coord2d, Coord2d)]) -> MultiLineString {
    MultiLineString::new(
        segs.iter()
            .map(|&(a, b)| LineString::new(CoordSeq::new(vec![a, b])).unwrap())
            .collect(),
    )
}

#[test]
fn the_sorted_segments_and_numbering_equal_the_ordered_maps_bit_for_bit() {
    let mut rng = Rng(0x0DD5_EED5);
    for case in 0..400 {
        let raw = lattice_segments(&mut rng, 1 + case % 40);
        let lines = as_borrowed_lines(&raw);
        let got = polygonize::distinct_segments(crate::predicates::borrowed(&lines));
        let want = reference_distinct_segments(&raw);
        assert_eq!(bits(&got), bits(&want), "case {case}");
        let (gc, ge) = polygonize::number_vertices(&got);
        let (wc, we) = reference_number_vertices(&want);
        assert_eq!(ge, we, "case {case}");
        let cb = |v: &[Coord2d]| {
            v.iter()
                .map(|p| (p.x.to_bits(), p.y.to_bits()))
                .collect::<Vec<_>>()
        };
        assert_eq!(cb(&gc), cb(&wc), "case {case}");
    }
}

/// The loop split as an ordered map computes it.
fn reference_split_loops(walk: &[usize]) -> Vec<Vec<usize>> {
    let mut loops = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    let mut at = std::collections::BTreeMap::new();
    for &v in walk {
        if let Some(&p) = at.get(&v) {
            let lp: Vec<usize> = stack.drain(p..).collect();
            for u in &lp {
                at.remove(u);
            }
            loops.push(lp);
        }
        at.insert(v, stack.len());
        stack.push(v);
    }
    if !stack.is_empty() {
        loops.push(stack);
    }
    loops
}

#[test]
fn the_loop_splitter_equals_the_ordered_map_walk_after_walk() {
    // Random walks over few vertices, so vertices repeat often; one
    // splitter serves every walk, so stale scratch would show.
    let mut rng = Rng(0x100F_5EED);
    let mut splitter = graph::LoopSplitter::new(12);
    for case in 0..3000 {
        let len = 1 + (rng.next() % 30) as usize;
        let walk: Vec<usize> = (0..len).map(|_| (rng.next() % 12) as usize).collect();
        let mut got: Vec<Vec<usize>> = Vec::new();
        splitter.split(walk.iter().copied(), |lp| got.push(lp.to_vec()));
        assert_eq!(got, reference_split_loops(&walk), "case {case}: {walk:?}");
    }
}

// ===========================================================================
// Rings whose area is beyond the range of f64
// ===========================================================================

/// Two valid squares whose union has a signed area beyond the range of
/// `f64`: each square's area is finite, and their union's is not.
fn overflowing_pair() -> Vec<Polygon> {
    let s = 9.0e153;
    vec![
        rect(-s, -s / 2.0, 0.0, s / 2.0),
        rect(0.0, -s / 2.0, s, s / 2.0),
    ]
}

#[test]
fn a_union_whose_area_overflows_is_refused_not_a_panic() {
    let cov = overflowing_pair();
    assert!(
        coverage_validate(&cov).is_empty(),
        "the input is a valid coverage"
    );
    assert_eq!(coverage_union(&cov), Err(GeometryError::NonFiniteResult));
}

#[test]
fn a_face_whose_area_overflows_is_reported_not_dropped() {
    let s = 9.0e153;
    let lines = MultiLineString::new(vec![line(&[
        (-s, -s / 2.0),
        (s, -s / 2.0),
        (s, s / 2.0),
        (-s, s / 2.0),
        (-s, -s / 2.0),
    ])]);
    let out = polygonize_noded(&lines).unwrap();
    assert!(out.polygons.is_empty());
    assert_eq!(
        out.invalid_rings.line_strings().len(),
        1,
        "the face is reported"
    );
}

#[test]
fn orientation_at_the_smallest_vertex_matches_the_area_sign() {
    // Star-shaped rings (simple), with collinear vertices on axis-parallel
    // runs, in both directions, rotated to start at their smallest vertex as
    // `ring_of` returns them.
    let mut rng = Rng(0x0C1E_7A7E);
    for case in 0..2000 {
        let n = 4 + (rng.next() % 30) as usize;
        let mut pts: Vec<Coord2d> = (0..n)
            .map(|k| {
                let t = k as f64 / n as f64 * std::f64::consts::TAU;
                let r = 3.0 + (rng.next() % 6) as f64;
                c((r * t.cos()).round(), (r * t.sin()).round())
            })
            .collect();
        pts.dedup();
        while pts.len() > 1 && pts.first() == pts.last() {
            pts.pop();
        }
        if case % 2 == 1 {
            pts.reverse();
        }
        let Ok(ring) = LinearRing::new(CoordSeq::new(pts)) else {
            continue;
        };
        if !ring.is_simple() || ring.signed_area() == 0.0 {
            continue;
        }
        let start = (0..ring.num_vertices())
            .min_by(|&i, &j| coord_exact_cmp(ring.coords()[i], ring.coords()[j]))
            .unwrap();
        let cs = ring.coords().as_slice();
        let rotated: Vec<Coord2d> = cs[start..].iter().chain(&cs[..start]).copied().collect();
        let ring = LinearRing::new(CoordSeq::new(rotated)).unwrap();
        assert_eq!(
            super::graph::runs_counterclockwise(&ring),
            ring.signed_area() > 0.0,
            "case {case}: {ring:?}"
        );
    }
}
