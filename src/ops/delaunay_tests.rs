// src/ops/delaunay_tests.rs

//! The Delaunay triangulation: hand cases for the tie rule, brute-force
//! empty-circle checks, canonical-order and permutation invariance,
//! degenerate input, and the oracle dump writer.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::unreadable_literal
)]

use arithmetica::Sign;

use super::delaunay::{Triangulation, delaunay};
use crate::codec::encode;
use crate::types::{Coord2d, MultiPoint, Point, UnvouchedReason};

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn multipoint(pts: &[(f64, f64)]) -> MultiPoint {
    MultiPoint::new(
        pts.iter()
            .map(|&(x, y)| Point::from_finite(c(x, y)))
            .collect(),
    )
}

fn tri_coords(t: &Triangulation) -> Vec<[(f64, f64); 3]> {
    t.triangles()
        .iter()
        .map(|tr| tr.map(|v| (t.points()[v as usize].x, t.points()[v as usize].y)))
        .collect()
}

fn edge_coords(t: &Triangulation) -> Vec<[(f64, f64); 2]> {
    // Read back from the MultiLineString, so the public form is what is
    // compared.
    t.edges()
        .line_strings()
        .iter()
        .map(|l| {
            let s = l.coords().as_slice();
            [(s[0].x, s[0].y), (s[1].x, s[1].y)]
        })
        .collect()
}

/// A deterministic generator (SplitMix64), so a failure reproduces.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    /// A float in [0, 1) with 53 random bits.
    fn f(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            v.swap(i, self.below(i + 1));
        }
    }
}

// ======================================================================
// The tie rule, pinned by hand
// ======================================================================

/// Every unit square of an integer grid is a cocircular face (its circle
/// has radius √2/2, and the nearest other lattice point is at distance
/// √10/2), so the canonical result cuts each square by the fan from its
/// lexicographically smallest corner.
#[test]
fn a_grid_is_cut_by_the_fan_rule_exactly() {
    let k = 4usize;
    let mut pts = Vec::new();
    for i in 0..=k {
        for j in 0..=k {
            pts.push((i as f64, j as f64));
        }
    }
    let t = delaunay(&multipoint(&pts)).unwrap();
    let mut expected = Vec::new();
    for i in 0..k {
        for j in 0..k {
            let (x, y) = (i as f64, j as f64);
            expected.push([(x, y), (x + 1.0, y), (x + 1.0, y + 1.0)]);
            expected.push([(x, y), (x + 1.0, y + 1.0), (x, y + 1.0)]);
        }
    }
    expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(tri_coords(&t), expected);
}

/// All lattice points on the circle of radius 65 (there are 36) are one
/// cocircular face: the canonical triangulation is the fan from the
/// lexicographically smallest, (-65, 0).
#[test]
fn the_circle_lattice_fans_from_its_smallest_vertex() {
    let pts = circle_lattice(65);
    assert_eq!(pts.len(), 36);
    let t = delaunay(&multipoint(&pts)).unwrap();
    assert_eq!(t.triangles().len(), 34);
    for tr in tri_coords(&t) {
        assert!(tr.contains(&(-65.0, 0.0)), "{tr:?} is not a fan triangle");
    }
}

/// Lattice points on the circle x² + y² = r², by scan.
fn circle_lattice(r: i64) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for x in -r..=r {
        let y2 = r * r - x * x;
        let y = (y2 as f64).sqrt().round() as i64;
        if y * y == y2 {
            out.push((x as f64, y as f64));
            if y != 0 {
                out.push((x as f64, -y as f64));
            }
        }
    }
    out
}

// ======================================================================
// Exactness, brute force
// ======================================================================

/// The defining property, checked against every input point with the same
/// exact predicate the library uses: no point strictly inside any
/// circumcircle. Random sets, clustered sets, and near-cocircular sets
/// (grids perturbed by a few ulps).
#[test]
fn every_circumcircle_is_empty_on_generated_sets() {
    let mut rng = Rng(0xDE1A_0001);
    for case in 0..160 {
        let pts = generated_set(&mut rng, case);
        let t = delaunay(&multipoint(&pts)).unwrap();
        check_empty_circles(&t);
        check_structure(&t);
    }
}

fn generated_set(rng: &mut Rng, case: usize) -> Vec<(f64, f64)> {
    let n = 3 + rng.below(26);
    let mut pts = Vec::with_capacity(n);
    match case % 4 {
        // Random in a box.
        0 => {
            for _ in 0..n {
                pts.push((rng.f() * 100.0, rng.f() * 100.0));
            }
        }
        // Clustered: a few tight clusters far apart.
        1 => {
            for _ in 0..n {
                let (cx, cy) = match rng.below(3) {
                    0 => (0.0, 0.0),
                    1 => (1e6, 3.0),
                    _ => (-7.0, 1e-3),
                };
                pts.push((cx + rng.f(), cy + rng.f()));
            }
        }
        // A grid perturbed by a few ulps: near-cocircular everywhere.
        2 => {
            let k = 2 + rng.below(4);
            for i in 0..=k {
                for j in 0..=k {
                    let mut x = i as f64;
                    let mut y = j as f64;
                    for _ in 0..rng.below(3) {
                        x = if rng.below(2) == 0 {
                            x.next_up()
                        } else {
                            x.next_down()
                        };
                    }
                    for _ in 0..rng.below(3) {
                        y = if rng.below(2) == 0 {
                            y.next_up()
                        } else {
                            y.next_down()
                        };
                    }
                    pts.push((x, y));
                }
            }
        }
        // Mixed magnitudes.
        _ => {
            for _ in 0..n {
                let e = [-300i32, -30, -3, 0, 3, 30, 300][rng.below(7)];
                let m = [-300i32, -30, -3, 0, 3, 30, 300][rng.below(7)];
                pts.push((
                    (rng.f() - 0.5) * 10f64.powi(e),
                    (rng.f() - 0.5) * 10f64.powi(m),
                ));
            }
        }
    }
    pts
}

fn check_empty_circles(t: &Triangulation) {
    for tr in t.triangles() {
        let [a, b, c] = tr.map(|v| t.points()[v as usize]);
        for p in t.points() {
            let s = arithmetica::incircle([a.x, a.y], [b.x, b.y], [c.x, c.y], [p.x, p.y]).unwrap();
            assert_ne!(s, Sign::Positive, "{p:?} is inside a circumcircle");
        }
    }
}

/// Structure: every triangle counterclockwise; every point a vertex; the
/// canonical order as documented; and Euler's relation `e = n + t − 1` for
/// a triangulation of the hull using every point.
fn check_structure(t: &Triangulation) {
    let n = t.points().len();
    let lex = |a: Coord2d, b: Coord2d| {
        a.x.partial_cmp(&b.x)
            .unwrap()
            .then(a.y.partial_cmp(&b.y).unwrap())
    };
    let mut seen = vec![false; n];
    for tr in t.triangles() {
        let [a, b, c] = tr.map(|v| t.points()[v as usize]);
        assert_eq!(
            arithmetica::orient2d([a.x, a.y], [b.x, b.y], [c.x, c.y]).unwrap(),
            Sign::Positive,
            "triangle not counterclockwise"
        );
        assert!(
            lex(a, b).is_lt() && lex(a, c).is_lt(),
            "triangle does not start at its smallest vertex"
        );
        for &v in tr {
            seen[v as usize] = true;
        }
    }
    if !t.triangles().is_empty() {
        assert!(seen.iter().all(|&s| s), "a point is not a vertex");
        assert_eq!(
            t.edges().line_strings().len(),
            n + t.triangles().len() - 1,
            "Euler count"
        );
        let trs = tri_coords(t);
        let mut sorted = trs.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(trs, sorted, "triangles not in canonical order");
    }
    let es = edge_coords(t);
    let mut sorted = es.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(es, sorted, "edges not in canonical order");
    for e in &es {
        assert!(e[0] < e[1], "edge not from its smaller endpoint");
    }
}

/// Points exactly on hull edges and on interior edges go through the
/// open-segment arm of the ghost test and the chord argument: all must end
/// as vertices with the empty-circle property intact.
#[test]
fn points_on_edges_are_vertices() {
    let cases: &[&[(f64, f64)]] = &[
        // On a hull edge.
        &[(0.0, 0.0), (4.0, 0.0), (2.0, 3.0), (2.0, 0.0)],
        // On an interior edge (the diagonal).
        &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (1.0, 1.0)],
        // Collinear run on the hull.
        &[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (1.5, 2.0)],
        // Collinear beyond a hull edge's span.
        &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (3.0, 0.0)],
    ];
    for pts in cases {
        let t = delaunay(&multipoint(pts)).unwrap();
        check_empty_circles(&t);
        check_structure(&t);
        assert_eq!(t.points().len(), pts.len());
    }
}

// ======================================================================
// Canonicality
// ======================================================================

/// Permutations of one point set give byte-identical outputs: the
/// triangulation depends only on the set.
#[test]
fn permutations_give_byte_identical_outputs() {
    let mut rng = Rng(0xDE1A_0002);
    let mut grid = Vec::new();
    for i in 0..6 {
        for j in 0..5 {
            grid.push((f64::from(i), f64::from(j)));
        }
    }
    let mut big: Vec<(f64, f64)> = (0..2000).map(|_| (rng.f() * 1e4, rng.f() * 1e4)).collect();
    // Make the big case carry ties too: a copy of a 5-by-5 block of integers.
    for i in 0..5 {
        for j in 0..5 {
            big.push((f64::from(i) + 20000.0, f64::from(j)));
        }
    }
    let random: Vec<(f64, f64)> = (0..57).map(|_| (rng.f(), rng.f())).collect();
    let cases: Vec<Vec<(f64, f64)>> = vec![
        grid,
        circle_lattice(65),
        circle_lattice(25),
        random,
        vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (5.0, 0.0)],
        big,
    ];
    for pts in cases {
        let base = delaunay(&multipoint(&pts)).unwrap();
        let (bt, be) = (tri_coords(&base), edge_coords(&base));
        let (bp, bl) = (encode(&base.to_polygons()), encode(&base.edges()));
        let mut shuffled = pts.clone();
        for _ in 0..20 {
            rng.shuffle(&mut shuffled);
            let t = delaunay(&multipoint(&shuffled)).unwrap();
            assert_eq!(tri_coords(&t), bt, "triangles differ under permutation");
            assert_eq!(edge_coords(&t), be, "edges differ under permutation");
            assert_eq!(encode(&t.to_polygons()), bp, "polygon bytes differ");
            assert_eq!(encode(&t.edges()), bl, "edge bytes differ");
        }
    }
}

// ======================================================================
// Degenerate input, duplicates, the error surface
// ======================================================================

#[test]
fn collinear_input_gives_the_chain() {
    // Shuffled collinear points: no triangles; the edges are between
    // lexicographically consecutive points.
    let t = delaunay(&multipoint(&[
        (4.0, 4.0),
        (0.0, 0.0),
        (3.0, 3.0),
        (1.0, 1.0),
    ]))
    .unwrap();
    assert!(t.triangles().is_empty());
    assert_eq!(
        edge_coords(&t),
        vec![
            [(0.0, 0.0), (1.0, 1.0)],
            [(1.0, 1.0), (3.0, 3.0)],
            [(3.0, 3.0), (4.0, 4.0)],
        ]
    );
    assert!(t.to_polygons().geometries().is_empty());

    let two = delaunay(&multipoint(&[(5.0, 1.0), (2.0, 1.0)])).unwrap();
    assert_eq!(edge_coords(&two), vec![[(2.0, 1.0), (5.0, 1.0)]]);

    let one = delaunay(&multipoint(&[(2.0, 1.0)])).unwrap();
    assert!(one.triangles().is_empty() && one.edges().line_strings().is_empty());
    assert_eq!(one.points().len(), 1);

    let none = delaunay(&MultiPoint::empty()).unwrap();
    assert!(none.points().is_empty() && none.edges().line_strings().is_empty());
}

#[test]
fn duplicates_and_signed_zeros_are_one_point() {
    let with_dupes = multipoint(&[
        (0.0, 0.0),
        (4.0, 0.0),
        (-0.0, 0.0),
        (2.0, 3.0),
        (4.0, 0.0),
        (0.0, -0.0),
        (7.0, 1.0),
    ]);
    let plain = multipoint(&[(0.0, 0.0), (4.0, 0.0), (2.0, 3.0), (7.0, 1.0)]);
    let a = delaunay(&with_dupes).unwrap();
    let b = delaunay(&plain).unwrap();
    assert_eq!(a.points().len(), 4);
    assert_eq!(tri_coords(&a), tri_coords(&b));
    assert_eq!(encode(&a.edges()), encode(&b.edges()));
}

/// The defensive refusal's class: a triangulation that fails certification
/// is the wrong answer, not a damaged one.
#[test]
fn the_defensive_reason_is_structurally_wrong() {
    assert!(UnvouchedReason::NotDelaunay.is_structurally_wrong());
    assert!(!UnvouchedReason::NotDelaunay.is_locally_invalid());
}

/// `delaunay` takes any `Relatable`: a polygon's distinct vertices
/// triangulate like the same points given as a multipoint.
#[test]
fn any_relatable_is_accepted() {
    use crate::types::{CoordSeq, LinearRing, Polygon};
    let ring = LinearRing::new(CoordSeq::new(vec![
        c(0.0, 0.0),
        c(4.0, 0.0),
        c(4.0, 4.0),
        c(0.0, 4.0),
    ]))
    .unwrap();
    let poly = Polygon::new_simple(ring).unwrap();
    let a = delaunay(&poly).unwrap();
    let b = delaunay(&multipoint(&[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
    ]))
    .unwrap();
    assert_eq!(tri_coords(&a), tri_coords(&b));
    assert_eq!(a.triangles().len(), 2);
}

// ======================================================================
// The oracle dump
// ======================================================================

/// Writes `delaunay_cases.csv` in the system temporary directory, for an
/// independent exact check: the case families above, each with the points
/// (as `f64` bits, so the file is exact) and the computed triangles and
/// edges.
///
/// ```text
/// cargo test -p geometres --release --lib ops::delaunay_tests::dump_delaunay_oracle_cases -- --ignored
/// ```
#[test]
#[ignore = "writes the case file delaunay_cases.csv to the temporary directory"]
fn dump_delaunay_oracle_cases() {
    use std::fmt::Write as _;
    let mut rng = Rng(0xDE1A_0003);
    let mut out = String::new();
    let mut cases = 0usize;

    let emit = |name: &str, pts: &[(f64, f64)], out: &mut String| {
        let t = delaunay(&multipoint(pts)).unwrap();
        writeln!(
            out,
            "CASE {name} {} {} {}",
            t.points().len(),
            t.triangles().len(),
            t.edges().line_strings().len()
        )
        .unwrap();
        for p in t.points() {
            writeln!(out, "P {} {}", p.x.to_bits(), p.y.to_bits()).unwrap();
        }
        for tr in t.triangles() {
            writeln!(out, "T {} {} {}", tr[0], tr[1], tr[2]).unwrap();
        }
        for l in t.edges().line_strings() {
            let s = l.coords().as_slice();
            let ix = |q: Coord2d| {
                t.points()
                    .iter()
                    .position(|&p| p.x.to_bits() == q.x.to_bits() && p.y.to_bits() == q.y.to_bits())
                    .unwrap()
            };
            writeln!(out, "E {} {}", ix(s[0]), ix(s[1])).unwrap();
        }
    };

    // 5,000+ random sets, most small (the oracle's tie check is O(n^2) per
    // triangle), some larger.
    for k in 0..4600 {
        let n = 3 + rng.below(14);
        let pts: Vec<(f64, f64)> = (0..n)
            .map(|_| {
                if k % 3 == 0 {
                    // Integer coordinates: cocircular ties are common.
                    (rng.below(12) as f64, rng.below(12) as f64)
                } else {
                    (rng.f() * 100.0, rng.f() * 100.0)
                }
            })
            .collect();
        emit(&format!("random{k}"), &pts, &mut out);
        cases += 1;
    }
    for k in 0..500 {
        let pts = generated_set(&mut rng, k);
        emit(&format!("family{k}"), &pts, &mut out);
        cases += 1;
    }
    // Grids: every unit square a cocircular face.
    for k in 2..=6i32 {
        let mut pts = Vec::new();
        for i in 0..=k {
            for j in 0..=k {
                pts.push((f64::from(i), f64::from(j)));
            }
        }
        emit(&format!("grid{k}"), &pts, &mut out);
        cases += 1;
    }
    // Many cocircular lattice points, with and without interior points.
    for r in [25i64, 65] {
        let pts = circle_lattice(r);
        emit(&format!("circle{r}"), &pts, &mut out);
        let mut with_center = pts.clone();
        with_center.push((0.0, 0.0));
        with_center.push((1.0, 2.0));
        emit(&format!("circle{r}c"), &with_center, &mut out);
        cases += 2;
    }
    // Collinear, short chains, mixed magnitudes.
    emit(
        "collinear",
        &[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (9.0, 9.0)],
        &mut out,
    );
    emit("pair", &[(3.0, 1.0), (0.0, 1.0)], &mut out);
    emit(
        "magnitudes",
        &[
            (1e-300, 0.0),
            (1e300, 1.0),
            (-1e300, 1e-300),
            (0.0, -1e300),
            (3.5, 7.25),
            (-2.0, 1e30),
        ],
        &mut out,
    );
    cases += 3;
    assert!(cases >= 5000, "{cases} cases");

    let path = std::env::temp_dir().join("delaunay_cases.csv");
    std::fs::write(&path, out).unwrap();
    println!("oracle dump: {cases} cases -> {}", path.display());
}
