// src/ops/basic_tests.rs

//! Exact properties of the basic operations and convex hull, on 10,000+
//! generated cases each (collinear- and duplicate-heavy, on a small integer
//! grid), with hand cases. One property is checked mechanically for every
//! operation: every output coordinate is an input coordinate, bit for bit.

// Test-code lint posture: seeds are written as hex tags (byte grouping is
// not meaningful); brute-force oracles favor the plain loop and the explicit
// comparison over the idiomatic form, so that they read as the definition
// they check.
#![allow(
    clippy::unusual_byte_groupings,
    clippy::collapsible_if,
    clippy::match_same_arms,
    clippy::type_complexity,
    clippy::nonminimal_bool,
    clippy::if_not_else,
    clippy::needless_range_loop,
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

use std::collections::{HashMap, HashSet};

use super::basic::{coords_of, lex, lex_seq};
use super::*;
use crate::GeometryError;
use crate::codec::encode;
use crate::exact::kernel::{on_segment as on_seg, orient};
use crate::predicates::Borrowed;
use crate::topology::de9im::IntersectionMatrix;
use crate::topology::label::{Dimension, Location};
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString,
    MultiPoint, MultiPolygon, Point, Polygon,
};
use arithmetica::Sign;

// ===========================================================================
// Generators
// ===========================================================================

/// SplitMix64: deterministic, no dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn coin(&mut self, k: u64) -> bool {
        self.below(k) == 0
    }
    fn grid(&mut self, n: u64) -> f64 {
        self.below(n) as f64
    }
    fn pt(&mut self, n: u64) -> Coord2d {
        Coord2d::new(self.grid(n), self.grid(n))
    }
}

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn line(p: &[(f64, f64)]) -> LineString {
    LineString::new(CoordSeq::new(p.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn ring(p: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(p.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new_simple(ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])).unwrap()
}

/// A random line string on an `n`-grid, 2–7 vertices, sometimes closed.
fn gen_line(r: &mut Rng, n: u64) -> Option<LineString> {
    let k = 2 + r.below(6) as usize;
    let mut v: Vec<Coord2d> = Vec::new();
    while v.len() < k {
        let p = r.pt(n);
        if v.last() != Some(&p) {
            v.push(p);
        }
    }
    if r.coin(5) && v.len() >= 3 && v[0] != v[v.len() - 1] {
        v.push(v[0]);
    }
    LineString::new(CoordSeq::new(v)).ok()
}

fn gen_lines(r: &mut Rng, n: u64, max: u64) -> MultiLineString {
    let k = 1 + r.below(max);
    MultiLineString::new((0..k).filter_map(|_| gen_line(r, n)).collect())
}

/// A polygon in cell `(i, j)` of side 10: a rectangle, a triangle with a
/// collinear vertex on its base, or a rectangle with a rectangular hole.
fn gen_polygon(r: &mut Rng, i: f64, j: f64) -> Polygon {
    let (x, y) = (10.0 * i, 10.0 * j);
    let w = 2.0 + r.grid(6);
    let h = 2.0 + r.grid(6);
    match r.below(3) {
        0 => rect(x, y, x + w, y + h),
        1 => Polygon::new_simple(ring(&[
            (x, y),
            (x + 1.0, y),
            (x + w, y),
            (x + w / 2.0, y + h),
        ]))
        .unwrap(),
        _ => {
            let outer = ring(&[(x, y), (x + 8.0, y), (x + 8.0, y + 8.0), (x, y + 8.0)]);
            let hole = ring(&[
                (x + 2.0, y + 2.0),
                (x + 2.0, y + 5.0),
                (x + 5.0, y + 5.0),
                (x + 5.0, y + 2.0),
            ]);
            Polygon::new(outer, vec![hole]).unwrap()
        }
    }
}

fn gen_multipolygon(r: &mut Rng) -> MultiPolygon {
    let k = 1 + r.below(3);
    let mut cells: Vec<(f64, f64)> = Vec::new();
    while (cells.len() as u64) < k {
        let cell = (r.grid(3), r.grid(3));
        if !cells.contains(&cell) {
            cells.push(cell);
        }
    }
    MultiPolygon::new(cells.iter().map(|&(i, j)| gen_polygon(r, i, j)).collect()).unwrap()
}

fn gen_points(r: &mut Rng, n: u64) -> MultiPoint {
    let k = r.below(6);
    MultiPoint::new((0..k).map(|_| Point::new(r.pt(n)).unwrap()).collect())
}

/// A random geometry of any type; collections hold polygonal leaves in
/// distinct cells (so `relate` accepts them), lines and points anywhere.
fn gen_geometry(r: &mut Rng) -> Geometry {
    match r.below(8) {
        0 => Geometry::Point(Point::new(r.pt(30)).unwrap()),
        1 => Geometry::MultiPoint(gen_points(r, 30)),
        2 => gen_line(r, 30).map_or(
            Geometry::MultiPoint(MultiPoint::empty()),
            Geometry::LineString,
        ),
        3 => Geometry::MultiLineString(gen_lines(r, 30, 4)),
        4 => {
            let (i, j) = (r.grid(3), r.grid(3));
            Geometry::Polygon(gen_polygon(r, i, j))
        }
        5 => Geometry::MultiPolygon(gen_multipolygon(r)),
        _ => {
            let mut gs = Vec::new();
            if r.coin(2) {
                gs.push(Geometry::MultiPolygon(gen_multipolygon(r)));
            }
            if r.coin(2) {
                gs.push(Geometry::MultiLineString(gen_lines(r, 30, 3)));
            }
            if r.coin(2) {
                if let Some(l) = gen_line(r, 30) {
                    gs.push(Geometry::LineString(l));
                }
            }
            if r.coin(2) {
                gs.push(Geometry::MultiPoint(gen_points(r, 30)));
            }
            let inner = GeometryCollection::new(gs).unwrap();
            if r.coin(4) {
                Geometry::GeometryCollection(
                    GeometryCollection::new(vec![
                        Geometry::GeometryCollection(inner),
                        Geometry::Point(Point::new(r.pt(30)).unwrap()),
                    ])
                    .unwrap(),
                )
            } else {
                Geometry::GeometryCollection(inner)
            }
        }
    }
}

// ===========================================================================
// Every output coordinate is an input coordinate
// ===========================================================================

fn bits(p: Coord2d) -> (u64, u64) {
    (p.x.to_bits(), p.y.to_bits())
}

/// Every coordinate of `out` is, bit for bit, a coordinate of `inputs` or
/// one of `extra`.
fn assert_test5(out: &Geometry, inputs: &[Coord2d], extra: &[Coord2d], what: &str) {
    let have: HashSet<(u64, u64)> = inputs.iter().chain(extra).copied().map(bits).collect();
    for p in coords_of(Borrowed::from_geometry(out)) {
        assert!(
            have.contains(&bits(p)),
            "{what}: output coordinate {p:?} is not an input coordinate"
        );
    }
}

fn mls(m: MultiLineString) -> Geometry {
    Geometry::MultiLineString(m)
}

fn all_coords(g: &Geometry) -> Vec<Coord2d> {
    coords_of(Borrowed::from_geometry(g))
}

const CASES: usize = 10_000;

// ===========================================================================
// convex_hull
// ===========================================================================

fn hull_points(r: &mut Rng) -> Vec<Coord2d> {
    match r.below(4) {
        // Scattered, duplicate-heavy on a small grid.
        0 | 1 => (0..r.below(20)).map(|_| r.pt(6)).collect(),
        // Collinear: on y = a x + b.
        2 => {
            let (a, b) = (r.grid(5) - 2.0, r.grid(5));
            (0..r.below(8))
                .map(|_| {
                    let x = r.grid(10);
                    c(x, a * x + b)
                })
                .collect()
        }
        // Vertical collinear.
        _ => {
            let x = r.grid(5);
            (0..r.below(8)).map(|_| c(x, r.grid(10))).collect()
        }
    }
}

fn check_hull(pts: &[Coord2d], h: &Geometry) {
    assert_test5(h, pts, &[], "convex_hull");
    let distinct: HashSet<(u64, u64)> = pts
        .iter()
        .map(|p| bits(Coord2d::new(p.x + 0.0, p.y + 0.0)))
        .collect();
    match h {
        Geometry::MultiPoint(m) => assert!(m.is_empty() && pts.is_empty()),
        Geometry::Point(p) => {
            assert_eq!(distinct.len(), 1);
            assert!(pts.iter().all(|q| *q == p.coord()));
        }
        Geometry::LineString(l) => {
            let (a, b) = (l.start(), l.end());
            assert_eq!(l.num_points(), 2);
            let mut sorted = pts.to_vec();
            sorted.sort_by(|p, q| lex(*p, *q));
            assert_eq!(a, sorted[0], "the smallest extreme first");
            assert_eq!(b, sorted[sorted.len() - 1], "the largest extreme last");
            for &p in pts {
                assert_eq!(orient(a, b, p), Sign::Zero, "all collinear");
                assert!(on_seg(p, a, b), "{p:?} within the extremes");
            }
        }
        Geometry::Polygon(poly) => {
            assert!(poly.holes().is_empty());
            let v = poly.exterior().coords().as_slice();
            let n = v.len();
            for i in 0..n {
                assert_eq!(
                    orient(v[i], v[(i + 1) % n], v[(i + 2) % n]),
                    Sign::Positive,
                    "every turn strictly left"
                );
            }
            for &p in pts {
                assert_ne!(
                    crate::exact::locate(p, poly).unwrap(),
                    Location::Exterior,
                    "{p:?} inside or on the hull"
                );
            }
        }
        other => panic!("unexpected hull {other:?}"),
    }
}

#[test]
fn hull_properties_on_generated_point_sets() {
    let mut r = Rng(0x48_554C_4C);
    let mut kinds: HashMap<&str, usize> = HashMap::new();
    for _ in 0..CASES {
        let pts = hull_points(&mut r);
        let mp = MultiPoint::new(pts.iter().map(|&p| Point::new(p).unwrap()).collect());
        let h = convex_hull(&mp).unwrap();
        check_hull(&pts, &h);
        *kinds
            .entry(match h {
                Geometry::MultiPoint(_) => "empty",
                Geometry::Point(_) => "point",
                Geometry::LineString(_) => "line",
                _ => "polygon",
            })
            .or_default() += 1;
    }
    for k in ["empty", "point", "line", "polygon"] {
        assert!(kinds.get(k).copied().unwrap_or(0) > 100, "{k}: {kinds:?}");
    }
}

#[test]
fn hull_of_every_geometry_type_uses_every_coordinate() {
    let mut r = Rng(0x48_554C_4D);
    for _ in 0..CASES {
        let g = gen_geometry(&mut r);
        let h = convex_hull(&g).unwrap();
        check_hull(&all_coords(&g), &h);
    }
}

#[test]
fn hull_hand_cases() {
    // A square with points on its edges and inside: the four corners, from
    // the smallest, counterclockwise.
    let pts = [
        (0., 0.),
        (2., 0.),
        (4., 0.),
        (4., 2.),
        (4., 4.),
        (2., 4.),
        (0., 4.),
        (0., 2.),
        (2., 2.),
    ];
    let mp = MultiPoint::new(
        pts.iter()
            .map(|&(x, y)| Point::new(c(x, y)).unwrap())
            .collect(),
    );
    let Geometry::Polygon(p) = convex_hull(&mp).unwrap() else {
        panic!()
    };
    assert_eq!(
        p.exterior().coords().as_slice(),
        &[c(0., 0.), c(4., 0.), c(4., 4.), c(0., 4.)]
    );
    // Collinear: the extremes.
    let l = convex_hull(&line(&[(3., 3.), (1., 1.), (2., 2.)])).unwrap();
    assert_eq!(l, Geometry::LineString(line(&[(1., 1.), (3., 3.)])));
    // One distinct point, given three times.
    let mp = MultiPoint::new(vec![Point::new(c(1., 2.)).unwrap(); 3]);
    assert_eq!(
        convex_hull(&mp).unwrap(),
        Geometry::Point(Point::new(c(1., 2.)).unwrap())
    );
    // Nothing.
    assert_eq!(
        convex_hull(&MultiPoint::empty()).unwrap(),
        Geometry::MultiPoint(MultiPoint::empty())
    );
}

/// Exactness at the limit: three points that differ only in their last
/// bits. A naive cross product rounds the turn to zero; `orient2d` does
/// not, so the near-collinear middle point is a hull vertex.
#[test]
fn hull_decides_a_turn_one_ulp_off_collinear() {
    let u = f64::EPSILON;
    let pts = [c(0., 0.), c(1., 1. + u), c(2., 2.)];
    let mp = MultiPoint::new(pts.iter().map(|&p| Point::new(p).unwrap()).collect());
    let Geometry::Polygon(p) = convex_hull(&mp).unwrap() else {
        panic!("a strict left turn makes a triangle")
    };
    assert_eq!(p.exterior().coords().len(), 3);
}

// ===========================================================================
// normalize
// ===========================================================================

/// `g` rewritten without changing it: rings rotated, holes and parts
/// permuted, lines reversed, points permuted.
fn scramble(r: &mut Rng, g: &Geometry) -> Geometry {
    fn rot_ring(r: &mut Rng, ring: &LinearRing) -> LinearRing {
        let v = ring.coords().as_slice();
        let k = r.below(v.len() as u64) as usize;
        LinearRing::new(CoordSeq::new(
            v[k..].iter().chain(&v[..k]).copied().collect(),
        ))
        .unwrap()
    }
    fn shuffle<T>(r: &mut Rng, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = r.below(i as u64 + 1) as usize;
            v.swap(i, j);
        }
    }
    fn poly(r: &mut Rng, p: &Polygon) -> Polygon {
        let mut holes: Vec<LinearRing> = p.holes().iter().map(|h| rot_ring(r, h)).collect();
        shuffle(r, &mut holes);
        Polygon::new(rot_ring(r, p.exterior()), holes).unwrap()
    }
    fn ln(r: &mut Rng, l: &LineString) -> LineString {
        let l = if r.coin(2) { l.reversed() } else { l.clone() };
        let v = l.coords().as_slice();
        if v.len() >= 4 && v[0] == v[v.len() - 1] {
            // A closed line written from another start.
            let body = &v[..v.len() - 1];
            let k = r.below(body.len() as u64) as usize;
            let mut w: Vec<Coord2d> = body[k..].iter().chain(&body[..k]).copied().collect();
            w.push(w[0]);
            return LineString::new(CoordSeq::new(w)).unwrap();
        }
        l
    }
    match g {
        Geometry::Point(p) => Geometry::Point(*p),
        Geometry::MultiPoint(m) => {
            let mut v = m.points().to_vec();
            shuffle(r, &mut v);
            Geometry::MultiPoint(MultiPoint::new(v))
        }
        Geometry::LineString(l) => Geometry::LineString(ln(r, l)),
        Geometry::MultiLineString(m) => {
            let mut v: Vec<LineString> = m.line_strings().iter().map(|l| ln(r, l)).collect();
            shuffle(r, &mut v);
            Geometry::MultiLineString(MultiLineString::new(v))
        }
        Geometry::Polygon(p) => Geometry::Polygon(poly(r, p)),
        Geometry::MultiPolygon(m) => {
            let mut v: Vec<Polygon> = m.polygons().iter().map(|p| poly(r, p)).collect();
            shuffle(r, &mut v);
            Geometry::MultiPolygon(MultiPolygon::new(v).unwrap())
        }
        Geometry::GeometryCollection(gc) => {
            let mut v: Vec<Geometry> = gc.geometries().iter().map(|g| scramble(r, g)).collect();
            shuffle(r, &mut v);
            Geometry::GeometryCollection(GeometryCollection::new(v).unwrap())
        }
    }
}

#[test]
fn normalize_is_idempotent_and_ignores_how_a_geometry_is_written() {
    let mut r = Rng(0x4E_4F52_4D);
    let mut changed = 0usize;
    for _ in 0..CASES {
        let g = gen_geometry(&mut r);
        let n = normalize(&g);
        assert_test5(&n, &all_coords(&g), &[], "normalize");
        let bytes = encode(&n);
        assert_eq!(encode(&normalize(&n)), bytes, "idempotent: {g:?}");
        let s = scramble(&mut r, &g);
        changed += usize::from(encode(&s) != encode(&g));
        assert_eq!(encode(&normalize(&s)), bytes, "invariant: {g:?} vs {s:?}");
        // The type's invariants hold: rebuilt through the checked
        // constructors, the value is accepted unchanged.
        if let Geometry::Polygon(p) = &n {
            let again = Polygon::new(p.exterior().clone(), p.holes().to_vec()).unwrap();
            assert_eq!(&again, p);
        }
    }
    assert!(
        changed > CASES / 2,
        "the scrambles change the writing: {changed}"
    );
}

#[test]
fn normalize_hand_cases() {
    // A ring starts at its smallest vertex and stays counterclockwise.
    let p = Polygon::new_simple(ring(&[(2., 0.), (2., 2.), (0., 2.), (0., 0.)])).unwrap();
    let Geometry::Polygon(n) = normalize(&Geometry::Polygon(p)) else {
        panic!()
    };
    assert_eq!(
        n.exterior().coords().as_slice(),
        &[c(0., 0.), c(2., 0.), c(2., 2.), c(0., 2.)]
    );
    // A line is reversed if it compares greater than its reverse.
    let l = normalize(&Geometry::LineString(line(&[(3., 0.), (1., 0.)])));
    assert_eq!(l, Geometry::LineString(line(&[(1., 0.), (3., 0.)])));
    // A closed line starts at its smallest vertex, toward the smaller
    // neighbor.
    let l = normalize(&Geometry::LineString(line(&[
        (1., 1.),
        (1., 0.),
        (0., 0.),
        (1., 1.),
    ])));
    assert_eq!(
        l,
        Geometry::LineString(line(&[(0., 0.), (1., 0.), (1., 1.), (0., 0.)]))
    );
    // A collection: by type (codec tag order: point, line, polygon, ...),
    // then by coordinates.
    let gc = GeometryCollection::new(vec![
        Geometry::Polygon(rect(0., 0., 1., 1.)),
        Geometry::Point(Point::new(c(5., 5.)).unwrap()),
        Geometry::LineString(line(&[(2., 2.), (3., 3.)])),
        Geometry::Point(Point::new(c(4., 4.)).unwrap()),
    ])
    .unwrap();
    let Geometry::GeometryCollection(n) = normalize(&Geometry::GeometryCollection(gc)) else {
        panic!()
    };
    let kinds: Vec<&str> = n
        .geometries()
        .iter()
        .map(|g| match g {
            Geometry::Point(p) if p.coord() == c(4., 4.) => "p44",
            Geometry::Point(_) => "p55",
            Geometry::LineString(_) => "line",
            _ => "poly",
        })
        .collect();
    assert_eq!(kinds, ["p44", "p55", "line", "poly"]);
}

// ===========================================================================
// boundary
// ===========================================================================

fn dim_of(g: &Geometry) -> Dimension {
    match g {
        Geometry::MultiPoint(m) if m.is_empty() => Dimension::Empty,
        Geometry::MultiPoint(_) => Dimension::Point,
        Geometry::MultiLineString(m) if m.is_empty() => Dimension::Empty,
        Geometry::MultiLineString(_) => Dimension::Line,
        Geometry::GeometryCollection(_) => Dimension::Line,
        other => panic!("unexpected boundary {other:?}"),
    }
}

#[test]
fn boundary_agrees_with_relate_on_generated_geometries() {
    let mut r = Rng(0x42_4E44_52);
    let mut kinds: HashMap<i8, usize> = HashMap::new();
    for _ in 0..CASES {
        let g = gen_geometry(&mut r);
        let b = boundary(&g).unwrap();
        assert_test5(&b, &all_coords(&g), &[], "boundary");
        let m = IntersectionMatrix::compute_geometry(&g, &g)
            .unwrap()
            .cells();
        // B(g) ∩ B(g) is the boundary itself.
        assert_eq!(dim_of(&b), m[4], "BB of relate(g, g) for {g:?}");
        *kinds.entry(dim_of(&b) as i8).or_default() += 1;
        // Every boundary point is located on the boundary by relate.
        for p in all_coords(&b) {
            let pt = Geometry::Point(Point::new(p).unwrap());
            let row = IntersectionMatrix::compute_geometry(&pt, &g)
                .unwrap()
                .cells();
            assert_eq!(row[1], Dimension::Point, "{p:?} in B({g:?})");
        }
    }
    for d in [Dimension::Empty, Dimension::Point, Dimension::Line] {
        assert!(
            kinds.get(&(d as i8)).copied().unwrap_or(0) > 300,
            "{kinds:?}"
        );
    }
}

#[test]
fn boundary_hand_cases() {
    // An open line: its endpoints, smallest first. A closed one: nothing.
    assert_eq!(
        boundary(&Geometry::LineString(line(&[(3., 0.), (0., 0.)]))).unwrap(),
        Geometry::MultiPoint(MultiPoint::new(vec![
            Point::new(c(0., 0.)).unwrap(),
            Point::new(c(3., 0.)).unwrap()
        ]))
    );
    assert_eq!(
        boundary(&Geometry::LineString(line(&[
            (0., 0.),
            (1., 0.),
            (0., 1.),
            (0., 0.)
        ])))
        .unwrap(),
        Geometry::MultiPoint(MultiPoint::empty())
    );
    // Mod 2: two lines meeting end to end share an interior point; three
    // lines meeting at one end make it a boundary point again.
    let two = MultiLineString::new(vec![
        line(&[(0., 0.), (1., 0.)]),
        line(&[(1., 0.), (2., 0.)]),
    ]);
    assert_eq!(
        boundary(&mls(two)).unwrap(),
        Geometry::MultiPoint(MultiPoint::new(vec![
            Point::new(c(0., 0.)).unwrap(),
            Point::new(c(2., 0.)).unwrap()
        ]))
    );
    let three = MultiLineString::new(vec![
        line(&[(0., 0.), (1., 0.)]),
        line(&[(1., 0.), (2., 0.)]),
        line(&[(1., 0.), (1., 1.)]),
    ]);
    let Geometry::MultiPoint(b) = boundary(&mls(three)).unwrap() else {
        panic!()
    };
    assert_eq!(b.num_points(), 4);
    // A polygon with a hole: two closed rings, orientation kept.
    let holed = Polygon::new(
        ring(&[(0., 0.), (4., 0.), (4., 4.), (0., 4.)]),
        vec![ring(&[(1., 1.), (1., 2.), (2., 2.), (2., 1.)])],
    )
    .unwrap();
    assert_eq!(
        boundary(&Geometry::Polygon(holed)).unwrap(),
        mls(MultiLineString::new(vec![
            line(&[(0., 0.), (4., 0.), (4., 4.), (0., 4.), (0., 0.)]),
            line(&[(1., 1.), (1., 2.), (2., 2.), (2., 1.), (1., 1.)]),
        ]))
    );
    // A collection: the square's ring, and the line's endpoint outside it;
    // the endpoint inside the square is absorbed.
    let gc = GeometryCollection::new(vec![
        Geometry::Polygon(rect(0., 0., 2., 2.)),
        Geometry::LineString(line(&[(1., 1.), (5., 1.)])),
        Geometry::Point(Point::new(c(9., 9.)).unwrap()),
    ])
    .unwrap();
    assert_eq!(
        boundary(&Geometry::GeometryCollection(gc)).unwrap(),
        Geometry::GeometryCollection(
            GeometryCollection::new(vec![
                mls(MultiLineString::new(vec![line(&[
                    (0., 0.),
                    (2., 0.),
                    (2., 2.),
                    (0., 2.),
                    (0., 0.)
                ])])),
                Geometry::MultiPoint(MultiPoint::new(vec![Point::new(c(5., 1.)).unwrap()])),
            ])
            .unwrap()
        )
    );
    // Overlapping polygonal leaves are refused, as relate refuses them.
    let bad = GeometryCollection::new(vec![
        Geometry::Polygon(rect(0., 0., 2., 2.)),
        Geometry::Polygon(rect(1., 1., 3., 3.)),
    ])
    .unwrap();
    assert!(matches!(
        boundary(&Geometry::GeometryCollection(bad)),
        Err(GeometryError::CollectionElementsOverlap { .. })
    ));
}

// ===========================================================================
// reverse, unique_points
// ===========================================================================

#[test]
fn reverse_and_unique_points() {
    let l = line(&[(0., 0.), (1., 0.), (1., 1.)]);
    assert_eq!(reverse(&l), line(&[(1., 1.), (1., 0.), (0., 0.)]));
    assert_eq!(reverse(&reverse(&l)), l);
    let m = MultiLineString::new(vec![l.clone(), line(&[(5., 5.), (6., 6.)])]);
    assert_eq!(
        reverse(&m),
        MultiLineString::new(vec![reverse(&l), line(&[(6., 6.), (5., 5.)])])
    );
    let p = Point::new(c(1., 2.)).unwrap();
    assert_eq!(reverse(&p), p);

    // First occurrence, and -0.0 counts as 0.0 with the first bits kept.
    let mp = MultiPoint::new(vec![
        Point::new(c(-0.0, 1.)).unwrap(),
        Point::new(c(2., 2.)).unwrap(),
        Point::new(c(0.0, 1.)).unwrap(),
        Point::new(c(2., 2.)).unwrap(),
    ]);
    let u = unique_points(&mp);
    assert_eq!(u.num_points(), 2);
    assert!(u.points()[0].x().is_sign_negative());
    let mut r = Rng(0x55_4E51);
    for _ in 0..CASES {
        let g = gen_geometry(&mut r);
        let u = unique_points(&g);
        let all = all_coords(&g);
        assert_test5(&Geometry::MultiPoint(u.clone()), &all, &[], "unique_points");
        let distinct: HashSet<(u64, u64)> = all.iter().map(|&p| super::basic::key(p)).collect();
        assert_eq!(u.num_points(), distinct.len());
    }
}

// ===========================================================================
// line_merge
// ===========================================================================

/// An undirected segment, its ends in lexicographic order.
fn useg(a: Coord2d, b: Coord2d) -> ((u64, u64), (u64, u64)) {
    let (p, q) = if lex(a, b) == std::cmp::Ordering::Greater {
        (b, a)
    } else {
        (a, b)
    };
    (super::basic::key(p), super::basic::key(q))
}

fn seg_multiset(lines: &[LineString], directed: bool) -> HashMap<((u64, u64), (u64, u64)), usize> {
    let mut m = HashMap::new();
    for l in lines {
        for w in l.coords().as_slice().windows(2) {
            let k = if directed {
                (super::basic::key(w[0]), super::basic::key(w[1]))
            } else {
                useg(w[0], w[1])
            };
            *m.entry(k).or_insert(0) += 1;
        }
    }
    m
}

fn endpoint_degrees(lines: &[LineString]) -> HashMap<(u64, u64), usize> {
    let mut d = HashMap::new();
    for l in lines {
        *d.entry(super::basic::key(l.start())).or_insert(0) += 1;
        *d.entry(super::basic::key(l.end())).or_insert(0) += 1;
    }
    d
}

#[test]
fn line_merge_properties_on_generated_line_work() {
    let mut r = Rng(0x4C_4D52_47);
    let (mut joined, mut loops) = (0usize, 0usize);
    for _ in 0..CASES {
        let input = gen_lines(&mut r, 4, 7);
        let lines = input.line_strings();
        let out = line_merge(&input);
        let all: Vec<Coord2d> = all_coords(&mls(input.clone()));
        assert_test5(&mls(out.clone()), &all, &[], "line_merge");
        assert_eq!(
            seg_multiset(out.line_strings(), false),
            seg_multiset(lines, false),
            "every input segment once: {input:?}"
        );
        let deg = endpoint_degrees(lines);
        for l in out.line_strings() {
            if l.start() != l.end() {
                for e in [l.start(), l.end()] {
                    assert_ne!(
                        deg[&super::basic::key(e)],
                        2,
                        "joinable further at {e:?}: {input:?}"
                    );
                }
            } else {
                loops += 1;
            }
        }
        joined += usize::from(out.num_line_strings() < input.num_line_strings());

        // Directed: the same coverage, every segment in its input direction.
        let d = line_merge_directed(&input);
        assert_eq!(
            seg_multiset(d.line_strings(), true),
            seg_multiset(lines, true)
        );
        assert_test5(&mls(d), &all, &[], "line_merge_directed");
    }
    assert!(
        joined > CASES / 10 && loops > 100,
        "joined {joined}, loops {loops}"
    );
}

#[test]
fn line_merge_hand_cases() {
    // Three lines in a row, the middle one backward: joined into one, in the
    // majority direction (two of three run left to right).
    let input = MultiLineString::new(vec![
        line(&[(0., 0.), (1., 0.)]),
        line(&[(2., 0.), (1., 0.)]),
        line(&[(2., 0.), (3., 0.)]),
    ]);
    assert_eq!(
        line_merge(&input),
        MultiLineString::new(vec![line(&[(0., 0.), (1., 0.), (2., 0.), (3., 0.)])])
    );
    // Directed merging stops where the directions disagree.
    assert_eq!(line_merge_directed(&input).num_line_strings(), 3);
    // A tie (one forward, one backward) starts at the smaller end.
    let tie = MultiLineString::new(vec![
        line(&[(3., 0.), (2., 0.)]),
        line(&[(1., 0.), (2., 0.)]),
    ]);
    assert_eq!(
        line_merge(&tie),
        MultiLineString::new(vec![line(&[(1., 0.), (2., 0.), (3., 0.)])])
    );
    // A degree-3 node ends every line at it.
    let star = MultiLineString::new(vec![
        line(&[(0., 0.), (1., 1.)]),
        line(&[(1., 1.), (2., 0.)]),
        line(&[(1., 1.), (1., 2.)]),
    ]);
    assert_eq!(line_merge(&star).num_line_strings(), 3);
    // A loop of two lines: one closed line from its smallest node.
    let lp = MultiLineString::new(vec![
        line(&[(1., 1.), (0., 1.), (0., 0.)]),
        line(&[(0., 0.), (1., 0.), (1., 1.)]),
    ]);
    let m = line_merge(&lp);
    assert_eq!(m.num_line_strings(), 1);
    assert_eq!(m.line_strings()[0].start(), c(0., 0.));
    assert_eq!(m.line_strings()[0].end(), c(0., 0.));
    // Polygons and points in a collection are ignored.
    let gc = GeometryCollection::new(vec![
        Geometry::Polygon(rect(0., 0., 1., 1.)),
        Geometry::LineString(line(&[(5., 5.), (6., 6.)])),
        Geometry::Point(Point::new(c(9., 9.)).unwrap()),
    ])
    .unwrap();
    assert_eq!(
        line_merge(&gc),
        MultiLineString::new(vec![line(&[(5., 5.), (6., 6.)])])
    );
}

// ===========================================================================
// line_sequence
// ===========================================================================

fn components_and_odd(lines: &[LineString]) -> (usize, bool) {
    let mut id: HashMap<(u64, u64), usize> = HashMap::new();
    let mut parent: Vec<usize> = Vec::new();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut node = |c: Coord2d, parent: &mut Vec<usize>| {
        *id.entry(super::basic::key(c)).or_insert_with(|| {
            parent.push(parent.len());
            parent.len() - 1
        })
    };
    let mut deg: HashMap<usize, usize> = HashMap::new();
    for l in lines {
        let (a, b) = (node(l.start(), &mut parent), node(l.end(), &mut parent));
        *deg.entry(a).or_default() += 1;
        *deg.entry(b).or_default() += 1;
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra] = rb;
    }
    let n = parent.len();
    let roots: HashSet<usize> = (0..n).map(|v| find(&mut parent, v)).collect();
    let mut odd: HashMap<usize, usize> = HashMap::new();
    for (&v, &d) in &deg {
        if d % 2 == 1 {
            *odd.entry(find(&mut parent, v)).or_default() += 1;
        }
    }
    (roots.len(), odd.values().all(|&k| k <= 2))
}

#[test]
fn line_sequence_properties_on_generated_line_work() {
    let mut r = Rng(0x53_4551);
    let (mut ok, mut refused) = (0usize, 0usize);
    for _ in 0..CASES {
        let input = gen_lines(&mut r, 4, 7);
        let lines = input.line_strings();
        let (ncomp, sequenceable) = components_and_odd(lines);
        match line_sequence(&input) {
            Ok(out) => {
                assert!(sequenceable, "sequenced an unsequenceable input: {input:?}");
                ok += 1;
                assert_eq!(out.num_line_strings(), lines.len());
                assert_eq!(
                    seg_multiset(out.line_strings(), false),
                    seg_multiset(lines, false)
                );
                // Each output line is an input line, possibly reversed.
                let mut pool: Vec<Vec<Coord2d>> = lines
                    .iter()
                    .map(|l| l.coords().as_slice().to_vec())
                    .collect();
                for l in out.line_strings() {
                    let v = l.coords().as_slice().to_vec();
                    let rv: Vec<Coord2d> = v.iter().rev().copied().collect();
                    let k = pool.iter().position(|p| *p == v || *p == rv).unwrap();
                    pool.remove(k);
                }
                // Runs of lines that start where the previous one ended:
                // exactly one per component.
                let ls = out.line_strings();
                let runs = 1
                    + (1..ls.len())
                        .filter(|&i| ls[i].start() != ls[i - 1].end())
                        .count();
                assert_eq!(runs, ncomp, "{input:?} -> {out:?}");
                assert_test5(
                    &mls(out),
                    &all_coords(&mls(input.clone())),
                    &[],
                    "line_sequence",
                );
            }
            Err(GeometryError::NotSequenceable) => {
                assert!(!sequenceable, "refused a sequenceable input: {input:?}");
                refused += 1;
            }
            Err(e) => panic!("{e:?}"),
        }
    }
    assert!(ok > 1000 && refused > 300, "ok {ok}, refused {refused}");
}

#[test]
fn line_sequence_hand_cases() {
    // Out of order and one reversed: sequenced from the degree-1 node where
    // a line starts.
    let input = MultiLineString::new(vec![
        line(&[(2., 0.), (3., 0.)]),
        line(&[(2., 0.), (1., 0.)]),
        line(&[(0., 0.), (1., 0.)]),
    ]);
    assert_eq!(
        line_sequence(&input).unwrap(),
        MultiLineString::new(vec![
            line(&[(0., 0.), (1., 0.)]),
            line(&[(1., 0.), (2., 0.)]),
            line(&[(2., 0.), (3., 0.)]),
        ])
    );
    // A star of three lines has three odd nodes plus the centre... four odd
    // nodes: not sequenceable.
    let star = MultiLineString::new(vec![
        line(&[(0., 0.), (1., 1.)]),
        line(&[(2., 0.), (1., 1.)]),
        line(&[(1., 2.), (1., 1.)]),
    ]);
    assert_eq!(line_sequence(&star), Err(GeometryError::NotSequenceable));
}

// ===========================================================================
// shared_paths
// ===========================================================================

/// The brute-force shared pieces: every segment of `a` split at every input
/// vertex strictly inside it, each piece classed by the `b` segments that
/// cover it.
fn brute_shared(
    a: &MultiLineString,
    b: &MultiLineString,
) -> [HashMap<((u64, u64), (u64, u64)), usize>; 2] {
    let verts: Vec<Coord2d> = all_coords(&mls(a.clone()))
        .into_iter()
        .chain(all_coords(&mls(b.clone())))
        .collect();
    let mut out = [HashMap::new(), HashMap::new()];
    for l in a.line_strings() {
        for w in l.coords().as_slice().windows(2) {
            let (p, q) = (w[0], w[1]);
            let k: fn(Coord2d) -> f64 = if p.x == q.x { |c| c.y } else { |c| c.x };
            let up = k(q) > k(p);
            let mut cuts: Vec<Coord2d> = verts
                .iter()
                .copied()
                .filter(|&v| on_seg(v, p, q) && v != p && v != q)
                .collect();
            cuts.push(p);
            cuts.push(q);
            cuts.sort_by(|x, y| k(*x).total_cmp(&k(*y)));
            cuts.dedup_by(|x, y| k(*x) == k(*y));
            if !up {
                cuts.reverse();
            }
            for piece in cuts.windows(2) {
                let (s, e) = (piece[0], piece[1]);
                for class in 0..2 {
                    let covered = b.line_strings().iter().any(|bl| {
                        bl.coords().as_slice().windows(2).any(|bw| {
                            let same = (k(bw[1]) > k(bw[0])) == up;
                            usize::from(!same) == class
                                && orient(p, q, bw[0]) == Sign::Zero
                                && orient(p, q, bw[1]) == Sign::Zero
                                && on_seg(s, bw[0], bw[1])
                                && on_seg(e, bw[0], bw[1])
                        })
                    });
                    if covered {
                        *out[class]
                            .entry((super::basic::key(s), super::basic::key(e)))
                            .or_insert(0) += 1;
                    }
                }
            }
        }
    }
    out
}

/// The output's pieces, split the same way.
fn split_output(
    m: &MultiLineString,
    verts: &[Coord2d],
) -> HashMap<((u64, u64), (u64, u64)), usize> {
    let mut out = HashMap::new();
    for l in m.line_strings() {
        for w in l.coords().as_slice().windows(2) {
            let (p, q) = (w[0], w[1]);
            let k: fn(Coord2d) -> f64 = if p.x == q.x { |c| c.y } else { |c| c.x };
            let up = k(q) > k(p);
            let mut cuts: Vec<Coord2d> = verts
                .iter()
                .copied()
                .filter(|&v| on_seg(v, p, q) && v != p && v != q)
                .collect();
            cuts.push(p);
            cuts.push(q);
            cuts.sort_by(|x, y| k(*x).total_cmp(&k(*y)));
            cuts.dedup_by(|x, y| k(*x) == k(*y));
            if !up {
                cuts.reverse();
            }
            for piece in cuts.windows(2) {
                *out.entry((super::basic::key(piece[0]), super::basic::key(piece[1])))
                    .or_insert(0) += 1;
            }
        }
    }
    out
}

#[test]
fn shared_paths_equal_brute_force_on_generated_line_work() {
    let mut r = Rng(0x53_4850);
    let mut found = [0usize; 2];
    for _ in 0..CASES {
        let a = gen_lines(&mut r, 4, 3);
        let b = gen_lines(&mut r, 4, 3);
        let sp = shared_paths(&a, &b);
        let verts: Vec<Coord2d> = all_coords(&mls(a.clone()))
            .into_iter()
            .chain(all_coords(&mls(b.clone())))
            .collect();
        let brute = brute_shared(&a, &b);
        for (class, got) in [&sp.same_direction, &sp.opposite_direction]
            .iter()
            .enumerate()
        {
            assert_test5(&mls((*got).clone()), &verts, &[], "shared_paths");
            assert_eq!(
                split_output(got, &verts),
                brute[class],
                "class {class}: {a:?} / {b:?}"
            );
            found[class] += usize::from(!got.is_empty());
        }
    }
    assert!(found[0] > 500 && found[1] > 500, "{found:?}");
}

#[test]
fn shared_paths_hand_cases() {
    // A vertical shared stretch, and pieces joined across a vertex of a.
    let a = line(&[(0., 0.), (0., 2.), (0., 4.)]);
    let b = line(&[(0., 1.), (0., 3.)]);
    let sp = shared_paths(&a, &b);
    assert_eq!(
        sp.same_direction,
        MultiLineString::new(vec![line(&[(0., 1.), (0., 2.), (0., 3.)])])
    );
    assert!(sp.opposite_direction.is_empty());
    // Opposite directions; a single shared point is not a path.
    let a = line(&[(0., 0.), (4., 0.)]);
    let b = MultiLineString::new(vec![
        line(&[(3., 0.), (1., 0.)]),
        line(&[(4., 0.), (5., 5.)]),
    ]);
    let sp = shared_paths(&a, &b);
    assert_eq!(
        sp.opposite_direction,
        MultiLineString::new(vec![line(&[(1., 0.), (3., 0.)])])
    );
    assert!(sp.same_direction.is_empty());
}

// ===========================================================================
// split_by_point
// ===========================================================================

#[test]
fn split_by_point_properties_on_generated_lines() {
    let mut r = Rng(0x53_504C);
    let (mut split, mut refused) = (0usize, 0usize);
    for _ in 0..CASES {
        let Some(l) = gen_line(&mut r, 6) else {
            continue;
        };
        let v = l.coords().as_slice();
        let p = match r.below(4) {
            0 => v[1 + r.below(v.len() as u64 - 1) as usize - usize::from(v.len() == 2)],
            1 => {
                let i = r.below(v.len() as u64 - 1) as usize;
                c((v[i].x + v[i + 1].x) * 0.5, (v[i].y + v[i + 1].y) * 0.5)
            }
            2 => v[0],
            _ => r.pt(6),
        };
        match split_by_point(&l, p) {
            Ok(parts) => {
                split += 1;
                assert!(v[0] != p && v[v.len() - 1] != p);
                let ps = parts.line_strings();
                // The parts meet at p, and concatenated give the line with
                // p written at each passage.
                let mut joined: Vec<Coord2d> = ps[0].coords().as_slice().to_vec();
                for w in ps.windows(2) {
                    assert_eq!(w[0].end(), p);
                    assert_eq!(w[1].start(), p);
                    joined.extend_from_slice(&w[1].coords().as_slice()[1..]);
                }
                let expected: Vec<Coord2d> = {
                    let mut e = vec![v[0]];
                    for w in v.windows(2) {
                        if w[0] != p && w[1] != p && on_seg(p, w[0], w[1]) {
                            e.push(p);
                        }
                        e.push(w[1]);
                    }
                    e
                };
                assert_eq!(joined, expected);
                assert_test5(&mls(parts), v, &[p], "split_by_point");
            }
            Err(GeometryError::PointNotOnLineInterior) => {
                refused += 1;
                let on = v.windows(2).any(|w| on_seg(p, w[0], w[1]));
                assert!(!on || v[0] == p || v[v.len() - 1] == p, "{l:?} at {p:?}");
            }
            Err(e) => panic!("{e:?}"),
        }
    }
    assert!(
        split > 3000 && refused > 2000,
        "split {split}, refused {refused}"
    );
}

#[test]
fn split_by_point_hand_cases() {
    let l = line(&[(0., 0.), (4., 0.), (4., 4.)]);
    assert_eq!(
        split_by_point(&l, c(2., 0.)).unwrap(),
        MultiLineString::new(vec![
            line(&[(0., 0.), (2., 0.)]),
            line(&[(2., 0.), (4., 0.), (4., 4.)])
        ])
    );
    assert_eq!(
        split_by_point(&l, c(4., 0.)).unwrap(),
        MultiLineString::new(vec![
            line(&[(0., 0.), (4., 0.)]),
            line(&[(4., 0.), (4., 4.)])
        ])
    );
    assert_eq!(
        split_by_point(&l, c(0., 0.)),
        Err(GeometryError::PointNotOnLineInterior)
    );
    // Near, not on: not snapped.
    assert_eq!(
        split_by_point(&l, c(2., f64::EPSILON)),
        Err(GeometryError::PointNotOnLineInterior)
    );
    // A line through p twice is split at both passages.
    let twice = line(&[(0., 0.), (2., 2.), (2., 0.), (0., 2.)]);
    assert_eq!(
        split_by_point(&twice, c(1., 1.))
            .unwrap()
            .num_line_strings(),
        3
    );
    // A closed line at its closing vertex: an endpoint.
    let closed = line(&[(0., 0.), (1., 0.), (0., 1.), (0., 0.)]);
    assert_eq!(
        split_by_point(&closed, c(0., 0.)),
        Err(GeometryError::PointNotOnLineInterior)
    );
    assert!(matches!(
        split_by_point(&l, c(f64::NAN, 0.)),
        Err(GeometryError::NonFiniteCoordinate { .. })
    ));
}

/// The order helper the operations share orders exactly, `-0.0` with `0.0`.
#[test]
fn the_lexicographic_order_is_exact() {
    assert_eq!(lex(c(-0.0, 1.), c(0.0, 1.)), std::cmp::Ordering::Equal);
    assert_eq!(
        lex(c(1., 1.), c(1., 1. + f64::EPSILON)),
        std::cmp::Ordering::Less
    );
    assert_eq!(
        lex_seq(&[c(0., 0.)], &[c(0., 0.), c(1., 1.)]),
        std::cmp::Ordering::Less
    );
}
