// src/predicates/fast_path_tests.rs

//! Tests for the paths that decide a predicate without the general
//! computation: the early-exit `intersects` of two areal operands, and the
//! windowed form of a large unindexed operand. Each path must give exactly
//! the answer of the path it short-cuts.

#![allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation
)]

use super::*;
use crate::ops::is_valid;
use crate::prepared::PreparedPolygonal;
use crate::types::{Coord2d, CoordSeq, LineString, LinearRing, MultiLineString, Point};

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new_simple(ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])).unwrap()
}

/// SplitMix64.
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
}

/// A star-shaped polygon on the integer grid around `(cx, cy)`, with an
/// optional square hole at its center; valid, or `None`.
fn star(r: &mut Rng, cx: f64, cy: f64, n: usize, rmax: u64, hole: bool) -> Option<Polygon> {
    let mut pts = Vec::new();
    for k in 0..n {
        let t = std::f64::consts::TAU * (k as f64 + 0.3 * (r.below(100) as f64 / 100.0)) / n as f64;
        let rad = (3 + r.below(rmax)) as f64;
        pts.push(c(
            (cx + rad * t.cos()).round(),
            (cy + rad * t.sin()).round(),
        ));
    }
    pts.dedup();
    let ext = LinearRing::new(CoordSeq::new(pts)).ok()?;
    let holes = if hole {
        vec![ring(&[
            (cx - 1.0, cy - 1.0),
            (cx - 1.0, cy + 1.0),
            (cx + 1.0, cy + 1.0),
            (cx + 1.0, cy - 1.0),
        ])]
    } else {
        Vec::new()
    };
    let p = Polygon::new(ext, holes).ok()?;
    is_valid(&p).then_some(p)
}

/// `intersects` and `disjoint` agree with the matrix on generated pairs
/// that touch at vertices, share edges, cross, nest and miss.
#[test]
fn intersects_agrees_with_relate_on_generated_areas() {
    let mut r = Rng(0x5EED_8800);
    let mut checked = [0usize; 2];
    while checked[0] + checked[1] < 4000 {
        let (na, ha) = (3 + r.below(10) as usize, r.below(3) == 0);
        let a = star(&mut r, 0.0, 0.0, na, 8, ha);
        let (dx, dy) = (r.below(25) as f64 - 12.0, r.below(25) as f64 - 12.0);
        let (nb, hb) = (3 + r.below(10) as usize, r.below(3) == 0);
        let b = star(&mut r, dx, dy, nb, 8, hb);
        let (Some(a), Some(b)) = (a, b) else { continue };
        let want = relate(&a, &b).unwrap().intersects();
        assert_eq!(intersects(&a, &b), Ok(want), "{a:?} {b:?}");
        assert_eq!(intersects(&b, &a), Ok(want), "{b:?} {a:?}");
        assert_eq!(disjoint(&a, &b), Ok(!want));
        checked[usize::from(want)] += 1;
    }
    assert!(checked[0] > 400 && checked[1] > 400, "{checked:?}");
}

/// Two bars forming a plus: no vertex of either lies in the other, and
/// only the edge-crossing search finds the common region.
#[test]
fn intersects_finds_a_crossing_with_no_vertex_inside() {
    let h = rect(0.0, 4.0, 10.0, 6.0);
    let v = rect(4.0, 0.0, 6.0, 10.0);
    assert!(relate(&h, &v).unwrap().intersects());
    assert_eq!(intersects(&h, &v), Ok(true));
    assert_eq!(intersects(&v, &h), Ok(true));
}

/// Containment, a hole, a touch along a hole's edge and a touch at a
/// corner, each against the matrix.
#[test]
fn intersects_sees_nesting_holes_and_touches() {
    let holed = Polygon::new(
        ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        vec![ring(&[(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0)])],
    )
    .unwrap();
    let cases = [
        (rect(1.0, 1.0, 2.0, 2.0), true),     // inside the shell
        (rect(4.0, 4.0, 6.0, 6.0), false),    // inside the hole
        (rect(4.0, 4.0, 7.0, 6.0), true),     // in the hole, touching its edge
        (rect(-5.0, -5.0, 20.0, 20.0), true), // covering everything
        (rect(10.0, 10.0, 12.0, 12.0), true), // a corner touch
        (rect(11.0, 0.0, 12.0, 1.0), false),  // apart
    ];
    for (p, want) in cases {
        assert_eq!(relate(&holed, &p).unwrap().intersects(), want, "{p:?}");
        assert_eq!(intersects(&holed, &p), Ok(want), "{p:?}");
        assert_eq!(intersects(&p, &holed), Ok(want), "{p:?}");
    }
}

/// Squares sharing an edge, and squares meeting at one corner: the edges
/// that meet lie on the boundary of each other's box, and must be found.
#[test]
fn edges_on_the_boundary_of_a_box_are_found() {
    let a = rect(0.0, 0.0, 2.0, 2.0);
    let side = rect(2.0, 0.0, 4.0, 2.0);
    let corner = rect(2.0, 2.0, 4.0, 4.0);
    assert_eq!(relate(&a, &side).unwrap().matrix().to_string(), "FF2F11212");
    assert_eq!(
        relate(&a, &corner).unwrap().matrix().to_string(),
        "FF2F01212"
    );
    assert!(touches(&a, &side).unwrap() && touches(&a, &corner).unwrap());
}

/// A polygon of about 600 edges, with a hole, far beyond the windowing
/// threshold.
fn large_holed() -> Polygon {
    let n = 600;
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * f64::from(k) / f64::from(n);
            let rad = if k % 2 == 0 { 100.0 } else { 97.0 };
            ((rad * t.cos()).round(), (rad * t.sin()).round())
        })
        .collect();
    let mut pts = pts;
    pts.dedup();
    let hole = ring(&[(-20.0, -20.0), (-20.0, 20.0), (20.0, 20.0), (20.0, -20.0)]);
    let p = Polygon::new(ring(&pts), vec![hole]).unwrap();
    assert!(is_valid(&p));
    p
}

/// Small operands placed against a large one, where `relate`,
/// `intersects` and the prepared mask window the large operand: every
/// matrix equals the one computed over the plain unindexed operand with no
/// window (`polygons_matrix` and `lines_polygon_matrix` on `Areal`), in
/// both orders, for polygons and lines, prepared and not.
///
/// The placements put the small operand inside the shell (so `locate`'s
/// ray must count edges far to the right of its window), in the hole,
/// across the outer boundary, across the hole's edge, and outside.
#[test]
fn a_windowed_operand_relates_as_the_plain_one() {
    let big = large_holed();
    let prepared = PreparedPolygonal::new(&big).unwrap();
    let mut n = 0;
    for x in (-110..=110).step_by(11) {
        for y in (-110..=110).step_by(13) {
            let (x, y) = (f64::from(x), f64::from(y));
            let small = rect(x, y, x + 6.0, y + 4.0);
            let want = plain_polygons(&big, &small);
            assert_eq!(
                prepared.relate(&Geometry::Polygon(small.clone())).unwrap(),
                want,
                "{x} {y}"
            );
            assert_eq!(relate(&big, &small).unwrap().matrix(), want, "{x} {y}");
            assert_eq!(
                relate(&small, &big).unwrap().matrix(),
                want.swap_ab(),
                "{x} {y}"
            );
            assert_eq!(intersects(&small, &big), Ok(want.intersects()), "{x} {y}");
            let line = LineString::new(CoordSeq::new(vec![c(x, y), c(x + 9.0, y + 5.0)])).unwrap();
            let want = plain_lines(&line, &big).swap_ab();
            assert_eq!(
                prepared
                    .relate(&Geometry::LineString(line.clone()))
                    .unwrap(),
                want,
                "line {x} {y}"
            );
            assert_eq!(relate(&big, &line).unwrap().matrix(), want, "line {x} {y}");
            n += 1;
        }
    }
    assert!(n > 300);
}

/// A point against lines, polygons and prepared polygons: `intersects`,
/// and the prepared `intersects`, `covers` and `contains`, equal the
/// matrix's on points at vertices, on edges, inside, in a hole and outside,
/// in both orders.
#[test]
fn a_point_meets_as_the_matrix_says() {
    let holed = Polygon::new(
        ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        vec![ring(&[(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0)])],
    )
    .unwrap();
    let line =
        LineString::new(CoordSeq::new(vec![c(-2.0, 5.0), c(4.0, 5.0), c(4.0, 12.0)])).unwrap();
    let multi = MultiLineString::new(vec![
        line.clone(),
        LineString::new(CoordSeq::new(vec![c(8.0, -1.0), c(8.0, 11.0)])).unwrap(),
    ]);
    let prepared = PreparedPolygonal::new(&holed).unwrap();
    let mut n = 0;
    for x in -4..=24 {
        for y in -4..=24 {
            let p = Point::new(c(f64::from(x) * 0.5, f64::from(y) * 0.5)).unwrap();
            let gp = Geometry::Point(p);
            let m = relate(&p, &holed).unwrap();
            assert_eq!(intersects(&p, &holed), Ok(m.intersects()), "{p:?}");
            assert_eq!(intersects(&holed, &p), Ok(m.intersects()), "{p:?}");
            assert_eq!(prepared.intersects(&gp), Ok(m.intersects()), "{p:?}");
            assert_eq!(prepared.covers(&gp), Ok(m.covered_by()), "{p:?}");
            assert_eq!(prepared.contains(&gp), Ok(m.within()), "{p:?}");
            for l in [
                &Geometry::LineString(line.clone()),
                &Geometry::MultiLineString(multi.clone()),
            ] {
                let ml = relate(&p, l).unwrap().intersects();
                assert_eq!(intersects(&p, l), Ok(ml), "{p:?}");
                assert_eq!(intersects(l, &p), Ok(ml), "{p:?}");
            }
            n += 1;
        }
    }
    assert!(n > 800);
}

/// The matrix of two polygons over the plain unindexed operands, with no
/// window: the reference the windowed paths must reproduce.
fn plain_polygons(a: &Polygon, b: &Polygon) -> IntersectionMatrix {
    use crate::topology::areal::Areal;
    crate::topology::polygon_analysis::polygons_matrix(&Areal::single(a), &Areal::single(b))
}

/// The matrix of a line against a polygon over the plain unindexed operand.
fn plain_lines(l: &LineString, p: &Polygon) -> IntersectionMatrix {
    use crate::topology::areal::Areal;
    use crate::topology::line_analysis::{LineSet, lines_polygon_matrix};
    lines_polygon_matrix(&LineSet::single(l), &Areal::single(p))
}
