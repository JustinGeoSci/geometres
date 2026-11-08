// src/ops/hausdorff_tests.rs

//! Tests for `ops::hausdorff`: the predicate equals an exact `i128` brute
//! force, the boundary cases answer exactly, and the indexed paths return
//! exactly what the exhaustive ones do.
//!
//! The brute force here shares nothing with the module: it extracts
//! vertices and linework from the public types itself and decides every
//! point-to-element comparison in integer arithmetic (an interior
//! projection compares `cross² ≤ D·|u|²`, so no division is needed).

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

use super::distance::Element;
use super::hausdorff::{directed_brute, directed_indexed, directed_within, targets, vertices};
use super::{hausdorff_distance, hausdorff_within};
use crate::GeometryError;
use crate::predicates::borrowed;
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString,
    MultiPoint, Point, Polygon,
};

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
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        lo + self.below((hi - lo + 1) as u64) as i64
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn pt(x: f64, y: f64) -> Point {
    Point::new(c(x, y)).unwrap()
}

fn line(v: &[(f64, f64)]) -> LineString {
    LineString::new(CoordSeq::new(v.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn ring(v: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(v.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new(ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]), vec![]).unwrap()
}

/// A rectangle with a rectangular hole strictly inside it.
fn holed(x0: f64, y0: f64, x1: f64, y1: f64, h: (f64, f64, f64, f64)) -> Polygon {
    let (hx0, hy0, hx1, hy1) = h;
    Polygon::new(
        ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]),
        vec![ring(&[(hx0, hy0), (hx0, hy1), (hx1, hy1), (hx1, hy0)])],
    )
    .unwrap()
}

/// A random geometry with integer coordinates in `[ox, ox + s] × [oy, oy + s]`.
fn random_geometry(rng: &mut Rng, ox: i64, oy: i64, s: i64) -> Geometry {
    let mut p = || (rng.int(ox, ox + s) as f64, rng.int(oy, oy + s) as f64);
    let kind = p().0 as i64 % 6;
    match kind.rem_euclid(6) {
        0 => Geometry::Point(pt(p().0, p().1)),
        1 => {
            let n = 1 + (p().0 as i64).rem_euclid(5);
            Geometry::MultiPoint(MultiPoint::new(
                (0..n)
                    .map(|_| {
                        let q = p();
                        pt(q.0, q.1)
                    })
                    .collect(),
            ))
        }
        2 | 3 => {
            let n = 2 + (p().0 as i64).rem_euclid(6);
            let mut v: Vec<(f64, f64)> = Vec::new();
            for _ in 0..n {
                let q = p();
                if v.last() != Some(&q) {
                    v.push(q);
                }
            }
            if v.len() < 2 {
                v = vec![(ox as f64, oy as f64), ((ox + s.max(1)) as f64, oy as f64)];
            }
            if kind.rem_euclid(6) == 3 {
                let w = vec![
                    (ox as f64, (oy + s) as f64),
                    ((ox + s.max(1)) as f64, oy as f64),
                ];
                Geometry::MultiLineString(MultiLineString::new(vec![line(&v), line(&w)]))
            } else {
                Geometry::LineString(line(&v))
            }
        }
        4 => {
            let (a, b) = (p(), p());
            let (x0, x1) = (a.0.min(b.0), a.0.max(b.0));
            let (y0, y1) = (a.1.min(b.1), a.1.max(b.1));
            if x1 - x0 < 1.0 || y1 - y0 < 1.0 {
                Geometry::Polygon(rect(
                    ox as f64,
                    oy as f64,
                    (ox + s.max(1)) as f64,
                    (oy + s.max(1)) as f64,
                ))
            } else {
                Geometry::Polygon(rect(x0, y0, x1, y1))
            }
        }
        _ => {
            let (x0, y0) = (ox as f64, oy as f64);
            let w = s.max(4) as f64;
            let q = (w / 4.0).floor();
            Geometry::Polygon(holed(
                x0,
                y0,
                x0 + w,
                y0 + w,
                (x0 + q, y0 + q, x0 + w - q, y0 + w - q),
            ))
        }
    }
}

// ===========================================================================
// The independent brute force
// ===========================================================================

type I2 = (i128, i128);

fn ivertices(g: &Geometry, out: &mut Vec<I2>) {
    let ic = |q: &Coord2d| (q.x as i128, q.y as i128);
    match g {
        Geometry::Point(p) => out.push(ic(&p.coord())),
        Geometry::MultiPoint(m) => out.extend(m.points().iter().map(|p| ic(&p.coord()))),
        Geometry::LineString(l) => out.extend(l.coords().iter().map(ic)),
        Geometry::MultiLineString(m) => {
            for l in m.line_strings() {
                out.extend(l.coords().iter().map(ic));
            }
        }
        Geometry::Polygon(p) => {
            for r in core::iter::once(p.exterior()).chain(p.holes().iter()) {
                out.extend(r.coords().iter().map(ic));
            }
        }
        Geometry::MultiPolygon(m) => {
            for p in m.polygons() {
                ivertices(&Geometry::Polygon(p.clone()), out);
            }
        }
        Geometry::GeometryCollection(gc) => {
            for g in gc.geometries() {
                ivertices(g, out);
            }
        }
    }
}

/// Linework as `(a, b)` pairs; a point is `(p, p)`.
fn ilinework(g: &Geometry, out: &mut Vec<(I2, I2)>) {
    let ic = |q: Coord2d| (q.x as i128, q.y as i128);
    match g {
        Geometry::Point(p) => out.push((ic(p.coord()), ic(p.coord()))),
        Geometry::MultiPoint(m) => {
            out.extend(m.points().iter().map(|p| (ic(p.coord()), ic(p.coord()))));
        }
        Geometry::LineString(l) => {
            out.extend(
                l.coords()
                    .as_slice()
                    .windows(2)
                    .map(|w| (ic(w[0]), ic(w[1]))),
            );
        }
        Geometry::MultiLineString(m) => {
            for l in m.line_strings() {
                ilinework(&Geometry::LineString(l.clone()), out);
            }
        }
        Geometry::Polygon(p) => {
            for r in core::iter::once(p.exterior()).chain(p.holes().iter()) {
                let v = r.coords().as_slice();
                for k in 0..v.len() {
                    out.push((ic(v[k]), ic(v[(k + 1) % v.len()])));
                }
            }
        }
        Geometry::MultiPolygon(m) => {
            for p in m.polygons() {
                ilinework(&Geometry::Polygon(p.clone()), out);
            }
        }
        Geometry::GeometryCollection(gc) => {
            for g in gc.geometries() {
                ilinework(g, out);
            }
        }
    }
}

/// Exactly: is `v` within `√dd` of the closed segment `(a, b)`?
fn iwithin(v: I2, (a, b): (I2, I2), dd: i128) -> bool {
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

fn brute_within(a: &Geometry, b: &Geometry, dd: i128) -> bool {
    let dir = |x: &Geometry, y: &Geometry| {
        let (mut vs, mut es) = (Vec::new(), Vec::new());
        ivertices(x, &mut vs);
        ilinework(y, &mut es);
        vs.iter().all(|&v| es.iter().any(|&e| iwithin(v, e, dd)))
    };
    dir(a, b) && dir(b, a)
}

/// An `f64` estimate of the true value, to aim `d` near it. Not part of
/// any check.
fn estimate(a: &Geometry, b: &Geometry) -> f64 {
    let seg = |v: I2, (s, t): (I2, I2)| {
        let (px, py) = (v.0 as f64, v.1 as f64);
        let (ax, ay, bx, by) = (s.0 as f64, s.1 as f64, t.0 as f64, t.1 as f64);
        let (ux, uy) = (bx - ax, by - ay);
        let l = ux * ux + uy * uy;
        let tt = if l == 0.0 {
            0.0
        } else {
            (((px - ax) * ux + (py - ay) * uy) / l).clamp(0.0, 1.0)
        };
        (px - ax - tt * ux).hypot(py - ay - tt * uy)
    };
    let dir = |x: &Geometry, y: &Geometry| {
        let (mut vs, mut es) = (Vec::new(), Vec::new());
        ivertices(x, &mut vs);
        ilinework(y, &mut es);
        vs.iter()
            .map(|&v| es.iter().map(|&e| seg(v, e)).fold(f64::INFINITY, f64::min))
            .fold(0.0, f64::max)
    };
    dir(a, b).max(dir(b, a))
}

// ===========================================================================
// The predicate is exact
// ===========================================================================

#[test]
fn the_predicate_equals_an_i128_brute_force() {
    let mut rng = Rng(0x48A5_D0F1);
    let mut cases = 0_u32;
    let mut trues = 0_u32;
    let mut equal_hits = 0_u32;
    let mut indexed = 0_u32;
    let top = 1_i64 << 26;
    while cases < 10_000 {
        // Families: a tiny shared window (exact hits common), a medium one,
        // the full range, and long lines that take the indexed path.
        let family = rng.below(4);
        let (a, b) = if family == 3 {
            let mut mk = |n: usize, ox: i64| {
                let mut v: Vec<(f64, f64)> = Vec::new();
                let (mut x, mut y) = (ox, rng.int(0, 1 << 20));
                for _ in 0..n {
                    x += rng.int(1, 3_000);
                    y += rng.int(-2_000, 2_000);
                    v.push((x as f64, y as f64));
                }
                Geometry::LineString(line(&v))
            };
            indexed += 1;
            (mk(60, 0), mk(60, 500))
        } else {
            let s = [16, 4_096, top - 1][family as usize];
            let ox = rng.int(0, top - 1 - s);
            let oy = rng.int(0, top - 1 - s);
            let sb = if rng.below(2) == 0 { s } else { (s / 2).max(1) };
            let a = random_geometry(&mut rng, ox, oy, s);
            let shift = rng.int(0, s / 4);
            let b = random_geometry(&mut rng, ox + shift, oy, sb);
            (a, b)
        };
        let e = estimate(&a, &b).floor() as i64;
        for k in [e - 1, e, e + 1, 0] {
            if k < 0 {
                continue;
            }
            let d = k as f64;
            let dd = i128::from(k) * i128::from(k);
            let want = brute_within(&a, &b, dd);
            let got = hausdorff_within(&a, &b, d).unwrap();
            assert_eq!(got, want, "d = {d}: {a:?} vs {b:?}");
            assert_eq!(hausdorff_within(&b, &a, d).unwrap(), want, "symmetric");
            trues += u32::from(want);
            // A boundary hit: the true squared value lies in (k² − 1, k²],
            // so `d = k` is the first integer that passes.
            if want && k > 0 && !brute_within(&a, &b, dd - 1) {
                equal_hits += 1;
            }
        }
        cases += 1;
    }
    assert!(
        trues > 2_000 && trues < 30_000,
        "both answers well represented: {trues}"
    );
    assert!(equal_hits > 50, "exact boundary hits: {equal_hits}");
    assert!(indexed > 1_000, "indexed cases: {indexed}");
}

// ===========================================================================
// The boundary cases
// ===========================================================================

#[test]
fn a_point_inside_a_polygon_is_measured_to_its_rings() {
    // A point at the center of a square.
    let square = rect(0.0, 0.0, 10.0, 10.0);
    let marks = MultiPoint::new(vec![
        pt(0.0, 0.0),
        pt(10.0, 0.0),
        pt(10.0, 10.0),
        pt(0.0, 10.0),
        pt(5.0, 5.0),
    ]);
    assert_eq!(hausdorff_distance(&square, &marks).unwrap(), 5.0);
    assert_eq!(hausdorff_distance(&marks, &square).unwrap(), 5.0);
    assert!(hausdorff_within(&square, &marks, 5.0).unwrap());
    assert!(!hausdorff_within(&square, &marks, 5.0_f64.next_down()).unwrap());
}

#[test]
fn a_vertex_nearest_a_segment_interior_is_measured_to_the_segment() {
    // (5, 1) is 1 from the segment's interior and √26 from its endpoints.
    let a = line(&[(0.0, 1.0), (5.0, 1.0), (10.0, 1.0)]);
    let b = line(&[(0.0, 0.0), (10.0, 0.0)]);
    assert!(hausdorff_within(&a, &b, 1.0).unwrap());
    assert!(!hausdorff_within(&a, &b, 1.0_f64.next_down()).unwrap());
    assert_eq!(hausdorff_distance(&a, &b).unwrap(), 1.0);
}

#[test]
fn an_asymmetric_pair_is_measured_both_ways() {
    // h(b → a) = 0 (the point is on the line); h(a → b) = 1000.
    let a = line(&[(0.0, 0.0), (1_000.0, 0.0)]);
    let b = pt(0.0, 0.0);
    for (x, y) in [
        (Geometry::LineString(a.clone()), Geometry::Point(b)),
        (Geometry::Point(b), Geometry::LineString(a.clone())),
    ] {
        assert!(!hausdorff_within(&x, &y, 10.0).unwrap());
        assert!(!hausdorff_within(&x, &y, 999.0).unwrap());
        assert!(hausdorff_within(&x, &y, 1_000.0).unwrap());
        assert_eq!(hausdorff_distance(&x, &y).unwrap(), 1_000.0);
    }
    // The directed values themselves.
    let va = vertices(borrowed(&a));
    let vb = vertices(borrowed(&b));
    let (ea, eb) = (targets(borrowed(&a)), targets(borrowed(&b)));
    assert_eq!(directed_brute(&vb, &ea), 0.0);
    assert_eq!(directed_brute(&va, &eb), 1_000.0);
}

#[test]
fn the_value_is_zero_exactly_when_every_vertex_is_on_the_other() {
    let square = rect(0.0, 0.0, 10.0, 10.0);
    let outline = line(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ]);
    assert_eq!(hausdorff_distance(&square, &outline).unwrap(), 0.0);
    assert!(hausdorff_within(&square, &outline, 0.0).unwrap());
    // A collinear vertex inserted: the same linework.
    let a = line(&[(0.0, 0.0), (10.0, 0.0)]);
    let b = line(&[(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]);
    assert_eq!(hausdorff_distance(&a, &b).unwrap(), 0.0);
    // One vertex a hair off: not zero, and positive.
    let off = line(&[(0.0, 0.0), (5.0, 1e-300), (10.0, 0.0)]);
    let h = hausdorff_distance(&a, &off).unwrap();
    assert!(h > 0.0);
    assert!(!hausdorff_within(&a, &off, 0.0).unwrap());
    // A vertex just off a long, far segment: the computed distance may
    // round to zero, and the guard keeps the value positive.
    let far = line(&[(1e15, 1e15), (3e15 + 1.0, 3e15)]);
    let near = MultiPoint::new(vec![pt(1e15, 1e15), pt(2e15, 2e15), pt(3e15 + 1.0, 3e15)]);
    let zero = hausdorff_within(&far, &near, 0.0).unwrap();
    let v = hausdorff_distance(&far, &near).unwrap();
    assert_eq!(v == 0.0, zero, "zero exactly when the zero test holds");
}

#[test]
fn empty_operands_and_invalid_tolerances() {
    let empty = MultiPoint::new(vec![]);
    let p = pt(1.0, 2.0);
    assert!(!hausdorff_within(&empty, &p, 10.0).unwrap());
    assert!(!hausdorff_within(&p, &empty, 10.0).unwrap());
    assert!(matches!(
        hausdorff_distance(&empty, &p),
        Err(GeometryError::EmptyOperand)
    ));
    assert!(matches!(
        hausdorff_distance(&p, &empty),
        Err(GeometryError::EmptyOperand)
    ));
    for d in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            hausdorff_within(&p, &p, d),
            Err(GeometryError::InvalidTolerance { .. })
        ));
    }
}

#[test]
fn an_overlapping_collection_is_accepted() {
    // No topology is consulted: the linework is well defined.
    let gc = GeometryCollection::new(vec![
        Geometry::Polygon(rect(0.0, 0.0, 10.0, 10.0)),
        Geometry::Polygon(rect(5.0, 5.0, 15.0, 15.0)),
    ])
    .unwrap();
    assert_eq!(hausdorff_distance(&gc, &gc).unwrap(), 0.0);
    let p = pt(20.0, 15.0);
    assert!(hausdorff_within(&gc, &p, 25.0).unwrap());
}

// ===========================================================================
// Indexing never changes the answer
// ===========================================================================

#[test]
fn the_indexed_paths_equal_the_exhaustive_ones() {
    let mut rng = Rng(0x1D3C_55AA);
    let mut checked = 0_u32;
    for case in 0..60 {
        let offset = [0.0, 1e6, 1e9][case % 3];
        let n = 80 + rng.below(120) as usize;
        let m = 80 + rng.below(120) as usize;
        let mk = |rng: &mut Rng, k: usize, shift: f64| {
            let (mut x, mut y) = (offset + shift, offset + rng.unit() * 100.0);
            let mut v = Vec::with_capacity(k);
            for _ in 0..k {
                x += 0.5 + rng.unit() * 10.0;
                y += (rng.unit() - 0.5) * 20.0;
                v.push(c(x, y));
            }
            v
        };
        let va = mk(&mut rng, n, 0.0);
        let vb = mk(&mut rng, m, 3.0);
        let eb: Vec<Element> = vb
            .windows(2)
            .map(|w| Element::Segment(w[0], w[1]))
            .collect();
        let delta = 64.0 * (f64::EPSILON / 2.0) * (offset + 3_000.0);
        assert!(va.len() * eb.len() > 2_048);
        let brute = directed_brute(&va, &eb);
        let idx = directed_indexed(&va, &eb, delta).unwrap();
        assert_eq!(brute.to_bits(), idx.to_bits(), "bit for bit");
        for d in [brute, brute.next_down(), brute * 0.5, brute * 2.0] {
            let want = va.iter().all(|&v| {
                eb.iter()
                    .any(|&e| super::distance::pair_within(Element::Point(v), e, d))
            });
            assert_eq!(directed_within(&va, &eb, d).unwrap(), want);
            checked += 1;
        }
    }
    assert!(checked >= 240);
}
