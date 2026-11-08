// src/ops/distance/tests.rs

//! `dwithin` is exact, `distance` is zero exactly on contact, the
//! indexed and prepared forms equal the direct ones, and the boundary cases
//! answer exactly.
//!
//! The oracle here is independent of the code under test: it enumerates
//! boundary elements straight from the geometry types and compares squared
//! distances in `i128`, which is exact for integer coordinates up to `2²⁶`
//! (the largest intermediate, `cross²`, is below `2¹¹⁰`). Only the
//! intersection branch reuses the crate, through `predicates::intersects`,
//! which is exact and has its own suites.
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
// Under the `3d` feature, `Geometry` has variants the test generator never
// builds; without it, the wildcard arm is unreachable.
#![allow(unreachable_patterns)]

use super::*;
use crate::predicates::intersects;
use crate::prepared::PreparedPolygonal;
use crate::types::{
    CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString, MultiPoint,
    MultiPolygon, Point, Polygon,
};

// ===========================================================================
// Fixtures
// ===========================================================================

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn pt(x: f64, y: f64) -> Point {
    Point::new(c(x, y)).unwrap()
}

fn line(pts: &[(f64, f64)]) -> LineString {
    LineString::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

/// An axis-aligned rectangle, counterclockwise.
fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new_simple(ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])).unwrap()
}

/// A deterministic xorshift64*, as the other suites use.
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

    /// An integer in `[lo, hi]`.
    fn int(&mut self, lo: i64, hi: i64) -> i64 {
        let span = (hi - lo + 1) as u64;
        lo + (self.next() % span) as i64
    }
}

// ===========================================================================
// The independent oracle
// ===========================================================================

/// An oracle element with integer coordinates, read straight from the types.
#[derive(Clone, Copy)]
enum OElem {
    P(i64, i64),
    S(i64, i64, i64, i64),
}

fn as_i(v: f64) -> i64 {
    assert_eq!(v, v.trunc(), "oracle inputs are integer-valued");
    v as i64
}

fn oracle_elements(g: &Geometry, out: &mut Vec<OElem>) {
    let seg = |a: Coord2d, b: Coord2d| OElem::S(as_i(a.x), as_i(a.y), as_i(b.x), as_i(b.y));
    match g {
        Geometry::Point(p) => out.push(OElem::P(as_i(p.coord().x), as_i(p.coord().y))),
        Geometry::MultiPoint(m) => {
            for p in m.points() {
                out.push(OElem::P(as_i(p.coord().x), as_i(p.coord().y)));
            }
        }
        Geometry::LineString(l) => {
            for w in l.coords().as_slice().windows(2) {
                out.push(seg(w[0], w[1]));
            }
        }
        Geometry::MultiLineString(m) => {
            for l in m.line_strings() {
                for w in l.coords().as_slice().windows(2) {
                    out.push(seg(w[0], w[1]));
                }
            }
        }
        Geometry::Polygon(p) => {
            for r in core::iter::once(p.exterior()).chain(p.holes().iter()) {
                for (a, b) in r.edges() {
                    out.push(seg(a, b));
                }
            }
        }
        Geometry::MultiPolygon(m) => {
            for p in m.polygons() {
                oracle_elements(&Geometry::Polygon(p.clone()), out);
            }
        }
        Geometry::GeometryCollection(gc) => {
            for g in gc.geometries() {
                oracle_elements(g, out);
            }
        }
        _ => unreachable!("the generator makes no other type"),
    }
}

/// `|pq|² ≤ d²`, exactly, in `i128`.
fn pp_within(px: i64, py: i64, qx: i64, qy: i64, d2: i128) -> bool {
    let (dx, dy) = ((px - qx) as i128, (py - qy) as i128);
    dx * dx + dy * dy <= d2
}

/// `dist(p, seg ab)² ≤ d²`, exactly: the clamped-projection case split with
/// the interior case cleared of its division, as `cross² ≤ d²·|u|²`.
fn ps_within(px: i64, py: i64, ax: i64, ay: i64, bx: i64, by: i64, d2: i128) -> bool {
    let (ux, uy) = ((bx - ax) as i128, (by - ay) as i128);
    let (wx, wy) = ((px - ax) as i128, (py - ay) as i128);
    let dot = wx * ux + wy * uy;
    let len2 = ux * ux + uy * uy;
    if dot <= 0 || len2 == 0 {
        return wx * wx + wy * wy <= d2;
    }
    if dot >= len2 {
        return pp_within(px, py, bx, by, d2);
    }
    let cross = ux * wy - uy * wx;
    cross * cross <= d2 * len2
}

/// Whether two oracle elements of **non-intersecting** geometries are within
/// `d`: the endpoint reduction, sound because every crossing was caught by
/// the intersection branch first.
fn oracle_pair_within(a: OElem, b: OElem, d2: i128) -> bool {
    match (a, b) {
        (OElem::P(px, py), OElem::P(qx, qy)) => pp_within(px, py, qx, qy, d2),
        (OElem::P(px, py), OElem::S(ax, ay, bx, by))
        | (OElem::S(ax, ay, bx, by), OElem::P(px, py)) => ps_within(px, py, ax, ay, bx, by, d2),
        (OElem::S(a1x, a1y, a2x, a2y), OElem::S(b1x, b1y, b2x, b2y)) => {
            ps_within(a1x, a1y, b1x, b1y, b2x, b2y, d2)
                || ps_within(a2x, a2y, b1x, b1y, b2x, b2y, d2)
                || ps_within(b1x, b1y, a1x, a1y, a2x, a2y, d2)
                || ps_within(b2x, b2y, a1x, a1y, a2x, a2y, d2)
        }
    }
}

fn oracle_dwithin(a: &Geometry, b: &Geometry, d: i64) -> bool {
    if intersects(a, b).unwrap() {
        return true;
    }
    let (mut ea, mut eb) = (Vec::new(), Vec::new());
    oracle_elements(a, &mut ea);
    oracle_elements(b, &mut eb);
    let d2 = (d as i128) * (d as i128);
    ea.iter()
        .any(|&x| eb.iter().any(|&y| oracle_pair_within(x, y, d2)))
}

// ===========================================================================
// The random generator
// ===========================================================================

const GRID: i64 = 1 << 26;

/// A random geometry with integer coordinates near `anchor`, every type.
fn random_geometry(rng: &mut Rng, ax: i64, ay: i64) -> Geometry {
    let o = |rng: &mut Rng| rng.int(0, 48);
    match rng.next() % 7 {
        0 => Geometry::Point(pt((ax + o(rng)) as f64, (ay + o(rng)) as f64)),
        1 => {
            let n = rng.int(1, 4);
            Geometry::MultiPoint(MultiPoint::new(
                (0..n)
                    .map(|_| pt((ax + o(rng)) as f64, (ay + o(rng)) as f64))
                    .collect(),
            ))
        }
        2 => {
            // Strictly increasing x keeps the vertices distinct.
            let x0 = ax + o(rng);
            Geometry::LineString(line(&[
                (x0 as f64, (ay + o(rng)) as f64),
                ((x0 + 1 + o(rng)) as f64, (ay + o(rng)) as f64),
                ((x0 + 60 + o(rng)) as f64, (ay + o(rng)) as f64),
            ]))
        }
        3 => {
            let x0 = ax + o(rng);
            let l1 = line(&[
                (x0 as f64, (ay + o(rng)) as f64),
                ((x0 + 1 + o(rng)) as f64, (ay + o(rng)) as f64),
            ]);
            let y1 = ay + 60 + o(rng);
            let l2 = line(&[
                ((ax + o(rng)) as f64, y1 as f64),
                ((ax + 1 + o(rng)) as f64, (y1 + 1 + o(rng)) as f64),
            ]);
            Geometry::MultiLineString(MultiLineString::new(vec![l1, l2]))
        }
        4 => {
            let (x0, y0) = (ax + o(rng), ay + o(rng));
            Geometry::Polygon(rect(
                x0 as f64,
                y0 as f64,
                (x0 + 1 + o(rng)) as f64,
                (y0 + 1 + o(rng)) as f64,
            ))
        }
        5 => {
            let (x0, y0) = (ax + o(rng), ay + o(rng));
            let w = 1 + o(rng);
            let a = rect(x0 as f64, y0 as f64, (x0 + w) as f64, (y0 + w) as f64);
            let x1 = x0 + w + 2 + o(rng);
            let b = rect(
                x1 as f64,
                y0 as f64,
                (x1 + 1 + o(rng)) as f64,
                (y0 + 1 + o(rng)) as f64,
            );
            Geometry::MultiPolygon(MultiPolygon::new(vec![a, b]).unwrap())
        }
        _ => {
            let p = pt((ax + o(rng)) as f64, (ay + o(rng)) as f64);
            let x0 = ax + 60 + o(rng);
            let l = line(&[
                (x0 as f64, (ay + o(rng)) as f64),
                ((x0 + 1 + o(rng)) as f64, (ay + 1 + o(rng)) as f64),
            ]);
            Geometry::GeometryCollection(
                GeometryCollection::new(vec![Geometry::Point(p), Geometry::LineString(l)]).unwrap(),
            )
        }
    }
}

/// A random pair: anchors drawn across the full `±2²⁶` grid, and, half the
/// time, `b`'s anchor pulled next to `a`'s so near-tie distances occur.
fn random_pair(rng: &mut Rng) -> (Geometry, Geometry) {
    let (ax, ay) = (rng.int(-GRID, GRID - 200), rng.int(-GRID, GRID - 200));
    let (bx, by) = if rng.next().is_multiple_of(2) {
        (ax + rng.int(-120, 120), ay + rng.int(-120, 120))
    } else {
        (rng.int(-GRID, GRID - 200), rng.int(-GRID, GRID - 200))
    };
    (random_geometry(rng, ax, ay), random_geometry(rng, bx, by))
}

const DS: [i64; 8] = [0, 1, 2, 3, 5, 10, 50, 141];

// ===========================================================================
// dwithin equals the exact oracle
// ===========================================================================

#[test]
fn dwithin_matches_the_exact_integer_oracle_on_random_pairs() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut checked = 0_u32;
    let mut true_count = 0_u32;
    for _ in 0..1_300 {
        let (a, b) = random_pair(&mut rng);
        for d in DS {
            let got = dwithin(&a, &b, d as f64).unwrap();
            let want = oracle_dwithin(&a, &b, d);
            assert_eq!(got, want, "dwithin disagrees with the oracle at d = {d}");
            checked += 1;
            true_count += u32::from(got);
        }
    }
    assert!(checked >= 10_000, "only {checked} cases");
    // The corpus exercises both answers.
    assert!(true_count > 500 && true_count < checked - 500);
}

#[test]
fn distance_is_zero_iff_intersects_and_consistent_with_dwithin() {
    let mut rng = Rng(0x0DD0_5EED_C0FF_EE01);
    let mut nonzero = 0_u32;
    for _ in 0..2_500 {
        let (a, b) = random_pair(&mut rng);
        let dist = distance(&a, &b).unwrap();
        let meets = intersects(&a, &b).unwrap();
        assert_eq!(dist == 0.0, meets, "distance zero iff intersects");
        nonzero += u32::from(dist > 0.0);
        let mut es = elements(crate::predicates::borrowed(&a));
        let mut eb = elements(crate::predicates::borrowed(&b));
        es.append(&mut eb);
        let delta = DIST_ERR_C * U * coordinate_magnitude(&es);
        for d in DS {
            let d = d as f64;
            if dwithin(&a, &b, d).unwrap() {
                assert!(dist < d + delta, "within d but distance {dist} ≥ {d} + δ");
            } else {
                assert!(
                    dist > d - delta,
                    "not within d but distance {dist} ≤ {d} − δ"
                );
            }
        }
    }
    assert!(nonzero > 500);
}

// ===========================================================================
// The boundary cases
// ===========================================================================

#[test]
fn the_boundary_cases_answer_exactly() {
    // A 3-4-5 pair: d exactly the distance is true; one ulp below is false.
    let (p, q) = (pt(0.0, 0.0), pt(3.0, 4.0));
    assert!(dwithin(&p, &q, 5.0).unwrap());
    assert!(!dwithin(&p, &q, 5.0_f64.next_down()).unwrap());
    assert_eq!(distance(&p, &q).unwrap(), 5.0);

    // Crossing segments: within 0, by the intersection branch, though every
    // endpoint is far from the other segment.
    let s1 = line(&[(-10.0, 0.0), (10.0, 0.0)]);
    let s2 = line(&[(0.0, -10.0), (0.0, 10.0)]);
    assert!(dwithin(&s1, &s2, 0.0).unwrap());
    assert_eq!(distance(&s1, &s2).unwrap(), 0.0);

    // A point inside a polygon: within every d, distance exactly 0, with no
    // boundary element of the polygon within 3 of it.
    let big = rect(-10.0, -10.0, 10.0, 10.0);
    let inside = pt(0.0, 0.0);
    assert!(dwithin(&big, &inside, 0.0).unwrap());
    assert!(dwithin(&big, &inside, 3.0).unwrap());
    assert_eq!(distance(&big, &inside).unwrap(), 0.0);

    // Parallel segments exactly 5 apart.
    let flat = line(&[(0.0, 0.0), (10.0, 0.0)]);
    let above = line(&[(2.0, 5.0), (8.0, 5.0)]);
    assert!(dwithin(&flat, &above, 5.0).unwrap());
    assert!(!dwithin(&flat, &above, 5.0_f64.next_down()).unwrap());
    assert_eq!(distance(&flat, &above).unwrap(), 5.0);

    // Closest approach at an interior point of the segment: the point
    // projects inside, so neither endpoint attains the distance.
    let p = pt(5.0, 3.0);
    assert!(dwithin(&flat, &p, 3.0).unwrap());
    assert!(!dwithin(&flat, &p, 3.0_f64.next_down()).unwrap());
    assert_eq!(distance(&flat, &p).unwrap(), 3.0);

    // Segment to segment nearest at interior points (both project inside).
    let short = line(&[(3.0, 4.0), (7.0, 4.0)]);
    assert!(dwithin(&flat, &short, 4.0).unwrap());
    assert!(!dwithin(&flat, &short, 4.0_f64.next_down()).unwrap());
    assert_eq!(distance(&flat, &short).unwrap(), 4.0);

    // An exact non-integer tie the predicate still decides exactly:
    // distance √2 between (0,0) and (1,1); d = fl(√2) is one side of it.
    let diag = pt(1.0, 1.0);
    let origin = pt(0.0, 0.0);
    let fl_sqrt2 = 2.0_f64.sqrt();
    // fl(√2) > √2 (fl(√2)² > 2, checked exactly in rationals), so fl(√2)
    // is within and its lower neighbor is not: the predicate places the
    // boundary between consecutive floats, at the irrational distance.
    assert!(dwithin(&origin, &diag, fl_sqrt2).unwrap());
    assert!(!dwithin(&origin, &diag, fl_sqrt2.next_down()).unwrap());
}

#[test]
fn errors_and_empties_are_as_documented() {
    let p = pt(0.0, 0.0);
    for bad in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            dwithin(&p, &p, bad),
            Err(GeometryError::InvalidTolerance { .. })
        ));
    }
    let empty = MultiPoint::empty();
    assert!(!dwithin(&p, &empty, 1.0e300).unwrap());
    assert!(!dwithin(&empty, &empty, 1.0e300).unwrap());
    assert!(matches!(
        distance(&p, &empty),
        Err(GeometryError::EmptyOperand)
    ));
    assert!(matches!(
        distance(&empty, &p),
        Err(GeometryError::EmptyOperand)
    ));

    // A collection whose polygonal leaves overlap is refused, as relate
    // refuses it. Both functions refuse it, ahead of the empty answer.
    let overlapping = Geometry::GeometryCollection(
        GeometryCollection::new(vec![
            Geometry::Polygon(rect(0.0, 0.0, 4.0, 4.0)),
            Geometry::Polygon(rect(1.0, 1.0, 5.0, 5.0)),
        ])
        .unwrap(),
    );
    assert!(matches!(
        dwithin(&overlapping, &p, 1.0),
        Err(GeometryError::CollectionElementsOverlap { .. })
    ));
    assert!(matches!(
        distance(&overlapping, &p),
        Err(GeometryError::CollectionElementsOverlap { .. })
    ));
    assert!(matches!(
        dwithin(&overlapping, &empty, 1.0),
        Err(GeometryError::CollectionElementsOverlap { .. })
    ));
}

// ===========================================================================
// The indexed forms equal the exhaustive scan
// ===========================================================================

/// Two multipoints big enough that `dwithin` and `distance` take the indexed
/// path. The test calls both paths directly and requires equality, bit for
/// bit in the case of `distance`.
#[test]
fn the_indexed_search_equals_the_exhaustive_scan() {
    let mut rng = Rng(0xFEED_FACE_CAFE_BEEF);
    for round in 0..12 {
        let n = 70; // 70 × 70 = 4 900 pairs > BRUTE_PAIRS
        let make = |rng: &mut Rng| -> Vec<Element> {
            (0..n)
                .map(|_| {
                    Element::Point(c(
                        rng.int(-1_000, 1_000) as f64,
                        rng.int(-1_000, 1_000) as f64,
                    ))
                })
                .collect()
        };
        let (ea, eb) = (make(&mut rng), make(&mut rng));
        // A disjoint corpus for distance: shift b away when they collide.
        for d in [0.0, 1.0, 7.0, 55.5, 2_000.0, 5_000.0] {
            assert_eq!(
                within_brute(&ea, &eb, d),
                within_indexed(&ea, &eb, d).unwrap(),
                "round {round}, d = {d}"
            );
            assert_eq!(
                within_brute(&ea, &eb, d),
                within_filtered(&ea, &eb, d),
                "filtered, round {round}, d = {d}"
            );
        }
        let delta = DIST_ERR_C * U * 1_000.0;
        let brute = distance_brute(&ea, &eb);
        let indexed = distance_indexed(&ea, &eb, delta).unwrap();
        assert_eq!(brute.to_bits(), indexed.to_bits(), "round {round}");
        let pruned = distance_pruned(&ea, &eb, delta);
        assert_eq!(brute.to_bits(), pruned.to_bits(), "pruned, round {round}");
    }

    // The same equality on mixed segment elements, asymmetric sizes (the
    // probe/indexed choice goes both ways).
    for round in 0..6 {
        let seg = |rng: &mut Rng| {
            let (x, y) = (rng.int(-2_000, 2_000), rng.int(-2_000, 2_000));
            Element::Segment(
                c(x as f64, y as f64),
                c((x + rng.int(1, 9)) as f64, (y + rng.int(1, 9)) as f64),
            )
        };
        let ea: Vec<Element> = (0..150).map(|_| seg(&mut rng)).collect();
        let eb: Vec<Element> = (0..40).map(|_| seg(&mut rng)).collect();
        for d in [0.0, 2.0, 10.0, 300.0] {
            assert_eq!(
                within_brute(&ea, &eb, d),
                within_indexed(&ea, &eb, d).unwrap(),
                "segment round {round}, d = {d}"
            );
            assert_eq!(
                within_brute(&eb, &ea, d),
                within_indexed(&eb, &ea, d).unwrap(),
                "segment round {round} swapped, d = {d}"
            );
            for (x, y) in [(&ea, &eb), (&eb, &ea)] {
                assert_eq!(
                    within_brute(x, y, d),
                    within_filtered(x, y, d),
                    "filtered segment round {round}, d = {d}"
                );
            }
        }
        let delta = DIST_ERR_C * U * 2_009.0;
        assert_eq!(
            distance_brute(&ea, &eb).to_bits(),
            distance_indexed(&ea, &eb, delta).unwrap().to_bits(),
            "segment round {round}"
        );
        for (x, y) in [(&ea, &eb), (&eb, &ea)] {
            assert_eq!(
                distance_brute(x, y).to_bits(),
                distance_pruned(x, y, delta).to_bits(),
                "pruned segment round {round}"
            );
        }
    }
}

// ===========================================================================
// The prepared forms equal the direct forms
// ===========================================================================

#[test]
fn the_prepared_forms_equal_the_direct_forms() {
    let holed = Polygon::new(
        ring(&[(0.0, 0.0), (30.0, 0.0), (30.0, 30.0), (0.0, 30.0)]),
        vec![ring(&[
            (10.0, 10.0),
            (10.0, 20.0),
            (20.0, 20.0),
            (20.0, 10.0),
        ])],
    )
    .unwrap();
    let mask = MultiPolygon::new(vec![holed, rect(50.0, 0.0, 60.0, 10.0)]).unwrap();
    let prepared = PreparedPolygonal::new(&mask).unwrap();
    assert!(!prepared.is_delegating());

    let others: Vec<Geometry> = vec![
        Geometry::Point(pt(15.0, 15.0)), // in the hole: distance 5 to its ring
        Geometry::Point(pt(5.0, 5.0)),   // inside: 0
        Geometry::Point(pt(40.0, 5.0)),  // between the parts
        Geometry::Point(pt(-3.0, -4.0)), // 5 from the corner
        Geometry::LineString(line(&[(35.0, 5.0), (45.0, 5.0)])),
        Geometry::LineString(line(&[(-5.0, 15.0), (65.0, 15.0)])), // crosses both
        Geometry::Polygon(rect(32.0, 12.0, 48.0, 28.0)),
        Geometry::MultiPoint(MultiPoint::new(vec![pt(100.0, 100.0), pt(31.0, 15.0)])),
    ];
    for g in &others {
        for d in [0.0, 1.0, 2.0, 5.0, 5.0_f64.next_down(), 25.0] {
            assert_eq!(
                prepared.dwithin(g, d).unwrap(),
                dwithin(&mask, g, d).unwrap(),
                "dwithin({g:?}, {d})"
            );
        }
        assert_eq!(
            prepared.distance(g).unwrap().to_bits(),
            distance(&mask, g).unwrap().to_bits(),
            "distance({g:?})"
        );
    }

    // A delegating mask (a self-crossing ring, accepted by `Polygon::new`,
    // ADR-27): every edge is a candidate, and the answers still match.
    let crossing =
        Polygon::new_simple(ring(&[(0.0, 0.0), (5.0, 4.0), (0.0, 3.0), (4.0, 0.0)])).unwrap();
    let dmask = MultiPolygon::new(vec![crossing]).unwrap();
    let dprep = PreparedPolygonal::new(&dmask).unwrap();
    assert!(dprep.is_delegating());
    let far = Geometry::Point(pt(10.0, 2.0));
    for d in [0.0, 5.0, 6.0, 10.0] {
        assert_eq!(
            dprep.dwithin(&far, d).unwrap(),
            dwithin(&dmask, &far, d).unwrap()
        );
    }
    assert_eq!(
        dprep.distance(&far).unwrap().to_bits(),
        distance(&dmask, &far).unwrap().to_bits()
    );

    // An empty mask: dwithin false, distance refused.
    let eprep = PreparedPolygonal::new(&MultiPolygon::empty()).unwrap();
    assert!(!eprep.dwithin(&far, 1.0e300).unwrap());
    assert!(matches!(
        eprep.distance(&far),
        Err(GeometryError::EmptyOperand)
    ));
}

/// The pruning radius must be at least `best`: a radius of `(best − δ)`
/// prunes a pair whose distance lies inside the slack, and the indexed
/// result parts from the exhaustive one. A mutation that flips the slack's
/// sign fails this case.
///
/// The construction: coordinates at magnitude `10⁶`, so
/// `δ = 64·u·10⁶ ≈ 7.1·10⁻⁹`. The seed pair sits at exactly `10⁶`; a second
/// indexed element sits at `10⁶ − 4·10⁻⁹`, inside the slack below `best`.
/// A radius of `(best − δ).next_down()` prunes it (its box gap exceeds
/// that radius by far more than the window's rounding) while the correct
/// radius `(best + δ).next_up()` keeps it.
#[test]
fn the_pruning_slack_is_load_bearing() {
    let m = 1.0e6_f64;
    let delta = DIST_ERR_C * U * m;
    let near = m - 4.0e-9;
    assert!(
        near < m && near > m - delta,
        "the fixture sits inside the slack"
    );
    let probe = vec![Element::Point(c(0.0, 0.0))];
    let indexed = vec![Element::Point(c(m, 0.0)), Element::Point(c(near, 0.0))];
    let brute = distance_brute(&probe, &indexed);
    let fast = distance_indexed(&probe, &indexed, delta).unwrap();
    assert!(brute < m, "the nearer pair is the minimum");
    assert_eq!(brute.to_bits(), fast.to_bits());
}

// ===========================================================================
// The helpers' own guarantees
// ===========================================================================

#[test]
fn grow_outward_contains_the_exact_grown_box() {
    // Exact check, by an error-free transformation: Knuth's TwoSum gives
    // `fl(v ± r) = s` and the exact rounding error `err` with
    // `v ± r = s + err` over the reals (no overflow in this range). The
    // returned lower bound `lo` is at or below `v − r` exactly when
    // `(s − lo) + err ≥ 0`; `s − lo` is exact (adjacent floats, Sterbenz),
    // and comparing it with `−err` is a comparison of two exact values.
    fn two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        let bb = s - a;
        (s, (a - (s - bb)) + (b - bb))
    }
    let mut rng = Rng(0xB0A7_1D3A_5EED_0001);
    let mut rounded_inward = 0_u32; // cases where fl(v − r) > v − r: the directed step mattered
    let mut checked = 0_u32;
    for _ in 0..40_000 {
        // Magnitudes 2^-60..2^60 for both, so v ± r is far from overflow.
        let e1 = (rng.next() % 121) as i32 - 60;
        let e2 = (rng.next() % 121) as i32 - 60;
        let frac = |rng: &mut Rng| 1.0 + (rng.next() >> 12) as f64 / (1_u64 << 52) as f64;
        let v = frac(&mut rng) * 2.0_f64.powi(e1);
        let v = if rng.next().is_multiple_of(2) { v } else { -v };
        let r = frac(&mut rng) * 2.0_f64.powi(e2);
        let e = grow_outward((c(v, v), c(v, v)), r);

        let (s, err) = two_sum(v, -r);
        let lo = e.min().x;
        assert!(
            s - lo >= -err,
            "lower bound above v − r: v = {v:e}, r = {r:e}"
        );
        rounded_inward += u32::from(err < 0.0);

        let (s, err) = two_sum(v, r);
        let hi = e.max().x;
        assert!(
            hi - s >= err,
            "upper bound below v + r: v = {v:e}, r = {r:e}"
        );
        checked += 1;
    }
    assert!(checked >= 40_000);
    // The corpus reaches the case the directed step exists for.
    assert!(
        rounded_inward > 1_000,
        "only {rounded_inward} inward roundings"
    );

    // The clamp: growth from the edge of the range stays finite.
    let e = grow_outward((c(-f64::MAX, 0.0), c(f64::MAX, 0.0)), f64::MAX);
    assert!(e.min().x.is_finite() && e.max().x.is_finite());
    assert_eq!(e.min().x, f64::MIN);
    assert_eq!(e.max().x, f64::MAX);
}

#[test]
fn norm_is_exact_on_pythagorean_inputs_and_survives_the_prescale() {
    assert_eq!(norm(3.0, 4.0), 5.0);
    assert_eq!(norm(-3.0, 4.0), 5.0);
    assert_eq!(norm(0.0, 0.0), 0.0);
    assert_eq!(norm(7.0, 0.0), 7.0);
    // Past 2^510 the squares would overflow; the prescale is exact
    // (powers of two), so the 3-4-5 triple stays exact at 2^520.
    let s = 2.0_f64.powi(520);
    assert_eq!(norm(3.0 * s, 4.0 * s), 5.0 * s);
    // And below 2^-511 the squares would flush to zero without it.
    let t = 2.0_f64.powi(-520);
    assert_eq!(norm(3.0 * t, 4.0 * t), 5.0 * t);
    // A non-finite difference propagates as +∞, never a finite wrong value.
    assert_eq!(norm(f64::INFINITY, 1.0), f64::INFINITY);
}

#[test]
fn point_segment_cases_are_exact_where_the_arithmetic_is() {
    let (a, b) = (c(0.0, 0.0), c(10.0, 0.0));
    assert_eq!(point_segment(c(5.0, 3.0), a, b), 3.0); // interior projection
    assert_eq!(point_segment(c(-3.0, 4.0), a, b), 5.0); // before a
    assert_eq!(point_segment(c(13.0, 4.0), a, b), 5.0); // past b
    assert_eq!(point_segment(c(7.0, 0.0), a, b), 0.0); // on the segment
    // A degenerate segment is the point distance.
    assert_eq!(point_segment(c(3.0, 4.0), a, a), 5.0);
    // Endpoints ~1.5e308 apart: the projection's products overflow, and the
    // fallback returns a finite endpoint distance, never a finite wrong one.
    let far = point_segment(c(0.0, 0.0), c(-f64::MAX, 3.0), c(f64::MAX, 4.0));
    assert!(
        far >= 3.0,
        "the fallback is an endpoint distance, ≥ the true one"
    );
}

// ===========================================================================
// The zero guard: a rounding never claims contact
// ===========================================================================

/// A point placed on a segment in `f64` and nudged by an ulp is often one
/// rounding off the segment: the exact test says the two do not meet, while
/// the clamped projection's residual cancels to `0.0`. Before the zero
/// guard, `distance` returned that `0.0`, which claims a contact the exact
/// test had just denied. Every probe here must keep `distance == 0.0`
/// exactly when the operands intersect, through the free function and the
/// prepared one, and the family must actually reach the guard.
#[test]
fn a_point_a_rounding_off_a_segment_is_not_at_distance_zero() {
    // The case first found: exact says apart; the raw arithmetic says 0.
    let seg = line(&[(834_029.0, 134_536.0), (1_346_245.0, 957_342.0)]);
    let p = pt(906_251.456, 250_551.646);
    assert!(!dwithin(&seg, &p, 0.0).unwrap());
    let raw = distance_brute(&elements(borrowed(&seg)), &elements(borrowed(&p)));
    assert_eq!(raw, 0.0, "the arithmetic rounds this pair to zero");
    assert_eq!(
        distance(&seg, &p).unwrap().to_bits(),
        1,
        "the guard's value"
    );

    // A sweep: segments at three magnitudes, points placed on them in f64
    // and nudged by an ulp half the time.
    let mut rng = Rng(0x2E80_6A4D);
    let (mut guarded, mut apart, mut prepared_checked) = (0_u32, 0_u32, 0_u32);
    for k in 0..20_000_u32 {
        let scale = [1.0, 1e3, 1e6][(k % 3) as usize];
        let (ax, ay) = (
            (rng.next() % 1_000_000) as f64 * scale,
            (rng.next() % 1_000_000) as f64 * scale,
        );
        let (bx, by) = (
            ax + ((rng.next() % 1_000_000) as f64 + 1.0) * scale,
            ay + (rng.next() % 1_000_000) as f64 * scale,
        );
        let t = (rng.next() % 1_000) as f64 / 1_000.0;
        let (px, mut py) = (ax + t * (bx - ax), ay + t * (by - ay));
        if rng.next().is_multiple_of(2) {
            py = py.next_up();
        }
        let s = line(&[(ax, ay), (bx, by)]);
        let q = pt(px, py);
        let meets = dwithin(&s, &q, 0.0).unwrap();
        let d = distance(&s, &q).unwrap();
        assert_eq!(d == 0.0, meets, "zero exactly on contact: {s:?} {q:?}");
        if !meets {
            apart += 1;
            let raw = distance_brute(&elements(borrowed(&s)), &elements(borrowed(&q)));
            guarded += u32::from(raw == 0.0);
        }
        // The prepared form, on a triangle with this segment as an edge,
        // its third corner on the far side from the point when apart.
        if k % 10 == 0 && !meets {
            for third in [(ax, by), (bx, ay)] {
                let tri = [(ax, ay), (bx, by), third];
                let area2 = (bx - ax) * (third.1 - ay) - (by - ay) * (third.0 - ax);
                if area2 == 0.0 {
                    continue;
                }
                let ccw = if area2 > 0.0 {
                    tri
                } else {
                    [tri[0], tri[2], tri[1]]
                };
                let poly = Polygon::new_simple(ring(&ccw)).unwrap();
                let mask = PreparedPolygonal::new(&poly).unwrap();
                let g = Geometry::Point(q);
                if mask.intersects(&g).unwrap() {
                    continue;
                }
                let free = distance(&poly, &q).unwrap();
                assert_eq!(mask.distance(&g).unwrap().to_bits(), free.to_bits());
                assert!(free > 0.0);
                prepared_checked += 1;
            }
        }
    }
    assert!(apart > 5_000, "apart cases: {apart}");
    assert!(guarded > 100, "the guard fired: {guarded}");
    assert!(prepared_checked > 200, "prepared cases: {prepared_checked}");
}

// ===========================================================================
// Contact without the matrix
// ===========================================================================

/// `dwithin` and `distance`'s zero against the definition they replace:
/// `intersects`, then the element scan. Generated pairs of polygons (with
/// holes), multipolygons, lines and points on a small grid, so that they
/// touch, cross with no endpoint near, nest, sit in holes and miss; and the
/// two crossing bars, which only the crossing test sees.
#[test]
fn matrix_free_contact_equals_intersects_then_elements() {
    fn star(rng: &mut Rng, cx: i64, cy: i64, n: usize, hole: bool) -> Option<Polygon> {
        let mut pts = Vec::new();
        for k in 0..n {
            let t = std::f64::consts::TAU * (k as f64 + 0.4) / n as f64;
            let r = rng.int(3, 9) as f64;
            pts.push(c(
                (cx as f64 + r * t.cos()).round(),
                (cy as f64 + r * t.sin()).round(),
            ));
        }
        pts.dedup();
        let ext = LinearRing::new(CoordSeq::new(pts)).ok()?;
        let holes = if hole {
            let h = [
                (cx - 1, cy - 1),
                (cx - 1, cy + 1),
                (cx + 1, cy + 1),
                (cx + 1, cy - 1),
            ];
            vec![
                LinearRing::new(CoordSeq::new(
                    h.iter().map(|&(x, y)| c(x as f64, y as f64)).collect(),
                ))
                .ok()?,
            ]
        } else {
            Vec::new()
        };
        let p = Polygon::new(ext, holes).ok()?;
        crate::ops::is_valid(&p).then_some(p)
    }
    let mut rng = Rng(0x0C0A_7AC7_0000_0001);
    let generate = |rng: &mut Rng| -> Option<Geometry> {
        let (x, y) = (rng.int(-12, 12), rng.int(-12, 12));
        Some(match rng.int(0, 5) {
            0 => Geometry::Point(Point::new(c(x as f64, y as f64)).ok()?),
            1 => Geometry::LineString(
                LineString::new(CoordSeq::new(vec![
                    c(x as f64, y as f64),
                    c((x + rng.int(-9, 9)) as f64, (y + rng.int(-9, 9)) as f64),
                    c((x + rng.int(-9, 9)) as f64, (y + rng.int(-9, 9)) as f64),
                ]))
                .ok()?,
            ),
            2 | 3 => {
                let n = rng.int(3, 12) as usize;
                let hole = rng.int(0, 2) == 0;
                Geometry::Polygon(star(rng, x, y, n, hole)?)
            }
            4 => {
                let a = star(rng, x, y, 6, false)?;
                let b = star(rng, x + 25, y, 5, false)?;
                Geometry::MultiPolygon(MultiPolygon::new(vec![a, b]).ok()?)
            }
            _ => Geometry::MultiPoint(MultiPoint::new(vec![
                Point::new(c(x as f64, y as f64)).ok()?,
                Point::new(c((x + 3) as f64, (y - 2) as f64)).ok()?,
            ])),
        })
    };
    let mut checked = 0;
    let mut met = 0;
    while checked < 3000 {
        let (Some(a), Some(b)) = (generate(&mut rng), generate(&mut rng)) else {
            continue;
        };
        let meet = intersects(&a, &b).unwrap();
        let (ea, eb) = (elements(borrowed(&a)), elements(borrowed(&b)));
        for d in [0.0, 0.5, 1.0, 3.0] {
            let want = meet || within_brute(&ea, &eb, d);
            assert_eq!(dwithin(&a, &b, d), Ok(want), "{a:?} {b:?} {d}");
        }
        assert_eq!(distance(&a, &b).unwrap() == 0.0, meet, "{a:?} {b:?}");
        checked += 1;
        met += usize::from(meet);
    }
    assert!(met > 300 && met < 2700, "{met}");
    // Two bars crossing with every endpoint far from the other bar.
    let bar = |x0: f64, y0: f64, x1: f64, y1: f64| {
        Polygon::new_simple(
            LinearRing::new(CoordSeq::new(vec![
                c(x0, y0),
                c(x1, y0),
                c(x1, y1),
                c(x0, y1),
            ]))
            .unwrap(),
        )
        .unwrap()
    };
    let (h, v) = (bar(0.0, 4.0, 10.0, 6.0), bar(4.0, 0.0, 6.0, 10.0));
    assert_eq!(dwithin(&h, &v, 0.0), Ok(true));
    assert_eq!(distance(&h, &v), Ok(0.0));
    let (lh, lv) = (
        LineString::new(CoordSeq::new(vec![c(0.0, 5.0), c(10.0, 5.0)])).unwrap(),
        LineString::new(CoordSeq::new(vec![c(5.0, 0.0), c(5.0, 10.0)])).unwrap(),
    );
    assert_eq!(dwithin(&lh, &lv, 1.0), Ok(true));
    assert_eq!(distance(&lh, &lv), Ok(0.0));
}
