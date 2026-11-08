// src/ops/concave_hull_tests.rs

//! The concave hull: the exact contract on generated families, the
//! ratio-to-length equivalence, canonicality, degenerate input, parameter
//! validation, and the rules the erosion must never break.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation
)]

use std::collections::{HashMap, HashSet};

use super::concave_hull::{ConcaveHullParameter, concave_hull};
use super::delaunay::{Triangulation, delaunay};
use crate::codec::encode;
use crate::exact::locate;
use crate::ops::convex_hull;
use crate::predicates::relate;
use crate::topology::label::Location;
use crate::types::{Coord2d, Geometry, MultiPoint, Point};

use ConcaveHullParameter::{MaxEdgeLength, MaxEdgeLengthRatio};

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
    fn f(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            v.swap(i, self.below(i + 1));
        }
    }
}

/// The delaunay test families, in brief: random, clustered, integer
/// lattices (cocircular ties), and a grid.
fn generated_set(rng: &mut Rng, case: usize) -> Vec<(f64, f64)> {
    let n = 4 + rng.below(23);
    let mut pts = Vec::with_capacity(n);
    match case % 4 {
        0 => {
            for _ in 0..n {
                pts.push((rng.f() * 100.0, rng.f() * 100.0));
            }
        }
        1 => {
            for _ in 0..n {
                let (cx, cy) = match rng.below(3) {
                    0 => (0.0, 0.0),
                    1 => (200.0, 3.0),
                    _ => (-70.0, 90.0),
                };
                pts.push((cx + 4.0 * rng.f(), cy + 4.0 * rng.f()));
            }
        }
        2 => {
            for _ in 0..n {
                pts.push((rng.below(10) as f64, rng.below(10) as f64));
            }
        }
        _ => {
            let k = 2 + rng.below(3);
            for i in 0..=k {
                for j in 0..=k {
                    pts.push((i as f64, j as f64));
                }
            }
        }
    }
    pts
}

/// The ring of a polygon result, as coordinates.
fn ring_of(g: &Geometry) -> Vec<Coord2d> {
    match g {
        Geometry::Polygon(p) => p.exterior().coords().as_slice().to_vec(),
        other => panic!("not a polygon: {other:?}"),
    }
}

/// The hull's whole contract, checked independently of the operation:
/// containment by [`locate`] (the unprepared path), every ring edge a
/// Delaunay edge, equality with the convex hull at ratio 1, and the fixed
/// point at threshold 0.
fn check_contract(pts: &[(f64, f64)]) {
    let g = multipoint(pts);
    let t = delaunay(&g).unwrap();
    if t.triangles().is_empty() {
        return;
    }
    let dedges: HashSet<(u64, u64, u64, u64)> = t
        .triangles()
        .iter()
        .flat_map(|tr| {
            (0..3).map(|k| {
                edge_key(
                    t.points()[tr[k] as usize],
                    t.points()[tr[(k + 1) % 3] as usize],
                )
            })
        })
        .collect();
    let convex = convex_hull(&g).unwrap();
    for parameter in [
        MaxEdgeLength(0.0),
        MaxEdgeLength(3.0),
        MaxEdgeLength(1e9),
        MaxEdgeLengthRatio(0.0),
        MaxEdgeLengthRatio(0.4),
        MaxEdgeLengthRatio(1.0),
    ] {
        let hull = concave_hull(&g, parameter).unwrap();
        let ring = ring_of(&hull);
        let Geometry::Polygon(ref poly) = hull else {
            unreachable!()
        };
        // Containment, exact, through the unprepared locate.
        for &(x, y) in pts {
            assert_ne!(
                locate(c(x, y), poly).unwrap(),
                Location::Exterior,
                "{parameter:?}: ({x}, {y}) outside"
            );
        }
        // Every ring edge is a Delaunay edge.
        for k in 0..ring.len() {
            let e = edge_key(ring[k], ring[(k + 1) % ring.len()]);
            assert!(
                dedges.contains(&e),
                "{parameter:?}: a ring edge is not a Delaunay edge"
            );
        }
        // Exact topological equality with the convex hull where everything
        // is kept; byte equality when the hull has no collinear border
        // points (then both rings are corner-only, lex-smallest first).
        if matches!(parameter, MaxEdgeLengthRatio(r) if r == 1.0) {
            assert!(
                relate(&hull, &convex).unwrap().equals(),
                "ratio 1 is not the convex hull"
            );
            if ring.len() == ring_of(&convex).len() {
                assert_eq!(
                    encode(&hull),
                    encode(&convex),
                    "ratio 1 differs from convex_hull"
                );
            }
        }
        // The fixed point: at threshold 0 no allowed removal remains.
        if matches!(parameter, MaxEdgeLength(l) if l == 0.0) {
            assert_fixed_point(&t, &ring);
        }
    }
}

fn edge_key(p: Coord2d, q: Coord2d) -> (u64, u64, u64, u64) {
    let a = (p.x.to_bits(), p.y.to_bits());
    let b = (q.x.to_bits(), q.y.to_bits());
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    (lo.0, lo.1, hi.0, hi.1)
}

/// No allowed removal remains: for each ring edge, the living triangle on
/// its left (the unique Delaunay triangle that traverses it in ring
/// direction) is blocked: more than one of its edges is on the ring, or
/// its opposite vertex is a ring vertex.
fn assert_fixed_point(t: &Triangulation, ring: &[Coord2d]) {
    let bits = |p: Coord2d| (p.x.to_bits(), p.y.to_bits());
    let mut by_directed: HashMap<_, usize> = HashMap::new();
    for (i, tr) in t.triangles().iter().enumerate() {
        for k in 0..3 {
            by_directed.insert(
                (
                    bits(t.points()[tr[k] as usize]),
                    bits(t.points()[tr[(k + 1) % 3] as usize]),
                ),
                i,
            );
        }
    }
    let ring_vertices: HashSet<_> = ring.iter().map(|&p| bits(p)).collect();
    let ring_edges: HashSet<_> = (0..ring.len())
        .map(|k| edge_key(ring[k], ring[(k + 1) % ring.len()]))
        .collect();
    for k in 0..ring.len() {
        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
        let tri = t.triangles()[by_directed[&(bits(a), bits(b))]];
        let coords = tr_coords(t, tri);
        let on_ring = (0..3)
            .filter(|&k| ring_edges.contains(&edge_key(coords[k], coords[(k + 1) % 3])))
            .count();
        if on_ring == 1 {
            let apex = (0..3)
                .map(|k| coords[k])
                .find(|&p| bits(p) != bits(a) && bits(p) != bits(b))
                .unwrap();
            assert!(
                ring_vertices.contains(&bits(apex)),
                "an allowed removal remains at threshold 0"
            );
        }
    }
}

fn tr_coords(t: &Triangulation, tr: [u32; 3]) -> [Coord2d; 3] {
    tr.map(|v| t.points()[v as usize])
}

#[test]
fn the_contract_holds_on_generated_sets() {
    let mut rng = Rng(0xC09C_0001);
    for case in 0..120 {
        let pts = generated_set(&mut rng, case);
        check_contract(&pts);
    }
}

/// The ratio form gives exactly the length form's result at the length
/// `shortest + r · (longest − shortest)` over the triangulation's edge
/// lengths, recomputed here independently.
#[test]
fn the_ratio_equals_the_length_it_maps_to() {
    let mut rng = Rng(0xC09C_0002);
    let mut compared = 0usize;
    for case in 0..300 {
        let pts = generated_set(&mut rng, case);
        let g = multipoint(&pts);
        let t = delaunay(&g).unwrap();
        if t.triangles().is_empty() {
            continue;
        }
        let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
        for l in t.edges().line_strings() {
            let s = l.coords().as_slice();
            let (dx, dy) = (s[1].x - s[0].x, s[1].y - s[0].y);
            let len = dx.mul_add(dx, dy * dy).sqrt();
            min = min.min(len);
            max = max.max(len);
        }
        for r in [0.1, 0.35, 0.6, 0.9] {
            let target = r * (max - min) + min;
            let by_ratio = concave_hull(&g, MaxEdgeLengthRatio(r)).unwrap();
            let by_length = concave_hull(&g, MaxEdgeLength(target)).unwrap();
            assert_eq!(
                encode(&by_ratio),
                encode(&by_length),
                "case {case}, ratio {r}"
            );
            compared += 1;
        }
    }
    assert!(compared >= 1000, "{compared} comparisons");
}

#[test]
fn permutations_give_byte_identical_results() {
    let mut rng = Rng(0xC09C_0003);
    for case in 0..24 {
        let pts = generated_set(&mut rng, case);
        let base = concave_hull(&multipoint(&pts), MaxEdgeLengthRatio(0.3)).unwrap();
        let mut shuffled = pts.clone();
        for _ in 0..10 {
            rng.shuffle(&mut shuffled);
            let again = concave_hull(&multipoint(&shuffled), MaxEdgeLengthRatio(0.3)).unwrap();
            assert_eq!(encode(&base), encode(&again), "case {case}");
        }
    }
}

#[test]
fn degenerate_input_follows_convex_hull() {
    let cases: &[&[(f64, f64)]] = &[
        &[],
        &[(2.0, 1.0)],
        &[(2.0, 1.0), (5.0, 7.0)],
        &[(0.0, 0.0), (1.0, 1.0), (3.0, 3.0), (9.0, 9.0)],
    ];
    for pts in cases {
        let g = multipoint(pts);
        let hull = concave_hull(&g, MaxEdgeLength(1.0)).unwrap();
        assert_eq!(encode(&hull), encode(&convex_hull(&g).unwrap()));
    }
}

#[test]
fn invalid_parameters_are_refused() {
    use crate::GeometryError;
    let g = multipoint(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)]);
    for parameter in [
        MaxEdgeLength(-1.0),
        MaxEdgeLength(f64::NAN),
        MaxEdgeLength(f64::INFINITY),
        MaxEdgeLengthRatio(-0.1),
        MaxEdgeLengthRatio(1.5),
        MaxEdgeLengthRatio(f64::NAN),
    ] {
        assert!(
            matches!(
                concave_hull(&g, parameter),
                Err(GeometryError::InvalidTolerance { .. })
            ),
            "{parameter:?} was not refused"
        );
    }
}

/// A far point's outermost triangle is an ear (two boundary edges), which
/// the removal rule never erodes, since removing it would strand the point
/// outside the hull. At every threshold the spike stays in.
#[test]
fn a_far_point_stays_in_the_hull() {
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            pts.push((f64::from(i), f64::from(j)));
        }
    }
    pts.push((1000.0, 2.0));
    let g = multipoint(&pts);
    for parameter in [
        MaxEdgeLength(0.0),
        MaxEdgeLength(2.0),
        MaxEdgeLengthRatio(0.0),
    ] {
        let hull = concave_hull(&g, parameter).unwrap();
        let Geometry::Polygon(ref poly) = hull else {
            panic!("not a polygon")
        };
        assert_ne!(locate(c(1000.0, 2.0), poly).unwrap(), Location::Exterior);
    }
}

/// A concave point set erodes strictly tighter than its convex hull: the
/// mouth of the C is longer than every edge along the shape, so a suitable
/// threshold removes it.
#[test]
fn a_concave_set_erodes_tighter_than_the_convex_hull() {
    // A filled C: the 7-by-7 grid with a notch cut from the right. (A
    // rim-only C stays convex: with every point on the boundary, each
    // candidate's apex is a boundary vertex and the pinch rule blocks every
    // removal, a case the contract test already covers.)
    let mut pts = Vec::new();
    for i in 0..=6i32 {
        for j in 0..=6i32 {
            if !(i >= 2 && (2..=4).contains(&j)) {
                pts.push((f64::from(i), f64::from(j)));
            }
        }
    }
    let g = multipoint(&pts);
    let concave = concave_hull(&g, MaxEdgeLength(1.5)).unwrap();
    let convex = convex_hull(&g).unwrap();
    let (Geometry::Polygon(a), Geometry::Polygon(b)) = (&concave, &convex) else {
        panic!("not polygons");
    };
    let (ca, cb) = (crate::ops::area(a), crate::ops::area(b));
    assert!(ca < 0.8 * cb, "no concavity: {ca} vs {cb}");
    assert!(!relate(&concave, &convex).unwrap().equals());
}
