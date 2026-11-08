// tests/operations.rs

//! The operations against properties their contracts state exactly, over
//! generated geometry.
//!
//! - Coverage: a grid of cells is a valid coverage; its union is one
//!   polygon of the grid's area, with a hole where a cell is missing; the
//!   noded grid lines polygonize into exactly the cells.
//! - Hulls: the convex hull covers every input point and its vertices are
//!   input points; a concave hull, when returned, covers every input point.
//! - Delaunay: the triangulation depends only on the point set, so a
//!   shuffled input gives the same triangles.
//! - Simplification: every result returned is valid; a coverage simplified
//!   as a coverage is still a valid coverage.
//! - The basic operations: `reverse` twice is the identity and keeps the
//!   point set, `normalize` is idempotent and keeps the point set, and
//!   `unique_points` is the set of stored coordinates.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss
)]

mod common;

use common::{Rng, coords, geometries, hole, jagged, line, points, rect};
use geometres::ops::{
    ConcaveHullParameter, area, concave_hull, convex_hull, coverage_simplify, coverage_union,
    coverage_validate, delaunay, is_valid, normalize, polygonize_noded, reverse,
    simplify_preserve_topology, simplify_vw, unique_points,
};
use geometres::predicates::{covers, equals};
use geometres::types::{Coord2d, MultiLineString, MultiPoint, Polygon};
use std::collections::BTreeSet;

/// The cells of an `n × n` grid of unit squares at `(ox, oy)`, skipping
/// the cells in `missing`.
fn grid(n: u32, ox: f64, oy: f64, missing: &[(u32, u32)]) -> Vec<Polygon> {
    let mut out = Vec::new();
    for j in 0..n {
        for i in 0..n {
            if !missing.contains(&(i, j)) {
                let (x, y) = (ox + f64::from(i), oy + f64::from(j));
                out.push(rect(x, y, x + 1.0, y + 1.0));
            }
        }
    }
    out
}

#[test]
fn a_grid_is_a_valid_coverage_and_unions_to_one_polygon() {
    let cells = grid(8, 100.0, -50.0, &[]);
    assert!(coverage_validate(&cells).is_empty());
    let u = coverage_union(&cells).unwrap();
    assert_eq!(u.polygons().len(), 1);
    assert_eq!(area(&u.polygons()[0]), 64.0);
    assert!(u.polygons()[0].holes().is_empty());
}

#[test]
fn a_missing_cell_becomes_a_hole_of_the_union() {
    let cells = grid(5, 0.0, 0.0, &[(2, 2)]);
    assert!(coverage_validate(&cells).is_empty());
    let u = coverage_union(&cells).unwrap();
    assert_eq!(u.polygons().len(), 1);
    let p = &u.polygons()[0];
    assert_eq!(p.holes().len(), 1);
    assert_eq!(area(p), 24.0);
    let want = Polygon::new(
        rect(0.0, 0.0, 5.0, 5.0).exterior().clone(),
        vec![hole(2.0, 2.0, 3.0, 3.0)],
    )
    .unwrap();
    assert!(equals(p, &want).unwrap());
}

#[test]
fn overlapping_cells_are_not_a_valid_coverage() {
    let mut cells = grid(3, 0.0, 0.0, &[]);
    cells.push(rect(0.5, 0.5, 1.5, 1.5));
    assert!(!coverage_validate(&cells).is_empty());
}

#[test]
fn noded_grid_lines_polygonize_into_the_cells() {
    // Every unit segment of a 6 × 6 grid, so the lines are noded.
    let n = 6;
    let mut segs = Vec::new();
    for k in 0..=n {
        for i in 0..n {
            let (k, i) = (f64::from(k), f64::from(i));
            segs.push(line(&[(i, k), (i + 1.0, k)]));
            segs.push(line(&[(k, i), (k, i + 1.0)]));
        }
    }
    let out = polygonize_noded(&MultiLineString::new(segs)).unwrap();
    assert_eq!(out.polygons.len(), 36);
    assert!(out.dangles.line_strings().is_empty());
    assert!(out.cut_edges.line_strings().is_empty());
    let mut got: Vec<f64> = out.polygons.iter().map(area).collect();
    got.dedup();
    assert_eq!(got, [1.0]);
}

fn input_points(mp: &MultiPoint) -> BTreeSet<(u64, u64)> {
    mp.points()
        .iter()
        .map(|p| (p.coord().x.to_bits(), p.coord().y.to_bits()))
        .collect()
}

#[test]
fn the_convex_hull_covers_the_points_and_uses_only_input_points() {
    let mut rng = Rng::new(0x5EED_0401);
    for k in 0..300 {
        let mp = points(&mut rng, 3 + k % 40, -30, 30);
        let hull = convex_hull(&mp).unwrap();
        assert!(covers(&hull, &mp).unwrap(), "{mp:?}");
        let input = input_points(&mp);
        for v in coords(&hull) {
            assert!(
                input.contains(&(v.x.to_bits(), v.y.to_bits())),
                "{v:?} is not an input point"
            );
        }
    }
}

#[test]
fn a_concave_hull_covers_the_points() {
    let mut rng = Rng::new(0x5EED_0402);
    let mut returned = 0;
    for k in 0..200 {
        let mp = points(&mut rng, 4 + k % 60, -40, 40);
        for r in [0.0, 0.3, 1.0] {
            if let Ok(h) = concave_hull(&mp, ConcaveHullParameter::MaxEdgeLengthRatio(r)) {
                assert!(covers(&h, &mp).unwrap(), "ratio {r}: {mp:?}");
                returned += 1;
            }
        }
    }
    assert!(returned > 300, "{returned}");
}

/// The triangles as sets of corner coordinates, independent of order.
fn triangle_set(pts: &[Coord2d], tris: &[[u32; 3]]) -> BTreeSet<[(u64, u64); 3]> {
    tris.iter()
        .map(|t| {
            let mut corners = t.map(|i| {
                let p = pts[i as usize];
                (p.x.to_bits(), p.y.to_bits())
            });
            corners.sort_unstable();
            corners
        })
        .collect()
}

#[test]
fn the_delaunay_triangulation_depends_only_on_the_point_set() {
    let mut rng = Rng::new(0x5EED_0403);
    for k in 0..200 {
        // Lattice points, so cocircular quadruples are common and the tie
        // rule is exercised.
        let mp = points(&mut rng, 4 + k % 50, 0, 8);
        let t1 = delaunay(&mp).unwrap();
        let mut shuffled = mp.points().to_vec();
        rng.shuffle(&mut shuffled);
        let t2 = delaunay(&MultiPoint::new(shuffled)).unwrap();
        assert_eq!(
            triangle_set(t1.points(), t1.triangles()),
            triangle_set(t2.points(), t2.triangles()),
            "{mp:?}"
        );
    }
}

#[test]
fn simplification_returns_only_valid_results() {
    let mut rng = Rng::new(0x5EED_0404);
    let mut returned = 0;
    for _ in 0..300 {
        let (x, y, n) = (
            rng.int(-100, 100),
            rng.int(-100, 100),
            3 + rng.below(40) as usize,
        );
        let p = jagged(&mut rng, x, y, n, 30);
        for t in [0.5, 2.0, 8.0] {
            if let Ok(q) = p.simplify(t) {
                assert!(
                    is_valid(&q),
                    "simplify({t}) returned an invalid polygon: {p:?}"
                );
                returned += 1;
            }
            if let Ok(q) = simplify_vw(&p, t * t) {
                assert!(
                    is_valid(&q),
                    "simplify_vw({}) returned an invalid polygon: {p:?}",
                    t * t
                );
                returned += 1;
            }
            if let Ok(q) = simplify_preserve_topology(&p, t) {
                assert!(is_valid(&q), "simplify_preserve_topology({t}): {p:?}");
                returned += 1;
            }
        }
    }
    assert!(returned > 1_500, "{returned}");
}

#[test]
fn a_simplified_coverage_is_still_a_coverage() {
    // A row of lots whose shared edges zigzag, so simplification has
    // something to remove and must remove it from both sides at once.
    let mut rng = Rng::new(0x5EED_0405);
    let xs: Vec<f64> = (0..=6).map(|i| f64::from(i) * 10.0).collect();
    let tops: Vec<Vec<f64>> = (0..xs.len())
        .map(|_| (0..=10).map(|_| rng.int(0, 2) * 0.5).collect())
        .collect();
    let mut lots = Vec::new();
    for w in 0..xs.len() - 1 {
        // The left edge of lot w is the right edge of lot w − 1: the same
        // zigzag vertices, in the other direction.
        let edge = |i: usize| -> Vec<(f64, f64)> {
            (0..=10)
                .map(|k| (xs[i] + tops[i][k], f64::from(k as u32) * 3.0))
                .collect()
        };
        let mut pts = edge(w + 1);
        let mut left = edge(w);
        left.reverse();
        pts.extend(left);
        lots.push(Polygon::new_simple(common::ring(&pts)).unwrap());
    }
    assert!(
        coverage_validate(&lots).is_empty(),
        "the generated lots are a coverage"
    );
    let mut returned = 0;
    for t in [0.25, 1.0, 4.0] {
        if let Ok(simplified) = coverage_simplify(&lots, t) {
            assert!(coverage_validate(&simplified).is_empty(), "threshold {t}");
            assert_eq!(simplified.len(), lots.len());
            let before: usize = lots.iter().map(|p| p.exterior().num_vertices()).sum();
            let after: usize = simplified.iter().map(|p| p.exterior().num_vertices()).sum();
            assert!(after <= before);
            returned += 1;
        }
    }
    assert!(returned > 0, "no threshold was simplified");
}

#[test]
fn normalize_is_idempotent_and_keeps_the_point_set() {
    for g in geometries(0x5EED_0406, 140) {
        let n = normalize(&g);
        assert_eq!(normalize(&n), n);
        assert!(equals(&g, &n).unwrap(), "{g:?}");
    }
}

#[test]
fn reverse_twice_is_the_identity_and_keeps_the_point_set() {
    let mut rng = Rng::new(0x5EED_0407);
    for _ in 0..200 {
        let n = 2 + rng.below(10) as usize;
        let l = common::walk(&mut rng, 0.0, 0.0, n);
        let r = reverse(&l);
        assert_eq!(reverse(&r), l);
        assert_eq!(
            r.coords().as_slice()[0],
            l.coords().as_slice()[l.coords().len() - 1]
        );
        assert!(equals(&l, &r).unwrap(), "{l:?}");
        let m = MultiLineString::new(vec![l.clone(), r.clone()]);
        assert_eq!(reverse(&reverse(&m)), m);
    }
}

#[test]
fn unique_points_are_the_stored_coordinates() {
    for g in geometries(0x5EED_0408, 140) {
        let want: BTreeSet<(u64, u64)> = coords(&g)
            .iter()
            .map(|p| (p.x.to_bits(), p.y.to_bits()))
            .collect();
        let got = unique_points(&g);
        let got_set: BTreeSet<(u64, u64)> = got
            .points()
            .iter()
            .map(|p| (p.coord().x.to_bits(), p.coord().y.to_bits()))
            .collect();
        assert_eq!(
            got.points().len(),
            got_set.len(),
            "a point is repeated: {g:?}"
        );
        assert_eq!(got_set, want, "{g:?}");
    }
}
