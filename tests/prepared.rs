// tests/prepared.rs

//! A prepared mask answers every query exactly as the unprepared operation
//! does.
//!
//! Preparation changes the cost of a query, never its result, so each
//! prepared answer is compared with the unprepared one for equality: the
//! matrix cell for cell, every predicate, the location of a point, a row of
//! points by `locate_on_line`, `dwithin` and `distance`. The masks are
//! generated polygons and multipolygons with holes, small (under the
//! indexing threshold) and large, and the queries are generated geometries
//! of every type that meet the masks in every way: crossing, along shared
//! edges, at single vertices and not at all.

#![allow(clippy::unwrap_used, clippy::float_cmp, clippy::cast_precision_loss)]

mod common;

use common::{Rng, c, geometries, hole, jagged, point, ring};
use geometres::ops::{distance, dwithin};
use geometres::predicates::relate;
use geometres::prepared::PreparedPolygonal;
use geometres::topology::label::{Dimension, Location};
use geometres::types::{Geometry, MultiPolygon, Polygon};

/// The masks: a small polygon, a large jagged polygon with holes, and a
/// multipolygon whose parts touch at a vertex.
fn masks() -> Vec<Geometry> {
    let mut rng = Rng::new(0x5EED_0101);
    let small = jagged(&mut rng, 2.0, 3.0, 6, 9);
    // 200 top vertices: well past the size at which the mask is indexed.
    let mut top: Vec<(f64, f64)> = vec![(0.0, 0.0), (25.0, 0.0)];
    for i in (0..=200).rev() {
        top.push((f64::from(i) * 0.125, 12.0 + rng.int(0, 12) * 0.5));
    }
    let big = Polygon::new(
        ring(&top),
        vec![
            hole(2.0, 2.0, 6.0, 5.0),
            hole(10.0, 1.0, 11.0, 9.0),
            hole(15.0, 4.0, 22.0, 6.0),
        ],
    )
    .unwrap();
    let a = jagged(&mut rng, 0.0, 0.0, 8, 10);
    // The second part's lower left corner is the first part's lower right.
    let b = Polygon::new_simple(ring(&[(8.0, 0.0), (20.0, -6.0), (20.0, 6.0)])).unwrap();
    let multi = MultiPolygon::new(vec![a, b]).unwrap();
    vec![
        Geometry::Polygon(small),
        Geometry::Polygon(big),
        Geometry::MultiPolygon(multi),
    ]
}

fn prepare(mask: &Geometry) -> PreparedPolygonal {
    match mask {
        Geometry::Polygon(p) => PreparedPolygonal::new(p).unwrap(),
        Geometry::MultiPolygon(m) => PreparedPolygonal::new(m).unwrap(),
        other => panic!("not a mask: {other:?}"),
    }
}

/// The location of `p` read from the unprepared matrix of the point
/// against the mask.
fn unprepared_location(mask: &Geometry, x: f64, y: f64) -> Location {
    let m = relate(&Geometry::Point(point(x, y)), mask)
        .unwrap()
        .matrix();
    let cells = m.cells();
    if cells[0] != Dimension::Empty {
        Location::Interior
    } else if cells[1] != Dimension::Empty {
        Location::Boundary
    } else {
        Location::Exterior
    }
}

#[test]
fn prepared_relate_and_predicates_equal_the_unprepared() {
    let queries = geometries(0x5EED_0102, 140);
    for mask in masks() {
        let pm = prepare(&mask);
        for g in &queries {
            let r = relate(&mask, g).unwrap();
            let ctx = || format!("mask {mask:?}\nquery {g:?}");
            assert_eq!(
                pm.relate(g).unwrap().cells(),
                r.matrix().cells(),
                "{}",
                ctx()
            );
            assert_eq!(pm.intersects(g).unwrap(), r.intersects(), "{}", ctx());
            assert_eq!(pm.disjoint(g).unwrap(), r.disjoint(), "{}", ctx());
            assert_eq!(pm.covers(g).unwrap(), r.covers(), "{}", ctx());
            assert_eq!(pm.contains(g).unwrap(), r.contains(), "{}", ctx());
            assert_eq!(pm.within(g).unwrap(), r.within(), "{}", ctx());
        }
    }
}

#[test]
fn prepared_distance_and_dwithin_equal_the_unprepared() {
    let queries = geometries(0x5EED_0103, 70);
    for mask in masks() {
        let pm = prepare(&mask);
        for g in &queries {
            let d = distance(&mask, g).unwrap();
            assert_eq!(
                pm.distance(g).unwrap().to_bits(),
                d.to_bits(),
                "{mask:?}\n{g:?}"
            );
            for r in [0.0, 0.5, 1.0, 3.0, d] {
                assert_eq!(
                    pm.dwithin(g, r).unwrap(),
                    dwithin(&mask, g, r).unwrap(),
                    "{r}: {g:?}"
                );
            }
        }
    }
}

#[test]
fn prepared_locate_equals_the_unprepared_location() {
    for mask in masks() {
        let pm = prepare(&mask);
        // Every half-integer and integer point of the masks' region: on
        // vertices, along edges and inside holes.
        for i in -4..=60 {
            for j in -16..=60 {
                let (x, y) = (f64::from(i) * 0.5, f64::from(j) * 0.5);
                assert_eq!(
                    pm.locate(c(x, y)),
                    unprepared_location(&mask, x, y),
                    "({x}, {y})"
                );
            }
        }
    }
}

#[test]
fn locate_on_line_equals_locate_point_by_point() {
    let mut rng = Rng::new(0x5EED_0104);
    for mask in masks() {
        let pm = prepare(&mask);
        for j in -16..=60 {
            let y = f64::from(j) * 0.5;
            // Every quarter-integer x, then the same values shuffled with
            // repeats: the answers come back in input order.
            let mut xs: Vec<f64> = (-8..=120).map(|i| f64::from(i) * 0.25).collect();
            let row = pm.locate_on_line(y, &xs).unwrap();
            for (k, &x) in xs.iter().enumerate() {
                assert_eq!(row[k], pm.locate(c(x, y)), "({x}, {y})");
            }
            xs.extend_from_within(..10);
            rng.shuffle(&mut xs);
            let row = pm.locate_on_line(y, &xs).unwrap();
            for (k, &x) in xs.iter().enumerate() {
                assert_eq!(row[k], pm.locate(c(x, y)), "({x}, {y}) shuffled");
            }
        }
    }
}

#[test]
fn locate_on_line_refuses_non_finite_input() {
    let pm = prepare(&masks()[0]);
    assert!(pm.locate_on_line(f64::NAN, &[0.0]).is_err());
    assert!(pm.locate_on_line(0.0, &[0.0, f64::INFINITY]).is_err());
    assert!(pm.locate_on_line(0.0, &[]).unwrap().is_empty());
}
