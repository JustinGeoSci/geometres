// tests/relate.rs

//! Properties of `relate` and the predicates that every exact
//! implementation must have, checked over generated geometry.
//!
//! The matrix is exact, so these hold without exception, not on average:
//!
//! - `relate(b, a)` is the transpose of `relate(a, b)`;
//! - every predicate, as a free function, as a `Relation` method and as an
//!   `IntersectionMatrix` method, is the same reading of one matrix;
//! - the matrix does not change when every coordinate is translated by a
//!   power of two or scaled by one, because on these integer coordinates
//!   both maps are exact and preserve every orientation;
//! - a point lies on a segment exactly when `arithmetica::orient2d` says it
//!   is collinear with the segment and it is inside the segment's box,
//!   including points one unit in the last place away from it.

#![allow(clippy::unwrap_used, clippy::float_cmp)]

mod common;

use arithmetica::{Sign, orient2d};
use common::{Rng, c, geometries, line, map, point};
use geometres::predicates::{
    contains, covered_by, covers, crosses, dimension, disjoint, equals, intersects, overlaps,
    relate, touches, within,
};
use geometres::topology::de9im::IntersectionMatrix;
use geometres::topology::label::Dimension;
use geometres::types::{Coord2d, Geometry};

/// The cells of the transpose of `m`, in `cells()` order (row-major,
/// interior, boundary, exterior).
fn transposed(m: &IntersectionMatrix) -> [Dimension; 9] {
    let a = m.cells();
    [a[0], a[3], a[6], a[1], a[4], a[7], a[2], a[5], a[8]]
}

/// Every ordered pair of `gs`.
fn pairs(gs: &[Geometry]) -> impl Iterator<Item = (&Geometry, &Geometry)> {
    gs.iter().flat_map(move |a| gs.iter().map(move |b| (a, b)))
}

#[test]
fn relate_of_the_swapped_pair_is_the_transpose() {
    let gs = geometries(0x5EED_0001, 42);
    let mut n = 0;
    for (a, b) in pairs(&gs) {
        let ab = relate(a, b).unwrap().matrix();
        let ba = relate(b, a).unwrap().matrix();
        assert_eq!(transposed(&ab), ba.cells(), "{a:?}\n{b:?}");
        n += 1;
    }
    assert_eq!(n, 42 * 42);
}

#[test]
fn every_predicate_reads_the_same_matrix() {
    let gs = geometries(0x5EED_0002, 42);
    for (a, b) in pairs(&gs) {
        let r = relate(a, b).unwrap();
        let m = r.matrix();
        let (da, db) = (dimension(a), dimension(b));
        let ctx = || format!("{a:?}\n{b:?}");
        assert_eq!(equals(a, b).unwrap(), r.equals(), "{}", ctx());
        assert_eq!(disjoint(a, b).unwrap(), r.disjoint(), "{}", ctx());
        assert_eq!(intersects(a, b).unwrap(), r.intersects(), "{}", ctx());
        assert_eq!(touches(a, b).unwrap(), r.touches(), "{}", ctx());
        assert_eq!(crosses(a, b).unwrap(), r.crosses(), "{}", ctx());
        assert_eq!(within(a, b).unwrap(), r.within(), "{}", ctx());
        assert_eq!(contains(a, b).unwrap(), r.contains(), "{}", ctx());
        assert_eq!(overlaps(a, b).unwrap(), r.overlaps(), "{}", ctx());
        assert_eq!(covers(a, b).unwrap(), r.covers(), "{}", ctx());
        assert_eq!(covered_by(a, b).unwrap(), r.covered_by(), "{}", ctx());

        assert_eq!(r.equals(), m.equals());
        assert_eq!(r.disjoint(), m.disjoint());
        assert_eq!(r.intersects(), m.intersects());
        assert_eq!(r.touches(), m.touches());
        assert_eq!(r.within(), m.within());
        assert_eq!(r.contains(), m.contains());
        assert_eq!(r.covers(), m.covers());
        assert_eq!(r.covered_by(), m.covered_by());
        assert_eq!(r.crosses(), m.crosses(da, db));
        assert_eq!(r.overlaps(), m.overlaps(da, db));

        // The relations between the predicates themselves.
        assert_eq!(r.intersects(), !r.disjoint(), "{}", ctx());
        assert_eq!(within(a, b).unwrap(), contains(b, a).unwrap(), "{}", ctx());
        assert_eq!(
            covered_by(a, b).unwrap(),
            covers(b, a).unwrap(),
            "{}",
            ctx()
        );
        assert_eq!(equals(a, b).unwrap(), equals(b, a).unwrap(), "{}", ctx());
        if r.within() {
            assert!(r.covered_by(), "{}", ctx());
        }
        if r.equals() {
            assert!(r.within() && r.contains(), "{}", ctx());
        }
    }
}

#[test]
fn every_geometry_equals_itself_and_intersects_itself() {
    for g in geometries(0x5EED_0003, 70) {
        let r = relate(&g, &g).unwrap();
        assert!(r.equals() && r.intersects() && !r.touches(), "{g:?}");
        assert!(r.covers() && r.covered_by(), "{g:?}");
    }
}

#[test]
fn the_matrix_survives_exact_translation_and_scaling() {
    let gs = geometries(0x5EED_0004, 35);
    // Integers below 2⁶ scaled by 4 and offset by 2³⁰ stay integers below
    // 2⁵³: every image coordinate is exact.
    let f = |p: Coord2d| c(p.x * 4.0 + 1_073_741_824.0, p.y * 4.0 - 536_870_912.0);
    let hs: Vec<Geometry> = gs.iter().map(|g| map(g, &f)).collect();
    for (i, a) in gs.iter().enumerate() {
        for (j, b) in gs.iter().enumerate() {
            let before = relate(a, b).unwrap().matrix();
            let after = relate(&hs[i], &hs[j]).unwrap().matrix();
            assert_eq!(before.cells(), after.cells(), "{a:?}\n{b:?}");
        }
    }
}

/// The exact answer for "is `p` on the closed segment `ab`".
fn on_segment(p: Coord2d, a: Coord2d, b: Coord2d) -> bool {
    orient2d([a.x, a.y], [b.x, b.y], [p.x, p.y]).unwrap() == Sign::Zero
        && a.x.min(b.x) <= p.x
        && p.x <= a.x.max(b.x)
        && a.y.min(b.y) <= p.y
        && p.y <= a.y.max(b.y)
}

#[test]
fn a_point_is_on_a_segment_exactly_when_orient2d_says_so() {
    let mut rng = Rng::new(0x5EED_0005);
    let (mut on, mut off) = (0, 0);
    for _ in 0..2_000 {
        // Large, non-integer endpoints: the region where an `f64`
        // cross product loses the sign.
        let a = c(rng.int(0, 1 << 20) + 0.5, rng.int(0, 1 << 20) * 0.25);
        let b = c(
            a.x + rng.int(1, 1 << 12) * 0.125,
            a.y + rng.int(-(1 << 12), 1 << 12) * 0.375,
        );
        let seg = Geometry::LineString(line(&[(a.x, a.y), (b.x, b.y)]));
        // A point by interpolation (usually a hair off the line), its
        // neighbors one unit in the last place away in `y`, and the exact
        // midpoint (on the line: the sums are exact at these magnitudes).
        let t = rng.int(1, 1023) / 1024.0;
        let p = c(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
        let mid = c(0.5 * (a.x + b.x), 0.5 * (a.y + b.y));
        let up = c(p.x, f64::from_bits(p.y.to_bits() + 1));
        for q in [p, up, mid, a, b] {
            let want = on_segment(q, a, b);
            let got = intersects(&Geometry::Point(point(q.x, q.y)), &seg).unwrap();
            assert_eq!(got, want, "point {q:?}, segment {a:?} → {b:?}");
            if want {
                on += 1;
            } else {
                off += 1;
            }
        }
    }
    // Both answers are exercised.
    assert!(on > 4_000 && off > 1_000, "on {on}, off {off}");
}
