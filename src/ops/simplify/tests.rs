// src/ops/simplify/tests.rs

//! Unit tests for `ops::simplify`: the RDP core on bare coordinate lists.
//! The vouch is tested through the public API in `tests/simplify.rs`.

#![allow(clippy::unwrap_used, clippy::float_cmp)]

use super::{segment_distance, simplify_chain, simplify_ring};
use crate::types::Coord2d;

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

#[test]
fn segment_distance_is_perpendicular_inside_and_endpoint_outside() {
    let (a, b) = (c(0.0, 0.0), c(10.0, 0.0));
    assert_eq!(segment_distance(c(5.0, 3.0), a, b), Some(3.0));
    assert_eq!(segment_distance(c(5.0, 0.0), a, b), Some(0.0));
    // Beyond the ends: distance to the nearer endpoint, a 3-4-5 triangle.
    assert_eq!(segment_distance(c(-3.0, 4.0), a, b), Some(5.0));
    assert_eq!(segment_distance(c(13.0, -4.0), a, b), Some(5.0));
    // Degenerate chord: distance to the point.
    assert_eq!(segment_distance(c(3.0, 4.0), a, a), Some(5.0));
    // Overflow in a product is reported, not read as a finite distance.
    assert_eq!(
        segment_distance(c(0.0, 1e200), c(-1e200, 0.0), c(1e200, 0.0)),
        None
    );
}

#[test]
fn collinear_interior_vertices_go_at_any_positive_tolerance() {
    let pts = [c(0.0, 0.0), c(1.0, 0.0), c(2.0, 0.0), c(3.0, 0.0)];
    assert_eq!(simplify_chain(&pts, 1e-300), vec![c(0.0, 0.0), c(3.0, 0.0)]);
}

#[test]
fn the_farthest_vertex_is_kept_and_the_chain_split_there() {
    let pts = [
        c(0.0, 0.0),
        c(1.0, 0.1),
        c(2.0, 2.0),
        c(3.0, 0.1),
        c(4.0, 0.0),
    ];
    // (1, 0.1) is 1.8/√8 ≈ 0.636 from the chord (0,0)-(2,2).
    assert_eq!(
        simplify_chain(&pts, 0.7),
        vec![c(0.0, 0.0), c(2.0, 2.0), c(4.0, 0.0)]
    );
    assert_eq!(simplify_chain(&pts, 0.6), pts.to_vec());
    assert_eq!(simplify_chain(&pts, 10.0), vec![c(0.0, 0.0), c(4.0, 0.0)]);
}

#[test]
fn a_vertex_exactly_at_the_tolerance_is_dropped() {
    // "Not greater than" drops it: the contract is "moves nothing by more
    // than the tolerance".
    let pts = [c(0.0, 0.0), c(1.0, 1.0), c(2.0, 0.0)];
    assert_eq!(simplify_chain(&pts, 1.0), vec![c(0.0, 0.0), c(2.0, 0.0)]);
    assert_eq!(simplify_chain(&pts, 0.999_999_999), pts.to_vec());
}

#[test]
fn ties_go_to_the_first_vertex_attaining_the_maximum() {
    // Two vertices at the same distance; the first splits, then the second
    // is on the chord of its own half and is dropped.
    let pts = [c(0.0, 0.0), c(1.0, 1.0), c(2.0, 1.0), c(3.0, 0.0)];
    // (2,1) is 1/√5 ≈ 0.447 from the chord (1,1)-(3,0) of its half.
    assert_eq!(
        simplify_chain(&pts, 0.4),
        vec![c(0.0, 0.0), c(1.0, 1.0), c(2.0, 1.0), c(3.0, 0.0)]
    );
    assert_eq!(
        simplify_chain(&pts, 0.5),
        vec![c(0.0, 0.0), c(1.0, 1.0), c(3.0, 0.0)]
    );
    let pts = [
        c(0.0, 0.0),
        c(1.0, 1.0),
        c(2.0, 1.0),
        c(3.0, 1.0),
        c(4.0, 0.0),
    ];
    let got = simplify_chain(&pts, 0.5);
    assert_eq!(
        got,
        vec![c(0.0, 0.0), c(1.0, 1.0), c(3.0, 1.0), c(4.0, 0.0)]
    );
}

#[test]
fn a_ring_anchors_at_vertex_zero_and_its_antipode() {
    // A square with a vertex in the middle of each edge: the midpoints go,
    // the corners stay, whichever vertex the ring starts at.
    let sq = [
        c(0.0, 0.0),
        c(1.0, 0.0),
        c(2.0, 0.0),
        c(2.0, 1.0),
        c(2.0, 2.0),
        c(1.0, 2.0),
        c(0.0, 2.0),
        c(0.0, 1.0),
    ];
    let corners = vec![c(0.0, 0.0), c(2.0, 0.0), c(2.0, 2.0), c(0.0, 2.0)];
    assert_eq!(simplify_ring(&sq, 0.1), corners);
    // Starting at a midpoint: the anchor is kept as a vertex of the ring.
    // It lies exactly on the chord between its neighbors, but RDP never
    // drops its anchor. That is the one asymmetry of anchoring, and the
    // vertex kept is the input's own rather than a constructed one.
    let rotated: Vec<Coord2d> = sq[1..].iter().chain(&sq[..1]).copied().collect();
    let got = simplify_ring(&rotated, 0.1);
    assert_eq!(got.len(), 5);
    assert_eq!(got[0], c(1.0, 0.0));
    for corner in &corners {
        assert!(got.contains(corner));
    }
}

#[test]
fn the_closing_edge_is_simplified_like_any_other() {
    // Vertices 0 and n-1 bracket a nearly straight closing stretch: both are
    // dropped in favor of the far corners. An open-chain port would have
    // kept them because they are "the ends".
    let ring = [
        c(1.0, 0.0),
        c(2.0, 0.0),
        c(2.0, 2.0),
        c(0.0, 2.0),
        c(0.0, 0.0),
        c(0.5, 0.0),
    ];
    let got = simplify_ring(&ring, 0.1);
    // Anchor (1,0) is kept by construction; (0.5,0) is between (0,0) and the
    // anchor on the closing arc and goes.
    assert!(!got.contains(&c(0.5, 0.0)));
    assert!(got.contains(&c(0.0, 0.0)) && got.contains(&c(2.0, 0.0)));
}

#[test]
fn a_ring_can_collapse_to_two_vertices() {
    // A thin triangle at a tolerance above its height: only the anchor and
    // its antipode survive. The vouch, not this function, refuses that.
    let tri = [c(0.0, 0.0), c(1.0, 0.1), c(2.0, 0.0)];
    assert_eq!(simplify_ring(&tri, 1.0), vec![c(0.0, 0.0), c(2.0, 0.0)]);
}

#[test]
fn a_long_monotone_chain_does_not_overflow_the_stack() {
    // Textbook recursive RDP recurses once per vertex on a monotone spiral;
    // the explicit stack does not care.
    let pts: Vec<Coord2d> = (0..200_000).map(f64::from).map(|t| c(t, t * t)).collect();
    let got = simplify_chain(&pts, 1e-9);
    assert_eq!(got.len(), pts.len());
}
