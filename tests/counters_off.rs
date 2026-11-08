// tests/counters_off.rs

//! Tests that the counters compile out.
//!
//! Built only **without** the `unstable-counters` feature. Three things are
//! checked:
//!
//! - this build has the feature off (so the run below is the one that
//!   matters);
//! - operations that carry `count!` sites run normally, so the empty
//!   expansion is sound;
//! - `geometres::counters` is unreachable: the `compile_fail` doctest on the
//!   hidden `geometres::__counters_absent_check` item names it and must fail
//!   to compile (`cargo test -p geometres --doc __counters_absent_check`).
//!
//! With the feature on this file compiles to nothing.

#![cfg(not(feature = "unstable-counters"))]
#![allow(clippy::unwrap_used)]

use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon};

fn square(x: f64) -> Polygon {
    let c = |a: f64, b: f64| Coord2d::new(a, b);
    Polygon::new_simple(
        LinearRing::new(CoordSeq::new(vec![
            c(x, 0.0),
            c(x + 1.0, 0.0),
            c(x + 1.0, 1.0),
            c(x, 1.0),
        ]))
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn this_build_has_the_counters_feature_off() {
    const { assert!(!cfg!(feature = "unstable-counters")) };
}

/// The areal DE-9IM analysis carries the `ring_edge_classifications` site.
#[test]
fn operations_with_counter_sites_run_without_the_counters() {
    let m = geometres::topology::de9im::IntersectionMatrix::compute(&square(0.0), &square(0.5));
    assert!(m.intersects());
}
