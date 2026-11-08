// tests/counters_thread_local.rs

//! The counters are per thread: two threads each run one counted
//! operation at the same time, and each thread's snapshot shows only its own
//! counts.
//!
//! Each operation is first run alone, on the test's own thread, to record its
//! count. Then two threads reset, meet at a barrier, run their operation,
//! meet again (so both operations have finished before either thread
//! reads), and snapshot. With process-global counters each snapshot would
//! hold the sum of both counts; with per-thread counters it holds exactly
//! its own.
//!
//! Built only with the `unstable-counters` feature. Runs under the default
//! parallel test runner.

#![cfg(feature = "unstable-counters")]
#![allow(clippy::unwrap_used)]

use std::sync::{Arc, Barrier};
use std::thread;

use geometres::counters::{self, CounterSnapshot};
use geometres::predicates::relate;
use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon};

fn polygon(pts: &[(f64, f64)]) -> Polygon {
    let ring = LinearRing::new(CoordSeq::new(
        pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
    ))
    .unwrap();
    Polygon::new_simple(ring).unwrap()
}

/// A regular-ish `n`-gon about `(cx, cy)` on a 1/64 grid.
fn ngon(cx: f64, cy: f64, r: f64, n: usize) -> Polygon {
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            (
                ((cx + r * t.cos()) * 64.0).round() / 64.0,
                ((cy + r * t.sin()) * 64.0).round() / 64.0,
            )
        })
        .collect();
    polygon(&pts)
}

/// `relate` of two overlapping polygons with `n` and `n + 1` vertices: the
/// areal analysis classifies every ring edge of each against the other.
fn counted(n: usize) -> CounterSnapshot {
    let a = ngon(0.0, 0.0, 10.0, n);
    let b = ngon(3.0, 1.0, 10.0, n + 1);
    counters::reset();
    relate(&a, &b).unwrap();
    counters::snapshot()
}

#[test]
fn counters_are_per_thread() {
    let small = counted(5);
    let large = counted(40);
    assert!(small.ring_edge_classifications > 0);
    assert!(
        large.ring_edge_classifications > small.ring_edge_classifications,
        "the two operations must count differently for the test to tell them apart"
    );

    let barrier = Arc::new(Barrier::new(2));
    let spawn = |n: usize| {
        let barrier = Arc::clone(&barrier);
        thread::spawn(move || {
            let a = ngon(0.0, 0.0, 10.0, n);
            let b = ngon(3.0, 1.0, 10.0, n + 1);
            counters::reset();
            barrier.wait();
            relate(&a, &b).unwrap();
            barrier.wait();
            counters::snapshot()
        })
    };
    let (t_small, t_large) = (spawn(5), spawn(40));
    assert_eq!(t_small.join().unwrap(), small);
    assert_eq!(t_large.join().unwrap(), large);
}

/// `reset` on one thread leaves another thread's counts alone.
#[test]
fn counters_reset_acts_on_the_calling_thread() {
    let mine = counted(12);
    thread::spawn(counters::reset).join().unwrap();
    assert_eq!(counters::snapshot(), mine);
}
