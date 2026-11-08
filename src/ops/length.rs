// src/ops/length.rs

//! Arc length and perimeter.
//!
//! Each edge length is `f64::hypot`, which avoids the overflow and underflow
//! that `sqrt(dx*dx + dy*dy)` has for very large or very small differences,
//! and the edge lengths are added with a compensated sum. Neither step has a
//! proven error bound in this crate: the result is an approximation, and its
//! accuracy is measured in the test suite.

use crate::accum::Accum;
use crate::types::{Coord2d, LineString, LinearRing, Polygon};

/// Length of the edge `a → b`: `hypot` of the two coordinate differences.
///
/// **Approximate, with no certified property.** Two roundings stand between
/// the inputs and the result. Each difference `b.x − a.x` and `b.y − a.y` is
/// rounded once (exact only when the operands are within a factor of two of
/// each other, by Sterbenz's lemma). Then `f64::hypot` is applied, and its
/// accuracy is the platform math library's: Rust does not specify it, so no
/// ulp bound is claimed here. `hypot` is used because, in practice, it has
/// none of the spurious overflow or underflow of `sqrt(dx² + dy²)` for very large or very small
/// differences.
///
/// A non-finite input gives a non-finite result; the result is never
/// negative. [`length`], [`ring_perimeter`] and [`perimeter`] are compensated
/// sums of this function.
#[must_use]
pub fn edge_length(a: Coord2d, b: Coord2d) -> f64 {
    (b.x - a.x).hypot(b.y - a.y)
}

/// Total arc length of a `LineString`: the sum of the Euclidean distances
/// between consecutive vertices.
pub fn length(line: &LineString) -> f64 {
    let mut acc = Accum::new();
    for w in line.coords().as_slice().windows(2) {
        acc.add(edge_length(w[0], w[1]));
    }
    acc.value()
}

/// Perimeter of a `LinearRing`, including the closing edge.
pub fn ring_perimeter(ring: &LinearRing) -> f64 {
    let mut acc = Accum::new();
    for (a, b) in ring.edges() {
        acc.add(edge_length(a, b));
    }
    acc.value()
}

/// Perimeter of a `Polygon`'s exterior ring. Hole rings are not included;
/// add `ring_perimeter(hole)` for each of `polygon.holes()` for the total.
pub fn perimeter(polygon: &Polygon) -> f64 {
    ring_perimeter(polygon.exterior())
}
