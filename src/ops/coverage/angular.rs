// src/ops/coverage/angular.rs

//! The counterclockwise order of directions at a vertex, decided exactly and
//! constructing nothing. The argument is in the parent module's
//! documentation.

use core::cmp::Ordering;

use arithmetica::Sign;

use crate::exact::kernel::orient;
use crate::types::Coord2d;

/// Which half-turn the direction from `v` to `a` lies in: `0` for the upper
/// half `[0, π)`, `1` for the lower half `[π, 2π)`.
///
/// Only signs of computed differences are read, and a correctly rounded
/// difference has the sign of the exact one (zero exactly when the
/// coordinates are equal). Comparisons with `0.0` treat `-0.0` as zero, as
/// a coordinate difference of zero must be.
fn half(v: Coord2d, a: Coord2d) -> u8 {
    let dx = a.x - v.x;
    let dy = a.y - v.y;
    u8::from(!(dy > 0.0 || (dy == 0.0 && dx > 0.0)))
}

/// The counterclockwise order, starting from the `+x` direction, of the
/// directions from `v` to `a` and from `v` to `b` (`a, b ≠ v`).
///
/// Different half-turns order by half. Within one half the directions span
/// less than a half turn, so `orient2d(v, a, b) > 0` means `b` is
/// counterclockwise of `a`. `Equal` means the same direction (collinear and
/// on the same side of `v`).
pub(crate) fn ccw_cmp(v: Coord2d, a: Coord2d, b: Coord2d) -> Ordering {
    let (ha, hb) = (half(v, a), half(v, b));
    if ha != hb {
        return ha.cmp(&hb);
    }
    match orient(v, a, b) {
        Sign::Positive => Ordering::Less,
        Sign::Negative => Ordering::Greater,
        Sign::Zero => Ordering::Equal,
    }
}
