// src/ops/hull.rs

//! [`convex_hull`] by Andrew's monotone chain.
//!
//! # Why it is exact
//!
//! The algorithm makes two kinds of decision: the order of the points, and
//! whether three points turn left. The order is the exact lexicographic
//! order of stored coordinates, and the turn is `orient2d` on stored
//! coordinates, exact for every finite input. A point stays on a chain only
//! on a strict left turn (`orient2d > 0`), so a point on the segment between
//! two hull vertices is not a hull vertex. Every hull vertex is
//! an input vertex, bit for bit; nothing is constructed.

use core::cmp::Ordering;

use arithmetica::Sign;

use crate::GeometryError;
use crate::exact::kernel::orient;
use crate::ops::basic::{coords_of, lex};
use crate::predicates::{Relatable, borrowed};
use crate::types::{
    Coord2d, CoordSeq, Geometry, LineString, LinearRing, MultiPoint, Point, Polygon,
};

/// One chain of the monotone chain: the points that turn strictly left, in
/// the order given.
fn chain<'a>(pts: impl Iterator<Item = &'a Coord2d>) -> Vec<Coord2d> {
    let mut h: Vec<Coord2d> = Vec::new();
    for &p in pts {
        while h.len() >= 2 && orient(h[h.len() - 2], h[h.len() - 1], p) != Sign::Positive {
            h.pop();
        }
        h.push(p);
    }
    h
}

/// The convex hull of every coordinate of `g`.
///
/// - the empty `MultiPoint` for a geometry with no coordinates;
/// - a `Point` for one distinct coordinate;
/// - a `LineString` of the two extreme points (lexicographically smallest
///   first) when every coordinate is collinear;
/// - otherwise a counterclockwise `Polygon` without holes, whose vertices
///   are the hull's corners: input vertices where the hull turns strictly
///   left, starting at the lexicographically smallest.
///
/// Distinct means unequal under `==`. See the module documentation for why
/// every decision is exact.
///
/// # Errors
///
/// None for valid input. The polygon is built through `Polygon::new`, which
/// re-checks the ring and its winding: a hull of points with a strict left
/// turn is a non-degenerate convex polygon, and the constructor confirms it
/// rather than the code assuming it. Its error, which cannot occur, is
/// returned rather than unwrapped.
pub fn convex_hull(g: &impl Relatable) -> Result<Geometry, GeometryError> {
    let mut pts = coords_of(borrowed(g));
    pts.sort_by(|a, b| lex(*a, *b));
    pts.dedup_by(|a, b| lex(*a, *b) == Ordering::Equal);
    match pts.len() {
        0 => return Ok(Geometry::MultiPoint(MultiPoint::empty())),
        1 => return Ok(Geometry::Point(Point::from_finite(pts[0]))),
        _ => {}
    }
    let mut lower = chain(pts.iter());
    let mut upper = chain(pts.iter().rev());
    lower.pop();
    upper.pop();
    lower.extend(upper);
    if lower.len() < 3 {
        // Every point is collinear: the two chains are the two extremes.
        let (first, last) = (pts[0], pts[pts.len() - 1]);
        return Ok(Geometry::LineString(LineString::from_valid(CoordSeq::new(
            vec![first, last],
        ))));
    }
    let ring = LinearRing::new(CoordSeq::new(lower))?;
    Ok(Geometry::Polygon(Polygon::new_simple(ring)?))
}
