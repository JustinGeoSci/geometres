// src/ops/shared.rs

//! [`shared_paths`]: the stretches two line work inputs have in common,
//! split by whether they run the same way.
//!
//! # Why no point is constructed
//!
//! Two segments share a stretch only when they are collinear, decided
//! exactly: `orient2d` is zero for both endpoints of one against the other.
//! On their common line every point is fixed by one coordinate (`x` on a
//! non-vertical line, `y` on a vertical one), so the four endpoints are
//! ordered along the line by comparing that stored coordinate, exactly. The
//! shared stretch runs from the larger of the two segments' smaller ends to
//! the smaller of their larger ends: both are segment endpoints, so both are
//! input vertices. A shared path is built from such stretches and the input
//! vertices where consecutive stretches meet, so every output coordinate is
//! an input coordinate and the operation constructs nothing. That property
//! places it in `geometres`.

use crate::exact::kernel::orient;
use crate::types::{Coord2d, CoordSeq, LineString, Lineal, MultiLineString};
use arithmetica::Sign;

/// The shared stretches of two line work inputs.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedPaths {
    /// Stretches along which both inputs run the same way, oriented as the
    /// first input runs.
    pub same_direction: MultiLineString,
    /// Stretches along which they run opposite ways, oriented as the first
    /// input runs.
    pub opposite_direction: MultiLineString,
}

/// The coordinate that orders points along a line: `x`, unless the line is
/// vertical (`a.x == b.x`), then `y`.
// Stored coordinates compared for identity, not computed values within a
// tolerance: the comparison is the exact decision (module documentation).
#[allow(clippy::float_cmp)]
fn along(a: Coord2d, b: Coord2d) -> fn(Coord2d) -> f64 {
    if a.x == b.x { |c| c.y } else { |c| c.x }
}

/// One shared stretch of a segment of the first input.
#[derive(Clone, Copy)]
struct Piece {
    /// The segment's index within its line string.
    seg: usize,
    /// Its ends, in the first input's direction.
    from: Coord2d,
    to: Coord2d,
    /// Whether it starts at the segment's start, and ends at its end.
    at_start: bool,
    at_end: bool,
}

/// The shared stretches of segment `p → q` with the segments of `b`, per
/// class (same direction, opposite), unioned and in `p → q` order.
// As in `along`: stored coordinates compared exactly.
#[allow(clippy::float_cmp)]
fn segment_pieces(seg: usize, p: Coord2d, q: Coord2d, b: &[LineString]) -> [Vec<Piece>; 2] {
    let k = along(p, q);
    let a_up = k(q) > k(p);
    // Intervals along the line, low end first, with the low and high points.
    let mut iv: [Vec<(Coord2d, Coord2d)>; 2] = [Vec::new(), Vec::new()];
    let (a_lo, a_hi) = if a_up { (p, q) } else { (q, p) };
    for l in b {
        for w in l.coords().as_slice().windows(2) {
            let (r, s) = (w[0], w[1]);
            if orient(p, q, r) != Sign::Zero || orient(p, q, s) != Sign::Zero {
                continue;
            }
            let b_up = k(s) > k(r);
            let (b_lo, b_hi) = if b_up { (r, s) } else { (s, r) };
            // The larger low end and the smaller high end; on a tie, the
            // first input's point, so its bits are kept.
            let lo = if k(b_lo) > k(a_lo) { b_lo } else { a_lo };
            let hi = if k(b_hi) < k(a_hi) { b_hi } else { a_hi };
            if k(lo) < k(hi) {
                iv[usize::from(a_up != b_up)].push((lo, hi));
            }
        }
    }
    iv.map(|mut v| {
        v.sort_by(|x, y| k(x.0).total_cmp(&k(y.0)));
        let mut merged: Vec<(Coord2d, Coord2d)> = Vec::new();
        for (lo, hi) in v {
            match merged.last_mut() {
                Some(last) if k(lo) <= k(last.1) => {
                    if k(hi) > k(last.1) {
                        last.1 = hi;
                    }
                }
                _ => merged.push((lo, hi)),
            }
        }
        if !a_up {
            merged.reverse();
        }
        merged
            .into_iter()
            .map(|(lo, hi)| {
                let (from, to) = if a_up { (lo, hi) } else { (hi, lo) };
                Piece {
                    seg,
                    from,
                    to,
                    at_start: k(from) == k(p),
                    at_end: k(to) == k(q),
                }
            })
            .collect()
    })
}

/// Joins one class's pieces of one line string where they meet at a vertex
/// of that line string, consecutive along it.
fn join(pieces: &[Piece], out: &mut Vec<LineString>) {
    let mut cur: Vec<Coord2d> = Vec::new();
    let mut last: Option<Piece> = None;
    for &pc in pieces {
        let continues = last.is_some_and(|l| l.at_end && pc.at_start && pc.seg == l.seg + 1);
        if !continues && !cur.is_empty() {
            out.push(LineString::from_valid(CoordSeq::new(core::mem::take(
                &mut cur,
            ))));
        }
        if cur.is_empty() {
            cur.push(pc.from);
        }
        cur.push(pc.to);
        last = Some(pc);
    }
    if !cur.is_empty() {
        // Every piece has distinct ends, and joined pieces meet at the
        // vertex the first ends on and the second starts on, written once.
        out.push(LineString::from_valid(CoordSeq::new(cur)));
    }
}

/// The paths `a` and `b` share, split by whether they run the same way.
///
/// For each segment of `a`, the stretches it shares with segments of `b`
/// (collinear and overlapping in more than a point, decided exactly, see the
/// module documentation) are found, unioned where `b` overlaps itself, and
/// classed by whether `b` runs the same way along them. Stretches of one
/// class that meet at a vertex of `a`, consecutive along `a`, are joined
/// into one line string. Output lines are oriented as `a` runs, and listed
/// in `a`'s order. A single shared point is not a path and is not reported.
///
/// The inputs are taken as line strings, not as point sets: where `a` runs
/// over the same stretch more than once, each pass is reported, with the
/// direction `b` has there.
///
/// Every output coordinate is a vertex of `a` or of `b`. The cost is one
/// exact collinearity test per pair of segments.
#[must_use]
pub fn shared_paths(a: &impl Lineal, b: &impl Lineal) -> SharedPaths {
    let bl = b.components();
    let mut same = Vec::new();
    let mut opposite = Vec::new();
    for l in a.components() {
        let mut classes: [Vec<Piece>; 2] = [Vec::new(), Vec::new()];
        for (i, w) in l.coords().as_slice().windows(2).enumerate() {
            let [s, o] = segment_pieces(i, w[0], w[1], bl);
            classes[0].extend(s);
            classes[1].extend(o);
        }
        join(&classes[0], &mut same);
        join(&classes[1], &mut opposite);
    }
    SharedPaths {
        same_direction: MultiLineString::new(same),
        opposite_direction: MultiLineString::new(opposite),
    }
}
