// src/ops/frechet.rs

//! The discrete Fréchet distance between two lines, `frechet_distance`, and
//! its exact predicate, `frechet_within`.
//!
//! # The definition
//!
//! Following Eiter & Mannila (*Computing Discrete Fréchet Distance*,
//! Technical Report CD-TR 94/64, Christian Doppler Laboratory for Expert
//! Systems, TU Wien, 1994): for two vertex sequences `p₁ … pₙ` and
//! `q₁ … qₘ`, a **coupling** is a sequence of pairs `(pᵢ, qⱼ)` that starts
//! at `(p₁, q₁)`, ends at `(pₙ, qₘ)`, and at each step advances `i`, `j`,
//! or both by one. The discrete Fréchet distance is the smallest, over all
//! couplings, of the largest point-to-point distance along the coupling.
//! It is the definition JTS's `DiscreteFrechetDistance` compute
//! (without densification), on the lines' vertices.
//!
//! It is defined here for two [`LineString`]s; other operand types are
//! refused by the signature. Direction matters: a line and its reverse are
//! at Fréchet distance at least `|p₁ − pₙ|`, not 0, because every coupling
//! starts by pairing the two first vertices.
//!
//! # The predicate is exact
//!
//! [`frechet_within`]`(a, b, d)` is Eiter and Mannila's dynamic program on
//! the `n × m` grid of pairs, with `arithmetica`'s exact
//! [`compare_point_distance`](arithmetica::compare_point_distance) deciding
//! each pair. A pair is **reachable** when its two points are within `d` and
//! it is the first pair or is reached from the pair before it in `i`, in
//! `j`, or in both, the three steps a coupling may take. Some coupling has
//! every pair within `d` exactly when the last pair is reachable, so the
//! answer is exact. One row is kept at a time, over the shorter line
//! (swapping the lines transposes the grid and maps the three steps onto
//! themselves): `O(n·m)` time, `O(min(n, m))` memory. A row with no
//! reachable pair ends the search, since every later row is reached only
//! through it. A pair's predicate is evaluated only when the pair could be
//! reached.
//!
//! # The quantity: an exact zero, and a bound
//!
//! [`frechet_distance`] is **exactly `0.0` if and only if the two lines
//! have the same vertices in the same order.** A valid [`LineString`] has no
//! two consecutive equal vertices, which makes the proof one sentence: a
//! coupling at distance 0 pairs equal points throughout, so a step that
//! advanced only `i` would pair `pᵢ` and `pᵢ₊₁` with the same `qⱼ` and make
//! them equal, which is impossible. So every step is diagonal, `n = m`, and
//! `pₖ = qₖ` for every `k`; the converse is the diagonal coupling.
//!
//! Otherwise the program runs again in `f64`, with `C(i, j) = max(|pᵢ − qⱼ|,
//! min(C(i−1, j), C(i, j−1), C(i−1, j−1)))`, and returns `C(n, m)`. Each
//! point-to-point distance is within `δ = c·u·M` of its true value (`c = 64`,
//! `u = 2⁻⁵³`, `M` the largest coordinate magnitude in either line; the bound
//! of [`distance`](super::distance()), whose module documentation derives it),
//! and **taking a maximum or a minimum of values each within `δ` of their true
//! values gives a value within `δ` of the true maximum or minimum** (the short
//! proof is in the module documentation of
//! [`hausdorff_distance`](super::hausdorff_distance): with `Xₖ` the largest
//! true value, `max x ≥ xₖ ≥ Xₖ − δ`, and `max x = xⱼ ≤ Xⱼ + δ ≤ max X + δ`).
//! The recurrence is built from maxima and minima only, so `C(n, m)` is within
//! `δ` of the true distance.
//!
//! The value is never a computed `0.0` for unequal lines: by the zero
//! argument every coupling of unequal lines pairs two distinct points, and
//! a computed distance between distinct points is positive (their
//! coordinate difference is nonzero, because subtraction of distinct
//! `f64`s never rounds to zero, and the scaled norm keeps it nonzero).
//!
//! **Consistency with the predicate:** if `frechet_within(a, b, d)` is
//! `false` then `frechet_distance(a, b) > d − δ`, and if `true` then
//! `frechet_distance(a, b) ≤ d + δ`. `frechet_distance(a, b) <= d` is not
//! the predicate.
//!
//! Neither function returns a coordinate.

use arithmetica::Sign;

use super::distance as dist;
use crate::GeometryError;
use crate::types::{Coord2d, LineString};

// ===========================================================================
// The public functions
// ===========================================================================

/// Whether the discrete Fréchet distance between the lines `a` and `b` is
/// at most `d`, decided **exactly**.
///
/// `true` when some coupling of the two vertex sequences (a monotone walk
/// from the first pair to the last; see the [module documentation](self))
/// keeps every coupled pair within `d`. Each pair is decided by
/// `arithmetica`'s exact point-to-point predicate. `O(n·m)` time,
/// `O(min(n, m))` memory.
///
/// ```
/// use geometres::ops::frechet_within;
/// use geometres::types::{Coord2d, CoordSeq, LineString};
///
/// let line = |v: &[(f64, f64)]| {
///     LineString::new(CoordSeq::new(v.iter().map(|&(x, y)| Coord2d::new(x, y)).collect()))
/// };
/// let a = line(&[(0.0, 0.0), (10.0, 0.0)])?;
/// let b = line(&[(10.0, 0.0), (0.0, 0.0)])?;
/// // Reversed: the first vertices are coupled, 10 apart.
/// assert!(frechet_within(&a, &b, 10.0)?);
/// assert!(!frechet_within(&a, &b, 10.0_f64.next_down())?);
/// # Ok::<(), geometres::GeometryError>(())
/// ```
///
/// # Errors
///
/// [`GeometryError::InvalidTolerance`] for a negative or non-finite `d`.
pub fn frechet_within(a: &LineString, b: &LineString, d: f64) -> Result<bool, GeometryError> {
    dist::check_distance(d)?;
    let (rows, cols) = by_length(a, b);
    let close = |p: Coord2d, q: Coord2d| {
        // `Err` is unreachable (validated coordinates, checked `d`); it is
        // read as "farther", never as a pass.
        arithmetica::compare_point_distance([p.x, p.y], [q.x, q.y], d)
            .is_ok_and(|s| s != Sign::Positive)
    };
    let mut row = vec![false; cols.len()];
    for (i, &p) in rows.iter().enumerate() {
        // `diag` is the previous row's entry one column to the left.
        let mut diag = false;
        let mut any = false;
        for (j, &q) in cols.iter().enumerate() {
            let up = row[j];
            let left = j > 0 && row[j - 1];
            let from = (i == 0 && j == 0) || up || left || diag;
            diag = up;
            row[j] = from && close(p, q);
            any |= row[j];
        }
        if !any {
            return Ok(false);
        }
    }
    Ok(row[cols.len() - 1])
}

/// The discrete Fréchet distance between the lines `a` and `b`: **exactly
/// `0.0` if and only if they have the same vertices in the same order**,
/// and otherwise an `f64` within `64·u·M` of the true value, with
/// `u = 2⁻⁵³` and `M` the largest coordinate magnitude in either line (the
/// definition and the derivation are in the [module documentation](self)).
///
/// For the question "is the Fréchet distance at most `d`?", use
/// [`frechet_within`], which is exact: if it is `false` then this value is
/// greater than `d − 64·u·M`, and if `true` then at most `d + 64·u·M`.
///
/// # Errors
///
/// None at present; the `Result` keeps the signature uniform with the
/// other distance functions.
pub fn frechet_distance(a: &LineString, b: &LineString) -> Result<f64, GeometryError> {
    if a.coords().as_slice() == b.coords().as_slice() {
        return Ok(0.0);
    }
    let (rows, cols) = by_length(a, b);
    let mut row = vec![f64::INFINITY; cols.len()];
    for (i, &p) in rows.iter().enumerate() {
        let mut diag = f64::INFINITY;
        for (j, &q) in cols.iter().enumerate() {
            let up = row[j];
            let left = if j > 0 { row[j - 1] } else { f64::INFINITY };
            let from = if i == 0 && j == 0 {
                0.0
            } else {
                up.min(left).min(diag)
            };
            diag = up;
            row[j] = point_distance(p, q).max(from);
        }
    }
    let value = row[cols.len() - 1];
    debug_assert!(
        value > 0.0,
        "unequal lines have a positive computed distance"
    );
    Ok(value)
}

// ===========================================================================
// Helpers
// ===========================================================================

/// The two vertex sequences, the longer first: the program keeps one row
/// over the shorter.
fn by_length<'g>(a: &'g LineString, b: &'g LineString) -> (&'g [Coord2d], &'g [Coord2d]) {
    let (p, q) = (a.coords().as_slice(), b.coords().as_slice());
    if p.len() >= q.len() { (p, q) } else { (q, p) }
}

/// The `f64` distance between two points: `distance`'s point-to-point
/// element distance, so the bound is the one derived there.
fn point_distance(p: Coord2d, q: Coord2d) -> f64 {
    dist::pair_distance(dist::Element::Point(p), dist::Element::Point(q))
}
