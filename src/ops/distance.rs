// src/ops/distance.rs

//! `dwithin`, the exact distance predicate, and `distance`, the measurement
//! with a documented bound and a certified zero.
//!
//! # The two answers, and which to use
//!
//! [`dwithin`]`(a, b, d)` decides "is the distance between `a` and `b` at
//! most `d`?" **exactly**, for every pair of geometries and every finite
//! `d ≥ 0`. It is the answer to use for that question.
//!
//! [`distance`]`(a, b)` returns the distance as an `f64`. It is **exactly
//! `0.0` if and only if `a` and `b` intersect**; otherwise it is an
//! approximation with the absolute bound below. This is the position
//! [`area`](super::area) and [`length`](super::length) hold under this
//! crate's charter: a quantity whose sign (here: zero or not) is certified
//! and whose magnitude has a documented bound. Neither function returns a
//! coordinate.
//!
//! # The decomposition, and why it is exact
//!
//! Both functions rest on one decomposition of "the distance between two
//! geometries":
//!
//! 1. **If `a` and `b` intersect, the distance is 0.** Intersection is
//!    decided by the exact [`relate`](crate::predicates::relate) machinery,
//!    so this covers every kind of contact (a point inside a polygon, a
//!    line crossing an area, boundaries touching) rather than only contacts
//!    between boundary elements.
//! 2. **Otherwise the distance is attained between boundary elements:** the
//!    segments of lines and of polygon rings, and the isolated points. For
//!    puntal and lineal geometries every point lies on an element already.
//!    For an areal `a`, suppose the infimum were attained at a point `p`
//!    interior to `a` nearest to some `q` in `b`. Every point of `b` is
//!    outside `a` (they do not intersect), so `q` is outside, and the
//!    segment `pq` leaves `a` through its boundary at some `p′` with
//!    `|p′q| < |pq|`, a contradiction. So the nearest point of `a` lies on
//!    its boundary, which is made of ring segments, and symmetrically
//!    for `b`.
//! 3. **For two segments that do not intersect, the minimum distance is
//!    attained at an endpoint of one of them against the other:** the
//!    distance between a point moving along one segment and the other
//!    segment is a convex function of the position, so its minimum over the
//!    closed segment is at an endpoint unless the interiors meet, and
//!    step 1 has already excluded every meeting. A composition that skipped
//!    step 1 would be wrong exactly there: two crossing segments have
//!    distance 0 with all four endpoint distances positive.
//! 4. **Point against segment, and point against point,** are
//!    `arithmetica`'s exact predicates
//!    ([`compare_segment_distance`](arithmetica::compare_segment_distance),
//!    [`compare_point_distance`](arithmetica::compare_point_distance))
//!    directly.
//!
//! So `dwithin(a, b, d)` is `true` exactly when step 1 holds, or some pair
//! of boundary elements passes its exact predicate with this `d`. Every
//! decision along the way is exact.
//!
//! # Contact without the matrix
//!
//! The intersection test of step 1 does not need the DE-9IM matrix. Call
//! two elements in contact when one's endpoint is within distance 0 of the
//! other (the step 4 predicates with `d = 0`) or they are segments that
//! cross properly.
//! Then `a` and `b` intersect exactly when some pair of their elements is in
//! contact, or some connected component of one lies wholly in the interior
//! of the other, which is then areal; the component's first vertex shows
//! it.
//!
//! - Contact or a located vertex is a common point, so either gives an
//!   intersection.
//! - Conversely, let `x` be common to `a` and `b`. If it lies on an element
//!   of each, the two elements meet. Two meeting segments either cross
//!   properly or have an endpoint of one on the other. A segment meeting a
//!   point has the point on it. So some pair is in contact.
//! - Otherwise `x` is interior to an areal operand, say `a`, and lies in a
//!   component `C` of `b` none of whose elements meets an element of `a`.
//!   If `C`'s elements meet no element of `a`, then `C`'s boundary, which
//!   is connected per ring and does not touch `a`'s boundary, lies wholly
//!   inside or wholly outside `a`. For a point or a line, `C` itself is
//!   inside, so its first vertex is too. For a polygon part, the argument
//!   of `predicates::intersects` for two areas applies (an outer ring in
//!   the other's exterior forces a contradiction through the holes), so
//!   either `C`'s first exterior vertex is inside `a`, or a component of
//!   `a` lies inside `C` and its first vertex shows that.
//!
//! So `dwithin` tests element pairs with the step 3–4 predicates plus the
//! crossing test, and locates one vertex per component; `distance` does the
//! same at `d = 0` for its exact zero. Both give exactly what step 1
//! followed by the element scan gives. A collection operand still goes
//! through `intersects`, which also carries its refusal.
//!
//! # The bound on `distance`
//!
//! Write `u = 2⁻⁵³` and let `M` be the largest coordinate magnitude in
//! either operand. For non-intersecting operands,
//!
//! > **|`distance(a, b)` − true distance| ≤ max(c·u·M, 2⁻¹⁰⁷⁴), with c = 64,**
//!
//! valid whenever every coordinate difference is finite (`M ≤ 2^1022`
//! suffices); past that a difference overflows and the affected pair
//! contributes `+∞`, never a finite wrong value. The derivation, with
//! `γk = kε/(1 − kε)` and every distance at most `2√2·M`:
//!
//! - **`norm`** (the only square root): the squares and their sum cost
//!   `(1 + θ3)`, the square root halves that and adds one rounding, and an
//!   exact power-of-two prescale keeps everything in range, so the result
//!   is within relative `γ3`, an absolute error of `γ3·2√2·M < 9·u·M`.
//! - **A computed difference** `fl(x − y)` is within `ε·|x − y| ≤ 2u·M` of
//!   the true difference, and the distance is 1-Lipschitz in each endpoint,
//!   so the four differences feeding a pair cost roughly `4·(2√2)u·M`;
//!   counted generously, `12·u·M` across both endpoints of both elements.
//! - **The clamped projection** (point to segment): with `w = p − a`,
//!   `t = (w·u)/(u·u)` clamped to `[0, 1]`. The dot product's absolute
//!   error is at most `γ2·(|wₓuₓ| + |w_yu_y|) ≤ γ2·4M·|u|`, and `|u|²`'s is
//!   relative `γ3`, so the clamped `t̂` is within `γ4 + 5γ3·M/|u|` of the
//!   true minimizer; the distance is `|u|`-Lipschitz in `t`, so the cost is
//!   at most `γ4·|u| + 5γ3·M ≤ 25·u·M`. Forming the residual
//!   `w − t̂·u` costs two roundings per coordinate, under `9·u·M` through
//!   the norm. A segment-to-segment distance is the minimum of four such
//!   values, and **the minimum of values each within `δ` of their true
//!   values is within `δ` of the true minimum**: for every pair,
//!   `computed ≥ true_min − δ` and the minimizing pair's computed value is
//!   `≤ its true + δ`.
//!
//! The terms sum below `55·u·M`; `c = 64` is the published constant, and
//! the test suite's exact oracle measures the realized ratio against it on
//! more than 50,000 cases.
//!
//! **The zero guard** accounts for the `2⁻¹⁰⁷⁴` term. The certified zero comes
//! from the exact intersection test, not from the arithmetic: for operands
//! that do not intersect, every computed pair distance could still round to
//! `0.0` (a point one rounding off a segment, whose residual cancels).
//! Returning that `0.0` would assert contact the exact test has just
//! denied. So a computed `0.0` for non-intersecting operands is returned as
//! the smallest positive `f64`, `2⁻¹⁰⁷⁴`. The true distance is then at most
//! `c·u·M` (the computed value `0` is within that of it), so the returned
//! value is within `max(c·u·M, 2⁻¹⁰⁷⁴)` of the truth.
//!
//! # Pruning never discards the minimum
//!
//! Large operands are searched through a static box index over one side's
//! elements (through the crate's private `situs` seam, so no `situs` type
//! appears in a public signature). For `dwithin` each query window is an
//! element's box **grown by `d`**; for `distance`, by `r = (best + δ)` rounded
//! upward, with `δ = c·u·M`. Growing is the one internal construction here, and
//! the crate-private `grow_outward` helper moves every bound outward with
//! directed rounding, so a grown box contains the exact grown box and a query
//! can only gain candidates. A pair the search skips has an axis gap, and so a
//! true distance, greater than `r ≥ best + δ`; its computed value is at least
//! its true value minus `δ`, hence greater than `best`. Skipping it therefore
//! changes neither `dwithin`'s answer nor the computed minimum, and the indexed
//! forms return **exactly** what the exhaustive scan returns.
//!
//! In fact, for a window grown by a *single* subtraction the directed step
//! is provably redundant: rounding to nearest is monotone, so from
//! `box_lo ≥ lo − d` (reals) follows `box_lo ≥ fl(lo − d)`, and a pair
//! within `d` is never excluded even by a to-nearest window. The directed
//! rounding is kept because its safety argument is local, requiring no
//! reasoning about which expression feeds the bound, and because compound
//! radii (`best + δ`) do need the upward step on the sum.
//!
//! # Errors and empties
//!
//! `d` must be finite and at least 0
//! ([`GeometryError::InvalidTolerance`]). `dwithin` with an empty operand
//! is `false`: no point of one geometry is within anything of the other.
//! `distance` with an empty operand is
//! [`GeometryError::EmptyOperand`]: Returning 0 there, which reads as "they
//! touch", is a plausible wrong answer this crate does not give. Collections
//! follow the `relate` semantics: one whose polygonal leaves overlap is
//! refused ([`GeometryError::CollectionElementsOverlap`]), by both
//! functions, through step 1.

use core::ops::ControlFlow;

use arithmetica::Sign;

use crate::GeometryError;
use crate::index::BoxIndex;
use crate::predicates::{Borrowed, Relatable, borrowed, intersects};
use crate::types::{Coord2d, Envelope};

#[cfg(test)]
mod tests;

/// Above this many element pairs, the indexed search replaces the
/// exhaustive scan. The answers are identical (module documentation); the
/// threshold is only a cost choice.
pub(crate) const BRUTE_PAIRS: usize = 2_048;

/// `u = 2⁻⁵³`, the unit roundoff.
pub(crate) const U: f64 = f64::EPSILON / 2.0;

/// The constant of the `distance` bound `c·u·M` (module documentation).
pub(crate) const DIST_ERR_C: f64 = 64.0;

// ===========================================================================
// The public functions
// ===========================================================================

/// Whether the distance between `a` and `b` is at most `d`, decided **exactly**.
///
/// `true` when some point of `a` and some point of `b` lie within `d` of
/// each other, endpoints and boundaries included; `d = 0` asks whether the
/// geometries intersect. Decided exactly for every input: between boundary
/// elements by `arithmetica`'s distance predicates and the exact crossing
/// test, and by locating one vertex of each component in the other operand
/// (see the [module documentation](self) for the decomposition and its
/// proof). No distance is ever computed.
///
/// ```
/// use geometres::ops::dwithin;
/// use geometres::types::{Coord2d, Point};
///
/// let p = Point::new(Coord2d::new(0.0, 0.0))?;
/// let q = Point::new(Coord2d::new(3.0, 4.0))?;
/// assert!(dwithin(&p, &q, 5.0)?);
/// assert!(!dwithin(&p, &q, 5.0_f64.next_down())?);
/// # Ok::<(), geometres::GeometryError>(())
/// ```
///
/// An empty operand gives `false`.
///
/// # Errors
///
/// - [`GeometryError::InvalidTolerance`] for a negative or non-finite `d`.
/// - [`GeometryError::CollectionElementsOverlap`] for a collection operand
///   whose polygonal leaves overlap, as `relate` refuses it.
/// - [`GeometryError::TooManyEntries`] for more than `u32::MAX` elements on
///   an indexed side.
pub fn dwithin(a: &impl Relatable, b: &impl Relatable, d: f64) -> Result<bool, GeometryError> {
    check_distance(d)?;
    let (ga, gb) = (borrowed(a), borrowed(b));
    if is_collection(ga) || is_collection(gb) {
        // The intersection test first: it also carries the collection
        // refusal, which empties must not mask.
        if intersects(a, b)? {
            return Ok(true);
        }
        let (ea, eb) = (elements(ga), elements(gb));
        if ea.is_empty() || eb.is_empty() {
            return Ok(false);
        }
        return within_elements(&ea, &eb, d);
    }
    let (ea, eb) = (elements(ga), elements(gb));
    if ea.is_empty() || eb.is_empty() {
        return Ok(false);
    }
    // The intersection test, split into its two kinds of contact (module
    // documentation, "Contact without the matrix"): an element contact is
    // found with the element pairs, and a component wholly inside the other
    // operand by locating one vertex of it.
    //
    // Contact at 0 is tried first: its pair test is the orientation kernel
    // and its box filter is tight, so touching operands, the common case for
    // neighbors, are decided cheaply. A point set has no segments to touch
    // with, so it skips the extra pass. The order of the alternatives does
    // not change the answer.
    let puntal = |g: Borrowed<'_>| matches!(g, Borrowed::Point(_) | Borrowed::MultiPoint(_));
    Ok(
        (d > 0.0 && !puntal(ga) && !puntal(gb) && within_elements(&ea, &eb, 0.0)?)
            || within_elements(&ea, &eb, d)?
            || component_inside(ga, gb)
            || component_inside(gb, ga),
    )
}

/// Whether some element pair of `ea` and `eb` passes its exact predicate
/// with `d`, or is a pair of properly crossing segments: by the filtered
/// scan for few pairs, through the box index for many.
fn within_elements(ea: &[Element], eb: &[Element], d: f64) -> Result<bool, GeometryError> {
    if ea.len().saturating_mul(eb.len()) <= BRUTE_PAIRS {
        Ok(within_filtered(ea, eb, d))
    } else {
        within_indexed(ea, eb, d)
    }
}

/// Whether `g` is a collection.
fn is_collection(g: Borrowed<'_>) -> bool {
    matches!(g, Borrowed::GeometryCollection(_))
}

/// Whether some connected component of `b` has a representative vertex
/// (the first vertex of each line, the first exterior vertex of each
/// polygon part, every point) not in the exterior of the areal `a`. `false`
/// when `a` is not areal: a line or a point set has no interior to hold a
/// component. Exact (`locate` on a stored coordinate).
fn component_inside(a: Borrowed<'_>, b: Borrowed<'_>) -> bool {
    use crate::topology::areal::{AreaOperand, Areal};
    use crate::topology::label::Location;
    let area = match a {
        Borrowed::Polygon(p) => Areal::single(p),
        Borrowed::MultiPolygon(m) if !m.is_empty() => Areal::multi(m),
        _ => return false,
    };
    let inside = |c: Coord2d| area.locate(c) != Location::Exterior;
    match b {
        Borrowed::Point(p) => inside(p.coord()),
        Borrowed::MultiPoint(m) => m.points().iter().any(|p| inside(p.coord())),
        Borrowed::LineString(l) => inside(l.start()),
        Borrowed::MultiLineString(m) => m.line_strings().iter().any(|l| inside(l.start())),
        Borrowed::Polygon(p) => inside(p.exterior().coords()[0]),
        Borrowed::MultiPolygon(m) => m
            .polygons()
            .iter()
            .any(|p| inside(p.exterior().coords()[0])),
        Borrowed::GeometryCollection(_) => false,
    }
}

/// The distance between `a` and `b`: **exactly `0.0` if and only if they
/// intersect**, and otherwise a positive `f64` within `max(64·u·M, 2⁻¹⁰⁷⁴)`
/// of the true distance, with `u = 2⁻⁵³` and `M` the largest coordinate
/// magnitude in either operand (the derivation, and the zero guard behind
/// the `2⁻¹⁰⁷⁴` term, are in the [module documentation](self)).
///
/// For the question "are they within `d` of each other?", use [`dwithin`],
/// which is exact. **`distance(a, b) <= d` is not that question's answer:**
/// within the bound of the true distance, the comparison can fall on either
/// side of `d`.
/// The two are consistent up to it: if `dwithin(a, b, d)` is `false` then
/// `distance(a, b) > d − 64·u·M`, and if `true` then
/// `distance(a, b) < d + 64·u·M`.
///
/// # Errors
///
/// - [`GeometryError::EmptyOperand`] when either operand has no points: the
///   distance to nothing is undefined, and 0 would read as contact.
/// - [`GeometryError::CollectionElementsOverlap`] and
///   [`GeometryError::TooManyEntries`], as for [`dwithin`].
pub fn distance(a: &impl Relatable, b: &impl Relatable) -> Result<f64, GeometryError> {
    let (ga, gb) = (borrowed(a), borrowed(b));
    let ea = elements(ga);
    let eb = elements(gb);
    let touching = if is_collection(ga) || is_collection(gb) {
        intersects(a, b)?
    } else if ea.is_empty() || eb.is_empty() {
        false
    } else {
        // Contact without the matrix, at `d = 0` (module documentation).
        within_elements(&ea, &eb, 0.0)? || component_inside(ga, gb) || component_inside(gb, ga)
    };
    if ea.is_empty() || eb.is_empty() {
        return Err(GeometryError::EmptyOperand);
    }
    if touching {
        return Ok(0.0);
    }
    let m = coordinate_magnitude(&ea).max(coordinate_magnitude(&eb));
    let delta = DIST_ERR_C * U * m;
    let d = if ea.len().saturating_mul(eb.len()) <= BRUTE_PAIRS {
        distance_pruned(&ea, &eb, delta)
    } else {
        distance_indexed(&ea, &eb, delta)?
    };
    Ok(certified_positive(d))
}

/// A computed distance between operands the exact test has certified
/// apart: `v` itself when positive, and otherwise `2⁻¹⁰⁷⁴`, the smallest
/// positive `f64`, so a rounded `0.0` never claims contact (the zero guard
/// of the [module documentation](self)).
pub(crate) fn certified_positive(v: f64) -> f64 {
    if v > 0.0 { v } else { f64::from_bits(1) }
}

/// `Ok(())` for a finite `d ≥ 0`.
pub(crate) fn check_distance(d: f64) -> Result<(), GeometryError> {
    if d.is_finite() && d >= 0.0 {
        Ok(())
    } else {
        Err(GeometryError::InvalidTolerance { tolerance: d })
    }
}

// ===========================================================================
// Elements: the boundary pieces distance is attained between
// ===========================================================================

/// One boundary element: an isolated point, or a segment of a line or ring.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Element {
    /// An isolated point.
    Point(Coord2d),
    /// A segment, by its stored endpoints.
    Segment(Coord2d, Coord2d),
}

impl Element {
    /// The element's box: its own coordinates, exactly.
    pub(crate) fn corners(self) -> (Coord2d, Coord2d) {
        match self {
            Self::Point(p) => (p, p),
            Self::Segment(p, q) => (
                Coord2d::new(p.x.min(q.x), p.y.min(q.y)),
                Coord2d::new(p.x.max(q.x), p.y.max(q.y)),
            ),
        }
    }
}

/// Every boundary element of `g`, in storage order: each polygon's exterior
/// then holes, ring edges in ring order; line segments in order; points.
/// Collections are walked depth first.
pub(crate) fn elements(g: Borrowed<'_>) -> Vec<Element> {
    let mut out = Vec::new();
    collect(g, &mut out);
    out
}

fn collect(g: Borrowed<'_>, out: &mut Vec<Element>) {
    match g {
        Borrowed::Point(p) => out.push(Element::Point(p.coord())),
        Borrowed::MultiPoint(m) => out.extend(m.points().iter().map(|p| Element::Point(p.coord()))),
        Borrowed::LineString(l) => {
            out.extend(
                l.coords()
                    .as_slice()
                    .windows(2)
                    .map(|w| Element::Segment(w[0], w[1])),
            );
        }
        Borrowed::MultiLineString(m) => {
            for l in m.line_strings() {
                collect(Borrowed::LineString(l), out);
            }
        }
        Borrowed::Polygon(p) => {
            for ring in core::iter::once(p.exterior()).chain(p.holes().iter()) {
                out.extend(ring.edges().map(|(a, b)| Element::Segment(a, b)));
            }
        }
        Borrowed::MultiPolygon(m) => {
            for p in m.polygons() {
                collect(Borrowed::Polygon(p), out);
            }
        }
        Borrowed::GeometryCollection(c) => {
            for g in c.geometries() {
                collect(Borrowed::from_geometry(g), out);
            }
        }
    }
}

/// The largest coordinate magnitude over the elements: the `M` of the
/// bound. An exact maximum of stored values.
pub(crate) fn coordinate_magnitude(es: &[Element]) -> f64 {
    let mut m = 0.0_f64;
    for e in es {
        let (p, q) = match *e {
            Element::Point(p) => (p, p),
            Element::Segment(p, q) => (p, q),
        };
        m = m
            .max(p.x.abs())
            .max(p.y.abs())
            .max(q.x.abs())
            .max(q.y.abs());
    }
    m
}

// ===========================================================================
// The exact pair predicate (dwithin)
// ===========================================================================

/// The sign `arithmetica` returns, with the `Err` arm written out.
///
/// The error is unreachable: every coordinate comes from a validated
/// geometry, so it is finite, and `d` has passed [`check_distance`]. As in
/// `exact::kernel::orient`, the fallback is stated rather than unwrapped;
/// `Zero` here means "exactly at distance `d`", which is within.
fn sign_or_zero(r: Result<Sign, arithmetica::NonFinite>) -> Sign {
    r.unwrap_or(Sign::Zero)
}

/// Whether two boundary elements of non-intersecting geometries are within
/// `d` of each other, decided exactly.
///
/// The decomposition's intersection test has already excluded every
/// contact, so two segments here never cross or touch and the endpoint
/// reduction (step 3) applies.
pub(crate) fn pair_within(a: Element, b: Element, d: f64) -> bool {
    let not_farther = |s: Sign| s != Sign::Positive;
    let pt = |p: Coord2d| [p.x, p.y];
    match (a, b) {
        (Element::Point(p), Element::Point(q)) => not_farther(sign_or_zero(
            arithmetica::compare_point_distance(pt(p), pt(q), d),
        )),
        (Element::Point(p), Element::Segment(s, t))
        | (Element::Segment(s, t), Element::Point(p)) => not_farther(sign_or_zero(
            arithmetica::compare_segment_distance(pt(p), pt(s), pt(t), d),
        )),
        (Element::Segment(a1, a2), Element::Segment(b1, b2)) => {
            // The four endpoint-to-segment comparisons of step 3, all of
            // them: a configuration nearest at an endpoint of `b` is seen
            // only from `b`'s side.
            [(a1, b1, b2), (a2, b1, b2), (b1, a1, a2), (b2, a1, a2)]
                .into_iter()
                .any(|(p, s, t)| {
                    not_farther(sign_or_zero(arithmetica::compare_segment_distance(
                        pt(p),
                        pt(s),
                        pt(t),
                        d,
                    )))
                })
        }
    }
}

/// Every pair, exhaustively.
#[cfg(test)]
pub(crate) fn within_brute(ea: &[Element], eb: &[Element], d: f64) -> bool {
    ea.iter().any(|&a| eb.iter().any(|&b| pair_meets(a, b, d)))
}

/// [`within_brute`]'s answer, with each pair first tested by box: a pair
/// whose second element's box misses the first's box grown outward by `d`
/// is skipped. That is the filter the indexed search applies (the grown
/// window rounds outward, so a pair within `d` always passes it), so the
/// answer is the same; only the exact comparisons of far pairs are not made.
pub(crate) fn within_filtered(ea: &[Element], eb: &[Element], d: f64) -> bool {
    let (Some(ba), Some(bb)) = (bounds(ea), bounds(eb)) else {
        return false;
    };
    // An element of one side that misses the other side's whole box grown
    // by `d` is in no passing pair: every element of the other side lies in
    // that box, and the growth rounds outward.
    let (wa, wb) = (grow_outward(ba, d), grow_outward(bb, d));
    let near_b: Vec<(Element, (Coord2d, Coord2d))> = eb
        .iter()
        .map(|&b| (b, b.corners()))
        .filter(|(_, c)| box_meets(*c, &wa))
        .collect();
    ea.iter().filter(|a| box_meets(a.corners(), &wb)).any(|&a| {
        let w = grow_outward(a.corners(), d);
        near_b
            .iter()
            .any(|&(b, c)| box_meets(c, &w) && pair_meets(a, b, d))
    })
}

/// [`pair_meets`] at `d = 0`, by the orientation kernel: a distance of at
/// most 0 is a distance of exactly 0, so a point is within 0 of a point
/// exactly when they are equal and of a segment exactly when it lies on the
/// closed segment (`on_segment`); two segments meet exactly when an endpoint
/// of one lies on the other or they cross properly. The same answers as the
/// distance predicates at 0, from predicates of lower degree.
fn pair_touches(a: Element, b: Element) -> bool {
    use crate::exact::kernel::{on_segment, segments_properly_cross};
    match (a, b) {
        (Element::Point(p), Element::Point(q)) => p == q,
        (Element::Point(p), Element::Segment(s, t))
        | (Element::Segment(s, t), Element::Point(p)) => on_segment(p, s, t),
        (Element::Segment(p, q), Element::Segment(r, s)) => {
            on_segment(p, r, s)
                || on_segment(q, r, s)
                || on_segment(r, p, q)
                || on_segment(s, p, q)
                || segments_properly_cross(p, q, r, s)
        }
    }
}

/// The box of every element, `None` for none.
fn bounds(es: &[Element]) -> Option<(Coord2d, Coord2d)> {
    es.iter().map(|e| e.corners()).reduce(|(l, h), (m, k)| {
        (
            Coord2d::new(l.x.min(m.x), l.y.min(m.y)),
            Coord2d::new(h.x.max(k.x), h.y.max(k.y)),
        )
    })
}

/// [`pair_within`], or two segments that cross properly: the pair test of
/// the matrix-free contact (module documentation). A crossing has distance
/// 0 while all four endpoint distances can exceed `d`, which is the one
/// contact the endpoint predicates do not see.
fn pair_meets(a: Element, b: Element, d: f64) -> bool {
    if d == 0.0 {
        return pair_touches(a, b);
    }
    if pair_within(a, b, d) {
        return true;
    }
    match (a, b) {
        (Element::Segment(p, q), Element::Segment(r, s)) => {
            crate::exact::kernel::segments_properly_cross(p, q, r, s)
        }
        _ => false,
    }
}

/// Whether the closed box `(lo, hi)` meets `w`: comparisons only.
fn box_meets((lo, hi): (Coord2d, Coord2d), w: &Envelope) -> bool {
    lo.x <= w.max().x && hi.x >= w.min().x && lo.y <= w.max().y && hi.y >= w.min().y
}

/// The same answer through the box index (module documentation): index one
/// side, query with each other-side box grown by `d`, stop at the first
/// passing pair.
pub(crate) fn within_indexed(
    ea: &[Element],
    eb: &[Element],
    d: f64,
) -> Result<bool, GeometryError> {
    // Index the larger side, probe with the smaller.
    let (probe, indexed) = if ea.len() <= eb.len() {
        (ea, eb)
    } else {
        (eb, ea)
    };
    let index = element_index(indexed)?;
    for &e in probe {
        let window = grow_outward(e.corners(), d);
        let found = index.try_candidates(&window, |id| {
            if pair_meets(e, indexed[id as usize], d) {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        });
        if found.is_break() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The elements' boxes as an index. More than `u32::MAX` elements is
/// [`GeometryError::TooManyEntries`].
pub(crate) fn element_index(es: &[Element]) -> Result<BoxIndex, GeometryError> {
    crate::index::check_count(es.len())?;
    let entries: Vec<((Coord2d, Coord2d), u32)> = es
        .iter()
        .enumerate()
        .map(|(i, e)| {
            let (lo, hi) = e.corners();
            // `check_count` passed, so `i` fits.
            ((lo, hi), u32::try_from(i).unwrap_or(u32::MAX))
        })
        .collect();
    Ok(BoxIndex::from_edges(&entries))
}

// ===========================================================================
// The grown window
// ===========================================================================

/// The box `(lo, hi)` grown outward by `r ≥ 0`: each lower bound moved to a
/// float **at or below** `lo − r`, each upper to one **at or above**
/// `hi + r`, clamped to the finite range.
///
/// `fl(lo − r)` can land above `lo − r`; one `next_down` then lands at or
/// below it, because consecutive floats bracket every real. So the returned
/// box contains the exact grown box, and a candidate filter over it can
/// only gain candidates, never lose one, which is the safe direction. The
/// clamp loses nothing: every stored coordinate is finite, so no candidate
/// box reaches past `±f64::MAX`. Crate-private; the coverage gap-width
/// check composes the same helper.
pub(crate) fn grow_outward((lo, hi): (Coord2d, Coord2d), r: f64) -> Envelope {
    debug_assert!(r >= 0.0, "a growth radius is non-negative");
    let down = |v: f64| {
        let g = (v - r).next_down();
        if g < f64::MIN { f64::MIN } else { g }
    };
    let up = |v: f64| {
        let g = (v + r).next_up();
        if g > f64::MAX { f64::MAX } else { g }
    };
    Envelope::new(
        Coord2d::new(down(lo.x), down(lo.y)),
        Coord2d::new(up(hi.x), up(hi.y)),
    )
    .expect("grown bounds are finite by the clamp and ordered by construction")
}

// ===========================================================================
// The f64 distances (the `distance` side)
// ===========================================================================

/// `√(dx² + dy²)` with an exact power-of-two prescale, so the bound needs
/// no assumption about the platform's `hypot`.
///
/// The scaling multiplies by `2^±600` exactly and keeps the squares normal
/// (magnitudes land in about `[2^-474, 2^424]` after scaling), the squares
/// and sum cost `(1 + θ3)`, and `sqrt` is correctly rounded, so the result
/// is within relative `γ3` of the true norm. When one component is more
/// than `2^600` below the other, its scaled square may be flushed, a
/// relative effect under `2^-1000` absorbed by the same `γ3`. A
/// non-finite input gives `+∞` (an overflowing difference; module
/// documentation), never a finite wrong value.
fn norm(dx: f64, dy: f64) -> f64 {
    let m = dx.abs().max(dy.abs());
    if m == 0.0 {
        return 0.0;
    }
    if !m.is_finite() {
        return f64::INFINITY;
    }
    let scale = if m >= 2.0_f64.powi(510) {
        2.0_f64.powi(-600)
    } else if m <= 2.0_f64.powi(-511) {
        2.0_f64.powi(600)
    } else {
        1.0
    };
    let (x, y) = (dx * scale, dy * scale);
    (x * x + y * y).sqrt() / scale
}

/// The `f64` distance from `p` to the closed segment `ab`: the clamped
/// projection, then the norm (bound in the module documentation). When the
/// projection's products overflow, falls back to the nearer endpoint,
/// which errs by at most the segment's (then astronomical) length.
fn point_segment(p: Coord2d, a: Coord2d, b: Coord2d) -> f64 {
    let (ux, uy) = (b.x - a.x, b.y - a.y);
    let (wx, wy) = (p.x - a.x, p.y - a.y);
    let dot = wx * ux + wy * uy;
    let len2 = ux * ux + uy * uy;
    if !dot.is_finite() || !len2.is_finite() {
        return norm(wx, wy).min(norm(p.x - b.x, p.y - b.y));
    }
    if dot <= 0.0 || len2 == 0.0 {
        return norm(wx, wy);
    }
    if dot >= len2 {
        return norm(p.x - b.x, p.y - b.y);
    }
    let t = dot / len2;
    norm(wx - t * ux, wy - t * uy)
}

/// The `f64` distance between two elements of non-intersecting geometries.
pub(crate) fn pair_distance(a: Element, b: Element) -> f64 {
    match (a, b) {
        (Element::Point(p), Element::Point(q)) => norm(p.x - q.x, p.y - q.y),
        (Element::Point(p), Element::Segment(s, t))
        | (Element::Segment(s, t), Element::Point(p)) => point_segment(p, s, t),
        (Element::Segment(a1, a2), Element::Segment(b1, b2)) => point_segment(a1, b1, b2)
            .min(point_segment(a2, b1, b2))
            .min(point_segment(b1, a1, a2))
            .min(point_segment(b2, a1, a2)),
    }
}

/// The minimum over every pair, exhaustively.
#[cfg(test)]
pub(crate) fn distance_brute(ea: &[Element], eb: &[Element]) -> f64 {
    let mut best = f64::INFINITY;
    for &a in ea {
        for &b in eb {
            best = best.min(pair_distance(a, b));
        }
    }
    best
}

/// [`distance_brute`]'s minimum, bit for bit, skipping each pair whose
/// second element's box misses the first's box grown by `(best + δ)`,
/// rounded upward: the pruning [`distance_indexed`] applies, by the same
/// argument (a skipped pair's computed value exceeds the running `best`, so
/// it cannot be the minimum). While no pair has been computed the radius is
/// infinite and nothing is skipped.
pub(crate) fn distance_pruned(ea: &[Element], eb: &[Element], delta: f64) -> f64 {
    // A first pass over the pairs whose boxes overlap, the likeliest to be
    // nearest, so the pruning radius is small from the start. Its values
    // are pair distances like any other, so the minimum is unchanged.
    let mut best = f64::INFINITY;
    for &a in ea {
        let (lo, hi) = a.corners();
        for &b in eb {
            let (m, k) = b.corners();
            if lo.x <= k.x && m.x <= hi.x && lo.y <= k.y && m.y <= hi.y {
                best = best.min(pair_distance(a, b));
            }
        }
    }
    for &a in ea {
        let r = (best + delta).next_up();
        if !r.is_finite() {
            for &b in eb {
                best = best.min(pair_distance(a, b));
            }
            continue;
        }
        let w = grow_outward(a.corners(), r);
        for &b in eb {
            if box_meets(b.corners(), &w) {
                best = best.min(pair_distance(a, b));
            }
        }
    }
    best
}

/// The same minimum through the box index, pruning with the shrinking
/// radius `(best + δ)` rounded upward. A pruned pair's computed value
/// exceeds `best` (module documentation), so the result **equals**
/// [`distance_brute`]'s, bit for bit.
pub(crate) fn distance_indexed(
    ea: &[Element],
    eb: &[Element],
    delta: f64,
) -> Result<f64, GeometryError> {
    let (probe, indexed) = if ea.len() <= eb.len() {
        (ea, eb)
    } else {
        (eb, ea)
    };
    let index = element_index(indexed)?;
    // A finite seed, so the first window is finite: the first pair.
    let mut best = pair_distance(probe[0], indexed[0]);
    for &e in probe {
        // Upward: `fl(best + δ)` can land below the exact sum; one
        // `next_up` lands at or above it, so the radius never understates
        // `best + δ` and the pruning argument holds.
        let r = (best + delta).next_up();
        if !r.is_finite() {
            // An infinite radius prunes nothing: scan this element.
            for &other in indexed {
                best = best.min(pair_distance(e, other));
            }
            continue;
        }
        let window = grow_outward(e.corners(), r);
        index.candidates(&window, |id| {
            best = best.min(pair_distance(e, indexed[id as usize]));
        });
    }
    Ok(best)
}
