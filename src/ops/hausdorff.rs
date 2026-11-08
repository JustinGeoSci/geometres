// src/ops/hausdorff.rs

//! The discrete Hausdorff distance, `hausdorff_distance`, and its exact
//! predicate, `hausdorff_within`.
//!
//! # The definition, and whose it is
//!
//! The **discrete Hausdorff distance** between `a` and `b` is the larger of
//! two directed values:
//!
//! - `h(a → b)`: the largest, over the vertices `v` of `a`, of the distance
//!   from `v` to the **linework** of `b`;
//! - `h(b → a)`: the same with the roles swapped.
//!
//! The linework of a geometry is its boundary elements, as in
//! [`distance`](super::distance()): the segments of its lines, the segments
//! of every polygon ring (exterior and holes), and its isolated points. **A
//! polygon is measured by its rings: its interior does not count.** A vertex
//! of `a` strictly inside a polygon of `b` contributes its distance to the
//! nearest ring of that polygon, not 0. The vertices of a geometry are its
//! stored coordinates: points, line vertices, and ring vertices.
//!
//! This is the definition of JTS's `DiscreteHausdorffDistance`
//! (`org.locationtech.jts.algorithm.distance`, Martin Davis). It visits every
//! coordinate of one geometry and measures it with `DistanceToPoint.computeDistance`,
//! which for a polygon measures to the exterior ring and to each hole, for a line
//! to each segment, and for a point to the point. Matching it means the number is
//! the one PostGIS's `ST_HausdorffDistance` and Shapely report for the same input
//! (without densification). A *region* reading (a polygon measured as a point set,
//! interior at distance 0) is a different quantity and is **not** what this function returns.
//! Densified variants are not offered.
//!
//! The discrete value never exceeds the continuous Hausdorff distance
//! between the same linework: it maximizes over a subset (the vertices) of
//! the points the continuous one maximizes over.
//!
//! # The predicate is exact
//!
//! [`hausdorff_within`]`(a, b, d)` is `true` exactly when every vertex of `a`
//! is within `d` of the linework of `b`, and every vertex of `b` is within `d`
//! of the linework of `a`. Each test is one of `arithmetica`'s exact
//! predicates, point to point
//! ([`compare_point_distance`](arithmetica::compare_point_distance)) or point
//! to closed segment
//! ([`compare_segment_distance`](arithmetica::compare_segment_distance)). These
//! are the predicates [`dwithin`](super::dwithin) decides boundary elements
//! with, so the answer is exact for every input. The first vertex that fails
//! ends the search.
//!
//! # The quantity: an exact zero, and a bound
//!
//! [`hausdorff_distance`] is **exactly `0.0` if and only if every vertex of
//! each geometry lies on the linework of the other**, decided by the same
//! exact predicates at `d = 0` (a point on a closed segment, or equal to a
//! point). Otherwise it is an `f64` within `δ = c·u·M` of the true value,
//! with `c = 64`, `u = 2⁻⁵³` and `M` the largest coordinate magnitude in
//! either operand. This is the bound of [`distance`](super::distance()),
//! whose module documentation derives it for one point-to-element distance.
//!
//! **The bound carries through maxima and minima.** Let `xᵢ` be computed
//! values with `|xᵢ − Xᵢ| ≤ δ` for each `i`. With `Xₖ` the largest true
//! value, `max x ≥ xₖ ≥ Xₖ − δ`; and `max x = xⱼ ≤ Xⱼ + δ ≤ max X + δ`. So
//! `|max x − max X| ≤ δ`, and the same argument with the inequalities
//! reversed gives `|min x − min X| ≤ δ`. The discrete Hausdorff distance is
//! a maximum (over both directions and their vertices) of minima (over the
//! other side's elements) of point-to-element distances, each within `δ`;
//! so the computed value is within `δ` of the true one.
//!
//! **One guard remains.** When the zero test fails, the true value is positive, but
//! every computed point-to-segment distance could round to `0.0` (a vertex
//! a hair off a segment). Returning `0.0` would read as "the zero test
//! holds"; the function returns the smallest positive `f64`, `2⁻¹⁰⁷⁴`,
//! instead. The true value is then at most `δ`, so the result is within
//! `max(δ, 2⁻¹⁰⁷⁴)` of it, the bound stated on the function.
//!
//! **Consistency with the predicate** follows: if
//! `hausdorff_within(a, b, d)` is `false`, the true value exceeds `d`, so
//! `hausdorff_distance(a, b) > d − δ`; if `true`, the true value is at most
//! `d`, so `hausdorff_distance(a, b) ≤ d + δ`. `hausdorff_distance(a, b) <=
//! d` is not the predicate: within `δ` of the true value, the
//! comparison can fall on either side.
//!
//! # Indexing never changes the answer
//!
//! Above a fixed number of vertex–element pairs (the threshold `distance`
//! uses), each direction indexes the target's elements in a static box index
//! (through the crate's private `situs` seam) and queries it once per vertex.
//! For the predicate the query window is the vertex grown outward by `d`; for
//! the quantity, by `(best + δ)` rounded upward, where `best` is the vertex's
//! smallest distance so far. These are the windows, and the pruning argument,
//! of [`distance`](super::distance()): an element outside the window has true
//! distance greater than `best + δ`, so its computed distance exceeds `best`
//! and cannot change the minimum. The minimum of `f64` values does not depend
//! on the order they are visited in, so the indexed result equals the
//! exhaustive one bit for bit.
//!
//! # Errors and empties
//!
//! `d` must be finite and at least 0 ([`GeometryError::InvalidTolerance`]).
//! An empty operand gives `false` from the predicate and
//! [`GeometryError::EmptyOperand`] from the quantity, as for `dwithin` and
//! `distance`. No topology is consulted, so a collection whose polygonal
//! leaves overlap is accepted: its linework is well defined.
//!
//! Neither function returns a coordinate.

use core::ops::ControlFlow;

use super::distance::{self as dist, Element};
use crate::GeometryError;
use crate::predicates::{Borrowed, Relatable, borrowed};
use crate::types::Coord2d;

// ===========================================================================
// The public functions
// ===========================================================================

/// Whether the discrete Hausdorff distance between `a` and `b` is at most
/// `d`, decided **exactly**.
///
/// `true` when every vertex of `a` is within `d` of the linework of `b`
/// (its segments, ring segments and points; a polygon's interior does not
/// count), and every vertex of `b` within `d` of the linework of `a`. Each
/// test is an exact `arithmetica` predicate; see the
/// [module documentation](self) for the definition.
///
/// ```
/// use geometres::ops::hausdorff_within;
/// use geometres::types::{Coord2d, CoordSeq, LineString};
///
/// let line = |v: &[(f64, f64)]| {
///     LineString::new(CoordSeq::new(v.iter().map(|&(x, y)| Coord2d::new(x, y)).collect()))
/// };
/// let a = line(&[(0.0, 0.0), (10.0, 0.0)])?;
/// let b = line(&[(0.0, 3.0), (10.0, 3.0)])?;
/// assert!(hausdorff_within(&a, &b, 3.0)?);
/// assert!(!hausdorff_within(&a, &b, 3.0_f64.next_down())?);
/// # Ok::<(), geometres::GeometryError>(())
/// ```
///
/// An empty operand gives `false`.
///
/// # Errors
///
/// - [`GeometryError::InvalidTolerance`] for a negative or non-finite `d`.
/// - [`GeometryError::TooManyEntries`] for more than `u32::MAX` elements on
///   an indexed side.
pub fn hausdorff_within(
    a: &impl Relatable,
    b: &impl Relatable,
    d: f64,
) -> Result<bool, GeometryError> {
    dist::check_distance(d)?;
    let (va, ea) = (vertices(borrowed(a)), targets(borrowed(a)));
    let (vb, eb) = (vertices(borrowed(b)), targets(borrowed(b)));
    if va.is_empty() || vb.is_empty() {
        return Ok(false);
    }
    Ok(directed_within(&va, &eb, d)? && directed_within(&vb, &ea, d)?)
}

/// The discrete Hausdorff distance between `a` and `b`: **exactly `0.0` if
/// and only if every vertex of each lies on the linework of the other**,
/// and otherwise an `f64` within `max(64·u·M, 2⁻¹⁰⁷⁴)` of the true value,
/// with `u = 2⁻⁵³` and `M` the largest coordinate magnitude in either
/// operand (the definition and the derivation are in the
/// [module documentation](self)).
///
/// For the question "is the Hausdorff distance at most `d`?", use
/// [`hausdorff_within`], which is exact: if it is `false` then this value
/// is greater than `d − 64·u·M`, and if `true` then at most `d + 64·u·M`.
///
/// # Errors
///
/// - [`GeometryError::EmptyOperand`] when either operand has no points.
/// - [`GeometryError::TooManyEntries`], as for [`hausdorff_within`].
pub fn hausdorff_distance(a: &impl Relatable, b: &impl Relatable) -> Result<f64, GeometryError> {
    let (va, ea) = (vertices(borrowed(a)), targets(borrowed(a)));
    let (vb, eb) = (vertices(borrowed(b)), targets(borrowed(b)));
    if va.is_empty() || vb.is_empty() {
        return Err(GeometryError::EmptyOperand);
    }
    if directed_within(&va, &eb, 0.0)? && directed_within(&vb, &ea, 0.0)? {
        return Ok(0.0);
    }
    let m = dist::coordinate_magnitude(&ea).max(dist::coordinate_magnitude(&eb));
    let delta = dist::DIST_ERR_C * dist::U * m;
    let h = directed(&va, &eb, delta)?.max(directed(&vb, &ea, delta)?);
    // The guard of the module documentation: the zero test failed, so the
    // true value is positive and a computed 0.0 must not stand for it.
    Ok(dist::certified_positive(h))
}

// ===========================================================================
// Vertices and targets
// ===========================================================================

/// Every vertex of `g`, in storage order: points, line vertices, and each
/// polygon ring's stored vertices (exterior then holes; rings are stored
/// without a closing vertex). Collections are walked depth first.
pub(super) fn vertices(g: Borrowed<'_>) -> Vec<Coord2d> {
    let mut out = Vec::new();
    collect_vertices(g, &mut out);
    out
}

fn collect_vertices(g: Borrowed<'_>, out: &mut Vec<Coord2d>) {
    match g {
        Borrowed::Point(p) => out.push(p.coord()),
        Borrowed::MultiPoint(m) => out.extend(m.points().iter().map(crate::types::Point::coord)),
        Borrowed::LineString(l) => out.extend_from_slice(l.coords().as_slice()),
        Borrowed::MultiLineString(m) => {
            for l in m.line_strings() {
                out.extend_from_slice(l.coords().as_slice());
            }
        }
        Borrowed::Polygon(p) => {
            for ring in core::iter::once(p.exterior()).chain(p.holes().iter()) {
                out.extend_from_slice(ring.coords().as_slice());
            }
        }
        Borrowed::MultiPolygon(m) => {
            for p in m.polygons() {
                collect_vertices(Borrowed::Polygon(p), out);
            }
        }
        Borrowed::GeometryCollection(c) => {
            for g in c.geometries() {
                collect_vertices(Borrowed::from_geometry(g), out);
            }
        }
    }
}

/// What a vertex is measured to: the linework of `g`, meaning its segments,
/// ring segments and points, as `distance` enumerates them.
pub(super) fn targets(g: Borrowed<'_>) -> Vec<Element> {
    dist::elements(g)
}

// ===========================================================================
// The exact directed predicate
// ===========================================================================

/// Whether every vertex in `vs` is within `d` of some element of `es`,
/// decided exactly. Stops at the first vertex that is not.
pub(super) fn directed_within(
    vs: &[Coord2d],
    es: &[Element],
    d: f64,
) -> Result<bool, GeometryError> {
    if vs.len().saturating_mul(es.len()) <= dist::BRUTE_PAIRS {
        return Ok(vs.iter().all(|&v| near_any(v, es, d)));
    }
    let index = dist::element_index(es)?;
    for &v in vs {
        // The vertex grown by `d`, outward: an element within `d` of `v`
        // has an axis gap of at most `d`, so its box meets the window.
        let window = dist::grow_outward((v, v), d);
        let found = index.try_candidates(&window, |id| {
            if dist::pair_within(Element::Point(v), es[id as usize], d) {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        });
        if found.is_continue() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether `v` is within `d` of some element of `es`, exhaustively.
fn near_any(v: Coord2d, es: &[Element], d: f64) -> bool {
    es.iter()
        .any(|&e| dist::pair_within(Element::Point(v), e, d))
}

// ===========================================================================
// The f64 directed value
// ===========================================================================

/// `h(vs → es)` in `f64`: the largest, over `vs`, of the smallest computed
/// distance to an element of `es`. Exhaustive below the pair threshold,
/// indexed above it, with the same result bit for bit (module
/// documentation).
pub(super) fn directed(vs: &[Coord2d], es: &[Element], delta: f64) -> Result<f64, GeometryError> {
    if vs.len().saturating_mul(es.len()) <= dist::BRUTE_PAIRS {
        Ok(directed_brute(vs, es))
    } else {
        directed_indexed(vs, es, delta)
    }
}

/// `h(vs → es)` over every pair.
pub(super) fn directed_brute(vs: &[Coord2d], es: &[Element]) -> f64 {
    let mut worst = 0.0_f64;
    for &v in vs {
        let mut best = f64::INFINITY;
        for &e in es {
            best = best.min(dist::pair_distance(Element::Point(v), e));
        }
        worst = worst.max(best);
    }
    worst
}

/// `h(vs → es)` through the box index. Each vertex's minimum is seeded
/// from the element nearest the previous vertex (consecutive vertices are
/// usually close, so the first window is already small), then searched in
/// the window grown by `(best + δ)` rounded upward.
pub(super) fn directed_indexed(
    vs: &[Coord2d],
    es: &[Element],
    delta: f64,
) -> Result<f64, GeometryError> {
    let index = dist::element_index(es)?;
    let mut worst = 0.0_f64;
    let mut seed = 0_usize;
    for &v in vs {
        let p = Element::Point(v);
        let mut best = dist::pair_distance(p, es[seed]);
        let mut nearest = seed;
        // Upward, as in `distance`: the radius never understates
        // `best + δ`, so the pruning argument holds.
        let r = (best + delta).next_up();
        if r.is_finite() {
            let window = dist::grow_outward((v, v), r);
            index.candidates(&window, |id| {
                let x = dist::pair_distance(p, es[id as usize]);
                if x < best {
                    best = x;
                    nearest = id as usize;
                }
            });
        } else {
            // An infinite radius prunes nothing: scan.
            for (i, &e) in es.iter().enumerate() {
                let x = dist::pair_distance(p, e);
                if x < best {
                    best = x;
                    nearest = i;
                }
            }
        }
        seed = nearest;
        worst = worst.max(best);
    }
    Ok(worst)
}
