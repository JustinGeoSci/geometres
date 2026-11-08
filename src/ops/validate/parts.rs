// src/ops/validate/parts.rs

//! Validity of multipolygon material: each part's issues, and the issues
//! between parts.

use crate::exact::kernel::{RingPlacement, place_ring};
use crate::types::{CoordSeq, LinearRing, Polygon, parts_meet_only_at_points};

use super::{ValidityIssue, validate};

/// Every issue that would keep `parts` from becoming a `MultiPolygon`: each
/// part's rings as [`validate`] takes them (ring 0 the exterior), with the
/// part's index.
///
/// The issues of each part come first, in part order, each tagged with its
/// part. Then the issues between parts, for every pair of parts that would
/// each construct (no issue a constructor refuses), in the order
/// `MultiPolygon::new` checks pairs, tagged with the issue's `part`:
///
/// - [`ValidityIssue::NestedParts`] when one part's exterior lies inside the
///   other's and meets none of its holes beyond isolated points, so its
///   interior is inside the other's body;
/// - [`ValidityIssue::PartsOverlap`] for any other overlap: crossing
///   exteriors, a shared stretch, or a part's exterior meeting the other's
///   holes.
///
/// A pair draws one of these if and only if `MultiPolygon::new` would refuse
/// it with `GeometryError::OverlappingParts`, by the same exact test
/// (interiors disjoint, no shared stretch); the first names the pair it
/// names. A pair involving a part that would not construct is not judged:
/// its regions are not defined.
#[must_use]
pub fn validate_multipolygon(parts: &[Vec<CoordSeq>]) -> Vec<(usize, ValidityIssue)> {
    let mut issues = Vec::new();
    let mut polygons: Vec<Option<Polygon>> = Vec::with_capacity(parts.len());
    for (part, rings) in parts.iter().enumerate() {
        let own = validate(rings);
        let constructs = !own.iter().any(ValidityIssue::constructor_refuses);
        issues.extend(own.into_iter().map(|i| (part, i)));
        polygons.push(if constructs { construct(rings) } else { None });
    }
    for (i, p) in polygons.iter().enumerate() {
        let Some(p) = p else { continue };
        for (j, q) in polygons.iter().enumerate().skip(i + 1) {
            let Some(q) = q else { continue };
            if let Some(issue) = between_parts((i, p), (j, q)) {
                let part = match issue {
                    ValidityIssue::NestedParts { part, .. }
                    | ValidityIssue::PartsOverlap { part, .. } => part,
                    _ => j,
                };
                issues.push((part, issue));
            }
        }
    }
    issues
}

/// The polygon of rings that drew no constructor-refused issue. Every
/// constructor check then passes, so this is `Some`; `None` is kept rather
/// than unwrapped.
fn construct(rings: &[CoordSeq]) -> Option<Polygon> {
    let mut lrs = rings.iter().map(|c| LinearRing::new(c.clone()).ok());
    let exterior = lrs.next()??;
    let holes: Option<Vec<LinearRing>> = lrs.collect();
    Polygon::new(exterior, holes?).ok()
}

/// `q`'s exterior meets none of `p`'s holes beyond isolated points: it is
/// `Outside` every hole. Only then does an exterior inside `p`'s exterior put
/// `q`'s interior inside `p`'s body (nested). A hole that `q`'s exterior lies
/// in, crosses, or shares a stretch with leaves part or all of `q`'s
/// interior outside `p`'s body, so the pair is an overlap, not a nesting.
fn clear_of_holes(q: &Polygon, p: &Polygon) -> bool {
    p.holes()
        .iter()
        .all(|h| place_ring(q.exterior(), h) == RingPlacement::Outside)
}

/// The issue between parts `i < j`, if they meet in more than finitely many
/// points.
pub(super) fn between_parts(
    (i, p): (usize, &Polygon),
    (j, q): (usize, &Polygon),
) -> Option<ValidityIssue> {
    if parts_meet_only_at_points(p, q) {
        return None;
    }
    Some(
        match (
            place_ring(q.exterior(), p.exterior()),
            place_ring(p.exterior(), q.exterior()),
        ) {
            (RingPlacement::Inside, RingPlacement::Outside) if clear_of_holes(q, p) => {
                ValidityIssue::NestedParts { part: j, inside: i }
            }
            (RingPlacement::Outside, RingPlacement::Inside) if clear_of_holes(p, q) => {
                ValidityIssue::NestedParts { part: i, inside: j }
            }
            _ => ValidityIssue::PartsOverlap { part: j, other: i },
        },
    )
}
