// src/ops/validate.rs

//! Validity enumeration: every reason a set of rings could not become a
//! `Polygon`, reported together rather than failing on the first.
//!
//! # What this is, and is not
//!
//! The constructors **reject**: `LinearRing::new` and `Polygon::new` return
//! the first defect they meet and nothing else. [`validate`] is the
//! complementary **non-rejecting** form for the same material: it runs every
//! per-ring check on every ring and returns all of the issues found, so a
//! caller fixing a geometry sees the whole list at once. An empty result
//! means the rings pass every check made here.
//!
//! The two are not duplicates. A `Polygon` that exists is already valid, so
//! `validate` cannot take one; it takes the raw material:
//!
//! - [`validate`] takes `&[CoordSeq]`, ring 0 the would-be exterior and
//!   rings 1.. the would-be holes, and applies the orientation rule by
//!   position.
//! - [`validate_line`] takes one `&CoordSeq` as a would-be `LineString`,
//!   where the only conditions are finite coordinates, two or more vertices
//!   and no consecutive repeat.
//!
//! `CoordSeq` is the one type in the crate that can hold anything a
//! constructor refuses (too few vertices, repeats, a self-crossing, either
//! winding, a non-finite coordinate) while still being coordinates and
//! nothing else: no record, schema, reference system or format.
//!
//! # The checks, per ring
//!
//! | Issue | The constructor's counterpart |
//! |---|---|
//! | [`ValidityIssue::NonFiniteCoordinate`] | `GeometryError::NonFiniteCoordinate` |
//! | [`ValidityIssue::TooFewPoints`] | `GeometryError::InsufficientVertices` |
//! | [`ValidityIssue::DuplicatePoint`] | `GeometryError::ConsecutiveDuplicateVertices` |
//! | [`ValidityIssue::NonFiniteArea`] | `GeometryError::NonFiniteResult` |
//! | [`ValidityIssue::DegenerateRing`] | `GeometryError::DegenerateRing` |
//! | [`ValidityIssue::RingOrientation`] | `GeometryError::IncorrectWindingOrder`, `IncorrectHoleWindingOrder` |
//! | [`ValidityIssue::SelfIntersection`] | none; the constructors do not test simplicity, and the overlay and simplify vouches do |
//!
//! # The checks, between rings
//!
//! | Issue | The constructor's counterpart |
//! |---|---|
//! | [`ValidityIssue::HoleOutsideExterior`] | `GeometryError::HoleNotInsideExterior` |
//! | [`ValidityIssue::RingsCross`] | `HoleNotInsideExterior` (against the exterior), `OverlappingHoles` (between holes) |
//! | [`ValidityIssue::DuplicateRings`] | `HoleNotInsideExterior`, `OverlappingHoles` |
//! | [`ValidityIssue::NestedHoles`] | `GeometryError::OverlappingHoles` |
//! | [`ValidityIssue::DisconnectedInterior`] | none; `Polygon::new` admits it (see below) |
//!
//! Each hole is placed against the exterior, then each pair of holes, in the
//! order `Polygon::new` checks them and by the same exact `place_ring`, so a
//! placement issue is reported if and only if `Polygon::new` refuses the
//! placement, and the first one names the hole or pair `Polygon::new` names.
//! Placement is skipped for a ring that is non-finite, has too few points,
//! has an area that is not a number, or is degenerate: it bounds no region
//! to place.
//!
//! **A disconnected interior is constructible.** `Polygon::new` accepts holes
//! that touch the exterior or each other at points, and a cycle of such
//! touches (a hole touching the exterior twice, or two holes touching each
//! other and each touching the exterior) cuts the interior in two. OGC calls
//! that polygon invalid; `Polygon::new` builds it, like a self-intersecting
//! ring, and [`validity_issues`] reports it. The test, its proof, and why it
//! counts touch points rather than ring pairs are in `ops/validate/between.rs`.
//! It runs only when every other check passes (an orientation issue aside):
//! its proof needs simple, correctly placed rings.
//!
//! # Multipolygons and constructed values
//!
//! [`validate_multipolygon`] runs [`validate`] on each part and adds the
//! issues between parts: [`ValidityIssue::PartsOverlap`] and
//! [`ValidityIssue::NestedParts`], the counterparts of
//! `GeometryError::OverlappingParts`.
//!
//! [`is_valid`] and [`validity_issues`] take a constructed geometry and
//! report what its constructor does not check: ring self-intersection and a
//! disconnected interior. Each issue carries the path of element and part
//! indices that leads to it ([`LocatedIssue`]).
//!
//! # Conventions, all the constructor's
//!
//! - A **closing vertex** (first repeated as last) is accepted and stripped,
//!   as `LinearRing::new` does. It is not a `DuplicatePoint`. Validation
//!   and construction follow one convention for the same ring.
//! - **Orientation** is the certified sign of the ring's area
//!   (`LinearRing::signed_area`): ring 0 must be positive, every other
//!   ring negative. An exactly zero area is a `DegenerateRing`, by the
//!   same exact test `Polygon::new` applies, and no orientation is reported
//!   for it.
//! - **Self-intersection** is the exact self-contact test of the vouches
//!   (`exact::kernel::self_contacts`): every pair of edges that cross, touch,
//!   pinch or fold back, by edge index. No intersection point is
//!   constructed, and none is offered.
//!
//! # Indices
//!
//! Every issue names its `ring` (0 the exterior) and, where it applies, a
//! vertex or edge index into **that ring as given, less the closing vertex
//! if it had one**. When a ring has consecutive repeats, orientation and
//! self-contact are measured on the ring with the repeats removed (the
//! repeats are already reported), and edge indices are mapped back: edge `e`
//! is the edge that leaves vertex `e` toward the next *distinct* vertex.
//!
//! # What is skipped
//!
//! A ring with a non-finite coordinate gets that issue for each offending
//! vertex and nothing else: area and orientation of such a ring are not
//! numbers. A ring with fewer than three distinct vertices gets
//! `TooFewPoints` (and any `DuplicatePoint`s) and nothing else: it has no
//! area and no pair of non-adjacent edges.

mod between;
mod constructed;
mod parts;

pub(crate) use between::disconnected_interior;
pub use constructed::{LocatedIssue, is_valid, validity_issues};
pub use parts::validate_multipolygon;

use crate::exact::kernel::self_contacts;
use crate::types::ring_is_degenerate;
use crate::types::{Coord2d, CoordSeq, LinearRing};
use between::Prepared;

// ===========================================================================
// The issue type
// ===========================================================================

/// One reason a ring could not become part of a valid `Polygon` or
/// `LineString`. Returned by [`validate`] and [`validate_line`], as many as
/// are found.
///
/// Every variant is a statement about coordinates and indices into the
/// caller's input. Nothing here names a record, a schema, a reference system
/// or a format; a consumer that wants per-feature reporting wraps a
/// `Vec<ValidityIssue>` with its own identifier.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidityIssue {
    /// A coordinate of the ring is NaN or infinite.
    ///
    /// Reported once per offending vertex. Nothing else is reported for a
    /// ring with this issue: its area and orientation are not defined.
    NonFiniteCoordinate {
        /// Index of the ring (0 the exterior) or, from [`validate_line`], 0.
        ring: usize,
        /// Index of the vertex in the ring as given.
        index: usize,
    },
    /// The ring has too few distinct vertices.
    ///
    /// Three for a polygon ring, two for a line; `found` counts the vertices
    /// after the closing vertex (if any) and consecutive repeats are removed.
    TooFewPoints {
        /// Index of the ring.
        ring: usize,
        /// Distinct vertices present.
        found: usize,
        /// Distinct vertices needed.
        required: usize,
    },
    /// Two consecutive vertices are identical.
    ///
    /// `index` is the first of the pair; the second is `index + 1`, or
    /// vertex 0 when `index` is the last vertex of a ring (the pair across
    /// the closing edge). The OGC closing vertex itself, the first vertex
    /// repeated as the last, is not a duplicate; it is stripped first, as
    /// `LinearRing::new` strips it.
    DuplicatePoint {
        /// Index of the ring.
        ring: usize,
        /// Index of the first vertex of the identical pair.
        index: usize,
    },
    /// The ring's area is not a finite number: its coordinates are so large
    /// that the shoelace products overflow. `Polygon::new` refuses such a
    /// ring with `GeometryError::NonFiniteResult`, before its degeneracy and
    /// winding checks, which are not defined for it; no orientation is
    /// reported, and the ring is not placed.
    NonFiniteArea {
        /// Index of the ring.
        ring: usize,
    },
    /// The ring is degenerate by the criterion `Polygon::new` applies: its
    /// exact signed area is zero (every vertex collinear, or lobes that
    /// cancel exactly, as in a symmetric bow tie). A permanent contract, with
    /// no tolerance: a ring with any nonzero area is not degenerate. No
    /// orientation is reported for a degenerate ring; a zero area has no
    /// winding.
    DegenerateRing {
        /// Index of the ring.
        ring: usize,
    },
    /// The ring is wound the wrong way for its position: ring 0 must be
    /// counterclockwise (positive certified area) and every other ring
    /// clockwise (negative), per OGC Simple Features and `Polygon::new`.
    RingOrientation {
        /// Index of the ring.
        ring: usize,
    },
    /// Two edges of the ring meet somewhere other than the vertex two
    /// consecutive edges share: a proper crossing, a vertex on a
    /// non-adjacent edge, a vertex visited twice, or a fold-back.
    ///
    /// `edges` are the indices of the two edges; edge `e` runs from vertex
    /// `e` to the next *distinct* vertex. For two consecutive edges that fold
    /// back the earlier comes first, so the pair across the closing edge is
    /// `(n − 1, 0)`; for any other pair the lower index is first. One
    /// contact is reported per pair of edges. The representation is
    /// exact and construction-free. No intersection point is offered: it
    /// would be a constructed coordinate, which the topology layer never
    /// produces, and the two edges locate the
    /// contact without it.
    SelfIntersection {
        /// Index of the ring.
        ring: usize,
        /// The two edges in contact, by index.
        edges: (usize, usize),
    },
    /// A hole that is not inside the exterior: disjoint from it, or
    /// surrounding it, touching it at points at most.
    HoleOutsideExterior {
        /// Index of the hole (1 or more).
        ring: usize,
    },
    /// A hole inside another hole.
    NestedHoles {
        /// Index of the hole that is inside.
        ring: usize,
        /// Index of the hole it is inside.
        inside: usize,
    },
    /// Two different rings cross, or overlap along a stretch of positive
    /// length.
    ///
    /// `ring` is the later of the two in input order and `other` the earlier
    /// (0 for the exterior). `edges` is a pair that witnesses the contact,
    /// `(edge of ring, edge of other)`, by input edge index: two edges lying
    /// along each other, two edges that properly cross, or, where a ring
    /// passes through the other at a point, the edge leaving that point and
    /// the first edge of the other ring holding it. No point is constructed.
    RingsCross {
        /// Index of the later ring.
        ring: usize,
        /// Index of the earlier ring.
        other: usize,
        /// `(edge of ring, edge of other)`.
        edges: (usize, usize),
    },
    /// Two rings bound the same point set (the same ring listed twice,
    /// perhaps from another start, in the other direction, or with extra
    /// collinear vertices). `ring` is the later in input order, `other`
    /// the earlier.
    DuplicateRings {
        /// Index of the later ring.
        ring: usize,
        /// Index of the earlier ring.
        other: usize,
    },
    /// The rings touch one another in a cycle, which cuts the interior into
    /// two or more pieces. Constructible: `Polygon::new` does not check it.
    ///
    /// The touch at vertex `vertex` of ring `ring` is one that closes a
    /// cycle (`ring` is the lowest ring with the touch point as a vertex).
    /// One is reported per polygon.
    DisconnectedInterior {
        /// Index of the ring holding the touch point as a vertex.
        ring: usize,
        /// Index of that vertex, in the ring as given.
        vertex: usize,
    },
    /// Two parts of a multipolygon overlap: their interiors meet, or their
    /// boundaries share a stretch. The counterpart of
    /// `GeometryError::OverlappingParts`. `part` is the later part, `other`
    /// the earlier. From [`validate_multipolygon`] and the multipolygon
    /// checks of [`validity_issues`]; [`ring`](Self::ring) is 0.
    PartsOverlap {
        /// Index of the later part.
        part: usize,
        /// Index of the earlier part.
        other: usize,
    },
    /// One part of a multipolygon lies inside another's body (its exterior
    /// inside the other's exterior, meeting none of its holes beyond
    /// isolated points). A part that lies in, crosses or shares a stretch
    /// with a hole is a [`PartsOverlap`](Self::PartsOverlap). Also a case of
    /// `GeometryError::OverlappingParts`.
    /// [`ring`](Self::ring) is 0.
    NestedParts {
        /// Index of the part that is inside.
        part: usize,
        /// Index of the part it is inside.
        inside: usize,
    },
}

impl ValidityIssue {
    /// The ring the issue is about (0 the exterior): for an issue between
    /// rings, the later ring; for an issue between parts, 0 (the exterior of
    /// the part named).
    #[must_use]
    pub fn ring(&self) -> usize {
        match *self {
            Self::NonFiniteCoordinate { ring, .. }
            | Self::TooFewPoints { ring, .. }
            | Self::DuplicatePoint { ring, .. }
            | Self::NonFiniteArea { ring }
            | Self::DegenerateRing { ring }
            | Self::RingOrientation { ring }
            | Self::SelfIntersection { ring, .. }
            | Self::HoleOutsideExterior { ring }
            | Self::NestedHoles { ring, .. }
            | Self::RingsCross { ring, .. }
            | Self::DuplicateRings { ring, .. }
            | Self::DisconnectedInterior { ring, .. } => ring,
            Self::PartsOverlap { .. } | Self::NestedParts { .. } => 0,
        }
    }

    /// `true` for the issues a constructor refuses: everything except
    /// self-intersection and a disconnected interior, which construction
    /// does not check (`Polygon::new` and `MultiPolygon::new` build such
    /// geometries; [`validity_issues`] reports them).
    #[must_use]
    pub fn constructor_refuses(&self) -> bool {
        !matches!(
            self,
            Self::SelfIntersection { .. } | Self::DisconnectedInterior { .. }
        )
    }
}

// ===========================================================================
// Entry points
// ===========================================================================

/// Every issue that would keep `rings` from becoming a `Polygon`: ring 0 is
/// the would-be exterior, rings 1.. the would-be holes.
///
/// Returns an empty `Vec` when the rings pass every check (the module doc
/// lists them). The per-ring issues come first, in ring order and within a
/// ring in the order of the module doc's first table; then the placement
/// issues, in `Polygon::new`'s order (each hole against the exterior, then
/// each pair of holes); then a disconnected interior, if every other check
/// passed (an orientation issue aside).
///
/// An empty slice has no rings and no issues: `Polygon::new` is not reached
/// with nothing, and there is nothing to report.
#[must_use]
pub fn validate(rings: &[CoordSeq]) -> Vec<ValidityIssue> {
    let mut issues = Vec::new();
    let prepared: Vec<Option<Prepared>> = rings
        .iter()
        .enumerate()
        .map(|(ring, coords)| validate_ring_into(ring, coords, &mut issues))
        .collect();
    between::placement_into(&prepared, &mut issues);
    let only_orientation = issues
        .iter()
        .all(|i| matches!(i, ValidityIssue::RingOrientation { .. }));
    if only_orientation && !prepared.is_empty() {
        let rings: Vec<&LinearRing> = prepared.iter().flatten().map(|p| &p.ring).collect();
        if let Some((ring, k)) = disconnected_interior(&rings) {
            issues.push(ValidityIssue::DisconnectedInterior {
                ring,
                vertex: prepared[ring].as_ref().map_or(k, |p| p.origin[k]),
            });
        }
    }
    issues
}

/// Every issue that would keep `line` from becoming a `LineString`: a
/// non-finite coordinate, fewer than two vertices, a consecutive repeat.
///
/// A line is not closed, so a last vertex equal to the first is an ordinary
/// pair of vertices; and a line may cross itself, so no self-contact is
/// reported. Issues name ring 0.
#[must_use]
pub fn validate_line(line: &CoordSeq) -> Vec<ValidityIssue> {
    let mut issues = Vec::new();
    let v = line.as_slice();
    if non_finite_into(0, v, &mut issues) {
        return issues;
    }
    for i in 0..v.len().saturating_sub(1) {
        if v[i] == v[i + 1] {
            issues.push(ValidityIssue::DuplicatePoint { ring: 0, index: i });
        }
    }
    let distinct = v.len() - issues.len();
    if distinct < 2 {
        issues.push(ValidityIssue::TooFewPoints {
            ring: 0,
            found: distinct,
            required: 2,
        });
    }
    issues
}

// ===========================================================================
// Per-ring checks
// ===========================================================================

/// Pushes every `NonFiniteCoordinate` of `v`; `true` if there was one.
fn non_finite_into(ring: usize, v: &[Coord2d], issues: &mut Vec<ValidityIssue>) -> bool {
    let before = issues.len();
    for (index, c) in v.iter().enumerate() {
        if !c.x.is_finite() || !c.y.is_finite() {
            issues.push(ValidityIssue::NonFiniteCoordinate { ring, index });
        }
    }
    issues.len() > before
}

/// The ring without the OGC closing vertex, as `LinearRing::new` sees it.
fn without_closing(v: &[Coord2d]) -> &[Coord2d] {
    match v {
        [first, .., last] if first == last => &v[..v.len() - 1],
        _ => v,
    }
}

/// Runs every per-ring check on `coords` as ring `ring`, appending to
/// `issues`, and returns the ring prepared for placement, or `None` if it
/// bounds no region to place (non-finite, too few points, an area that is
/// not a number, degenerate).
fn validate_ring_into(
    ring: usize,
    coords: &CoordSeq,
    issues: &mut Vec<ValidityIssue>,
) -> Option<Prepared> {
    let v = coords.as_slice();
    if non_finite_into(ring, v, issues) {
        return None;
    }
    let v = without_closing(v);
    let n = v.len();

    // Consecutive repeats, including the pair across the closing edge, and
    // the distinct ring with the original index of each kept vertex.
    let mut distinct: Vec<Coord2d> = Vec::with_capacity(n);
    let mut origin: Vec<usize> = Vec::with_capacity(n);
    for (i, &p) in v.iter().enumerate() {
        if n >= 2 && p == v[(i + 1) % n] {
            issues.push(ValidityIssue::DuplicatePoint { ring, index: i });
        }
        if distinct.last() != Some(&p) {
            distinct.push(p);
            origin.push(i);
        }
    }
    // The pair across the closing edge: the last kept vertex equal to the
    // first is the same repeat the loop reported at index n - 1.
    if distinct.len() > 1 && distinct.last() == distinct.first() {
        distinct.pop();
        origin.pop();
    }

    if distinct.len() < 3 {
        issues.push(ValidityIssue::TooFewPoints {
            ring,
            found: distinct.len(),
            required: 3,
        });
        return None;
    }

    // Every condition `LinearRing::new` tests now holds, so this succeeds;
    // the `Err` arm cannot be reached, and is not unwrapped.
    let Ok(lr) = LinearRing::new(CoordSeq::new(distinct)) else {
        return None;
    };
    // As `Polygon::new`: a non-finite area first, then degeneracy, then
    // winding. A ring whose area is not a number has no winding to judge and
    // no region to place.
    let area = lr.signed_area();
    let placeable = if !area.is_finite() {
        issues.push(ValidityIssue::NonFiniteArea { ring });
        false
    } else if ring_is_degenerate(&lr) {
        issues.push(ValidityIssue::DegenerateRing { ring });
        false
    } else {
        let wound_right = if ring == 0 { area > 0.0 } else { area < 0.0 };
        if !wound_right {
            issues.push(ValidityIssue::RingOrientation { ring });
        }
        true
    };
    for (i, j) in self_contacts(lr.coords().as_slice(), false) {
        issues.push(ValidityIssue::SelfIntersection {
            ring,
            edges: (origin[i], origin[j]),
        });
    }
    placeable.then_some(Prepared { ring: lr, origin })
}

#[cfg(test)]
mod between_tests;
#[cfg(test)]
mod tests;
