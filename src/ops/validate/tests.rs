// src/ops/validate/tests.rs

//! Unit tests for `ops::validate`: each check in isolation, the index
//! conventions, and the one property the constructors cannot offer: every
//! issue at once.

#![allow(clippy::unwrap_used)]

use super::{ValidityIssue, validate, validate_line};
use crate::types::{Coord2d, CoordSeq, LinearRing, Polygon};

fn seq(pts: &[(f64, f64)]) -> CoordSeq {
    CoordSeq::new(pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect())
}

/// A counterclockwise unit square, open.
fn ccw() -> CoordSeq {
    seq(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
}

#[test]
fn a_valid_polygon_has_no_issues_and_its_constructor_agrees() {
    let hole = seq(&[(2.0, 2.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0)]);
    let rings = [ccw(), hole.clone()];
    assert!(validate(&rings).is_empty());
    assert!(
        Polygon::new(
            LinearRing::new(ccw()).unwrap(),
            vec![LinearRing::new(hole).unwrap()]
        )
        .is_ok()
    );
    assert!(validate(&[]).is_empty());
}

#[test]
fn the_closing_vertex_is_accepted_as_the_constructor_accepts_it() {
    // Closed input is the constructor's convention too: no DuplicatePoint.
    let closed = seq(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ]);
    assert!(validate(core::slice::from_ref(&closed)).is_empty());
    assert!(LinearRing::new(closed).is_ok());
    // But the closing vertex repeated twice is a repeat across the closing
    // edge, at the last index after one strip. That is the constructor's
    // `ConsecutiveDuplicateVertices { index: 4 }`.
    let twice = seq(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
        (0.0, 0.0),
    ]);
    assert_eq!(
        validate(core::slice::from_ref(&twice)),
        vec![ValidityIssue::DuplicatePoint { ring: 0, index: 4 }]
    );
    assert!(matches!(
        LinearRing::new(twice),
        Err(crate::GeometryError::ConsecutiveDuplicateVertices { index: 4 })
    ));
}

#[test]
fn a_clockwise_exterior_and_a_counterclockwise_hole_are_orientation_issues() {
    let cw = seq(&[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)]);
    assert_eq!(
        validate(core::slice::from_ref(&cw)),
        vec![ValidityIssue::RingOrientation { ring: 0 }]
    );
    let ccw_hole = seq(&[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)]);
    assert_eq!(
        validate(&[ccw(), ccw_hole]),
        vec![ValidityIssue::RingOrientation { ring: 1 }]
    );
    // Both wrong at once: both reported, in ring order.
    assert_eq!(
        validate(&[cw, seq(&[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0)])]),
        vec![
            ValidityIssue::RingOrientation { ring: 0 },
            ValidityIssue::RingOrientation { ring: 1 },
        ]
    );
}

#[test]
fn a_bow_tie_reports_the_crossing_edges() {
    // (0,0) → (10,10) → (10,0) → (0,10): edge 0 crosses edge 2 at (5,5).
    let bow = seq(&[(0.0, 0.0), (10.0, 10.0), (10.0, 0.0), (0.0, 10.0)]);
    let issues = validate(&[bow]);
    // A bow tie has zero signed area (the two lobes cancel exactly), so it
    // is degenerate as well as self-crossing. Both are reported; neither
    // hides the other.
    assert_eq!(
        issues,
        vec![
            ValidityIssue::DegenerateRing { ring: 0 },
            ValidityIssue::SelfIntersection {
                ring: 0,
                edges: (0, 2)
            },
        ]
    );
}

#[test]
fn an_unbalanced_bow_tie_is_a_crossing_with_an_orientation_verdict() {
    // Lobes of different size: net area nonzero, so orientation is judged;
    // the crossing is edge 0 × edge 2.
    let bow = seq(&[(0.0, 0.0), (10.0, 10.0), (10.0, 0.0), (0.0, 4.0)]);
    let issues = validate(&[bow]);
    assert!(
        issues.contains(&ValidityIssue::SelfIntersection {
            ring: 0,
            edges: (0, 2)
        }),
        "{issues:?}"
    );
    assert_eq!(issues.len(), 2, "{issues:?}");
}

#[test]
fn every_self_contact_is_reported_not_only_the_first() {
    // A zigzag whose edge 0 crosses edges 2 and 4: two contacts on one ring.
    //   (0,0) → (20,2) → (20,0) → (10,4) → (0,0)? Build a comb instead:
    // (0,5) (30,5) (30,0) (20,10) (10,0) (5,10): edge 0 (y = 5) is
    // crossed by edges 2, 3 and 4.
    let comb = seq(&[
        (0.0, 5.0),
        (30.0, 5.0),
        (30.0, 0.0),
        (20.0, 10.0),
        (10.0, 0.0),
        (5.0, 10.0),
    ]);
    let crossings: Vec<(usize, usize)> = validate(&[comb])
        .into_iter()
        .filter_map(|i| match i {
            ValidityIssue::SelfIntersection { edges, .. } => Some(edges),
            _ => None,
        })
        .collect();
    assert_eq!(crossings, vec![(0, 2), (0, 3), (0, 4)]);
}

#[test]
fn a_pinch_and_a_t_touch_are_self_intersections_too() {
    // Pinch: (5,5) visited twice. Edges 1 and 4 share the vertex... the
    // contact is reported between the edges that meet at the repeated
    // vertex, by the exact test's own rule.
    let pinch = seq(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (5.0, 5.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (5.0, 5.0),
    ]);
    let issues = validate(&[pinch]);
    assert!(
        issues
            .iter()
            .any(|i| matches!(i, ValidityIssue::SelfIntersection { .. })),
        "{issues:?}"
    );
    // T-touch: vertex (5,0) of edge 2–3 lies on edge 0.
    let t = seq(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (5.0, 0.0),
        (0.0, 10.0),
    ]);
    let issues = validate(&[t]);
    assert!(
        issues
            .iter()
            .any(|i| matches!(i, ValidityIssue::SelfIntersection { .. })),
        "{issues:?}"
    );
}

#[test]
fn too_few_points_counts_distinct_vertices() {
    assert_eq!(
        validate(&[seq(&[(0.0, 0.0), (1.0, 1.0)])]),
        vec![ValidityIssue::TooFewPoints {
            ring: 0,
            found: 2,
            required: 3
        }]
    );
    assert_eq!(
        validate(&[seq(&[(0.0, 0.0)])]),
        vec![ValidityIssue::TooFewPoints {
            ring: 0,
            found: 1,
            required: 3
        }]
    );
    assert_eq!(
        validate(&[seq(&[])]),
        vec![ValidityIssue::TooFewPoints {
            ring: 0,
            found: 0,
            required: 3
        }]
    );
    // Three vertices, one repeated: the repeat and then two distinct.
    assert_eq!(
        validate(&[seq(&[(0.0, 0.0), (0.0, 0.0), (1.0, 1.0)])]),
        vec![
            ValidityIssue::DuplicatePoint { ring: 0, index: 0 },
            ValidityIssue::TooFewPoints {
                ring: 0,
                found: 2,
                required: 3
            },
        ]
    );
}

#[test]
fn duplicates_are_reported_and_the_rest_is_judged_on_the_distinct_ring() {
    // A downstream fixture: (0,0) (0,0) (10,0) (10,10) (0,0), closed, with a
    // repeat at the start. Distinct ring (0,0) (10,0) (10,10): a valid
    // counterclockwise triangle, so the repeat is the only issue.
    let d = seq(&[
        (0.0, 0.0),
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 0.0),
    ]);
    assert_eq!(
        validate(&[d]),
        vec![ValidityIssue::DuplicatePoint { ring: 0, index: 0 }]
    );
    // A repeat inside a clockwise bow tie: repeat, degeneracy and crossing,
    // with the crossing's edge indices in the input's numbering. Edge 3 is
    // the edge leaving input vertex 3, which skips the repeat at 2.
    let d = seq(&[
        (0.0, 0.0),
        (10.0, 10.0),
        (10.0, 10.0),
        (10.0, 0.0),
        (0.0, 10.0),
    ]);
    assert_eq!(
        validate(&[d]),
        vec![
            ValidityIssue::DuplicatePoint { ring: 0, index: 1 },
            ValidityIssue::DegenerateRing { ring: 0 },
            ValidityIssue::SelfIntersection {
                ring: 0,
                edges: (0, 3)
            },
        ]
    );
}

#[test]
fn a_collinear_ring_is_degenerate_and_gets_no_orientation_verdict() {
    // Three collinear vertices: zero area, and the closing edge (2) folds
    // back over edges 1 and 0: two fold-back contacts, the second being
    // the consecutive pair across the closing edge, reported as (2, 0).
    let flat = seq(&[(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]);
    assert_eq!(
        validate(&[flat]),
        vec![
            ValidityIssue::DegenerateRing { ring: 0 },
            ValidityIssue::SelfIntersection {
                ring: 0,
                edges: (1, 2)
            },
            ValidityIssue::SelfIntersection {
                ring: 0,
                edges: (2, 0)
            },
        ]
    );
    // A nearly collinear ring, a sliver about 1e-17 of its length thick,
    // has an exact signed area of −1e−15: nonzero, so it is not degenerate.
    // It is wound clockwise in the exterior position, so the issue is its
    // orientation. No contact. (A test with a tolerance would call it
    // degenerate.)
    let sliver = seq(&[(0.0, 0.0), (10.0, 1e-16), (20.0, 0.0)]);
    assert_eq!(
        validate(&[sliver]),
        vec![ValidityIssue::RingOrientation { ring: 0 }]
    );
}

#[test]
fn a_non_finite_coordinate_is_reported_per_vertex_and_nothing_else_is() {
    let bad = seq(&[
        (0.0, 0.0),
        (f64::NAN, 0.0),
        (10.0, f64::INFINITY),
        (0.0, 10.0),
    ]);
    assert_eq!(
        validate(&[bad]),
        vec![
            ValidityIssue::NonFiniteCoordinate { ring: 0, index: 1 },
            ValidityIssue::NonFiniteCoordinate { ring: 0, index: 2 },
        ]
    );
}

#[test]
fn issues_name_their_ring_and_come_in_ring_order() {
    let rings = [
        seq(&[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)]), // CW exterior
        seq(&[(2.0, 2.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0)]),     // fine hole
        seq(&[(6.0, 6.0), (7.0, 7.0)]),                             // too few
    ];
    let issues = validate(&rings);
    assert_eq!(
        issues.iter().map(ValidityIssue::ring).collect::<Vec<_>>(),
        vec![0, 2]
    );
    assert_eq!(
        issues[1],
        ValidityIssue::TooFewPoints {
            ring: 2,
            found: 2,
            required: 3
        }
    );
}

#[test]
fn lines_need_two_distinct_vertices_and_no_repeats_and_may_cross() {
    assert!(validate_line(&seq(&[(0.0, 0.0), (5.0, 5.0)])).is_empty());
    assert_eq!(
        validate_line(&seq(&[(0.0, 0.0)])),
        vec![ValidityIssue::TooFewPoints {
            ring: 0,
            found: 1,
            required: 2
        }]
    );
    assert_eq!(
        validate_line(&seq(&[(0.0, 0.0), (0.0, 0.0)])),
        vec![
            ValidityIssue::DuplicatePoint { ring: 0, index: 0 },
            ValidityIssue::TooFewPoints {
                ring: 0,
                found: 1,
                required: 2
            },
        ]
    );
    // A closed line is fine; a self-crossing line is fine.
    assert!(validate_line(&seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])).is_empty());
    assert!(validate_line(&seq(&[(0.0, 0.0), (10.0, 10.0), (10.0, 0.0), (0.0, 10.0)])).is_empty());
    assert_eq!(
        validate_line(&seq(&[(0.0, f64::NAN), (1.0, 0.0)])),
        vec![ValidityIssue::NonFiniteCoordinate { ring: 0, index: 0 }]
    );
}
