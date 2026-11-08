// tests/validity.rs

//! Validation: one hand case for each reason `validate` and
//! `validate_multipolygon` report, and the agreement of the constructors,
//! `validate`, `validity_issues` and `is_valid` on generated polygons.

#![allow(clippy::unwrap_used, clippy::cast_possible_truncation)]

mod common;

use common::{Rng, c, jagged, ring};
use geometres::ops::{ValidityIssue, is_valid, validate, validate_multipolygon, validity_issues};
use geometres::types::{CoordSeq, Polygon};

fn seq(pts: &[(f64, f64)]) -> CoordSeq {
    CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())
}

fn square(x0: f64, y0: f64, s: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)]
}

fn square_cw(x0: f64, y0: f64, s: f64) -> Vec<(f64, f64)> {
    let mut v = square(x0, y0, s);
    v.reverse();
    v
}

/// The issues of `rings`, which must include one matching `want`.
fn reports(rings: &[Vec<(f64, f64)>], want: impl Fn(&ValidityIssue) -> bool, what: &str) {
    let rings: Vec<CoordSeq> = rings.iter().map(|r| seq(r)).collect();
    let issues = validate(&rings);
    assert!(issues.iter().any(want), "{what}: got {issues:?}");
}

#[test]
fn a_valid_polygon_has_no_issues() {
    let rings = [seq(&square(0.0, 0.0, 10.0)), seq(&square_cw(2.0, 2.0, 2.0))];
    assert_eq!(validate(&rings), []);
}

#[test]
fn each_ring_issue_is_reported() {
    reports(
        &[vec![(0.0, 0.0), (f64::NAN, 0.0), (1.0, 1.0)]],
        |i| matches!(i, ValidityIssue::NonFiniteCoordinate { ring: 0, index: 1 }),
        "non-finite coordinate",
    );
    reports(
        &[vec![(0.0, 0.0), (1.0, 0.0)]],
        |i| {
            matches!(
                i,
                ValidityIssue::TooFewPoints {
                    ring: 0,
                    found: 2,
                    required: 3
                }
            )
        },
        "too few points",
    );
    reports(
        &[vec![(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 1.0)]],
        |i| matches!(i, ValidityIssue::DuplicatePoint { ring: 0, index: 1 }),
        "duplicate point",
    );
    reports(
        &[vec![(0.0, 0.0), (1e300, 0.0), (1e300, 1e300), (0.0, 1e300)]],
        |i| matches!(i, ValidityIssue::NonFiniteArea { ring: 0 }),
        "area overflows",
    );
    reports(
        &[vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)]],
        |i| matches!(i, ValidityIssue::DegenerateRing { ring: 0 }),
        "collinear ring",
    );
    reports(
        &[square_cw(0.0, 0.0, 4.0)],
        |i| matches!(i, ValidityIssue::RingOrientation { ring: 0 }),
        "clockwise exterior",
    );
    reports(
        &[square(0.0, 0.0, 10.0), square(2.0, 2.0, 2.0)],
        |i| matches!(i, ValidityIssue::RingOrientation { ring: 1 }),
        "counterclockwise hole",
    );
    // A bow tie with unequal lobes: nonzero area, crossing edges.
    reports(
        &[vec![(0.0, 0.0), (4.0, 0.0), (0.0, 8.0), (4.0, 6.0)]],
        |i| matches!(i, ValidityIssue::SelfIntersection { ring: 0, .. }),
        "self-intersection",
    );
}

#[test]
fn each_placement_issue_is_reported() {
    reports(
        &[square(0.0, 0.0, 4.0), square_cw(10.0, 10.0, 1.0)],
        |i| matches!(i, ValidityIssue::HoleOutsideExterior { ring: 1 }),
        "hole outside the exterior",
    );
    reports(
        &[
            square(0.0, 0.0, 10.0),
            square_cw(1.0, 1.0, 6.0),
            square_cw(2.0, 2.0, 1.0),
        ],
        |i| matches!(i, ValidityIssue::NestedHoles { ring: 2, inside: 1 }),
        "nested holes",
    );
    reports(
        &[
            square(0.0, 0.0, 10.0),
            square_cw(1.0, 1.0, 4.0),
            square_cw(3.0, 3.0, 4.0),
        ],
        |i| {
            matches!(
                i,
                ValidityIssue::RingsCross {
                    ring: 2,
                    other: 1,
                    ..
                }
            )
        },
        "crossing holes",
    );
    reports(
        &[
            square(0.0, 0.0, 10.0),
            square_cw(1.0, 1.0, 2.0),
            square_cw(1.0, 1.0, 2.0),
        ],
        |i| matches!(i, ValidityIssue::DuplicateRings { ring: 2, other: 1 }),
        "the same hole twice",
    );
    // Two diamond holes in a row, touching each other at a vertex and the
    // exterior at one vertex each: the touches close a cycle that cuts the
    // interior into a top and a bottom.
    reports(
        &[
            square(0.0, 0.0, 4.0),
            vec![(0.0, 2.0), (1.0, 3.0), (2.0, 2.0), (1.0, 1.0)],
            vec![(2.0, 2.0), (3.0, 3.0), (4.0, 2.0), (3.0, 1.0)],
        ],
        |i| matches!(i, ValidityIssue::DisconnectedInterior { .. }),
        "disconnected interior",
    );
}

#[test]
fn overlapping_parts_of_a_multipolygon_are_reported() {
    let parts = vec![
        vec![seq(&square(0.0, 0.0, 4.0))],
        vec![seq(&square(2.0, 2.0, 4.0))],
    ];
    let issues = validate_multipolygon(&parts);
    assert!(
        issues
            .iter()
            .any(|(_, i)| matches!(i, ValidityIssue::PartsOverlap { part: 1, other: 0 })),
        "{issues:?}"
    );
    let apart = vec![
        vec![seq(&square(0.0, 0.0, 1.0))],
        vec![seq(&square(1.0, 1.0, 1.0))],
    ];
    assert_eq!(validate_multipolygon(&apart), []);
}

#[test]
fn a_constructible_polygon_with_a_disconnected_interior_is_not_valid() {
    // `Polygon::new` admits a disconnected interior; `validity_issues` and
    // `is_valid` report it.
    let p = Polygon::new(
        ring(&square(0.0, 0.0, 4.0)),
        vec![
            ring(&[(0.0, 2.0), (1.0, 3.0), (2.0, 2.0), (1.0, 1.0)]),
            ring(&[(2.0, 2.0), (3.0, 3.0), (4.0, 2.0), (3.0, 1.0)]),
        ],
    )
    .unwrap();
    assert!(!is_valid(&p));
    assert!(
        validity_issues(&p)
            .iter()
            .any(|li| matches!(li.issue, ValidityIssue::DisconnectedInterior { .. }))
    );
}

#[test]
fn constructors_validate_and_validity_issues_agree_on_generated_polygons() {
    let mut rng = Rng::new(0x5EED_0301);
    for _ in 0..500 {
        let (x, y, n) = (
            rng.int(-50, 50),
            rng.int(-50, 50),
            2 + rng.below(30) as usize,
        );
        let p = jagged(&mut rng, x, y, n, 20);
        let rings: Vec<CoordSeq> = std::iter::once(p.exterior())
            .chain(p.holes())
            .map(|r| r.coords().clone())
            .collect();
        assert_eq!(validate(&rings), [], "{p:?}");
        assert!(is_valid(&p) && validity_issues(&p).is_empty());
        // Reversing the exterior makes it clockwise: the constructor and
        // `validate` refuse it for the same reason.
        let mut rev: Vec<_> = p.exterior().coords().iter().copied().collect();
        rev.reverse();
        let rev = CoordSeq::new(rev);
        assert!(
            validate(std::slice::from_ref(&rev))
                .contains(&ValidityIssue::RingOrientation { ring: 0 })
        );
        assert!(Polygon::new_simple(geometres::types::LinearRing::new(rev).unwrap()).is_err());
    }
}
