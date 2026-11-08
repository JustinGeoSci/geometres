// src/ops/validate/between_tests.rs

//! Tests for the issues between rings and between parts, the disconnected
//! interior, the constructed-value checks, and the agreement of `validate`
//! and `validate_multipolygon` with the constructors on random material.

#![allow(clippy::unwrap_used)]
// test-lint-posture: see the workspace lint comment in Cargo.toml
#![allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

use super::{
    LocatedIssue, ValidityIssue, is_valid, validate, validate_multipolygon, validity_issues,
};
use crate::GeometryError;
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiPolygon, Point,
    Polygon,
};

use ValidityIssue as V;

fn seq(pts: &[(f64, f64)]) -> CoordSeq {
    CoordSeq::new(pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect())
}

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(seq(pts)).unwrap()
}

/// The 8 × 8 square exterior, counterclockwise.
fn ext() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (0.0, 8.0)]
}

/// A clockwise axis-aligned square hole.
fn hole(x0: f64, y0: f64, s: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x0, y0 + s), (x0 + s, y0 + s), (x0 + s, y0)]
}

fn rings(rs: &[Vec<(f64, f64)>]) -> Vec<CoordSeq> {
    rs.iter().map(|r| seq(r)).collect()
}

fn polygon(rs: &[Vec<(f64, f64)>]) -> Result<Polygon, GeometryError> {
    let mut it = rs.iter().map(|r| ring(r));
    let e = it.next().unwrap();
    Polygon::new(e, it.collect())
}

// ===========================================================================
// Placement issues, by hand
// ===========================================================================

#[test]
fn a_hole_outside_or_surrounding_the_exterior_is_reported() {
    let outside = rings(&[ext(), hole(20.0, 20.0, 2.0)]);
    assert_eq!(validate(&outside), vec![V::HoleOutsideExterior { ring: 1 }]);
    // A clockwise ring around the exterior: outside it, as a region.
    let around = rings(&[ext(), hole(-2.0, -2.0, 12.0)]);
    assert_eq!(validate(&around), vec![V::HoleOutsideExterior { ring: 1 }]);
    // Touching the exterior from outside at a vertex is still outside.
    let touching = rings(&[ext(), hole(8.0, 8.0, 2.0)]);
    assert_eq!(
        validate(&touching),
        vec![V::HoleOutsideExterior { ring: 1 }]
    );
    assert_eq!(
        polygon(&[ext(), hole(20.0, 20.0, 2.0)]).unwrap_err(),
        GeometryError::HoleNotInsideExterior { hole_index: 0 }
    );
}

#[test]
fn a_hole_inside_another_hole_is_reported_with_both_indices() {
    let rs = rings(&[ext(), hole(1.0, 1.0, 6.0), hole(3.0, 3.0, 1.0)]);
    assert_eq!(validate(&rs), vec![V::NestedHoles { ring: 2, inside: 1 }]);
    let rs = rings(&[ext(), hole(3.0, 3.0, 1.0), hole(1.0, 1.0, 6.0)]);
    assert_eq!(validate(&rs), vec![V::NestedHoles { ring: 1, inside: 2 }]);
    assert_eq!(
        polygon(&[ext(), hole(3.0, 3.0, 1.0), hole(1.0, 1.0, 6.0)]).unwrap_err(),
        GeometryError::OverlappingHoles {
            first: 0,
            second: 1
        }
    );
}

#[test]
fn a_hole_crossing_the_exterior_names_a_properly_crossing_edge_pair() {
    // Hole edges: 0 up the left side x = 6, 1 across y = 10 (outside), 2 down
    // x = 10, 3 back along y = 2. Edge 1 of the exterior (x = 8, upward)
    // properly crosses the hole's edge 3... the first proper crossing in
    // hole-edge order is the hole's edge 1? No: edge 1 is at y = 10, above the
    // exterior. The first is the hole's edge 3 (y = 2) against exterior edge 1.
    let h = vec![(6.0, 2.0), (6.0, 10.0), (10.0, 10.0), (10.0, 2.0)];
    let issues = validate(&rings(&[ext(), h]));
    assert_eq!(
        issues,
        vec![V::RingsCross {
            ring: 1,
            other: 0,
            edges: (0, 2)
        }]
    );
}

#[test]
fn a_hole_passing_through_the_exterior_at_vertices_is_a_crossing() {
    // An L-shaped exterior, and a hole leaving it through two of its
    // boundary vertices, with no proper crossing.
    let l = vec![
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 2.0),
        (2.0, 2.0),
        (2.0, 4.0),
        (0.0, 4.0),
    ];
    // Clockwise triangle from the reflex corner's neighbors out into the
    // notch: (4,2) and (2,4) are exterior vertices, (4,4) is outside.
    let h = vec![(1.0, 1.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0)];
    let issues = validate(&rings(&[l.clone(), h.clone()]));
    assert_eq!(issues.len(), 1, "{issues:?}");
    let V::RingsCross { ring, other, edges } = issues[0] else {
        panic!("{issues:?}");
    };
    assert_eq!((ring, other), (1, 0));
    // The witness edge of the hole leaves an exterior boundary point, and
    // the exterior edge named holds that point.
    let hv = &h[edges.0];
    let (a, b) = (l[edges.1], l[(edges.1 + 1) % l.len()]);
    let on = (hv.0 - a.0) * (b.1 - a.1) == (hv.1 - a.1) * (b.0 - a.0);
    assert!(on, "{issues:?}");
    assert_eq!(
        polygon(&[l, h]).unwrap_err(),
        GeometryError::HoleNotInsideExterior { hole_index: 0 }
    );
}

#[test]
fn holes_sharing_a_stretch_name_the_edges_along_each_other() {
    // Two holes side by side sharing the segment x = 4, y in [2, 4].
    let a = hole(2.0, 2.0, 2.0); // edges: 0 left, 1 top, 2 right (x = 4), 3 bottom
    let b = hole(4.0, 2.0, 2.0); // edges: 0 left (x = 4), ...
    let issues = validate(&rings(&[ext(), a, b]));
    assert_eq!(
        issues,
        vec![V::RingsCross {
            ring: 2,
            other: 1,
            edges: (0, 2)
        }]
    );
}

#[test]
fn duplicate_rings_are_reported_however_they_are_listed() {
    // The same hole from another start vertex.
    let a = hole(2.0, 2.0, 2.0);
    let mut b = a.clone();
    b.rotate_left(2);
    assert_eq!(
        validate(&rings(&[ext(), a.clone(), b])),
        vec![V::DuplicateRings { ring: 2, other: 1 }]
    );
    // With an extra collinear vertex.
    let c = vec![(2.0, 2.0), (2.0, 3.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0)];
    assert_eq!(
        validate(&rings(&[ext(), a, c])),
        vec![V::DuplicateRings { ring: 2, other: 1 }]
    );
    // A hole equal to the exterior, reversed to be clockwise.
    let mut e = ext();
    e.reverse();
    assert_eq!(
        validate(&rings(&[ext(), e])),
        vec![V::DuplicateRings { ring: 1, other: 0 }]
    );
}

#[test]
fn between_ring_indices_are_input_indices() {
    // The crossing hole of the test above, given closed and with a repeat
    // at vertex 1: the repeat is reported, and the witness edges are
    // counted in the ring as given.
    let h = vec![
        (6.0, 2.0),
        (6.0, 10.0),
        (6.0, 10.0),
        (10.0, 10.0),
        (10.0, 2.0),
        (6.0, 2.0),
    ];
    let issues = validate(&rings(&[ext(), h]));
    assert_eq!(
        issues,
        vec![
            V::DuplicatePoint { ring: 1, index: 1 },
            V::RingsCross {
                ring: 1,
                other: 0,
                edges: (0, 2)
            }
        ]
    );
}

#[test]
fn an_overflowing_area_is_reported_as_the_constructor_refuses_it() {
    let huge = vec![(0.0, 0.0), (1e200, 0.0), (1e200, 1e200), (0.0, 1e200)];
    assert_eq!(
        validate(&rings(core::slice::from_ref(&huge))),
        vec![V::NonFiniteArea { ring: 0 }]
    );
    assert_eq!(
        polygon(core::slice::from_ref(&huge)).unwrap_err(),
        GeometryError::NonFiniteResult
    );
}

#[test]
fn placement_is_skipped_for_rings_that_bound_nothing() {
    // A collinear ring folds back on itself (its self-contacts are
    // reported) and is degenerate; it draws no placement issue, while the
    // ring after it still does.
    let collinear = vec![(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)];
    let issues = validate(&rings(&[ext(), collinear, hole(20.0, 20.0, 1.0)]));
    assert_eq!(issues[0], V::DegenerateRing { ring: 1 });
    assert_eq!(issues.last(), Some(&V::HoleOutsideExterior { ring: 2 }));
    assert!(
        issues[1..issues.len() - 1]
            .iter()
            .all(|i| matches!(i, V::SelfIntersection { ring: 1, .. })),
        "{issues:?}"
    );
}

// ===========================================================================
// A disconnected interior
// ===========================================================================

/// The rings, the expected disconnection (if any), and that the polygon
/// constructs.
fn check_connectivity(rs: &[Vec<(f64, f64)>], expected: Option<(usize, usize)>) {
    let p = polygon(rs).expect("constructible");
    let raw = validate(&rings(rs));
    let got = validity_issues(&p);
    match expected {
        None => {
            assert!(raw.is_empty(), "{raw:?}");
            assert!(got.is_empty(), "{got:?}");
            assert!(is_valid(&p));
        }
        Some((ring, vertex)) => {
            let issue = V::DisconnectedInterior { ring, vertex };
            assert_eq!(raw, vec![issue]);
            assert_eq!(
                got,
                vec![LocatedIssue {
                    path: vec![],
                    issue
                }]
            );
            assert!(!is_valid(&p));
        }
    }
}

#[test]
fn a_hole_touching_the_exterior_at_two_points_disconnects() {
    // A diamond touching the left and right sides at (0,4) and (8,4).
    let diamond = vec![(0.0, 4.0), (4.0, 6.0), (8.0, 4.0), (4.0, 2.0)];
    // Touch points in discovery order: hole vertices on the exterior, (0,4)
    // then (8,4). The second closes the cycle; it is the hole's vertex 2
    // (not an exterior vertex).
    check_connectivity(&[ext(), diamond], Some((1, 2)));
}

#[test]
fn two_holes_touching_each_other_and_the_exterior_disconnect() {
    let a = vec![(0.0, 4.0), (2.0, 6.0), (4.0, 4.0), (2.0, 2.0)];
    let b = vec![(4.0, 4.0), (6.0, 6.0), (8.0, 4.0), (6.0, 2.0)];
    // Incidences: (ext,(0,4)), (a,(0,4)); (ext,(8,4)), (b,(8,4)); then
    // (a,(4,4)), and (b,(4,4)) closes the cycle. (4,4) is vertex 2 of a.
    check_connectivity(&[ext(), a, b], Some((1, 2)));
}

#[test]
fn a_chain_of_three_holes_across_the_polygon_disconnects() {
    let a = vec![(0.0, 4.0), (1.0, 5.0), (2.0, 4.0), (1.0, 3.0)];
    let b = vec![(2.0, 4.0), (4.0, 5.0), (6.0, 4.0), (4.0, 3.0)];
    let c = vec![(6.0, 4.0), (7.0, 5.0), (8.0, 4.0), (7.0, 3.0)];
    // The last incidence to be added closes the cycle: (c, (6,4)), where
    // (6,4) is vertex 2 of b, the lowest ring holding it as a vertex.
    check_connectivity(&[ext(), a, b, c], Some((2, 2)));
}

#[test]
fn a_hole_touching_the_exterior_twice_inside_its_edges_disconnects() {
    // Both touches are hole vertices in the interiors of exterior edges.
    let h = vec![(0.0, 3.0), (3.0, 8.0), (5.0, 5.0)];
    // Clockwise? (0,3) -> (3,8) -> (5,5): orient = (3)(2) - (5)(5) < 0. Yes.
    check_connectivity(&[ext(), h], Some((1, 1)));
}

#[test]
fn touches_that_do_not_disconnect_are_not_reported() {
    // One touch with the exterior.
    check_connectivity(
        &[ext(), vec![(0.0, 4.0), (2.0, 6.0), (4.0, 4.0), (2.0, 2.0)]],
        None,
    );
    // Two holes touching each other only.
    check_connectivity(
        &[
            ext(),
            vec![(1.0, 4.0), (2.0, 6.0), (4.0, 4.0), (2.0, 2.0)],
            vec![(4.0, 4.0), (6.0, 6.0), (7.0, 4.0), (6.0, 2.0)],
        ],
        None,
    );
    // A chain of three holes touching the exterior at one end only.
    check_connectivity(
        &[
            ext(),
            vec![(0.0, 4.0), (1.0, 5.0), (2.0, 4.0), (1.0, 3.0)],
            vec![(2.0, 4.0), (4.0, 5.0), (6.0, 4.0), (4.0, 3.0)],
            vec![(6.0, 4.0), (6.5, 5.0), (7.0, 4.0), (6.5, 3.0)],
        ],
        None,
    );
}

#[test]
fn three_rings_meeting_at_one_point_do_not_disconnect() {
    // The exterior and two holes all through (0,4): a star in the touch
    // graph, not the triangle a ring-pair graph would see.
    let a = vec![(0.0, 4.0), (1.0, 7.0), (3.0, 6.0)];
    let b = vec![(0.0, 4.0), (3.0, 2.0), (1.0, 1.0)];
    check_connectivity(&[ext(), a.clone(), b.clone()], None);
    // A third hole through the same point is still a star.
    let c = vec![(0.0, 4.0), (4.0, 5.0), (4.0, 3.0)];
    check_connectivity(&[ext(), a, b, c], None);
}

#[test]
fn disconnection_is_judged_only_when_the_rings_are_otherwise_valid() {
    // A wrongly wound exterior still has its disconnection reported.
    let mut cw = ext();
    cw.reverse();
    let diamond = vec![(0.0, 4.0), (4.0, 6.0), (8.0, 4.0), (4.0, 2.0)];
    let issues = validate(&rings(&[cw, diamond.clone()]));
    assert!(
        issues.contains(&V::RingOrientation { ring: 0 }),
        "{issues:?}"
    );
    assert!(
        issues
            .iter()
            .any(|i| matches!(i, V::DisconnectedInterior { .. })),
        "{issues:?}"
    );
    // A self-intersecting hole: the proof's premise fails, so no verdict.
    let bow = vec![(1.0, 1.0), (3.0, 3.0), (3.0, 1.0), (1.0, 3.0)];
    let issues = validate(&rings(&[ext(), diamond, bow]));
    assert!(
        !issues
            .iter()
            .any(|i| matches!(i, V::DisconnectedInterior { .. })),
        "{issues:?}"
    );
}

#[test]
fn the_reported_vertex_is_counted_in_the_ring_as_given() {
    // The diamond, closed and with a repeated first vertex: its vertex 2 in
    // the stored ring is vertex 3 as given.
    let diamond = vec![
        (0.0, 4.0),
        (0.0, 4.0),
        (4.0, 6.0),
        (8.0, 4.0),
        (4.0, 2.0),
        (0.0, 4.0),
    ];
    let issues = validate(&rings(&[ext(), diamond]));
    assert_eq!(
        issues,
        vec![V::DuplicatePoint { ring: 1, index: 0 },],
        "a repeat is a constructor issue, so connectivity is not judged"
    );
    // Closed only: indices unchanged.
    let closed = vec![(0.0, 4.0), (4.0, 6.0), (8.0, 4.0), (4.0, 2.0), (0.0, 4.0)];
    assert_eq!(
        validate(&rings(&[ext(), closed])),
        vec![V::DisconnectedInterior { ring: 1, vertex: 2 }]
    );
}

// ===========================================================================
// Constructed values
// ===========================================================================

#[test]
fn constructed_checks_report_self_intersection_and_disconnection_with_paths() {
    // A polygon whose exterior pinches: constructible, but not valid.
    let pinched = polygon(&[vec![
        (0.0, 0.0),
        (4.0, 0.0),
        (2.0, 2.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (2.0, 2.0),
    ]])
    .unwrap();
    let got = validity_issues(&pinched);
    assert!(!got.is_empty());
    assert!(
        got.iter()
            .all(|i| i.path.is_empty() && matches!(i.issue, V::SelfIntersection { ring: 0, .. }))
    );

    let disconnected =
        polygon(&[ext(), vec![(0.0, 4.0), (4.0, 6.0), (8.0, 4.0), (4.0, 2.0)]]).unwrap();
    let square = polygon(&[vec![(20.0, 0.0), (22.0, 0.0), (22.0, 2.0), (20.0, 2.0)]]).unwrap();
    let mp = MultiPolygon::new(vec![square.clone(), disconnected.clone()]).unwrap();
    assert_eq!(
        validity_issues(&mp),
        vec![LocatedIssue {
            path: vec![1],
            issue: V::DisconnectedInterior { ring: 1, vertex: 2 }
        }]
    );

    let inner = GeometryCollection::new(vec![
        Geometry::Point(Point::new(Coord2d::new(0.0, 0.0)).unwrap()),
        Geometry::MultiPolygon(mp.clone()),
    ])
    .unwrap();
    let gc = GeometryCollection::new(vec![
        Geometry::Polygon(square),
        Geometry::GeometryCollection(inner),
        Geometry::Polygon(pinched),
    ])
    .unwrap();
    let got = validity_issues(&gc);
    assert_eq!(
        got[0],
        LocatedIssue {
            path: vec![1, 1, 1],
            issue: V::DisconnectedInterior { ring: 1, vertex: 2 }
        }
    );
    assert!(got[1..].iter().all(|i| i.path == vec![2]));
    assert!(!is_valid(&gc));
    assert!(!is_valid(&Geometry::GeometryCollection(gc)));
}

#[test]
fn points_lines_and_empty_values_are_always_valid() {
    let crossing = LineString::new(seq(&[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0)])).unwrap();
    assert!(is_valid(&crossing));
    assert!(is_valid(&Point::new(Coord2d::new(1.0, 1.0)).unwrap()));
    assert!(is_valid(&MultiPolygon::empty()));
    assert!(is_valid(&GeometryCollection::empty()));
    assert!(validity_issues(&polygon(&[ext(), hole(2.0, 2.0, 2.0)]).unwrap()).is_empty());
}

// ===========================================================================
// Multipolygon material
// ===========================================================================

fn square(x0: f64, y0: f64, s: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)]
}

#[test]
fn overlapping_and_nested_parts_are_told_apart() {
    let overlap = vec![
        rings(&[square(0.0, 0.0, 4.0)]),
        rings(&[square(2.0, 2.0, 4.0)]),
    ];
    assert_eq!(
        validate_multipolygon(&overlap),
        vec![(1, V::PartsOverlap { part: 1, other: 0 })]
    );
    let nested = vec![
        rings(&[square(0.0, 0.0, 8.0)]),
        rings(&[square(2.0, 2.0, 2.0)]),
    ];
    assert_eq!(
        validate_multipolygon(&nested),
        vec![(1, V::NestedParts { part: 1, inside: 0 })]
    );
    let nested_first = vec![
        rings(&[square(2.0, 2.0, 2.0)]),
        rings(&[square(0.0, 0.0, 8.0)]),
    ];
    assert_eq!(
        validate_multipolygon(&nested_first),
        vec![(0, V::NestedParts { part: 0, inside: 1 })]
    );
    // An island in a hole, touching it at a point: valid.
    let island = vec![
        rings(&[square(0.0, 0.0, 8.0), hole(2.0, 2.0, 4.0)]),
        rings(&[vec![(2.0, 2.0), (4.0, 3.0), (3.0, 4.0)]]),
    ];
    assert!(validate_multipolygon(&island).is_empty());
}

/// A part whose exterior lies inside another part's exterior but meets
/// one of its holes along a stretch, or crosses it, is not nested: its
/// interior is partly or wholly in the hole. It is an overlap.
#[test]
fn a_part_meeting_a_hole_is_an_overlap_not_nested() {
    // In the hole, sharing its left and right sides.
    let in_hole = vec![
        rings(&[square(10.0, 10.0, 80.0), hole(30.0, 30.0, 40.0)]),
        rings(&[vec![(30.0, 40.0), (70.0, 40.0), (70.0, 60.0), (30.0, 60.0)]]),
    ];
    assert_eq!(
        validate_multipolygon(&in_hole),
        vec![(1, V::PartsOverlap { part: 1, other: 0 })]
    );
    // Crossing the hole: vertices on the hole's vertices, edges across it.
    let crossing = vec![
        rings(&[
            vec![(20.0, 20.0), (420.0, 20.0), (420.0, 380.0), (20.0, 380.0)],
            vec![
                (220.0, 340.0),
                (300.0, 240.0),
                (340.0, 60.0),
                (140.0, 100.0),
                (60.0, 200.0),
                (180.0, 240.0),
            ],
        ]),
        rings(&[vec![(60.0, 200.0), (340.0, 60.0), (220.0, 340.0)]]),
    ];
    assert_eq!(
        validate_multipolygon(&crossing),
        vec![(1, V::PartsOverlap { part: 1, other: 0 })]
    );
}

#[test]
fn part_issues_carry_their_part_and_unbuildable_parts_are_not_paired() {
    let mut cw = square(0.0, 0.0, 4.0);
    cw.reverse();
    let parts = vec![
        rings(&[square(10.0, 0.0, 4.0)]),
        rings(&[cw]),
        rings(&[square(11.0, 1.0, 1.0)]),
    ];
    assert_eq!(
        validate_multipolygon(&parts),
        vec![
            (1, V::RingOrientation { ring: 0 }),
            (2, V::NestedParts { part: 2, inside: 0 })
        ]
    );
}

// ===========================================================================
// Agreement with the constructors, on random material
// ===========================================================================

/// SplitMix64: a fixed, dependency-free generator, so every run sees the
/// same cases.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn chance(&mut self, percent: u64) -> bool {
        self.below(100) < percent
    }
}

/// A star-shaped ring around `(cx, cy)`: vertices at increasing angles on
/// 16 sectors, radii in `rmin..=rmax`, rounded to the integer grid. The
/// grid makes touches, shared stretches, collinear vertices and repeats
/// common, which is where the checks are exercised.
fn star(rng: &mut Rng, cx: i64, cy: i64, rmin: u64, rmax: u64, clockwise: bool) -> Vec<(f64, f64)> {
    let k = 3 + rng.below(4) as usize;
    let mut sectors: Vec<u64> = (0..16).collect();
    for i in 0..16 {
        let j = i + rng.below(16 - i as u64) as usize;
        sectors.swap(i, j);
    }
    let mut chosen: Vec<u64> = sectors[..k].to_vec();
    chosen.sort_unstable();
    let mut pts: Vec<(f64, f64)> = chosen
        .iter()
        .map(|&s| {
            let theta = s as f64 * core::f64::consts::TAU / 16.0;
            let r = (rmin + rng.below(rmax - rmin + 1)) as f64;
            (
                (cx as f64 + (r * theta.cos()).round()),
                (cy as f64 + (r * theta.sin()).round()),
            )
        })
        .collect();
    if clockwise {
        pts.reverse();
    }
    pts
}

/// One random ring set: an exterior (a star, or the 16 × 16 square) and up
/// to three holes, with the occasional wrong winding, repeat, pinch or
/// closing vertex. Hole vertices are often snapped onto a vertex of an
/// earlier ring or onto the square's boundary, which makes touches, and so
/// disconnected interiors, common.
fn random_rings(rng: &mut Rng) -> Vec<CoordSeq> {
    let square_exterior = rng.chance(50);
    let mut rs = if square_exterior {
        vec![vec![(0.0, 0.0), (16.0, 0.0), (16.0, 16.0), (0.0, 16.0)]]
    } else {
        let flip = rng.chance(4);
        vec![star(rng, 8, 8, 4, 8, flip)]
    };
    if square_exterior && rng.chance(4) {
        rs[0].reverse();
    }
    for _ in 0..rng.below(4) {
        let (cx, cy) = (2 + rng.below(13) as i64, 2 + rng.below(13) as i64);
        let clockwise = !rng.chance(4);
        let mut h = star(rng, cx, cy, 1, 3, clockwise);
        for _ in 0..2 {
            if rng.chance(45) {
                let i = rng.below(h.len() as u64) as usize;
                let src = rng.below(rs.len() as u64) as usize;
                h[i] = if square_exterior && src == 0 && rng.chance(50) {
                    // A point on the square's boundary, usually mid-edge.
                    let t = rng.below(17) as f64;
                    match rng.below(4) {
                        0 => (t, 0.0),
                        1 => (16.0, t),
                        2 => (t, 16.0),
                        _ => (0.0, t),
                    }
                } else {
                    let r = &rs[src];
                    r[rng.below(r.len() as u64) as usize]
                };
            }
        }
        rs.push(h);
    }
    for r in &mut rs {
        if rng.chance(3) && !r.is_empty() {
            let i = rng.below(r.len() as u64) as usize;
            r.insert(i, r[i]);
        }
        if rng.chance(3) && r.len() >= 4 {
            // A pinch: a vertex visited again two steps later.
            let i = rng.below(r.len() as u64 - 2) as usize;
            r.insert(i + 2, r[i]);
        }
        if rng.chance(5) {
            r.push(r[0]);
        }
    }
    rs.iter().map(|r| seq(r)).collect()
}

/// The order in which the constructors meet each kind of refused issue:
/// `LinearRing::new` on every ring, then `Polygon::new`'s ring checks, then
/// its placement checks.
fn phase(i: &ValidityIssue) -> u8 {
    match i {
        V::NonFiniteCoordinate { .. } | V::TooFewPoints { .. } | V::DuplicatePoint { .. } => 0,
        V::NonFiniteArea { .. } | V::DegenerateRing { .. } | V::RingOrientation { .. } => 1,
        _ => 2,
    }
}

/// Builds the polygon as a caller would: each ring, in order, then the
/// polygon. The error, with the ring it came from for `LinearRing::new`.
fn construct(rs: &[CoordSeq]) -> Result<Polygon, (Option<usize>, GeometryError)> {
    let mut lrs = Vec::new();
    for (k, c) in rs.iter().enumerate() {
        lrs.push(LinearRing::new(c.clone()).map_err(|e| (Some(k), e))?);
    }
    let e = lrs.remove(0);
    Polygon::new(e, lrs).map_err(|e| (None, e))
}

/// `true` if `first` (validate's first refused issue in constructor order)
/// names the defect `err` names.
fn names_the_same_defect(err: &(Option<usize>, GeometryError), first: &ValidityIssue) -> bool {
    match (err, *first) {
        (
            (Some(k), GeometryError::NonFiniteCoordinate { .. }),
            V::NonFiniteCoordinate { ring, .. },
        )
        | ((Some(k), GeometryError::InsufficientVertices { .. }), V::TooFewPoints { ring, .. })
        | (
            (Some(k), GeometryError::ConsecutiveDuplicateVertices { .. }),
            V::DuplicatePoint { ring, .. },
        ) => *k == ring,
        ((None, GeometryError::DegenerateRing), V::DegenerateRing { .. })
        | ((None, GeometryError::NonFiniteResult), V::NonFiniteArea { .. }) => true,
        ((None, GeometryError::IncorrectWindingOrder), V::RingOrientation { ring }) => ring == 0,
        (
            (None, GeometryError::IncorrectHoleWindingOrder { hole_index }),
            V::RingOrientation { ring },
        ) => ring == hole_index + 1,
        ((None, GeometryError::HoleNotInsideExterior { hole_index }), issue) => {
            let h = hole_index + 1;
            matches!(issue,
                V::HoleOutsideExterior { ring } if ring == h)
                || matches!(issue, V::RingsCross { ring, other: 0, .. } if ring == h)
                || matches!(issue, V::DuplicateRings { ring, other: 0 } if ring == h)
        }
        ((None, GeometryError::OverlappingHoles { first, second }), issue) => {
            let (a, b) = (first + 1, second + 1);
            matches!(issue, V::RingsCross { ring, other, .. } if (ring, other) == (b, a))
                || matches!(issue, V::DuplicateRings { ring, other } if (ring, other) == (b, a))
                || matches!(issue, V::NestedHoles { ring, inside } if (ring, inside) == (a, b) || (ring, inside) == (b, a))
        }
        _ => false,
    }
}

#[test]
fn validate_agrees_with_the_constructors_on_random_ring_sets() {
    const CASES: usize = 12_000;
    let mut rng = Rng(0x00C0_FFEE_5EED_0003);
    let mut tally = std::collections::BTreeMap::<&str, usize>::new();
    for case in 0..CASES {
        let rs = random_rings(&mut rng);
        let issues = validate(&rs);
        let mut refused: Vec<&ValidityIssue> =
            issues.iter().filter(|i| i.constructor_refuses()).collect();
        refused.sort_by_key(|i| phase(i)); // stable: keeps validate's order within a phase
        match construct(&rs) {
            Ok(p) => {
                assert!(
                    refused.is_empty(),
                    "case {case}: constructed, but {issues:?}\n{rs:?}"
                );
                // What construction does not check, `validity_issues` on
                // the constructed polygon reports exactly as `validate` does.
                let rest: Vec<ValidityIssue> = issues.clone();
                let got: Vec<ValidityIssue> =
                    validity_issues(&p).into_iter().map(|l| l.issue).collect();
                assert_eq!(got, rest, "case {case}\n{rs:?}");
                let key = if rest
                    .iter()
                    .any(|i| matches!(i, V::DisconnectedInterior { .. }))
                {
                    "constructed, disconnected"
                } else if rest.is_empty() {
                    "constructed, valid"
                } else {
                    "constructed, self-intersecting"
                };
                *tally.entry(key).or_default() += 1;
            }
            Err(err) => {
                assert!(
                    !refused.is_empty(),
                    "case {case}: refused {err:?}, but validate found nothing refused: {issues:?}\n{rs:?}"
                );
                assert!(
                    names_the_same_defect(&err, refused[0]),
                    "case {case}: constructor {err:?}, validate first {:?}\nall {issues:?}\n{rs:?}",
                    refused[0]
                );
                let key = match err.1 {
                    GeometryError::HoleNotInsideExterior { .. } => "refused, hole placement",
                    GeometryError::OverlappingHoles { .. } => "refused, hole pair",
                    _ => "refused, ring",
                };
                *tally.entry(key).or_default() += 1;
            }
        }
    }
    println!("validate agreement over {CASES} ring sets: {tally:?}");
    for key in [
        "constructed, valid",
        "constructed, disconnected",
        "constructed, self-intersecting",
        "refused, hole placement",
        "refused, hole pair",
        "refused, ring",
    ] {
        assert!(
            tally.get(key).copied().unwrap_or(0) >= 50,
            "too few cases of {key}: {tally:?}"
        );
    }
}

#[test]
fn validate_multipolygon_agrees_with_multipolygon_new_on_random_parts() {
    const CASES: usize = 4_000;
    let mut rng = Rng(0x00C0_FFEE_5EED_0004);
    let (mut overlapping, mut disjoint) = (0, 0);
    for case in 0..CASES {
        let n = 2 + rng.below(3) as usize;
        let parts: Vec<Vec<CoordSeq>> = (0..n)
            .map(|_| {
                let (cx, cy) = (rng.below(16) as i64, rng.below(16) as i64);
                let mut rs = vec![star(&mut rng, cx, cy, 2, 5, false)];
                if rng.chance(30) {
                    rs.push(star(&mut rng, cx, cy, 1, 1, true));
                }
                rs.iter().map(|r| seq(r)).collect()
            })
            .collect();
        let polygons: Vec<Polygon> = match parts.iter().map(|p| construct(p).ok()).collect() {
            Some(ps) => ps,
            None => continue,
        };
        let issues = validate_multipolygon(&parts);
        let pair = issues
            .iter()
            .find(|(_, i)| matches!(i, V::PartsOverlap { .. } | V::NestedParts { .. }));
        match MultiPolygon::new(polygons) {
            Ok(_) => {
                assert!(pair.is_none(), "case {case}: {issues:?}");
                disjoint += 1;
            }
            Err(GeometryError::OverlappingParts { first, second }) => {
                let Some((_, issue)) = pair else {
                    panic!("case {case}: refused ({first}, {second}), but {issues:?}");
                };
                let ok = match *issue {
                    V::PartsOverlap { part, other } => (other, part) == (first, second),
                    V::NestedParts { part, inside } => {
                        (part, inside) == (second, first) || (part, inside) == (first, second)
                    }
                    _ => false,
                };
                assert!(ok, "case {case}: ({first}, {second}) vs {issue:?}");
                overlapping += 1;
            }
            Err(e) => panic!("case {case}: {e:?}"),
        }
    }
    println!("validate_multipolygon agreement: {overlapping} overlapping, {disjoint} disjoint");
    assert!(
        overlapping >= 200 && disjoint >= 200,
        "{overlapping} / {disjoint}"
    );
}
