// src/prepared/tests.rs

//! The constructed tests behind `prepared`: the far-part rule, delegation,
//! the predicates against the full matrix, and the empty cases. The
//! generated identity sweeps are in the integration test `prepared`.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use super::*;
use crate::topology::de9im::IntersectionMatrix;
use crate::types::{Coord2d, CoordSeq, LineString, LinearRing, MultiPoint, Point};

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

fn square(x: f64, y: f64, s: f64) -> Polygon {
    Polygon::new_simple(ring(&[(x, y), (x + s, y), (x + s, y + s), (x, y + s)])).unwrap()
}

fn geometry(mp: &MultiPolygon) -> Geometry {
    Geometry::MultiPolygon(mp.clone())
}

fn unprepared(mask: &MultiPolygon, g: &Geometry) -> IntersectionMatrix {
    IntersectionMatrix::compute_geometry(&geometry(mask), g).unwrap()
}

/// Every method of a prepared mask against the unprepared computation.
fn check(prepared: &PreparedPolygonal, mask: &MultiPolygon, g: &Geometry) {
    let m = unprepared(mask, g);
    assert_eq!(prepared.relate(g).unwrap(), m, "relate: {g:?}");
    assert_eq!(
        prepared.intersects(g).unwrap(),
        m.intersects(),
        "intersects"
    );
    assert_eq!(prepared.disjoint(g).unwrap(), m.disjoint(), "disjoint");
    assert_eq!(prepared.covers(g).unwrap(), m.covers(), "covers");
    assert_eq!(prepared.contains(g).unwrap(), m.contains(), "contains");
    assert_eq!(prepared.within(g).unwrap(), m.within(), "within");
}

/// A set of features around a mask: inside, outside, touching a vertex,
/// along an edge, inside a hole, crossing, covering, far away.
fn features() -> Vec<Geometry> {
    let line = |pts: &[(f64, f64)]| {
        Geometry::LineString(
            LineString::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap(),
        )
    };
    let point = |x, y| Geometry::Point(Point::new(c(x, y)).unwrap());
    vec![
        Geometry::Polygon(square(1.0, 1.0, 1.0)),
        Geometry::Polygon(square(4.0, 4.0, 1.0)),
        Geometry::Polygon(square(3.0, 0.0, 1.0)),
        Geometry::Polygon(square(-1.0, -1.0, 20.0)),
        Geometry::Polygon(square(2.5, 2.5, 0.2)),
        Geometry::Polygon(square(100.0, 100.0, 1.0)),
        Geometry::Polygon(square(3.0, 3.0, 2.0)),
        line(&[(0.5, 0.5), (2.5, 2.5)]),
        line(&[(3.0, 0.0), (3.0, 3.0)]),
        line(&[(-5.0, 1.0), (50.0, 1.0)]),
        line(&[(100.0, 0.0), (101.0, 1.0)]),
        line(&[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]),
        point(1.0, 1.0),
        point(3.0, 0.0),
        point(2.5, 2.5),
        point(50.0, 50.0),
        Geometry::MultiPoint(MultiPoint::new(vec![
            Point::new(c(0.0, 0.0)).unwrap(),
            Point::new(c(9.0, 9.0)).unwrap(),
        ])),
    ]
}

/// A mask of three parts: a holed square, a square touching it at a corner,
/// and a far square.
fn mask3() -> MultiPolygon {
    let holed = Polygon::new(
        ring(&[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)]),
        vec![ring(&[(2.0, 2.0), (2.0, 2.8), (2.8, 2.8), (2.8, 2.0)])],
    )
    .unwrap();
    MultiPolygon::new(vec![holed, square(3.0, 3.0, 1.0), square(8.0, 8.0, 1.0)]).unwrap()
}

#[test]
fn every_method_matches_the_unprepared_computation() {
    let mask = mask3();
    let prepared = PreparedPolygonal::new(&mask).unwrap();
    assert!(!prepared.is_delegating());
    for g in features() {
        check(&prepared, &mask, &g);
    }
}

#[test]
fn a_polygon_mask_prepares_like_its_one_part_multipolygon() {
    let p = square(0.0, 0.0, 3.0);
    let prepared = PreparedPolygonal::new(&p).unwrap();
    for g in features() {
        let m = IntersectionMatrix::compute_geometry(&Geometry::Polygon(p.clone()), &g).unwrap();
        assert_eq!(prepared.relate(&g).unwrap(), m, "{g:?}");
    }
}

/// The non-candidate exterior-cell rule, tested directly: a mask with a far
/// part, against features that meet only the near parts. The far part
/// raises `IE` and `BE` and changes no other cell.
#[test]
fn a_far_part_contributes_exactly_ie_area_and_be_line() {
    let near = square(0.0, 0.0, 3.0);
    let far = square(500.0, 500.0, 1.0);
    let with_far = MultiPolygon::new(vec![near.clone(), far]).unwrap();
    let without = MultiPolygon::from(near);
    let prepared = PreparedPolygonal::new(&with_far).unwrap();
    for g in features() {
        let full = unprepared(&with_far, &g);
        assert_eq!(prepared.relate(&g).unwrap(), full, "{g:?}");
        // Against the near part alone, every cell but IE and BE is the same.
        let mut near_only = unprepared(&without, &g);
        near_only.interior_exterior = Dimension::Area;
        near_only.boundary_exterior = Dimension::Line;
        assert_eq!(near_only, full, "the rule, on {g:?}");
    }
}

/// No part is a candidate: the feature is in the gap between parts.
#[test]
fn a_feature_between_parts_is_disjoint_by_one_far_part() {
    let mask = MultiPolygon::new(vec![square(0.0, 0.0, 1.0), square(10.0, 10.0, 1.0)]).unwrap();
    let prepared = PreparedPolygonal::new(&mask).unwrap();
    for g in [
        Geometry::Polygon(square(5.0, 5.0, 1.0)),
        Geometry::Point(Point::new(c(5.0, 0.5)).unwrap()),
        Geometry::LineString(
            LineString::new(CoordSeq::new(vec![c(5.0, 0.0), c(5.0, 10.5)])).unwrap(),
        ),
    ] {
        check(&prepared, &mask, &g);
    }
}

#[test]
fn the_envelope_early_out_gives_the_disjoint_matrix_for_every_type() {
    let mask = mask3();
    let prepared = PreparedPolygonal::new(&mask).unwrap();
    let far = [
        Geometry::Polygon(square(100.0, 100.0, 1.0)),
        Geometry::Point(Point::new(c(-100.0, 0.0)).unwrap()),
        Geometry::LineString(
            LineString::new(CoordSeq::new(vec![c(0.0, 100.0), c(1.0, 100.0)])).unwrap(),
        ),
        // A closed line: an empty mod-2 boundary.
        Geometry::LineString(
            LineString::new(CoordSeq::new(vec![
                c(0.0, 100.0),
                c(1.0, 100.0),
                c(0.0, 101.0),
                c(0.0, 100.0),
            ]))
            .unwrap(),
        ),
        Geometry::MultiPoint(MultiPoint::new(vec![Point::new(c(200.0, 200.0)).unwrap()])),
    ];
    for g in far {
        check(&prepared, &mask, &g);
    }
}

#[test]
fn empty_operands_delegate_and_agree() {
    let empty = MultiPolygon::empty();
    let prepared = PreparedPolygonal::new(&empty).unwrap();
    assert!(!prepared.is_delegating());
    assert_eq!(prepared.envelope(), None);
    for g in features() {
        check(&prepared, &empty, &g);
    }
    let mask = mask3();
    let prepared = PreparedPolygonal::new(&mask).unwrap();
    let empty_g = Geometry::MultiPoint(MultiPoint::new(Vec::new()));
    check(&prepared, &mask, &empty_g);
}

// ===========================================================================
// Delegation
// ===========================================================================

mod delegation {
    use super::*;

    /// A ring that crosses itself (a bowtie), as `Polygon::new` accepts it
    /// (ADR-27: the constructors do not check simplicity).
    pub(super) fn self_crossing() -> MultiPolygon {
        // Asymmetric, so its signed area is not zero: (0,0)-(5,4) crosses
        // (0,3)-(4,0) at about (1.94, 1.55).
        let bowtie =
            Polygon::new_simple(ring(&[(0.0, 0.0), (5.0, 4.0), (0.0, 3.0), (4.0, 0.0)])).unwrap();
        MultiPolygon::new(vec![bowtie]).unwrap()
    }

    /// A ring pinched at (1, 1): the inverted hole.
    pub(super) fn pinched() -> MultiPolygon {
        let p = Polygon::new_simple(ring(&[
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (1.0, 1.0),
            (2.0, 2.5),
            (1.0, 2.8),
            (1.0, 1.0),
            (0.0, 3.0),
        ]))
        .unwrap();
        MultiPolygon::new(vec![square(10.0, 10.0, 1.0), p]).unwrap()
    }

    #[test]
    fn a_self_crossing_mask_delegates_and_names_the_ring() {
        let mask = self_crossing();
        let prepared = PreparedPolygonal::new(&mask).unwrap();
        assert!(prepared.is_delegating());
        let r = prepared.delegation_reason().unwrap();
        assert_eq!((r.part, r.ring), (0, 0));
        assert_eq!(
            r.edges,
            first_self_contact(mask.polygons()[0].exterior()).unwrap()
        );
        for g in features() {
            check(&prepared, &mask, &g);
        }
    }

    /// A small mask takes its verdict on first need: points asked first
    /// get the unprepared answers without taking it, and the verdict taken
    /// afterwards is the one an eager check gives, on a simple mask and on
    /// a delegating one.
    #[test]
    fn a_small_mask_answers_points_before_its_verdict() {
        for (mask, delegates) in [(self_crossing(), true), (mask3(), false)] {
            let prepared = PreparedPolygonal::new(&mask).unwrap();
            let pending = |p: &PreparedPolygonal| match &p.kind {
                Kind::Small(v) => v.get().is_none(),
                _ => panic!("a mask of at most {SMALL_MASK_EDGES} edges is small"),
            };
            assert!(pending(&prepared));
            for g in features() {
                if let Geometry::Point(pt) = &g {
                    let m = unprepared(&mask, &g);
                    assert_eq!(
                        prepared.locate(pt.coord()) != Location::Exterior,
                        m.intersects()
                    );
                    assert_eq!(prepared.intersects(&g).unwrap(), m.intersects());
                    assert_eq!(prepared.covers(&g).unwrap(), m.covers());
                    assert_eq!(prepared.contains(&g).unwrap(), m.contains());
                }
            }
            assert!(pending(&prepared), "a point question took the verdict");
            let copy = prepared.clone();
            assert_eq!(prepared.is_delegating(), delegates);
            assert_eq!(
                prepared.delegation_reason(),
                super::super::first_non_simple_ring(mask.polygons())
            );
            assert!(!pending(&prepared));
            for g in features() {
                check(&copy, &mask, &g);
            }
            assert_eq!(copy.is_delegating(), delegates);
        }
    }

    #[test]
    fn a_pinched_mask_delegates_and_names_the_ring() {
        let mask = pinched();
        let prepared = PreparedPolygonal::new(&mask).unwrap();
        assert!(prepared.is_delegating());
        let r = prepared.delegation_reason().unwrap();
        assert_eq!((r.part, r.ring), (1, 0), "the second part's exterior");
        for g in features() {
            check(&prepared, &mask, &g);
        }
    }

    #[test]
    fn a_non_simple_hole_delegates_too() {
        let p = Polygon::new(
            ring(&[(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 5.0)]),
            vec![
                ring(&[
                    (1.0, 1.0),
                    (1.0, 3.0),
                    (3.0, 3.0),
                    (3.0, 1.0),
                    (2.0, 2.0),
                    (1.0, 3.0),
                    (1.0, 1.0),
                ])
                .clone(),
            ],
        );
        // Such a hole may be refused by `Polygon::new` for other reasons;
        // only an accepted one can be prepared.
        if let Ok(p) = p {
            let prepared = PreparedPolygonal::new(&p).unwrap();
            assert!(prepared.is_delegating());
            assert_eq!(prepared.delegation_reason().unwrap().ring, 1);
        }
    }

    /// The delegation counter shows the unprepared path ran. The counters
    /// are thread-local, so a delegation test running in parallel on another
    /// thread does not add to the count.
    #[cfg(feature = "unstable-counters")]
    #[test]
    fn delegation_runs_the_unprepared_path() {
        let mask = self_crossing();
        let prepared = PreparedPolygonal::new(&mask).unwrap();
        let g = Geometry::Polygon(square(0.5, 0.5, 0.2));
        let before = crate::counters::snapshot().prepared_delegations;
        prepared.relate(&g).unwrap();
        prepared.covers(&g).unwrap();
        let after = crate::counters::snapshot().prepared_delegations;
        assert_eq!(after - before, 2);

        let simple = PreparedPolygonal::new(&mask3()).unwrap();
        let before = crate::counters::snapshot().prepared_delegations;
        simple.relate(&g).unwrap();
        assert_eq!(crate::counters::snapshot().prepared_delegations, before);
    }
}

// ===========================================================================
// The exact-point oracle (deep-W)
// ===========================================================================

mod oracle {
    use super::*;
    use crate::topology::label::Location;
    use arithmetica::RationalCoord2d;

    /// Crossings of two grid segments that land on representable points, so
    /// the rational branch can be checked against the stored one.
    fn representable_crossings() -> Vec<(RationalCoord2d, Coord2d)> {
        let mut out = Vec::new();
        for (p, q, r, s, at) in [
            ((0.0, 0.0), (4.0, 4.0), (0.0, 4.0), (4.0, 0.0), (2.0, 2.0)),
            ((1.0, 0.0), (1.0, 6.0), (0.0, 3.0), (6.0, 3.0), (1.0, 3.0)),
            ((0.0, 0.0), (6.0, 3.0), (0.0, 3.0), (6.0, 0.0), (3.0, 1.5)),
            ((0.5, 0.5), (2.5, 2.5), (0.5, 2.5), (2.5, 0.5), (1.5, 1.5)),
            ((0.0, 1.0), (3.0, 1.0), (2.0, 0.0), (2.0, 5.0), (2.0, 1.0)),
        ] {
            let x = RationalCoord2d::from_segment_intersection(
                [p.0, p.1],
                [q.0, q.1],
                [r.0, r.1],
                [s.0, s.1],
            )
            .unwrap()
            .unwrap();
            out.push((x, c(at.0, at.1)));
        }
        out
    }

    /// Directions to probe from a point, stored.
    fn directions(at: Coord2d) -> Vec<Coord2d> {
        [
            (1.0, 0.0),
            (-1.0, 0.0),
            (0.0, 1.0),
            (0.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
            (1.0, -1.0),
            (-1.0, -1.0),
            (2.0, 1.0),
            (-2.0, 1.0),
        ]
        .iter()
        .map(|&(dx, dy)| c(at.x + dx, at.y + dy))
        .collect()
    }

    #[test]
    fn the_rational_branch_agrees_with_the_stored_branch() {
        for (mask, name) in [(mask3(), "mask3"), (self_touching_mask(), "touching")] {
            let prepared = PreparedPolygonal::new(&mask).unwrap();
            for (x, at) in representable_crossings() {
                assert_eq!(
                    prepared.locate(at),
                    exact_point::locate(&ExactPoint::Rational(&x), &Areal::multi(&mask)),
                    "{name}: locate at {at:?}"
                );
                for toward in directions(at) {
                    let stored = prepared.location_left_of(ExactPoint::Stored(at), toward);
                    let rational = prepared.location_left_of(ExactPoint::Rational(&x), toward);
                    assert_eq!(stored, rational, "{name}: left of {at:?} toward {toward:?}");
                    assert_ne!(stored, Location::Boundary);
                }
            }
        }
    }

    /// A mask whose rings pass through the probe points, so the wedge
    /// branch is exercised: a square with a hole touching the exterior at
    /// a vertex is not valid, so two squares touching at (2, 2) serve.
    fn self_touching_mask() -> MultiPolygon {
        MultiPolygon::new(vec![square(0.0, 0.0, 2.0), square(2.0, 2.0, 2.0)]).unwrap()
    }

    /// The oracle against a brute-force definition: the location of a
    /// stored point a small step along the direction and a small step to its
    /// left, where that point is off every ring. On a half-unit grid a step
    /// of 1/8 along and 1/8 left of an axis-aligned or diagonal direction
    /// never lands on a ring.
    #[test]
    fn left_of_matches_a_point_just_left_of_the_ray() {
        let mask = mask3();
        let prepared = PreparedPolygonal::new(&mask).unwrap();
        let plain = Areal::multi(&mask);
        let mut checked = 0;
        for xi in 0..=20 {
            for yi in 0..=20 {
                let at = c(f64::from(xi) * 0.5 - 1.0, f64::from(yi) * 0.5 - 1.0);
                for toward in directions(at) {
                    let (dx, dy) = (toward.x - at.x, toward.y - at.y);
                    let len = dx.hypot(dy);
                    let (ux, uy) = (dx / len, dy / len);
                    // Along by 1/8, then left (rotate by +90°) by 1/64: the
                    // sideways step is much smaller than the step along, so
                    // the probe is on the ray's side of any edge the ray
                    // leaves at an angle.
                    let probe = c(at.x + ux / 8.0 - uy / 64.0, at.y + uy / 8.0 + ux / 64.0);
                    if plain.locate(probe) == Location::Boundary {
                        continue;
                    }
                    assert_eq!(
                        prepared.location_left_of(ExactPoint::Stored(at), toward),
                        plain.locate(probe),
                        "at {at:?} toward {toward:?}"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 3_000, "{checked}");
    }

    #[test]
    fn locate_and_segment_events_match_the_unprepared_areal() {
        let mask = mask3();
        let prepared = PreparedPolygonal::new(&mask).unwrap();
        let plain = Areal::multi(&mask);
        for xi in 0..=12 {
            for yi in 0..=12 {
                let p = c(f64::from(xi) * 0.75 - 1.0, f64::from(yi) * 0.75 - 1.0);
                assert_eq!(prepared.locate(p), plain.locate(p), "{p:?}");
                let q = c(9.0 - p.x, p.y + 0.5);
                assert_eq!(
                    prepared.segment_events(p, q),
                    segment_events(p, q, &plain),
                    "{p:?}→{q:?}"
                );
            }
        }
    }

    #[test]
    fn edges_near_reports_exactly_the_edges_whose_box_meets_the_window() {
        let mask = mask3();
        let prepared = PreparedPolygonal::new(&mask).unwrap();
        let w = Envelope::new(c(2.5, 2.5), c(3.5, 3.5)).unwrap();
        let mut got = Vec::new();
        prepared.edges_near(&w, |r| got.push(r));
        let mut want = Vec::new();
        for (part, p) in mask.polygons().iter().enumerate() {
            let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
            for (ring, r) in rings.enumerate() {
                for (edge, (a, b)) in r.edges().enumerate() {
                    let meets = a.x.min(b.x) <= w.max().x
                        && w.min().x <= a.x.max(b.x)
                        && a.y.min(b.y) <= w.max().y
                        && w.min().y <= a.y.max(b.y);
                    if meets {
                        want.push(EdgeRef { part, ring, edge });
                    }
                }
            }
        }
        assert_eq!(got, want);
        for r in got {
            let (a, b) = prepared.edge(r);
            assert!(
                mask.polygons()[r.part]
                    .exterior()
                    .edges()
                    .any(|e| e == (a, b))
                    || r.ring > 0
            );
        }
    }
}

// ===========================================================================
// locate_on_line
// ===========================================================================

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::range_plus_one
)]
mod on_line {
    use super::*;
    use crate::topology::label::Location;

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
    }

    /// A star-shaped ring of `n` integer vertices around `(cx, cy)`, radii
    /// in `[rmin, rmax]`, counterclockwise; `None` if construction refuses.
    fn star(r: &mut Rng, cx: f64, cy: f64, n: usize, rmin: u64, rmax: u64) -> Option<LinearRing> {
        let pts: Vec<Coord2d> = (0..n)
            .map(|k| {
                let t = k as f64 / n as f64 * core::f64::consts::TAU;
                let rad = (rmin + r.below(rmax - rmin + 1)) as f64;
                c((cx + rad * t.cos()).round(), (cy + rad * t.sin()).round())
            })
            .collect();
        LinearRing::new(CoordSeq::new(pts)).ok()
    }

    /// A star with a square hole at its center, or the star alone.
    fn holed_star(r: &mut Rng, cx: f64, cy: f64, n: usize) -> Option<Polygon> {
        let ext = star(r, cx, cy, n, 8, 14)?;
        let hole = ring(&[
            (cx - 2.0, cy - 2.0),
            (cx - 2.0, cy + 2.0),
            (cx + 2.0, cy + 2.0),
            (cx + 2.0, cy - 2.0),
        ]);
        Polygon::new(ext.clone(), vec![hole])
            .or_else(|_| Polygon::new_simple(ext))
            .ok()
    }

    /// Masks of every kind: indexed (over 64 edges), small, delegating,
    /// multipart with holes.
    fn masks(r: &mut Rng) -> Vec<MultiPolygon> {
        let mut out = vec![mask3(), delegation::self_crossing(), delegation::pinched()];
        // Over 64 edges with a self-crossing part: a delegating mask.
        let big = loop {
            if let Some(p) = holed_star(r, 0.0, 0.0, 90) {
                break p;
            }
        };
        let bowtie =
            Polygon::new_simple(ring(&[(40.0, 0.0), (45.0, 4.0), (40.0, 3.0), (44.0, 0.0)]))
                .unwrap();
        out.push(MultiPolygon::new(vec![big, bowtie]).unwrap());
        for _ in 0..12 {
            let n = [6, 12, 40, 90, 160][r.below(5) as usize];
            let parts: Vec<Polygon> = (0..1 + r.below(3))
                .filter_map(|k| {
                    holed_star(
                        r,
                        40.0 * k as f64,
                        f64::from(u8::try_from(k).unwrap()) * 3.0,
                        n,
                    )
                })
                .collect();
            if let Ok(mp) = MultiPolygon::new(parts) {
                out.push(mp);
            }
        }
        out
    }

    #[test]
    fn every_point_on_a_line_is_located_as_locate_locates_it() {
        let mut r = Rng(0x0011_7E5E);
        let mut seen = [0usize; 3];
        let mut kinds = [0usize; 3];
        for mask in masks(&mut r) {
            let prepared = PreparedPolygonal::new(&mask).unwrap();
            kinds[match &prepared.kind {
                Kind::Simple { .. } => 0,
                Kind::Small(_) => 1,
                Kind::Delegating(_) => 2,
            }] += 1;
            let (lo, hi) = match prepared.envelope() {
                Some(e) => (e.min(), e.max()),
                None => continue,
            };
            let vertex_xs: Vec<f64> = mask
                .polygons()
                .iter()
                .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes()))
                .flat_map(|rg| rg.coords().as_slice().iter().map(|v| v.x))
                .collect();
            let mut y = lo.y - 1.0;
            while y <= hi.y + 1.0 {
                let mut xs: Vec<f64> = Vec::new();
                let mut x = lo.x - 1.0;
                while x <= hi.x + 1.0 {
                    xs.push(if x == 0.0 && r.below(2) == 0 { -0.0 } else { x });
                    x += 0.5;
                }
                xs.extend(vertex_xs.iter().filter(|_| r.below(3) == 0));
                xs.sort_by(f64::total_cmp);
                if r.below(3) == 0 {
                    // Unsorted, with repeats.
                    let extra: Vec<f64> = xs.iter().copied().filter(|_| r.below(4) == 0).collect();
                    xs.extend(extra);
                    for i in (1..xs.len()).rev() {
                        xs.swap(i, r.below(i as u64 + 1) as usize);
                    }
                }
                let got = prepared.locate_on_line(y, &xs).unwrap();
                assert_eq!(got.len(), xs.len());
                for (i, &x) in xs.iter().enumerate() {
                    let want = prepared.locate(c(x, y));
                    assert_eq!(got[i], want, "mask {mask:?}\n row {y}, x {x}");
                    seen[match want {
                        Location::Boundary => 0,
                        Location::Interior => 1,
                        Location::Exterior => 2,
                    }] += 1;
                }
                y += 0.5;
            }
        }
        assert!(
            kinds.iter().all(|&k| k > 0),
            "every kind of mask occurs: {kinds:?}"
        );
        assert!(
            seen.iter().all(|&k| k > 500),
            "every location occurs: {seen:?}"
        );
    }

    #[test]
    fn non_finite_input_is_refused_and_empty_cases_answer() {
        let prepared = PreparedPolygonal::new(&mask3()).unwrap();
        assert_eq!(
            prepared
                .locate_on_line(f64::NAN, &[f64::INFINITY])
                .unwrap_err(),
            GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 1
            }
        );
        assert_eq!(
            prepared
                .locate_on_line(1.0, &[0.0, 1.0, f64::NEG_INFINITY])
                .unwrap_err(),
            GeometryError::NonFiniteCoordinate {
                coord_index: 2,
                axis: 0
            }
        );
        assert_eq!(prepared.locate_on_line(1.0, &[]).unwrap(), vec![]);
        assert_eq!(
            prepared.locate_on_line(50.0, &[0.0, 1.0]).unwrap(),
            vec![Location::Exterior; 2],
            "a row outside the mask's box"
        );
        let empty = PreparedPolygonal::new(&MultiPolygon::empty()).unwrap();
        assert_eq!(
            empty.locate_on_line(0.0, &[0.0]).unwrap(),
            vec![Location::Exterior]
        );
    }
}
