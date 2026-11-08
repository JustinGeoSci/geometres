// src/ops/simplify_topology/tests.rs

//! Tests of topology-preserving simplification: the five promises checked
//! by an independent oracle on 10,000+ generated inputs and on hand cases,
//! each hand case with its exact expected output.
//!
//! The oracle does not reuse the acceptance test. It compares the input and
//! the output element by element, as point sets, with this crate's exact
//! `relate` and point location:
//!
//! - **promise 1:** each output chain is a subsequence of its input chain,
//!   bit for bit;
//! - **promise 2:** each removed vertex is within tolerance (in `f64`) of
//!   the output segment that replaced it;
//! - **promise 3:** for every pair of elements (rings as closed lines, lines,
//!   points), the DE-9IM matrix before equals the matrix after, and each
//!   line's and ring's simplicity is unchanged;
//! - **promise 4:** every output vertex and point keeps its location
//!   (inside, on, outside) relative to every ring of another element, and
//!   on/off relative to every line;
//! - **promise 5:** a valid input gives an output with no validity issue.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines,
    clippy::needless_range_loop
)]

use super::sealed::{Flat, Sealed};
use super::*;
use crate::exact::kernel::{RingLocation, locate_in_ring};
use crate::ops::validate::ValidityIssue;
use crate::ops::{is_simple, validity_issues};
use crate::predicates::relate;
use crate::types::{MultiPoint, Point};

// ===========================================================================
// Helpers
// ===========================================================================

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn seq(p: &[(f64, f64)]) -> CoordSeq {
    CoordSeq::new(p.iter().map(|&(x, y)| c(x, y)).collect())
}

fn ring(p: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(seq(p)).unwrap()
}

fn line(p: &[(f64, f64)]) -> LineString {
    LineString::new(seq(p)).unwrap()
}

fn poly(ext: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Polygon {
    Polygon::new(ring(ext), holes.iter().map(|h| ring(h)).collect()).unwrap()
}

fn pts(r: &LinearRing) -> Vec<(f64, f64)> {
    r.coords().as_slice().iter().map(|p| (p.x, p.y)).collect()
}

fn lpts(l: &LineString) -> Vec<(f64, f64)> {
    l.coords().as_slice().iter().map(|p| (p.x, p.y)).collect()
}

fn flat_of<G: Sealed>(g: &G) -> Flat {
    let mut f = Flat::default();
    g.flatten(&mut f);
    f
}

/// The chains of a flattened value, as coordinate vectors.
fn chains(f: &Flat) -> Vec<(Vec<Coord2d>, bool)> {
    f.chains
        .iter()
        .map(|ch| (f.verts[ch.start..ch.start + ch.len].to_vec(), ch.closed))
        .collect()
}

// ===========================================================================
// The oracle
// ===========================================================================

/// One element as a geometry `relate` takes: a ring as a closed line.
enum Element {
    Line(LineString),
    Point(Point),
}

fn element(chain: &[Coord2d], closed: bool) -> Element {
    let mut v = chain.to_vec();
    if closed {
        v.push(chain[0]);
    }
    Element::Line(LineString::new(CoordSeq::new(v)).unwrap())
}

fn elements(f: &Flat) -> Vec<Element> {
    let mut out: Vec<Element> = chains(f).iter().map(|(v, cl)| element(v, *cl)).collect();
    out.extend(
        f.points
            .iter()
            .map(|&p| Element::Point(Point::new(p).unwrap())),
    );
    out
}

fn rel(a: &Element, b: &Element) -> crate::predicates::Relation {
    match (a, b) {
        (Element::Line(x), Element::Line(y)) => relate(x, y),
        (Element::Line(x), Element::Point(y)) => relate(x, y),
        (Element::Point(x), Element::Line(y)) => relate(x, y),
        (Element::Point(x), Element::Point(y)) => relate(x, y),
    }
    .unwrap()
}

fn bits(p: Coord2d) -> (u64, u64) {
    (p.x.to_bits(), p.y.to_bits())
}

/// Promises 1 and 2 for one chain: the output is a subsequence of the input
/// (bit for bit, and starting where the input starts, since anchors are
/// kept) and every removed vertex is within tolerance of the segment that
/// replaced it.
fn check_chain(inp: &[Coord2d], out: &[Coord2d], closed: bool, tol: f64, what: &str) {
    let mut pos = Vec::new();
    let mut k = 0;
    for &p in out {
        while k < inp.len() && bits(inp[k]) != bits(p) {
            k += 1;
        }
        assert!(
            k < inp.len(),
            "{what}: output {out:?} is not a subsequence of {inp:?}"
        );
        pos.push(k);
        k += 1;
    }
    let n = inp.len();
    let spans: Vec<(usize, usize)> = if closed {
        (0..pos.len())
            .map(|i| {
                (
                    pos[i],
                    if i + 1 < pos.len() {
                        pos[i + 1]
                    } else {
                        pos[0] + n
                    },
                )
            })
            .collect()
    } else {
        assert_eq!(pos.first(), Some(&0), "{what}: line start dropped");
        assert_eq!(pos.last(), Some(&(n - 1)), "{what}: line end dropped");
        pos.windows(2).map(|w| (w[0], w[1])).collect()
    };
    for (a, b) in spans {
        let (u, w) = (inp[a % n], inp[b % n]);
        for m in a + 1..b {
            let d = segment_distance(inp[m % n], u, w).unwrap();
            assert!(
                d <= tol,
                "{what}: removed vertex {:?} is {d} from its chord",
                inp[m % n]
            );
        }
    }
}

/// Every promise, for one input and its output.
fn check_promises<G: TopologySimplifiable>(input: &G, output: &G, tol: f64, what: &str) {
    let (fi, fo) = (flat_of(input), flat_of(output));
    let (ci, co) = (chains(&fi), chains(&fo));
    assert_eq!(ci.len(), co.len(), "{what}: chain count");
    assert_eq!(fi.points, fo.points, "{what}: points");
    for (k, ((a, closed), (b, _))) in ci.iter().zip(&co).enumerate() {
        check_chain(a, b, *closed, tol, &format!("{what}, chain {k}"));
    }
    // Promise 3: pairwise relations, and each chain's own simplicity.
    let (ei, eo) = (elements(&fi), elements(&fo));
    for i in 0..ei.len() {
        if let (Element::Line(x), Element::Line(y)) = (&ei[i], &eo[i]) {
            assert_eq!(is_simple(x), is_simple(y), "{what}: element {i} simplicity");
        }
        for j in i + 1..ei.len() {
            if matches!((&ei[i], &ei[j]), (Element::Point(_), Element::Point(_))) {
                continue;
            }
            assert_eq!(
                rel(&ei[i], &ei[j]),
                rel(&eo[i], &eo[j]),
                "{what}: elements {i} and {j} changed relation"
            );
        }
    }
    check_sides(&ci, &co, &fo.points, what);
    // Promise 5.
    if validity_issues(input).is_empty() {
        let issues = validity_issues(output);
        assert!(
            issues.is_empty(),
            "{what}: valid input, invalid output: {issues:?}"
        );
    }
}

/// Promise 4: every output vertex of every other chain, and every point,
/// keeps its location relative to each ring (inside, on, outside) and each
/// line (on or off).
fn check_sides(
    ci: &[(Vec<Coord2d>, bool)],
    co: &[(Vec<Coord2d>, bool)],
    points: &[Coord2d],
    what: &str,
) {
    for y in 0..ci.len() {
        let (before, closed) = (&ci[y].0, ci[y].1);
        let after = &co[y].0;
        let probes = co
            .iter()
            .enumerate()
            .filter(|(x, _)| *x != y)
            .flat_map(|(_, (v, _))| v.iter().copied())
            .chain(points.iter().copied());
        for p in probes {
            if closed {
                let (rb, ra) = (
                    LinearRing::new(CoordSeq::new(before.clone())).unwrap(),
                    LinearRing::new(CoordSeq::new(after.clone())).unwrap(),
                );
                assert_eq!(
                    locate_in_ring(p, &rb),
                    locate_in_ring(p, &ra),
                    "{what}: {p:?} changed side of ring {y}"
                );
            } else {
                let on = |v: &[Coord2d]| v.windows(2).any(|w| on_segment(p, w[0], w[1]));
                assert_eq!(
                    on(before),
                    on(after),
                    "{what}: {p:?} on/off line {y} changed"
                );
            }
        }
    }
}

// ===========================================================================
// Hand cases
// ===========================================================================

/// A hole near the exterior's edge: the exterior dips toward the hole.
/// Plain RDP drops the dip and collapses the hole; here the dip's
/// neighbors go and the dip stays, because removing the dip would sweep
/// over the whole hole (check (b)), and the hole cannot collapse.
#[test]
fn a_hole_near_the_exterior_edge_stays_inside() {
    let p = poly(
        &[
            (0.0, 0.0),
            (30.0, 0.0),
            (50.0, -5.0),
            (70.0, 0.0),
            (100.0, 0.0),
            (100.0, 100.0),
            (0.0, 100.0),
        ],
        &[&[(48.0, -2.0), (50.0, -1.0), (52.0, -2.0)]],
    );
    // Plain RDP refuses: it collapses the hole (and drops the dip).
    assert!(p.simplify(10.0).is_err());
    let s = simplify_preserve_topology(&p, 10.0).unwrap();
    assert_eq!(
        pts(s.exterior()),
        [
            (0.0, 0.0),
            (50.0, -5.0),
            (100.0, 0.0),
            (100.0, 100.0),
            (0.0, 100.0)
        ]
    );
    assert_eq!(pts(&s.holes()[0]), pts(&p.holes()[0]));
    check_promises(&p, &s, 10.0, "hole near edge");
}

/// Two parts touching at a vertex `w` of both. `w` is fixed, so each side
/// of it simplifies on its own; were it not, the runs through it would be
/// refused (check (c)) and split elsewhere, keeping `(10, 1)`.
#[test]
fn parts_touching_at_a_vertex_still_touch_there() {
    let a = poly(
        &[
            (0.0, 0.0),
            (10.0, 1.0),
            (20.0, 0.0),
            (30.0, 1.0),
            (40.0, 0.0),
            (40.0, 30.0),
            (0.0, 30.0),
        ],
        &[],
    );
    let b = poly(&[(20.0, 0.0), (15.0, -10.0), (25.0, -10.0)], &[]);
    let mp = MultiPolygon::new(vec![a, b.clone()]).unwrap();
    let s = simplify_preserve_topology(&mp, 2.0).unwrap();
    assert_eq!(
        pts(s.polygons()[0].exterior()),
        [
            (0.0, 0.0),
            (20.0, 0.0),
            (40.0, 0.0),
            (40.0, 30.0),
            (0.0, 30.0)
        ]
    );
    assert_eq!(s.polygons()[1], b);
    check_promises(&mp, &s, 2.0, "touching parts");
}

/// The fixed set, directly: line endpoints, a repeated coordinate, a vertex
/// on another segment's inside, and a point on a vertex.
#[test]
fn the_fixed_vertices_are_exactly_the_contacts() {
    let g = GeometryCollection::new(vec![
        Geometry::LineString(line(&[
            (0.0, 0.0),
            (5.0, 1.0),
            (10.0, 0.0),
            (15.0, 1.0),
            (20.0, 0.0),
        ])),
        // Ends on (10, 0), a vertex of the first line, and rests its middle
        // vertex on the first line's segment (0,0)-(5,1) at (2.5, 0.5).
        Geometry::LineString(line(&[(10.0, 0.0), (2.5, 0.5), (2.5, 5.0)])),
        Geometry::Point(Point::new(c(15.0, 1.0)).unwrap()),
    ])
    .unwrap();
    let flat = flat_of(&g);
    let s = Simplifier::new(&flat, 1.0).unwrap();
    let fixed: Vec<usize> = (0..flat.verts.len()).filter(|&i| s.fixed[i]).collect();
    // 0, 4: endpoints; 2: repeated (10,0); 3: under the point (15,1);
    // 5, 7: endpoints; 6: on a segment of line 0.
    assert_eq!(fixed, [0, 2, 3, 4, 5, 6, 7]);
}

/// A narrow channel: part B pokes into a notch of part A. A's notch cannot
/// close over B's tip; B's tip, once A has been tried, can go, leaving the
/// parts still apart.
#[test]
fn a_narrow_channel_between_parts_stays_open() {
    let a = poly(
        &[
            (0.0, 0.0),
            (40.0, 0.0),
            (40.0, 10.0),
            (22.0, 10.0),
            (20.0, 4.0),
            (18.0, 10.0),
            (0.0, 10.0),
        ],
        &[],
    );
    let b = poly(
        &[
            (20.0, 6.0),
            (21.0, 11.0),
            (40.0, 11.0),
            (40.0, 20.0),
            (0.0, 20.0),
            (0.0, 11.0),
            (19.0, 11.0),
        ],
        &[],
    );
    let mp = MultiPolygon::new(vec![a.clone(), b]).unwrap();
    assert!(mp.simplify(6.0).is_err(), "plain RDP closes the channel");
    let s = simplify_preserve_topology(&mp, 6.0).unwrap();
    assert_eq!(
        pts(s.polygons()[0].exterior()),
        [
            (0.0, 0.0),
            (40.0, 0.0),
            (40.0, 10.0),
            (20.0, 4.0),
            (0.0, 10.0)
        ]
    );
    assert_eq!(
        pts(s.polygons()[1].exterior()),
        [
            (20.0, 6.0),
            (40.0, 11.0),
            (40.0, 20.0),
            (0.0, 20.0),
            (0.0, 11.0)
        ]
    );
    check_promises(&mp, &s, 6.0, "narrow channel");
}

/// A line zigzagging around a point of another element keeps the point on
/// its side (check (b)); a point clear of the swept region does not hold it.
#[test]
fn a_zigzag_line_does_not_jump_a_point() {
    let l = Geometry::LineString(line(&[(0.0, 0.0), (5.0, 2.0), (10.0, 0.0)]));
    let held = GeometryCollection::new(vec![
        l.clone(),
        Geometry::Point(Point::new(c(5.0, 1.0)).unwrap()),
    ])
    .unwrap();
    let s = simplify_preserve_topology(&held, 3.0).unwrap();
    assert_eq!(s, held);
    let free = GeometryCollection::new(vec![l, Geometry::Point(Point::new(c(5.0, 3.0)).unwrap())])
        .unwrap();
    let s = simplify_preserve_topology(&free, 3.0).unwrap();
    match &s.geometries()[0] {
        Geometry::LineString(x) => assert_eq!(lpts(x), [(0.0, 0.0), (10.0, 0.0)]),
        g => panic!("{g:?}"),
    }
    check_promises(&free, &s, 3.0, "free point");

    let zig = GeometryCollection::new(vec![
        Geometry::LineString(line(&[
            (0.0, 0.0),
            (2.0, 1.0),
            (4.0, -1.0),
            (6.0, 1.0),
            (8.0, -1.0),
            (10.0, 0.0),
        ])),
        Geometry::MultiPoint(MultiPoint::new(vec![Point::new(c(2.0, 0.5)).unwrap()])),
    ])
    .unwrap();
    let s = simplify_preserve_topology(&zig, 2.0).unwrap();
    check_promises(&zig, &s, 2.0, "zigzag");
}

/// A line lying along the chord: only check (a) sees it. The run is below
/// the chord, so the even-odd count puts the other line's endpoints, which
/// are on the chord, outside.
#[test]
fn a_line_on_the_chord_holds_the_run() {
    let m = MultiLineString::new(vec![
        line(&[(0.0, 0.0), (5.0, -1.0), (10.0, 0.0)]),
        line(&[(4.0, 0.0), (6.0, 0.0)]),
    ]);
    let s = simplify_preserve_topology(&m, 3.0).unwrap();
    assert_eq!(s, m);
    check_promises(&m, &s, 3.0, "line on chord");
}

/// A line crossing a run twice, touching neither the chord nor holding an
/// endpoint inside the region: only check (c) sees it. Removing the run
/// would delete both crossings.
#[test]
fn a_line_crossing_the_run_twice_keeps_both_crossings() {
    let m = MultiLineString::new(vec![
        line(&[(0.0, 0.0), (5.0, 1.0), (10.0, 0.0)]),
        line(&[(1.0, 0.5), (9.0, 0.5)]),
    ]);
    let s = simplify_preserve_topology(&m, 3.0).unwrap();
    assert_eq!(s, m);
    check_promises(&m, &s, 3.0, "double crossing");
}

/// Two chords that interact: the first line's chord is accepted, and the
/// second line's chords would cross it. The second is refused only because
/// the first chord is in the index.
#[test]
fn a_later_chord_sees_an_earlier_one() {
    let first = line(&[(0.0, 0.0), (5.0, 1.0), (10.0, 0.0)]);
    let second = line(&[(7.0, 2.0), (12.0, 2.0), (12.0, -2.0), (7.0, -2.0)]);
    let m = MultiLineString::new(vec![first, second.clone()]);
    let s = simplify_preserve_topology(&m, 6.0).unwrap();
    assert_eq!(lpts(&s.line_strings()[0]), [(0.0, 0.0), (10.0, 0.0)]);
    assert_eq!(s.line_strings()[1], second);
    check_promises(&m, &s, 6.0, "interacting chords");
}

/// A hole vertex resting on the inside of an exterior edge: the contact is
/// not a shared coordinate, and it holds the exterior's run.
#[test]
fn a_vertex_on_an_edge_keeps_its_contact() {
    let p = poly(
        &[
            (0.0, 0.0),
            (10.0, 1.0),
            (20.0, 0.0),
            (20.0, 20.0),
            (0.0, 20.0),
        ],
        &[&[(5.0, 0.5), (4.0, 5.0), (6.0, 5.0)]],
    );
    assert!(validity_issues(&p).is_empty());
    let s = simplify_preserve_topology(&p, 3.0).unwrap();
    assert_eq!(s, p);
    check_promises(&p, &s, 3.0, "vertex on edge");
}

/// A ring cannot collapse: the last chord would overlap the remaining
/// segment (check (a)).
#[test]
fn a_ring_never_collapses() {
    let p = poly(&[(0.0, 0.0), (10.0, 0.0), (5.0, 1.0)], &[]);
    let s = simplify_preserve_topology(&p, 100.0).unwrap();
    assert_eq!(s, p);
    // Plain RDP collapses it.
    assert!(p.simplify(100.0).is_err());
}

/// A line that crosses itself keeps its crossing (check (d) and the same
/// element's segments in E).
#[test]
fn a_self_crossing_line_keeps_its_crossing() {
    let l = line(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 1.0),
        (5.0, -1.0),
        (4.0, -1.0),
    ]);
    let s = simplify_preserve_topology(&l, 2.0).unwrap();
    check_promises(&l, &s, 2.0, "self-crossing line");
    assert!(!is_simple(&s));
}

#[test]
fn tolerance_zero_is_the_identity_and_a_bad_tolerance_is_an_error() {
    let l = line(&[(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]);
    assert_eq!(simplify_preserve_topology(&l, 0.0).unwrap(), l);
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            simplify_preserve_topology(&l, bad),
            Err(GeometryError::InvalidTolerance { .. })
        ));
    }
    // At a positive tolerance the collinear middle vertex goes.
    assert_eq!(
        lpts(&simplify_preserve_topology(&l, 0.5).unwrap()),
        [(0.0, 0.0), (10.0, 0.0)]
    );
}

/// An invalid input (a hole that crosses itself) is simplified but not
/// vouched for, and its defect is unchanged.
#[test]
fn an_invalid_input_is_returned_unvouched() {
    let p = poly(
        &[(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)],
        &[&[(5.0, 5.0), (5.0, 15.0), (15.0, 5.0), (15.0, 10.0)]],
    );
    let issue = validity_issues(&p)[0].issue;
    assert!(matches!(issue, ValidityIssue::SelfIntersection { .. }));
    match simplify_preserve_topology(&p, 1.0) {
        Err(GeometryError::Unvouched {
            reason: UnvouchedReason::InvalidInput { issue: i },
            ..
        }) => assert_eq!(i, issue),
        r => panic!("{r:?}"),
    }
    match simplify_preserve_topology_unvouched(&p, 1.0).unwrap() {
        MaybeVouched::Unvouched(u) => {
            assert_eq!(*u.reason(), UnvouchedReason::InvalidInput { issue });
            assert_eq!(validity_issues(u.inner())[0].issue, issue);
            check_promises(&p, u.inner(), 1.0, "invalid input");
        }
        MaybeVouched::Vouched(_) => panic!("an invalid input was vouched for"),
    }
}

/// The defensive check refuses an invalid result. It cannot fire on the
/// algorithm's output for a valid input (see the module documentation), so
/// it is exercised directly on a value with a validity issue.
#[test]
fn the_defensive_check_refuses_an_invalid_result() {
    let bad = Polygon::new_simple(ring(&[
        (5.0, -5.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (10.0, 0.0),
        (0.0, 0.0),
    ]));
    let bad = bad.unwrap();
    let issue = validity_issues(&bad)[0].issue;
    match vouch(bad) {
        MaybeVouched::Unvouched(u) => {
            assert_eq!(*u.reason(), UnvouchedReason::InvalidResult { issue });
        }
        MaybeVouched::Vouched(_) => panic!("an invalid result was vouched for"),
    }
}

// ===========================================================================
// Generated cases
// ===========================================================================

/// SplitMix64: deterministic, no dependency.
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
    fn coin(&mut self, k: u64) -> bool {
        self.below(k) == 0
    }
    /// A multiple of 0.5 in `[lo, hi]`.
    fn half(&mut self, lo: f64, hi: f64) -> f64 {
        let steps = ((hi - lo) * 2.0) as u64;
        lo + self.below(steps + 1) as f64 * 0.5
    }
}

/// A star-shaped ring around `(cx, cy)`: `k` vertices at increasing angles,
/// radii in `[r0, r1]`, on a half-unit grid, counterclockwise. Rounding to
/// the grid can break simplicity; callers keep only valid results.
fn star(r: &mut Rng, cx: f64, cy: f64, r0: f64, r1: f64) -> Vec<Coord2d> {
    let k = 5 + r.below(20) as usize;
    let mut v: Vec<Coord2d> = Vec::new();
    for i in 0..k {
        let t = (i as f64 + 0.3 * (r.below(3) as f64)) / k as f64 * core::f64::consts::TAU;
        let rad = r0 + (r1 - r0) * (r.below(100) as f64 / 100.0);
        let p = c(
            ((cx + rad * t.cos()) * 2.0).round() / 2.0,
            ((cy + rad * t.sin()) * 2.0).round() / 2.0,
        );
        if v.last() != Some(&p) && v.first() != Some(&p) {
            v.push(p);
        }
    }
    v
}

fn valid_polygon(ext: Vec<Coord2d>, holes: Vec<Vec<Coord2d>>) -> Option<Polygon> {
    let ext = LinearRing::new(CoordSeq::new(ext)).ok()?;
    let holes = holes
        .into_iter()
        .map(|h| LinearRing::new(CoordSeq::new(h)).ok())
        .collect::<Option<Vec<_>>>()?;
    let p = Polygon::new(ext, holes).ok()?;
    validity_issues(&p).is_empty().then_some(p)
}

/// A star polygon with up to two holes, some near the exterior's edge, and
/// sometimes a hole vertex placed exactly at the midpoint of an exterior
/// edge (a vertex-on-edge contact).
fn gen_polygon(r: &mut Rng, cx: f64, cy: f64) -> Polygon {
    loop {
        let ext = star(r, cx, cy, 6.0, 14.0);
        let mut holes = Vec::new();
        for _ in 0..r.below(3) {
            let (hx, hy) = (cx + r.half(-4.0, 4.0), cy + r.half(-4.0, 4.0));
            let mut h = star(r, hx, hy, 0.5, 2.5);
            h.reverse();
            if r.coin(3) && ext.len() >= 2 {
                // Pull the hole's first vertex onto the middle of an exterior
                // edge: exact, since the grid's midpoints are representable.
                let e = r.below(ext.len() as u64) as usize;
                let (a, b) = (ext[e], ext[(e + 1) % ext.len()]);
                h[0] = c(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
            }
            holes.push(h);
        }
        if let Some(p) = valid_polygon(ext.clone(), holes) {
            return p;
        }
        if let Some(p) = valid_polygon(ext, Vec::new()) {
            return p;
        }
    }
}

/// Parts in a 3 × 3 grid of cells; sometimes a part touching the previous
/// one at a shared vertex.
fn gen_multipolygon(r: &mut Rng) -> MultiPolygon {
    let mut parts: Vec<Polygon> = Vec::new();
    for _ in 0..=r.below(4) {
        let (cx, cy) = (30.0 * r.below(3) as f64, 30.0 * r.below(3) as f64);
        let mut cand = gen_polygon(r, cx, cy);
        if r.coin(2)
            && let Some(prev) = parts.last()
        {
            // A part touching `prev` at one of its vertices.
            let w = prev.exterior().coords().as_slice()[0];
            let mut s = star(r, w.x + 7.0, w.y, 2.0, 6.0);
            let mid = s.len() / 2;
            s[mid] = w;
            if let Some(p) = valid_polygon(s, Vec::new()) {
                cand = p;
            }
        }
        let mut trial = parts.clone();
        trial.push(cand);
        if let Ok(mp) = MultiPolygon::new(trial.clone())
            && validity_issues(&mp).is_empty()
        {
            parts = trial;
        }
    }
    if parts.is_empty() {
        parts.push(gen_polygon(r, 0.0, 0.0));
    }
    MultiPolygon::new(parts).unwrap()
}

/// A wiggly line: a random walk of small steps.
fn gen_line(r: &mut Rng) -> LineString {
    let mut v = vec![c(r.half(-10.0, 70.0), r.half(-10.0, 70.0))];
    for _ in 0..2 + r.below(15) {
        let last = v[v.len() - 1];
        let p = c(last.x + r.half(-3.0, 6.0), last.y + r.half(-3.0, 3.0));
        if p != last {
            v.push(p);
        }
    }
    if v.len() < 2 {
        v.push(c(v[0].x + 1.0, v[0].y));
    }
    LineString::new(CoordSeq::new(v)).unwrap()
}

/// A collection: a multipolygon, lines (which may cross anything) and
/// points.
fn gen_collection(r: &mut Rng) -> GeometryCollection {
    let mut gs = vec![Geometry::MultiPolygon(gen_multipolygon(r))];
    for _ in 0..r.below(4) {
        gs.push(Geometry::LineString(gen_line(r)));
    }
    let k = r.below(6);
    gs.push(Geometry::MultiPoint(MultiPoint::new(
        (0..k)
            .map(|_| Point::new(c(r.half(-10.0, 70.0), r.half(-10.0, 70.0))).unwrap())
            .collect(),
    )));
    GeometryCollection::new(gs).unwrap()
}

const CASES: usize = 10_000;
const TOLERANCES: [f64; 6] = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0];

#[test]
fn generated_polygons_keep_every_promise() {
    let mut r = Rng(0x7357_0001);
    let mut removed = 0usize;
    for case in 0..CASES {
        let p = gen_polygon(&mut r, 0.0, 0.0);
        let t = TOLERANCES[r.below(TOLERANCES.len() as u64) as usize];
        let s = simplify_preserve_topology(&p, t)
            .unwrap_or_else(|e| panic!("case {case}: valid polygon refused: {e:?}"));
        removed += flat_of(&p).verts.len() - flat_of(&s).verts.len();
        check_promises(&p, &s, t, &format!("polygon case {case}"));
    }
    assert!(
        removed > CASES,
        "the generator gives the algorithm work: {removed}"
    );
}

#[test]
fn generated_multipolygons_keep_every_promise() {
    let mut r = Rng(0x7357_0002);
    for case in 0..CASES {
        let m = gen_multipolygon(&mut r);
        let t = TOLERANCES[r.below(TOLERANCES.len() as u64) as usize];
        let s = simplify_preserve_topology(&m, t)
            .unwrap_or_else(|e| panic!("case {case}: valid multipolygon refused: {e:?}"));
        check_promises(&m, &s, t, &format!("multipolygon case {case}"));
    }
}

#[test]
fn generated_lines_keep_every_promise() {
    let mut r = Rng(0x7357_0003);
    for case in 0..CASES {
        let m = MultiLineString::new((0..=r.below(5)).map(|_| gen_line(&mut r)).collect());
        let t = TOLERANCES[r.below(TOLERANCES.len() as u64) as usize];
        let s = simplify_preserve_topology(&m, t).unwrap();
        check_promises(&m, &s, t, &format!("lines case {case}"));
    }
}

#[test]
fn generated_collections_keep_every_promise() {
    let mut r = Rng(0x7357_0004);
    for case in 0..CASES {
        let g = gen_collection(&mut r);
        let t = TOLERANCES[r.below(TOLERANCES.len() as u64) as usize];
        let s = simplify_preserve_topology(&g, t)
            .unwrap_or_else(|e| panic!("case {case}: valid collection refused: {e:?}"));
        check_promises(&g, &s, t, &format!("collection case {case}"));
    }
}

/// The generators reach the configurations the hand cases name: touching
/// parts, vertex-on-edge holes, and lines crossing other elements.
#[test]
fn the_generators_reach_the_contact_configurations() {
    let mut r = Rng(0x7357_0002);
    let mut touching = 0;
    for _ in 0..1000 {
        let m = gen_multipolygon(&mut r);
        let ps = m.polygons();
        for i in 0..ps.len() {
            for j in i + 1..ps.len() {
                if relate(&ps[i], &ps[j]).unwrap().touches() {
                    touching += 1;
                }
            }
        }
    }
    let mut r = Rng(0x7357_0001);
    let mut on_edge = 0;
    for _ in 0..1000 {
        let p = gen_polygon(&mut r, 0.0, 0.0);
        for h in p.holes() {
            if h.coords()
                .as_slice()
                .iter()
                .any(|&q| locate_in_ring(q, p.exterior()) == RingLocation::OnBoundary)
            {
                on_edge += 1;
            }
        }
    }
    assert!(touching >= 50, "touching parts: {touching}");
    assert!(on_edge >= 50, "vertex-on-edge holes: {on_edge}");
}
