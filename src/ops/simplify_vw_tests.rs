// src/ops/simplify_vw_tests.rs

//! Two checks on `ops::simplify_vw`. First, the elimination is the
//! specified one (an independent scan, the properties, the veto). Second,
//! the vouch refuses exactly the results `validate` rejects.
//!
//! The scan reference here shares nothing with the module but the area
//! expression, which the specification fixes: it recomputes every
//! removable vertex's area from its current neighbors on every step and
//! takes the smallest `(area, index)`. The queue-based core must select
//! exactly the same vertices in the same order.

#![allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

use super::simplify_vw::{ChainKind, Elimination, effective_area, eliminate};
use super::{simplify_vw, validate};
use crate::GeometryError;
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString,
    MultiPoint, MultiPolygon, Point, Polygon, UnvouchedReason,
};

// ===========================================================================
// Fixtures
// ===========================================================================

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
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn line(v: &[(f64, f64)]) -> LineString {
    LineString::new(CoordSeq::new(v.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

/// A random chain: grid (ties and collinear runs common), continuous, or
/// with repeats of earlier vertices; consecutive repeats removed.
fn random_chain(rng: &mut Rng, kind: ChainKind) -> Vec<Coord2d> {
    let n = 3 + rng.below(40) as usize;
    let style = rng.below(3);
    let off = [0.0, 1e3, 1e9][rng.below(3) as usize];
    let mut v: Vec<Coord2d> = Vec::new();
    while v.len() < n {
        let p = match style {
            0 => c(off + rng.below(8) as f64, off + rng.below(8) as f64),
            1 if !v.is_empty() && rng.below(6) == 0 => v[rng.below(v.len() as u64) as usize],
            _ => c(off + rng.unit() * 100.0, off + rng.unit() * 100.0),
        };
        if v.last() != Some(&p) {
            v.push(p);
        }
    }
    if kind == ChainKind::Ring {
        while v.len() >= 2 && v[0] == v[v.len() - 1] {
            v.pop();
        }
        while v.len() < 3 {
            let p = c(off + v.len() as f64 * 7.0 + 0.5, off + 3.0);
            v.push(p);
        }
    } else if rng.below(5) == 0 {
        // A closed line now and then.
        let first = v[0];
        if v[v.len() - 1] != first {
            v.push(first);
        }
    }
    v
}

/// A threshold aimed at the chain's own initial areas, so that some but
/// not all vertices go; now and then 0 or a huge value.
fn threshold_for(rng: &mut Rng, pts: &[Coord2d], kind: ChainKind) -> f64 {
    let n = pts.len();
    let mut areas: Vec<f64> = (0..n)
        .filter_map(|i| {
            let (a, b) = if kind == ChainKind::Ring {
                ((i + n - 1) % n, (i + 1) % n)
            } else if i == 0 || i + 1 == n {
                return None;
            } else {
                (i - 1, i + 1)
            };
            Some(effective_area(pts[a], pts[i], pts[b]))
        })
        .collect();
    areas.sort_by(f64::total_cmp);
    match rng.below(8) {
        0 => 0.0,
        1 => 1e300,
        2 if !areas.is_empty() => areas[rng.below(areas.len() as u64) as usize],
        _ if !areas.is_empty() => {
            let k = areas[rng.below(areas.len() as u64) as usize];
            k * (0.5 + rng.unit() * 2.0)
        }
        _ => 1.0,
    }
}

// ===========================================================================
// The independent scan
// ===========================================================================

/// The specification as a plain scan: on every step, every removable vertex
/// (alive, not fixed, neighbors distinct) is measured against its current
/// neighbors, and the smallest `(area, index)` is selected.
fn scan(
    pts: &[Coord2d],
    kind: ChainKind,
    threshold: f64,
    allow: &mut dyn FnMut(usize, usize, usize) -> bool,
) -> Elimination {
    let n = pts.len();
    let mut alive: Vec<usize> = (0..n).collect();
    let mut keep = vec![true; n];
    let mut fixed = vec![false; n];
    let mut removed = Vec::new();
    let closed_line = kind == ChainKind::Line && n >= 2 && pts[0] == pts[n - 1];
    let floor = if kind == ChainKind::Ring || closed_line {
        3
    } else {
        2
    };
    let mut last = 0.0_f64;
    loop {
        if alive.len() <= floor {
            break;
        }
        let m = alive.len();
        let mut best: Option<(f64, usize, usize)> = None; // area, vertex, slot
        for slot in 0..m {
            let i = alive[slot];
            let ends = kind == ChainKind::Line && (slot == 0 || slot + 1 == m);
            if ends || fixed[i] {
                continue;
            }
            let (a, b) = (alive[(slot + m - 1) % m], alive[(slot + 1) % m]);
            if pts[a] == pts[b] {
                continue;
            }
            let area = effective_area(pts[a], pts[i], pts[b]);
            if best.is_none_or(|(ba, bi, _)| area < ba || (area == ba && i < bi)) {
                best = Some((area, i, slot));
            }
        }
        let Some((area, i, slot)) = best else { break };
        if !(area < threshold || area == 0.0) {
            break;
        }
        let (a, b) = (alive[(slot + m - 1) % m], alive[(slot + 1) % m]);
        if !allow(i, a, b) {
            fixed[i] = true;
            continue;
        }
        keep[i] = false;
        alive.remove(slot);
        last = area.max(last);
        removed.push((i, last));
    }
    Elimination { keep, removed }
}

// ===========================================================================
// The elimination is the specified one
// ===========================================================================

#[test]
fn the_core_equals_a_plain_scan() {
    let mut rng = Rng(0x5157_7A11);
    let mut removals = 0_usize;
    for case in 0..6_000 {
        let kind = if case % 2 == 0 {
            ChainKind::Line
        } else {
            ChainKind::Ring
        };
        let pts = random_chain(&mut rng, kind);
        let t = threshold_for(&mut rng, &pts, kind);
        let got = eliminate(&pts, kind, t, &mut |_, _, _| true);
        let want = scan(&pts, kind, t, &mut |_, _, _| true);
        assert_eq!(got.keep, want.keep, "kept vertices, case {case}");
        assert_eq!(got.removed.len(), want.removed.len());
        for (g, w) in got.removed.iter().zip(&want.removed) {
            assert_eq!(g.0, w.0, "removal order, case {case}");
            assert_eq!(g.1.to_bits(), w.1.to_bits(), "recorded area, case {case}");
        }
        removals += got.removed.len();
    }
    assert!(removals > 20_000, "removals exercised: {removals}");
}

/// The elimination's properties, checked on one elimination.
fn check_properties(pts: &[Coord2d], kind: ChainKind, t: f64, e: &Elimination) {
    let n = pts.len();
    let kept: Vec<usize> = (0..n).filter(|&i| e.keep[i]).collect();
    if kind == ChainKind::Line {
        assert!(e.keep[0] && e.keep[n - 1], "line endpoints kept");
    } else {
        assert!(kept.len() >= 3.min(n), "a ring keeps at least three");
    }
    // Every recorded area below the threshold (or 0), and monotone.
    let mut last = 0.0_f64;
    for &(_, a) in &e.removed {
        assert!(a < t || a == 0.0, "removed at {a} with threshold {t}");
        assert!(a >= last, "recorded areas are monotone");
        last = a;
    }
    // On the result: every remaining removable vertex is at or above the
    // threshold and not zero, unless the floor stopped the elimination.
    let m = kept.len();
    let closed_line = kind == ChainKind::Line && n >= 2 && pts[0] == pts[n - 1];
    let floor = if kind == ChainKind::Ring || closed_line {
        3
    } else {
        2
    };
    if m > floor {
        for slot in 0..m {
            if kind == ChainKind::Line && (slot == 0 || slot + 1 == m) {
                continue;
            }
            let (a, i, b) = (kept[(slot + m - 1) % m], kept[slot], kept[(slot + 1) % m]);
            if pts[a] == pts[b] {
                continue;
            }
            let area = effective_area(pts[a], pts[i], pts[b]);
            assert!(area >= t && area != 0.0, "vertex {i} left at {area} < {t}");
        }
    }
    // No consecutive repeats were created.
    for s in 0..m {
        if kind == ChainKind::Line && s + 1 == m {
            break;
        }
        assert_ne!(pts[kept[s]], pts[kept[(s + 1) % m]], "consecutive repeat");
    }
}

#[test]
fn the_properties_hold_on_every_case() {
    let mut rng = Rng(0xB0B0_5EED);
    for case in 0..6_000 {
        let kind = if case % 3 == 0 {
            ChainKind::Ring
        } else {
            ChainKind::Line
        };
        let pts = random_chain(&mut rng, kind);
        let t = threshold_for(&mut rng, &pts, kind);
        let e = eliminate(&pts, kind, t, &mut |_, _, _| true);
        check_properties(&pts, kind, t, &e);
    }
}

#[test]
fn endpoints_and_the_ring_floor_are_kept() {
    // Every interior vertex is tiny; the endpoints are not removable.
    let l = line(&[(0.0, 0.0), (1.0, 0.01), (2.0, 0.0), (3.0, 0.01), (4.0, 0.0)]);
    let out = simplify_vw(&l, 1e9).unwrap();
    assert_eq!(out.coords().as_slice(), &[c(0.0, 0.0), c(4.0, 0.0)]);
    // A closed line keeps one vertex between its equal ends.
    let closed = line(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]);
    let out = simplify_vw(&closed, 1e9).unwrap();
    assert_eq!(out.num_points(), 3);
    // A ring stops at three.
    let ring: Vec<Coord2d> = (0..12)
        .map(|k| {
            let a = f64::from(k) * std::f64::consts::TAU / 12.0;
            c(a.cos() * 10.0, a.sin() * 10.0)
        })
        .collect();
    let e = eliminate(&ring, ChainKind::Ring, 1e9, &mut |_, _, _| true);
    assert_eq!(e.keep.iter().filter(|&&k| k).count(), 3);
}

#[test]
fn threshold_zero_removes_exactly_the_zero_areas() {
    // Collinear runs go; a bend stays.
    let l = line(&[
        (0.0, 0.0),
        (1.0, 0.0),
        (2.0, 0.0),
        (3.0, 1.0),
        (4.0, 2.0),
        (5.0, 2.0),
    ]);
    let out = simplify_vw(&l, 0.0).unwrap();
    assert_eq!(
        out.coords().as_slice(),
        &[c(0.0, 0.0), c(2.0, 0.0), c(4.0, 2.0), c(5.0, 2.0)]
    );
    // A vertex whose neighbors coincide is never removed (it would leave a
    // consecutive repeat), although its area is 0.
    let spike = line(&[(0.0, 0.0), (5.0, 5.0), (0.0, 0.0), (9.0, 0.0)]);
    let out = simplify_vw(&spike, 0.0).unwrap();
    assert_eq!(out.num_points(), 4);
}

#[test]
fn invalid_thresholds_are_refused() {
    let l = line(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]);
    for t in [-1.0, f64::NAN, f64::INFINITY, (-0.0_f64).next_down()] {
        assert!(matches!(
            simplify_vw(&l, t),
            Err(GeometryError::InvalidTolerance { .. })
        ));
    }
}

// ===========================================================================
// The veto (crate-private; coverage simplification relies on it)
// ===========================================================================

#[test]
fn the_veto_is_honored() {
    let mut rng = Rng(0x7E70_0001);
    let mut vetoed_total = 0_usize;
    for case in 0..3_000 {
        let kind = if case % 2 == 0 {
            ChainKind::Line
        } else {
            ChainKind::Ring
        };
        let pts = random_chain(&mut rng, kind);
        let t = threshold_for(&mut rng, &pts, kind);
        // Reject every vertex whose index is a multiple of 3.
        let mut asked = Vec::new();
        let got = eliminate(&pts, kind, t, &mut |i, _, _| {
            asked.push(i);
            i % 3 != 0
        });
        let want = scan(&pts, kind, t, &mut |i, _, _| i % 3 != 0);
        assert_eq!(got, want, "case {case}");
        for i in asked.iter().filter(|&&i| i % 3 == 0) {
            assert!(got.keep[*i], "a vetoed vertex is kept");
            vetoed_total += 1;
        }
    }
    assert!(vetoed_total > 1_000, "vetoes exercised: {vetoed_total}");
    // The callback sees the vertex's current neighbors.
    let pts = vec![
        c(0.0, 0.0),
        c(1.0, 0.1),
        c(2.0, 0.0),
        c(3.0, 0.2),
        c(4.0, 0.0),
    ];
    let mut seen = Vec::new();
    eliminate(&pts, ChainKind::Line, 10.0, &mut |i, a, b| {
        seen.push((i, a, b));
        true
    });
    assert_eq!(seen[0], (1, 0, 2));
    assert!(seen.iter().all(|&(i, a, b)| a < i && i < b));
}

// ===========================================================================
// The vouch refuses exactly what `validate` rejects
// ===========================================================================

/// A comb: a base rectangle with thin teeth along its top, randomly
/// jittered, so that removing small-area vertices can cut across a
/// neighboring tooth. Counterclockwise. Optionally a hole in the base.
fn comb(rng: &mut Rng, ox: f64, oy: f64, hole: bool) -> Polygon {
    let teeth = 3 + rng.below(6) as usize;
    let w = 10.0;
    let width = teeth as f64 * w;
    let mut v = vec![c(ox, oy), c(ox + width, oy), c(ox + width, oy + 5.0)];
    for k in (0..teeth).rev() {
        let x0 = ox + k as f64 * w;
        let h = 3.0 + rng.unit() * 20.0;
        let gap = 0.05 + rng.unit() * 2.0;
        let lean = (rng.unit() - 0.5) * 3.0;
        v.push(c(x0 + w - gap, oy + 5.0));
        v.push(c(x0 + w - gap + lean, oy + 5.0 + h));
        v.push(c(x0 + gap * 0.5 + lean, oy + 5.0 + h + rng.unit()));
        v.push(c(x0 + gap * 0.5, oy + 5.0 + rng.unit() * 0.5));
    }
    v.push(c(ox, oy + 5.0));
    v.dedup();
    let ext = LinearRing::new(CoordSeq::new(v)).unwrap();
    let holes = if hole {
        let (hx, hy) = (ox + 1.0 + rng.unit() * (width - 4.0), oy + 1.0);
        vec![
            LinearRing::new(CoordSeq::new(vec![
                c(hx, hy),
                c(hx, hy + 2.5 + rng.unit()),
                c(hx + 1.0 + rng.unit(), hy + 3.9),
                c(hx + 2.0, hy + 0.01),
            ]))
            .unwrap(),
        ]
    } else {
        Vec::new()
    };
    Polygon::new(ext, holes).unwrap_or_else(|_| {
        Polygon::new(
            LinearRing::new(CoordSeq::new(vec![
                c(ox, oy),
                c(ox + 1.0, oy),
                c(ox + 1.0, oy + 1.0),
            ]))
            .unwrap(),
            vec![],
        )
        .unwrap()
    })
}

/// The kept rings of `p`, with no constructor or vouch in between.
fn raw_rings(p: &Polygon, t: f64) -> Vec<CoordSeq> {
    std::iter::once(p.exterior())
        .chain(p.holes())
        .map(|r| {
            let pts = r.coords().as_slice();
            let e = eliminate(pts, ChainKind::Ring, t, &mut |_, _, _| true);
            CoordSeq::new(
                pts.iter()
                    .zip(&e.keep)
                    .filter_map(|(&p, &k)| k.then_some(p))
                    .collect(),
            )
        })
        .collect()
}

/// `(returned, refused, refused by SimplifiedInvalid)` over `polys` at `t`,
/// asserting on every polygon that a result is returned exactly when
/// `validate` accepts the raw kept rings.
fn vouch_matches_validate(polys: &[Polygon], t: f64) -> (usize, usize, usize) {
    let (mut ok, mut refused, mut invalid) = (0, 0, 0);
    for p in polys {
        let valid = validate(&raw_rings(p, t)).is_empty();
        match simplify_vw(p, t) {
            Ok(q) => {
                assert!(valid, "returned a result validate rejects: {p:?} at {t}");
                assert!(crate::ops::validity_issues(&q).is_empty());
                ok += 1;
            }
            Err(e) => {
                assert!(
                    !valid,
                    "refused a result validate accepts ({e}): {p:?} at {t}"
                );
                refused += 1;
                invalid += usize::from(matches!(
                    e,
                    GeometryError::Unvouched {
                        reason: UnvouchedReason::SimplifiedInvalid { .. },
                        ..
                    }
                ));
            }
        }
    }
    (ok, refused, invalid)
}

#[test]
fn the_vouch_refuses_exactly_what_validate_rejects() {
    let mut rng = Rng(0x0C0B_0B5E);
    let polys: Vec<Polygon> = (0..600)
        .map(|k| comb(&mut rng, f64::from(k % 7) * 1e3, 0.0, k % 2 == 0))
        .collect();
    let (mut ok, mut refused, mut invalid) = (0, 0, 0);
    for t in [0.0, 0.1, 0.5, 2.0, 8.0, 30.0] {
        let (o, r, i) = vouch_matches_validate(&polys, t);
        ok += o;
        refused += r;
        invalid += i;
    }
    assert!(ok > 500, "returned: {ok}");
    assert!(refused > 100, "refused: {refused}");
    assert!(
        invalid > 50,
        "refused by the validity check itself: {invalid}"
    );
}

#[test]
fn every_type_simplifies_and_collections_keep_their_points() {
    let l = line(&[(0.0, 0.0), (5.0, 0.1), (10.0, 0.0), (10.0, 10.0)]);
    assert_eq!(simplify_vw(&l, 1.0).unwrap().num_points(), 3);
    let ml = MultiLineString::new(vec![l.clone(), l.clone()]);
    let out = simplify_vw(&ml, 1.0).unwrap();
    assert!(out.line_strings().iter().all(|x| x.num_points() == 3));
    let sq = |x: f64| {
        Polygon::new(
            LinearRing::new(CoordSeq::new(vec![
                c(x, 0.0),
                c(x + 5.0, 0.01),
                c(x + 10.0, 0.0),
                c(x + 10.0, 10.0),
                c(x, 10.0),
            ]))
            .unwrap(),
            vec![],
        )
        .unwrap()
    };
    assert_eq!(
        simplify_vw(&sq(0.0), 1.0)
            .unwrap()
            .exterior()
            .num_vertices(),
        4
    );
    let mp = MultiPolygon::new(vec![sq(0.0), sq(20.0)]).unwrap();
    assert!(
        simplify_vw(&mp, 1.0)
            .unwrap()
            .polygons()
            .iter()
            .all(|p| p.exterior().num_vertices() == 4)
    );
    let gc = GeometryCollection::new(vec![
        Geometry::Point(Point::new(c(1.0, 1.0)).unwrap()),
        Geometry::MultiPoint(MultiPoint::new(vec![Point::new(c(2.0, 2.0)).unwrap()])),
        Geometry::LineString(l),
        Geometry::Polygon(sq(40.0)),
    ])
    .unwrap();
    let out = simplify_vw(&gc, 1.0).unwrap();
    assert_eq!(out.geometries()[0], gc.geometries()[0]);
    assert_eq!(out.geometries()[1], gc.geometries()[1]);
}

// ===========================================================================
// External polygon files: parsing, and the vouch check on them
// ===========================================================================

/// The text of every file named in the environment variable `var`, a list
/// of paths separated by the platform's path separator (`;` on Windows,
/// `:` elsewhere), concatenated; empty when the variable is unset. A named
/// file that cannot be read is an error, so a mistyped path is not a
/// silently empty run.
fn files_named_by(var: &str) -> String {
    let Some(list) = std::env::var_os(var) else {
        return String::new();
    };
    std::env::split_paths(&list)
        .map(|p| {
            std::fs::read_to_string(&p)
                .unwrap_or_else(|e| panic!("{var} names {}: {e}", p.display()))
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn numbers(line: &str, skip: usize) -> Vec<Coord2d> {
    let v: Vec<f64> = line
        .split_whitespace()
        .skip(skip)
        .filter_map(|t| t.parse().ok())
        .collect();
    v.chunks_exact(2).map(|p| c(p[0], p[1])).collect()
}

/// `POLY n` / `RING k x y …` records: each record's rings, closing vertex
/// dropped, consecutive repeats removed.
fn poly_records(text: &str) -> Vec<Vec<Vec<Coord2d>>> {
    let mut out = Vec::new();
    let mut lines = text.lines();
    while let Some(l) = lines.next() {
        let mut it = l.split_whitespace();
        if it.next() != Some("POLY") {
            continue;
        }
        let n: usize = it.next().and_then(|v| v.parse().ok()).unwrap_or(0);
        let mut rings = Vec::new();
        for _ in 0..n {
            let Some(rl) = lines.next() else { break };
            let mut r = numbers(rl, 2);
            r.dedup();
            while r.len() > 1 && r[0] == r[r.len() - 1] {
                r.pop();
            }
            rings.push(r);
        }
        out.push(rings);
    }
    out
}

/// `MLINE n` / `LINE k x y …` records: every line's vertices, consecutive
/// repeats removed.
fn line_records(text: &str) -> Vec<Vec<Coord2d>> {
    text.lines()
        .filter(|l| l.starts_with("LINE "))
        .map(|l| {
            let mut v = numbers(l, 2);
            v.dedup();
            v
        })
        .filter(|v| v.len() >= 2)
        .collect()
}

/// A record's rings as polygons: exterior counterclockwise, holes
/// clockwise, sorted by area sign as a clockwise-outside file requires
/// (a negative ring opens a part, a positive one is a hole of the last).
fn polygons_of(records: &[Vec<Vec<Coord2d>>]) -> Vec<Polygon> {
    let mut out = Vec::new();
    for rec in records {
        let mut parts: Vec<(LinearRing, Vec<LinearRing>)> = Vec::new();
        let mut ok = true;
        for r in rec {
            let Ok(ring) = LinearRing::new(CoordSeq::new(r.clone())) else {
                ok = false;
                break;
            };
            let mut rev = r.clone();
            rev.reverse();
            let reversed = LinearRing::new(CoordSeq::new(rev)).unwrap();
            let a = ring.signed_area();
            if a < 0.0 {
                parts.push((reversed, Vec::new()));
            } else if let Some(last) = parts.last_mut() {
                last.1.push(reversed);
            } else {
                parts.push((ring, Vec::new()));
            }
        }
        if ok {
            out.extend(
                parts
                    .into_iter()
                    .filter_map(|(e, h)| Polygon::new(e, h).ok()),
            );
        }
    }
    out
}

/// The polygons of the files named by `GEOMETRES_VW_POLYGONS`.
fn external_polygons() -> Vec<Polygon> {
    polygons_of(&poly_records(&files_named_by("GEOMETRES_VW_POLYGONS")))
}

/// The vouch check on the polygons of the files named by
/// `GEOMETRES_VW_POLYGONS`, with the refusal rates beside `simplify`'s
/// (Douglas–Peucker) at the comparable setting: tolerance `t` against
/// threshold `t²`. Writes `vw_vouch.txt` in the temporary directory.
///
/// The files hold `POLY n` / `RING k x y …` records, one polygon per
/// record, rings in either winding. The synthetic version of this check is
/// [`the_vouch_refuses_exactly_what_validate_rejects`].
#[test]
#[ignore = "needs GEOMETRES_VW_POLYGONS to name polygon files"]
fn the_vouch_matches_validate_on_external_polygons() {
    use std::fmt::Write as _;
    let polys = external_polygons();
    assert!(!polys.is_empty(), "GEOMETRES_VW_POLYGONS names no polygons");
    let mut out = format!("{} polygons\n", polys.len());
    for t in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0] {
        let area = t * t;
        let (ok, refused, invalid) = vouch_matches_validate(&polys, area);
        let rdp_refused = polys.iter().filter(|p| p.simplify(t).is_err()).count();
        writeln!(
            out,
            "t = {t}: simplify_vw(t^2) returned {ok}, refused {refused} ({invalid} by the validity \
             check), rate {:.4}; simplify(t) refused {rdp_refused}, rate {:.4}",
            refused as f64 / polys.len() as f64,
            rdp_refused as f64 / polys.len() as f64
        )
        .unwrap();
    }
    print!("{out}");
    let path = std::env::temp_dir().join("vw_vouch.txt");
    std::fs::write(&path, out).unwrap();
    println!("wrote {}", path.display());
}

// ===========================================================================
// The elimination case writer, for an independent exact check
// ===========================================================================

/// Writes `vw_cases.csv` in the temporary directory, one chain per line,
/// `kind;threshold;coords;kept;removed` with `kind` `L` (line) or `R`
/// (ring), floats as hex bit patterns, `coords` as `x,y` joined by `|`,
/// `kept` the kept indices joined by spaces, and `removed` the removal
/// order as `i:area` joined by `|`. Random chains, then every ring of the
/// files named by `GEOMETRES_VW_POLYGONS` and every line of the files named
/// by `GEOMETRES_VW_LINES` (`LINE k x y …` records), when they are set.
#[test]
#[ignore = "writes the elimination case file vw_cases.csv to the temporary directory"]
fn write_vw_cases_csv() {
    use std::fmt::Write as _;
    let mut rng = Rng(0x7657_2026);
    let mut out = String::new();
    let emit = |rng: &mut Rng, pts: &[Coord2d], kind: ChainKind, out: &mut String| {
        if pts.len() < 2 || (kind == ChainKind::Ring && pts.len() < 3) {
            return 0;
        }
        let t = threshold_for(rng, pts, kind);
        let e = eliminate(pts, kind, t, &mut |_, _, _| true);
        let coords: Vec<String> = pts
            .iter()
            .map(|p| format!("{:x},{:x}", p.x.to_bits(), p.y.to_bits()))
            .collect();
        let kept: Vec<String> = (0..pts.len())
            .filter(|&i| e.keep[i])
            .map(|i| i.to_string())
            .collect();
        let removed: Vec<String> = e
            .removed
            .iter()
            .map(|(i, a)| format!("{i}:{:x}", a.to_bits()))
            .collect();
        writeln!(
            out,
            "{};{:x};{};{};{}",
            if kind == ChainKind::Line { 'L' } else { 'R' },
            t.to_bits(),
            coords.join("|"),
            kept.join(" "),
            removed.join("|")
        )
        .unwrap();
        1
    };
    let mut written = 0_usize;
    for case in 0..20_500 {
        let kind = if case % 2 == 0 {
            ChainKind::Line
        } else {
            ChainKind::Ring
        };
        let pts = random_chain(&mut rng, kind);
        written += emit(&mut rng, &pts, kind, &mut out);
    }
    let mut external = 0_usize;
    for rec in poly_records(&files_named_by("GEOMETRES_VW_POLYGONS")) {
        for r in rec {
            external += emit(&mut rng, &r, ChainKind::Ring, &mut out);
        }
    }
    for l in line_records(&files_named_by("GEOMETRES_VW_LINES")) {
        external += emit(&mut rng, &l, ChainKind::Line, &mut out);
    }
    let path = std::env::temp_dir().join("vw_cases.csv");
    std::fs::write(&path, out).unwrap();
    println!(
        "wrote {written} random and {external} external chains to {}",
        path.display()
    );
}
