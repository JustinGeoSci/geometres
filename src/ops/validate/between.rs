// src/ops/validate/between.rs

//! Issues between the rings of one polygon: placement (a hole outside the
//! exterior, nested or crossing rings, duplicate rings) and the connectivity
//! of the interior.
//!
//! # Placement
//!
//! Each hole is placed against the exterior, then each pair of holes against
//! each other, in the order `Polygon::new` checks them, by the same exact
//! [`place_ring`] it uses. So a set of rings draws a placement issue here if
//! and only if `Polygon::new` refuses its placement, and the first issue
//! names the pair `Polygon::new` names.
//!
//! # Connectivity: the theorem and its proof
//!
//! **Setting:** Simple rings (an exterior and holes) such that every hole
//! lies inside the exterior, the holes are interior-disjoint, and any two
//! rings meet in at most finitely many points: exactly what `Polygon::new`
//! accepts, plus ring simplicity, which construction does not check.
//!
//! **The touch graph:** One node per ring, one node per distinct point where
//! two or more rings meet, and an edge joining a ring to a point wherever
//! the ring passes through the point.
//!
//! **Theorem:** The interior of the polygon is connected if and only if the
//! touch graph has no cycle. More precisely, the interior has `1 + β`
//! components, where `β` is the number of independent cycles of the graph.
//!
//! **Proof:** Work on the sphere `S²`, the plane plus a point at infinity.
//! By the Jordan–Schoenflies theorem each simple ring bounds a closed disk on
//! each side. Let `K` be the union of the closed disk outside the exterior
//! (the side containing infinity) and the closed disks of the holes. The
//! interior of the polygon is exactly `S² ∖ K`. By the setting, two of these
//! disks meet only at the finitely many points where their rings touch. `K`
//! is compact and locally contractible (a finite CW complex), so Alexander
//! duality (A. Hatcher, *Algebraic Topology*, Cambridge University Press,
//! 2002, Theorem 3.44) gives `H̃₀(S² ∖ K) ≅ H̃¹(K)`: the number of components
//! of the interior, minus one, is the rank of `H¹(K)`. Replace each disk by a
//! star joining an interior point to that disk's touch points. A closed
//! disk and a star are both contractible complexes containing the same
//! touch points, so the replacement is a homotopy equivalence relative to
//! those points, and gluing preserves it. `K` is therefore homotopy
//! equivalent to the touch graph (with each star's center standing for its
//! ring), and the rank of `H¹` of a graph is its number of independent
//! cycles. ∎
//!
//! **Why points, not pairs:** A graph with one edge per pair of touching
//! rings is wrong where three or more rings meet at one point: the exterior
//! and two holes meeting at one point make a triangle of pairs, yet the
//! interior is connected (the three disks share one point, a contractible
//! union). With point nodes the same configuration is a star, with no cycle.
//! Two distinct touch points between the same two rings are a cycle of
//! length four (ring, point, ring, point), and that does disconnect.
//!
//! **Exactness:** Every touch point is a vertex of one of the two rings (see
//! `ring_touch_points`), an input coordinate, found by exact point-on-segment
//! tests and compared by value. The cycle test is union–find on integers. No
//! coordinate is constructed and no tolerance is involved.

use crate::exact::kernel::{
    RingPlacement, contact_witness, place_ring, ring_touch_points, same_point_set,
};
use crate::types::{Coord2d, LinearRing};

use super::ValidityIssue;

/// A ring of the material, prepared for placement: the ring with repeats
/// removed, and the input index of each of its vertices.
pub(super) struct Prepared {
    pub(super) ring: LinearRing,
    pub(super) origin: Vec<usize>,
}

impl Prepared {
    /// An edge pair of the prepared rings, as input edge indices.
    fn edges(&self, other: &Prepared, (e, f): (usize, usize)) -> (usize, usize) {
        (self.origin[e], other.origin[f])
    }
}

/// Appends the placement issues of `rings` (index 0 the exterior; `None` for
/// a ring whose placement is undefined: non-finite, too few points, an area
/// that is not a number, or degenerate), in `Polygon::new`'s order: every hole against the exterior,
/// then every pair of holes.
pub(super) fn placement_into(rings: &[Option<Prepared>], issues: &mut Vec<ValidityIssue>) {
    let Some((first, holes)) = rings.split_first() else {
        return;
    };
    if let Some(ext) = first {
        for (k, h) in holes.iter().enumerate() {
            let ring = k + 1;
            if let Some(h) = h
                && let Some(issue) = against_exterior(ring, h, ext)
            {
                issues.push(issue);
            }
        }
    }
    for (i, a) in rings.iter().enumerate().skip(1) {
        let Some(a) = a else { continue };
        for (j, b) in rings.iter().enumerate().skip(i + 1) {
            let Some(b) = b else { continue };
            if let Some(issue) = between_holes((i, a), (j, b)) {
                issues.push(issue);
            }
        }
    }
}

/// The issue of hole `ring` against the exterior, if it is not inside.
fn against_exterior(ring: usize, h: &Prepared, ext: &Prepared) -> Option<ValidityIssue> {
    match place_ring(&h.ring, &ext.ring) {
        RingPlacement::Inside => None,
        RingPlacement::Outside => Some(ValidityIssue::HoleOutsideExterior { ring }),
        RingPlacement::SharedEdge if same_point_set(&h.ring, &ext.ring) => {
            Some(ValidityIssue::DuplicateRings { ring, other: 0 })
        }
        // `Crossing`, a partial `SharedEdge`, or a placement added later.
        _ => Some(crossing(ring, h, 0, ext)),
    }
}

/// The issue between holes `i < j`, if their interiors are not disjoint.
fn between_holes((i, a): (usize, &Prepared), (j, b): (usize, &Prepared)) -> Option<ValidityIssue> {
    let ab = place_ring(&a.ring, &b.ring);
    let ba = place_ring(&b.ring, &a.ring);
    let shared = ab == RingPlacement::SharedEdge || ba == RingPlacement::SharedEdge;
    let crossing_pair = ab == RingPlacement::Crossing || ba == RingPlacement::Crossing;
    if crossing_pair {
        return Some(crossing(j, b, i, a));
    }
    if shared {
        return Some(if same_point_set(&a.ring, &b.ring) {
            ValidityIssue::DuplicateRings { ring: j, other: i }
        } else {
            crossing(j, b, i, a)
        });
    }
    match (ab, ba) {
        (RingPlacement::Outside, RingPlacement::Outside) => None,
        (RingPlacement::Inside, _) => Some(ValidityIssue::NestedHoles { ring: i, inside: j }),
        (_, RingPlacement::Inside) => Some(ValidityIssue::NestedHoles { ring: j, inside: i }),
        // No other combination of two simple rings meeting at points; a
        // placement added later is reported as a crossing rather than
        // passed.
        _ => Some(crossing(j, b, i, a)),
    }
}

/// `RingsCross` for rings `ring` and `other`, with a witnessing edge pair.
fn crossing(ring: usize, r: &Prepared, other: usize, o: &Prepared) -> ValidityIssue {
    let edges = contact_witness(&r.ring, &o.ring)
        .or_else(|| contact_witness(&o.ring, &r.ring).map(|(f, e)| (e, f)))
        .map_or((r.origin[0], o.origin[0]), |pair| r.edges(o, pair));
    ValidityIssue::RingsCross { ring, other, edges }
}

/// The touch whose incidence first closes a cycle of the touch graph, as
/// `(ring, vertex)`: the point is vertex `vertex` of ring `ring` (the lowest
/// ring having it as a vertex), by index into that ring's stored vertices.
/// `None` if the graph is a forest, that is, if the interior is connected.
///
/// The rings must satisfy the module's setting (simple, placed as
/// `Polygon::new` requires); the result is otherwise meaningless, and
/// callers do not ask.
pub(crate) fn disconnected_interior(rings: &[&LinearRing]) -> Option<(usize, usize)> {
    let n = rings.len();
    let boxes: Vec<Option<(f64, f64, f64, f64)>> = rings.iter().map(|r| ring_box(r)).collect();
    // Nodes 0..n are rings; n.. are touch points, in order of discovery.
    let mut points: Vec<Coord2d> = Vec::new();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut incidences: Vec<(usize, usize)> = Vec::new();
    for r in 0..n {
        for s in (r + 1)..n {
            if !boxes_meet(boxes[r], boxes[s]) {
                continue;
            }
            for p in ring_touch_points(rings[r], rings[s]) {
                let id = if let Some(k) = points.iter().position(|&q| q == p) {
                    n + k
                } else {
                    points.push(p);
                    parent.push(parent.len());
                    n + points.len() - 1
                };
                for ring in [r, s] {
                    if incidences.contains(&(ring, id)) {
                        continue;
                    }
                    incidences.push((ring, id));
                    let (x, y) = (find(&mut parent, ring), find(&mut parent, id));
                    if x == y {
                        return vertex_of(rings, points[id - n]);
                    }
                    parent[x] = y;
                }
            }
        }
    }
    None
}

/// The lowest ring having `p` as a vertex, and the vertex's index.
fn vertex_of(rings: &[&LinearRing], p: Coord2d) -> Option<(usize, usize)> {
    rings.iter().enumerate().find_map(|(r, ring)| {
        ring.coords()
            .as_slice()
            .iter()
            .position(|&v| v == p)
            .map(|k| (r, k))
    })
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

fn ring_box(r: &LinearRing) -> Option<(f64, f64, f64, f64)> {
    r.envelope()
        .ok()
        .map(|e| (e.min().x, e.min().y, e.max().x, e.max().y))
}

/// Closed boxes meet (compared exactly). A missing box is treated as
/// meeting everything, so the screen never skips a pair it cannot judge.
fn boxes_meet(a: Option<(f64, f64, f64, f64)>, b: Option<(f64, f64, f64, f64)>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a.0 <= b.2 && b.0 <= a.2 && a.1 <= b.3 && b.1 <= a.3,
        _ => true,
    }
}
