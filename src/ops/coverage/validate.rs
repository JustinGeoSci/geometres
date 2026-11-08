// src/ops/coverage/validate.rs

//! The coverage issues: overlaps by exact `relate`, unmatched edges by exact
//! edge identity and `orient2d`, and narrow gaps by `ops::distance`'s exact
//! segment test. Candidate pairs come from box indexes through the private
//! seam, which never loses a pair whose boxes meet.

use std::collections::BTreeSet;

use arithmetica::Sign;

use crate::exact::kernel::{first_self_contact, on_segment, orient};
use crate::index::BoxIndex;
use crate::ops::distance::{Element, grow_outward, pair_within};
use crate::predicates::relate;
use crate::topology::label::Dimension;
use crate::types::{Coord2d, Envelope, Polygon};

use super::CoverageIssue;
use super::graph::{coord_exact_cmp, edge_key};

/// One ring edge of the input, with where it came from. (`edge` mirrors
/// the issue fields it fills.)
#[derive(Clone, Copy)]
#[allow(clippy::struct_field_names)]
pub(super) struct Edge {
    pub(super) polygon: usize,
    pub(super) ring: usize,
    pub(super) edge: usize,
    pub(super) a: Coord2d,
    pub(super) b: Coord2d,
}

/// Every ring edge of every polygon, in `(polygon, ring, edge)` order.
pub(super) fn edges_of(polys: &[Polygon]) -> Vec<Edge> {
    let mut out = Vec::new();
    for (polygon, p) in polys.iter().enumerate() {
        for (ring, r) in core::iter::once(p.exterior()).chain(p.holes()).enumerate() {
            for (edge, (a, b)) in r.edges().enumerate() {
                out.push(Edge {
                    polygon,
                    ring,
                    edge,
                    a,
                    b,
                });
            }
        }
    }
    out
}

/// A box index over the edges, ids their positions. Edge counts beyond
/// `u32::MAX` are not indexable; such an input is beyond any coverage this
/// crate can hold in memory, and the ids saturate rather than wrap.
fn edge_index(edges: &[Edge]) -> BoxIndex {
    let entries: Vec<((Coord2d, Coord2d), u32)> = edges
        .iter()
        .enumerate()
        .map(|(i, e)| ((e.a, e.b), u32::try_from(i).unwrap_or(u32::MAX)))
        .collect();
    BoxIndex::from_edges(&entries)
}

/// Every issue, sorted; with `gap = Some(width)`, narrow gaps too.
pub(super) fn issues(polys: &[Polygon], gap: Option<f64>) -> Vec<CoverageIssue> {
    let mut found: BTreeSet<CoverageIssue> = BTreeSet::new();
    for (polygon, p) in polys.iter().enumerate() {
        for (ring, r) in core::iter::once(p.exterior()).chain(p.holes()).enumerate() {
            if first_self_contact(r.coords().as_slice()).is_some() {
                found.insert(CoverageIssue::NonSimpleRing { polygon, ring });
            }
        }
    }
    overlaps(polys, &mut found);
    let edges = edges_of(polys);
    let index = edge_index(&edges);
    index.self_join(|i, j| {
        let (s, t) = (edges[i as usize], edges[j as usize]);
        if s.polygon != t.polygon {
            unmatched(s, t, &mut found);
        }
    });
    if let Some(width) = gap.filter(|&w| w > 0.0) {
        narrow_gaps(&edges, &index, width, &mut found);
    }
    found.into_iter().collect()
}

// ===========================================================================
// Overlaps
// ===========================================================================

/// `Overlap` for every pair of polygons whose interiors meet. Candidates
/// are pairs whose exterior envelopes meet; each is decided by the exact
/// matrix (`II ≠ F`).
fn overlaps(polys: &[Polygon], found: &mut BTreeSet<CoverageIssue>) {
    let boxes: Vec<(Envelope, u32)> = polys
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            let e = p.exterior().envelope().ok()?;
            Some((e, u32::try_from(i).unwrap_or(u32::MAX)))
        })
        .collect();
    BoxIndex::from_boxes(&boxes).self_join(|i, j| {
        let (a, b) = (i.min(j) as usize, i.max(j) as usize);
        let interiors_meet =
            relate(&polys[a], &polys[b]).is_ok_and(|r| r.matrix().cells()[0] != Dimension::Empty);
        if interiors_meet {
            found.insert(CoverageIssue::Overlap { a, b });
        }
    });
}

// ===========================================================================
// Unmatched edges
// ===========================================================================

fn unmatched_issue(e: Edge, other: usize) -> CoverageIssue {
    CoverageIssue::UnmatchedEdge {
        polygon: e.polygon,
        ring: e.ring,
        edge: e.edge,
        other,
    }
}

/// The unmatched-edge issues of one pair of edges from different polygons.
///
/// - **Identical edges** (either direction) are matched, with nothing to report.
/// - **A collinear overlap of positive length** between different edges leaves
///   both unmatched. Collinear is `orient2d = 0` for both endpoints of
///   one against the other; along a line the lexicographic order of
///   coordinates is monotone, so the overlap has positive length exactly
///   when the larger of the two lower ends is strictly below the smaller of
///   the two upper ends. Both tests compare stored values.
/// - **A vertex of one edge strictly inside the other** leaves the touched edge
///   unmatched (the coverage is not noded at that vertex).
fn unmatched(s: Edge, t: Edge, found: &mut BTreeSet<CoverageIssue>) {
    if edge_key(s.a, s.b) == edge_key(t.a, t.b) {
        return;
    }
    if collinear_overlap(s, t) {
        found.insert(unmatched_issue(s, t.polygon));
        found.insert(unmatched_issue(t, s.polygon));
    }
    for p in [t.a, t.b] {
        if strictly_inside(p, s) {
            found.insert(unmatched_issue(s, t.polygon));
        }
    }
    for p in [s.a, s.b] {
        if strictly_inside(p, t) {
            found.insert(unmatched_issue(t, s.polygon));
        }
    }
}

/// `p` lies on the closed edge and is neither endpoint, decided exactly.
fn strictly_inside(p: Coord2d, e: Edge) -> bool {
    p != e.a && p != e.b && on_segment(p, e.a, e.b)
}

/// The two edges are collinear and overlap in a stretch of positive length.
pub(super) fn collinear_overlap(s: Edge, t: Edge) -> bool {
    if orient(s.a, s.b, t.a) != Sign::Zero || orient(s.a, s.b, t.b) != Sign::Zero {
        return false;
    }
    let (s_lo, s_hi) = ordered(s.a, s.b);
    let (t_lo, t_hi) = ordered(t.a, t.b);
    let lo = if coord_exact_cmp(s_lo, t_lo).is_ge() {
        s_lo
    } else {
        t_lo
    };
    let hi = if coord_exact_cmp(s_hi, t_hi).is_le() {
        s_hi
    } else {
        t_hi
    };
    coord_exact_cmp(lo, hi).is_lt()
}

fn ordered(a: Coord2d, b: Coord2d) -> (Coord2d, Coord2d) {
    if coord_exact_cmp(a, b).is_le() {
        (a, b)
    } else {
        (b, a)
    }
}

// ===========================================================================
// Narrow gaps
// ===========================================================================

/// The candidate window of a gap test: the edge's box grown by `width`,
/// outward-rounded, so it contains the exact grown box (`ops::distance`'s
/// helper and its proof). A pair within `width` has boxes within `width` on each axis,
/// so the other edge's box meets this window.
pub(super) fn gap_window(a: Coord2d, b: Coord2d, width: f64) -> Envelope {
    let lo = Coord2d::new(a.x.min(b.x), a.y.min(b.y));
    let hi = Coord2d::new(a.x.max(b.x), a.y.max(b.y));
    grow_outward((lo, hi), width)
}

/// `NarrowGap` for every pair of edges of different polygons that are not
/// the same edge, do not touch, and lie within `width`, decided exactly.
fn narrow_gaps(edges: &[Edge], index: &BoxIndex, width: f64, found: &mut BTreeSet<CoverageIssue>) {
    for (i, s) in edges.iter().enumerate() {
        let window = gap_window(s.a, s.b, width);
        index.candidates(&window, |j| {
            let j = j as usize;
            if j <= i {
                return;
            }
            let t = edges[j];
            if s.polygon == t.polygon || edge_key(s.a, s.b) == edge_key(t.a, t.b) || touch(*s, t) {
                return;
            }
            if pair_within(
                Element::Segment(s.a, s.b),
                Element::Segment(t.a, t.b),
                width,
            ) {
                found.insert(CoverageIssue::NarrowGap {
                    polygon: s.polygon,
                    ring: s.ring,
                    edge: s.edge,
                    other: t.polygon,
                    other_ring: t.ring,
                    other_edge: t.edge,
                });
            }
        });
    }
}

/// The two closed edges share a point, decided exactly: a proper crossing,
/// or an endpoint of one on the other.
pub(super) fn touch(s: Edge, t: Edge) -> bool {
    let o1 = orient(s.a, s.b, t.a);
    let o2 = orient(s.a, s.b, t.b);
    let o3 = orient(t.a, t.b, s.a);
    let o4 = orient(t.a, t.b, s.b);
    let opposite = |x: Sign, y: Sign| {
        matches!(
            (x, y),
            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
        )
    };
    (opposite(o1, o2) && opposite(o3, o4))
        || on_segment(t.a, s.a, s.b)
        || on_segment(t.b, s.a, s.b)
        || on_segment(s.a, t.a, t.b)
        || on_segment(s.b, t.a, t.b)
}
