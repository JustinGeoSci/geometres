// src/ops/coverage_simplify.rs

//! Coverage simplification: Visvalingam–Whyatt simplification of a
//! polygonal coverage that keeps it a coverage.
//!
//! # The promise
//!
//! For a valid coverage `polys` (one that [`coverage_validate`] accepts),
//! `coverage_simplify(polys, threshold)` returns one polygon per input
//! polygon, in input order, such that:
//!
//! 1. **Every output vertex is an input vertex,** bit for bit. Vertices are
//!    only removed.
//! 2. **Every shared edge is simplified once,** and both polygons that share
//!    it use the same result, so neighbors still match vertex for vertex.
//! 3. **No edge crosses, touches or jumps over another** that it did not
//!    before: every removal is decided by the exact checks of
//!    [`simplify_preserve_topology`](super::simplify_preserve_topology).
//! 4. **Every node is kept,** and no polygon disappears: every ring keeps at
//!    least three vertices.
//! 5. **The result is a valid coverage,** checked exactly by
//!    [`coverage_validate`]; a result that fails is refused, never returned.
//!
//! The method follows JTS's `CoverageSimplifier`
//! (`org.locationtech.jts.coverage`, Martin Davis), which "simplifies the
//! boundaries of the polygons in a polygonal coverage while preserving the
//! original coverage topology" by a topology-preserving Visvalingam–Whyatt
//! on the coverage's edges, keeping node points unchanged.
//!
//! # The edge graph
//!
//! **Nodes:** A ring vertex is a node when three or more distinct segments
//! of the coverage meet at its coordinate, or when the ring passes there
//! from one sharing to another: its incoming and outgoing segments are not
//! used by the same rings (from a shared stretch to the coverage's outer
//! boundary, or from one neighbor to another). JTS's statement of the same
//! set: "inner vertices shared by three or more polygons, or boundary
//! vertices shared by two or more". Coordinates and segments are compared
//! by exact equality (`-0.0` equals `0.0`).
//!
//! **Edges:** Every ring is cut at its nodes into pieces, node to node. Two
//! pieces are the same edge when their vertex sequences are equal by exact
//! coordinate equality, in either direction; each edge is stored in the
//! direction whose sequence is lexicographically smaller. A ring with no
//! node is one closed edge starting at its lexicographically smallest
//! vertex. That starting vertex is treated as a node (never removed), so a
//! closed edge is identified the same way from either ring that uses it
//! (an island and the hole it fills). In a valid coverage an edge is used
//! by one ring (it lies on the coverage's outer boundary) or by two, in
//! opposite directions.
//!
//! # The simplification
//!
//! Edges are processed in a fixed order, by their stored vertex sequence
//! (so first by their first vertex), each once, with the elimination core
//! of [`simplify_vw`](super::simplify_vw()):
//!
//! - an open edge is a line, its two end nodes fixed; a closed edge is a
//!   ring of at least three vertices, its starting node fixed;
//! - the effective areas are the same `f64` selection rule, with the same
//!   threshold semantics (below the threshold, or exactly 0, is removed);
//! - **the veto:** before a vertex `v` between its current neighbors `a`
//!   and `c` is removed, the run `a, v, c` and its chord `a → c` are put to
//!   the four exact checks of the topology-preserving simplifier against
//!   every other current segment of the whole coverage: (a) the chord meets
//!   no segment other than at `a` or `c`, and overlaps none; (b) no endpoint
//!   of another segment lies inside the triangle `a, v, c`; (c) no segment
//!   meets the run other than at `a` or `c`; (d) the run does not touch
//!   itself. Candidates come from a dynamic R\*-tree of current segments
//!   (the crate's private `situs` seam), which filters with no false
//!   negatives; every decision is then exact (`orient2d` and comparisons of
//!   stored coordinates). A rejected vertex is kept, and the next one is
//!   tried. An accepted removal replaces the run's two segments by the
//!   chord in the index.
//!
//! # Why the result is a coverage
//!
//! - **Shared edges match.** Each shared edge is one chain, simplified once;
//!   both polygons rebuild their rings from its kept vertices. So no
//!   `UnmatchedEdge` can appear: every segment of the result is used by
//!   exactly the rings that used the edge it lies on.
//! - **Nothing crosses or jumps.** One accepted removal replaces the run
//!   `a, v, c` by its chord and changes nothing else. The topology-preserving
//!   simplifier's argument (its module documentation) shows that a
//!   replacement passing (a)–(d) changes no contact between the edge and any
//!   other segment, and moves no other segment from one side of the edge to
//!   the other: the region the edge sweeps is the closed triangle, which
//!   holds nothing else but `a` and `c`. Applied removal by removal, every
//!   contact and side relation of the input survives. So no ring loses
//!   simplicity, no `Overlap` appears, and no hole leaves its exterior.
//! - **Nothing collapses.** A ring cannot be left with two vertices (open
//!   edges keep their nodes, closed edges keep three), and no ring can lose
//!   all its area: that needs a chord lying along a remaining segment, which
//!   (a) refuses. A ring is degenerate exactly when its exact signed area
//!   is zero, so a ring made very thin, but not flat,
//!   is a valid ring.
//!
//! **The vouch.** Every polygon is built by [`Polygon::new`], and the result
//! is checked by [`coverage_validate`]. By the argument both always pass; a
//! failure would be a defect in this crate, and is refused (with
//! [`UnvouchedReason::NotACoverage`], or the constructor's error), never
//! returned.
//!
//! # The threshold
//!
//! The threshold is an **area**, with [`simplify_vw`](super::simplify_vw())'s
//! units and semantics. JTS's `CoverageSimplifier` takes a distance tolerance
//! and documents it as "the square root of the area tolerance used in the
//! Visvalingam-Whyatt algorithm". So a JTS tolerance `t` corresponds to a
//! threshold of `t * t` here.
//!
//! # Invalid input
//!
//! A set of polygons that is not a valid coverage is refused with
//! [`GeometryError::InvalidCoverage`] carrying the first issue, as
//! [`coverage_union`](super::coverage_union) refuses it: the shared-edge
//! argument has no meaning without matching edges.

use crate::GeometryError;
use crate::index::SegmentIndex;
use crate::ops::coverage::graph::{Key, coord_exact_cmp, key, seq_cmp};
use crate::ops::simplify::check_tolerance;
use crate::ops::simplify_topology::{Obstacle, Run, clear_of, run_touches_itself};
use crate::ops::simplify_vw::{ChainKind, eliminate};
use crate::ops::{CoverageIssue, coverage_validate};
use crate::types::{
    Coord2d, CoordSeq, LinearRing, MaybeVouched, Polygon, Unvouched, UnvouchedReason,
};
use std::collections::{BTreeMap, BTreeSet};

/// "No next vertex": the end of an open edge.
const END: u32 = u32::MAX;

// ===========================================================================
// The public entry points
// ===========================================================================

/// Visvalingam–Whyatt simplification of a polygonal coverage that keeps it a
/// coverage: every shared edge simplified once, no edge crossing or jumping
/// over another, every node kept. See the [module documentation](self) for
/// the edge graph, the exact checks and the proof.
///
/// Returns one polygon per input polygon, in input order; every output
/// vertex is an input vertex. `threshold` is an area, as for
/// [`simplify_vw`](super::simplify_vw()).
/// ```
/// use geometres::ops::coverage_simplify;
/// use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon};
///
/// let poly = |v: &[(f64, f64)]| {
///     Polygon::new_simple(LinearRing::new(CoordSeq::new(
///         v.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
///     ))?)
/// };
/// // Two squares sharing a slightly bent edge through (1, 0.5).
/// let left = poly(&[(0.0, 0.0), (1.0, 0.0), (1.01, 0.5), (1.0, 1.0), (0.0, 1.0)])?;
/// let right = poly(&[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.01, 0.5)])?;
/// let out = coverage_simplify(&[left, right], 0.01)?;
/// assert_eq!(out[0].exterior().num_vertices(), 4);
/// assert_eq!(out[1].exterior().num_vertices(), 4);
/// # Ok::<(), geometres::GeometryError>(())
/// ```
///
/// # Errors
///
/// - [`GeometryError::InvalidTolerance`] for a negative or non-finite
///   threshold.
/// - [`GeometryError::InvalidCoverage`] with the first issue
///   [`coverage_validate`] finds, if `polys` is not a valid coverage.
/// - [`GeometryError::TooManyEntries`] for more than `u32::MAX − 1` ring
///   vertices.
/// - [`GeometryError::Unvouched`] with [`UnvouchedReason::NotACoverage`],
///   or a constructor's error, if the closing checks fail: not expected
///   (module documentation).
pub fn coverage_simplify(polys: &[Polygon], threshold: f64) -> Result<Vec<Polygon>, GeometryError> {
    simplify(polys, threshold, Mode::All)
}

/// [`coverage_simplify`], simplifying only the **shared** edges: every edge
/// on the coverage's outer boundary (used by one polygon) is left unchanged,
/// byte for byte, so the result fits exactly inside the original coverage's
/// outline. JTS's `simplifyInner`.
///
/// # Errors
///
/// As [`coverage_simplify`].
pub fn coverage_simplify_inner(
    polys: &[Polygon],
    threshold: f64,
) -> Result<Vec<Polygon>, GeometryError> {
    simplify(polys, threshold, Mode::Inner)
}

/// Which edges are simplified.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Mode {
    /// Every edge.
    All,
    /// Only edges used by two rings.
    Inner,
}

fn simplify(polys: &[Polygon], threshold: f64, mode: Mode) -> Result<Vec<Polygon>, GeometryError> {
    check_tolerance(threshold)?;
    if let Some(&first) = coverage_validate(polys).first() {
        return Err(GeometryError::InvalidCoverage { first });
    }
    let graph = EdgeGraph::build(polys)?;
    let mut state = State::new(&graph);
    for &e in &graph.order {
        if mode == Mode::Inner && graph.edges[e].users < 2 {
            continue;
        }
        state.simplify_edge(&graph, e, threshold);
    }
    vouch(rebuild(polys, &graph, &state)?)
}

/// The closing check: `out` if it is a valid coverage, refused with its
/// first issue otherwise. By the module documentation's argument the
/// refusal is unreachable for a valid input; it is a defense.
pub(super) fn vouch(out: Vec<Polygon>) -> Result<Vec<Polygon>, GeometryError> {
    let first: Option<CoverageIssue> = coverage_validate(&out).first().copied();
    match first {
        None => Ok(out),
        Some(first) => {
            MaybeVouched::Unvouched(Unvouched::new(out, UnvouchedReason::NotACoverage { first }))
                .vouched()
        }
    }
}

// ===========================================================================
// The edge graph
// ===========================================================================

/// One edge: a node-to-node chain, stored in its canonical direction.
#[derive(Debug)]
pub(super) struct Edge {
    /// The vertices, canonical direction. A closed edge's first vertex is
    /// repeated at the end.
    pub(super) verts: Vec<Coord2d>,
    /// How many rings use it: 1 (outer boundary) or 2 (shared).
    pub(super) users: usize,
}

impl Edge {
    pub(super) fn closed(&self) -> bool {
        self.verts[0] == self.verts[self.verts.len() - 1]
    }
    /// The vertices the elimination runs on: a closed edge without its
    /// repeated last vertex.
    fn chain(&self) -> &[Coord2d] {
        if self.closed() {
            &self.verts[..self.verts.len() - 1]
        } else {
            &self.verts
        }
    }
}

/// Where a ring vertex lives: its edge, and its position in the edge's
/// chain.
type Slot = (usize, usize);

/// The edges of a coverage, and each ring vertex's place among them.
#[derive(Debug)]
pub(super) struct EdgeGraph {
    pub(super) edges: Vec<Edge>,
    /// Edge ids in processing order: by stored vertex sequence.
    pub(super) order: Vec<usize>,
    /// `slots[ring][k]`: where vertex `k` of ring `ring` lives. Rings are
    /// numbered polygon by polygon, exterior then holes.
    pub(super) slots: Vec<Vec<Slot>>,
}

/// Every ring of the coverage, polygon by polygon, exterior then holes.
fn rings_of(polys: &[Polygon]) -> Vec<&[Coord2d]> {
    polys
        .iter()
        .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes()))
        .map(|r| r.coords().as_slice())
        .collect()
}

/// An undirected segment's identity.
fn seg_key(a: Coord2d, b: Coord2d) -> (Key, Key) {
    let (ka, kb) = (key(a), key(b));
    if ka <= kb { (ka, kb) } else { (kb, ka) }
}

impl EdgeGraph {
    pub(super) fn build(polys: &[Polygon]) -> Result<Self, GeometryError> {
        let rings = rings_of(polys);
        let total: usize = rings.iter().map(|r| r.len()).sum();
        if u32::try_from(total).map_or(true, |n| n == u32::MAX) {
            return Err(GeometryError::TooManyEntries { count: total });
        }
        let nodes = node_flags(&rings);
        let mut by_seq: BTreeMap<Vec<Key>, usize> = BTreeMap::new();
        let mut edges: Vec<Edge> = Vec::new();
        let mut slots: Vec<Vec<Slot>> = Vec::with_capacity(rings.len());
        for (r, ring) in rings.iter().enumerate() {
            let n = ring.len();
            let mut cuts: Vec<usize> = (0..n).filter(|&k| nodes[r][k]).collect();
            if cuts.is_empty() {
                // No node: one closed edge from the smallest vertex.
                let start = (0..n)
                    .min_by(|&i, &j| coord_exact_cmp(ring[i], ring[j]).then(i.cmp(&j)))
                    .unwrap_or(0);
                cuts.push(start);
            }
            let mut ring_slots = vec![(0, 0); n];
            for (j, &from) in cuts.iter().enumerate() {
                let to = if j + 1 < cuts.len() {
                    cuts[j + 1]
                } else {
                    cuts[0] + n
                };
                let piece: Vec<Coord2d> = (from..=to).map(|k| ring[k % n]).collect();
                let rev: Vec<Coord2d> = piece.iter().rev().copied().collect();
                let forward = seq_cmp(&piece, &rev) != core::cmp::Ordering::Greater;
                let canon = if forward { piece } else { rev };
                let ckey: Vec<Key> = canon.iter().map(|&c| key(c)).collect();
                let len = canon.len();
                let id = *by_seq.entry(ckey).or_insert_with(|| {
                    edges.push(Edge {
                        verts: canon,
                        users: 0,
                    });
                    edges.len() - 1
                });
                edges[id].users += 1;
                // The chain of a closed edge drops its repeated last
                // vertex, which is position 0 again.
                let chain_len = if edges[id].closed() { len - 1 } else { len };
                for t in 0..len - 1 {
                    let pos = if forward { t } else { len - 1 - t };
                    ring_slots[(from + t) % n] = (id, pos % chain_len);
                }
            }
            slots.push(ring_slots);
        }
        let mut order: Vec<usize> = (0..edges.len()).collect();
        order.sort_by(|&a, &b| seq_cmp(&edges[a].verts, &edges[b].verts).then(a.cmp(&b)));
        Ok(Self {
            edges,
            order,
            slots,
        })
    }
}

/// `nodes[ring][k]`: vertex `k` of `ring` is a node (module documentation):
/// three or more distinct segments meet at its coordinate, or its incoming
/// and outgoing segments are used by different rings.
fn node_flags(rings: &[&[Coord2d]]) -> Vec<Vec<bool>> {
    let mut users: BTreeMap<(Key, Key), Vec<usize>> = BTreeMap::new();
    let mut incident: BTreeMap<Key, BTreeSet<(Key, Key)>> = BTreeMap::new();
    for (r, ring) in rings.iter().enumerate() {
        let n = ring.len();
        for k in 0..n {
            let (a, b) = (ring[k], ring[(k + 1) % n]);
            let s = seg_key(a, b);
            users.entry(s).or_default().push(r);
            incident.entry(key(a)).or_default().insert(s);
            incident.entry(key(b)).or_default().insert(s);
        }
    }
    rings
        .iter()
        .map(|ring| {
            let n = ring.len();
            (0..n)
                .map(|k| {
                    let v = ring[k];
                    let (p, q) = (ring[(k + n - 1) % n], ring[(k + 1) % n]);
                    let degree = incident.get(&key(v)).map_or(0, BTreeSet::len);
                    degree >= 3 || users.get(&seg_key(p, v)) != users.get(&seg_key(v, q))
                })
                .collect()
        })
        .collect()
}

// ===========================================================================
// The simplification state
// ===========================================================================

/// The current segments of every edge, indexed, and which vertices are
/// alive. Vertices are numbered globally, edge after edge (`gid`); a
/// segment is named by the gid of its start vertex and ends at that
/// vertex's `next`.
pub(super) struct State {
    verts: Vec<Coord2d>,
    /// The first gid of each edge's chain.
    pub(super) base: Vec<usize>,
    next: Vec<u32>,
    pub(super) alive: Vec<bool>,
    index: SegmentIndex,
}

impl State {
    pub(super) fn new(graph: &EdgeGraph) -> Self {
        let mut verts = Vec::new();
        let mut base = Vec::with_capacity(graph.edges.len());
        let mut next = Vec::new();
        for e in &graph.edges {
            let chain = e.chain();
            let start = verts.len();
            base.push(start);
            verts.extend_from_slice(chain);
            for k in 0..chain.len() {
                let succ = if k + 1 < chain.len() {
                    Some(start + k + 1)
                } else if e.closed() {
                    Some(start)
                } else {
                    None
                };
                next.push(succ.map_or(END, gid));
            }
        }
        let mut index = SegmentIndex::new();
        for (i, &n) in next.iter().enumerate() {
            if n != END {
                index.insert(verts[i], verts[n as usize], gid(i));
            }
        }
        crate::count!(
            coverage_simplify_index_updates,
            next.iter().filter(|&&n| n != END).count()
        );
        let alive = vec![true; verts.len()];
        Self {
            verts,
            base,
            next,
            alive,
            index,
        }
    }

    /// Simplifies edge `e` with the elimination core and the exact veto.
    pub(super) fn simplify_edge(&mut self, graph: &EdgeGraph, e: usize, threshold: f64) {
        let edge = &graph.edges[e];
        let chain = edge.chain();
        let base = self.base[e];
        let (kind, closed) = if edge.closed() {
            (ChainKind::Ring, true)
        } else {
            (ChainKind::Line, false)
        };
        eliminate(chain, kind, threshold, &mut |i, b, a| {
            // A closed edge's starting vertex is a node.
            if closed && i == 0 {
                return false;
            }
            self.try_remove(base + b, base + i, base + a)
        });
    }

    /// The veto: the run `b, i, a` passes the four exact checks against
    /// every other current segment, and if so the removal is committed to
    /// the index.
    fn try_remove(&mut self, b: usize, i: usize, a: usize) -> bool {
        crate::count!(coverage_simplify_veto_checks);
        let pts = vec![self.verts[b], self.verts[i], self.verts[a]];
        let (mut lo, mut hi) = (pts[0], pts[0]);
        for p in &pts {
            lo = Coord2d::new(lo.x.min(p.x), lo.y.min(p.y));
            hi = Coord2d::new(hi.x.max(p.x), hi.y.max(p.y));
        }
        let run = Run {
            ids: vec![gid(b), gid(i), gid(a)],
            pts,
            lo,
            hi,
        };
        if !self.accept(&run) {
            crate::count!(coverage_simplify_veto_rejections);
            return false;
        }
        let removed_first = self.index.remove(run.pts[0], run.pts[1], run.ids[0]);
        let removed_second = self.index.remove(run.pts[1], run.pts[2], run.ids[1]);
        debug_assert!(
            removed_first && removed_second,
            "current segments are indexed"
        );
        self.index.insert(run.pts[0], run.pts[2], run.ids[0]);
        crate::count!(coverage_simplify_index_updates, 3);
        self.next[b] = gid(a);
        self.alive[i] = false;
        true
    }

    /// The acceptance test of the topology-preserving simplifier on a
    /// three-vertex run: (d), then (a)–(c) against every candidate segment
    /// other than the run's own two.
    fn accept(&self, run: &Run) -> bool {
        if run_touches_itself(run) {
            return false;
        }
        crate::count!(coverage_simplify_index_queries);
        let own = [run.ids[0], run.ids[1]];
        let mut ok = true;
        self.index.candidates(run.lo, run.hi, |id| {
            if ok && !own.contains(&id) {
                let i = id as usize;
                let seg = Obstacle::Segment(self.verts[i], self.verts[self.next[i] as usize]);
                ok = clear_of(run, seg);
            }
        });
        ok
    }
}

/// A gid as `u32`; every gid is below the count checked in
/// [`EdgeGraph::build`].
fn gid(i: usize) -> u32 {
    u32::try_from(i).unwrap_or(END)
}

// ===========================================================================
// Rebuilding
// ===========================================================================

/// Each polygon's rings, in their original order and rotation, keeping the
/// vertices still alive in their edges; built through the constructors.
pub(super) fn rebuild(
    polys: &[Polygon],
    graph: &EdgeGraph,
    state: &State,
) -> Result<Vec<Polygon>, GeometryError> {
    let rings = rings_of(polys);
    let mut r = 0;
    let mut out = Vec::with_capacity(polys.len());
    for p in polys {
        let mut built: Vec<LinearRing> = Vec::with_capacity(1 + p.holes().len());
        for _ in 0..=p.holes().len() {
            let kept: Vec<Coord2d> = rings[r]
                .iter()
                .zip(&graph.slots[r])
                .filter(|&(_, &(e, pos))| state.alive[state.base[e] + pos])
                .map(|(&c, _)| c)
                .collect();
            built.push(LinearRing::new(CoordSeq::new(kept))?);
            r += 1;
        }
        let exterior = built.remove(0);
        out.push(Polygon::new(exterior, built)?);
    }
    Ok(out)
}
