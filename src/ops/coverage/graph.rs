// src/ops/coverage/graph.rs

//! The shared machinery of coverage union and polygonize: exact vertex and
//! edge identity, a half-edge graph traced with the face on the left, the
//! split of a boundary walk into simple loops, and hole assignment by exact
//! placement. Every coordinate here is an input coordinate.

use core::cmp::Ordering;

use crate::exact::{RingPlacement, place_ring};
use crate::types::{Coord2d, CoordSeq, LinearRing, Polygon};

use super::angular::ccw_cmp;

// ===========================================================================
// Identity
// ===========================================================================

/// A coordinate's identity: its two values' bit patterns with `-0.0`
/// canonicalized to `0.0` (`v + 0.0` does exactly that and changes nothing
/// else), so equal keys are exactly equal coordinates.
pub(in crate::ops) type Key = (u64, u64);

pub(in crate::ops) fn key(c: Coord2d) -> Key {
    ((c.x + 0.0).to_bits(), (c.y + 0.0).to_bits())
}

/// An undirected edge's identity: its two endpoint keys, smaller first.
pub(super) fn edge_key(a: Coord2d, b: Coord2d) -> (Key, Key) {
    let (ka, kb) = (key(a), key(b));
    if ka <= kb { (ka, kb) } else { (kb, ka) }
}

/// The lexicographic order of two coordinate sequences, element by element
/// with [`coord_exact_cmp`], then by length.
pub(in crate::ops) fn seq_cmp(a: &[Coord2d], b: &[Coord2d]) -> Ordering {
    for (x, y) in a.iter().zip(b) {
        let c = coord_exact_cmp(*x, *y);
        if c != Ordering::Equal {
            return c;
        }
    }
    a.len().cmp(&b.len())
}

/// The exact lexicographic order with `-0.0 == 0.0`: `partial_cmp` on finite
/// values is total and treats the two zeros as equal.
pub(in crate::ops) fn coord_exact_cmp(a: Coord2d, b: Coord2d) -> Ordering {
    a.x.partial_cmp(&b.x)
        .unwrap_or(Ordering::Equal)
        .then(a.y.partial_cmp(&b.y).unwrap_or(Ordering::Equal))
}

// ===========================================================================
// The half-edge graph
// ===========================================================================

/// An embedded graph over input coordinates: each undirected edge is two
/// half-edges, `2i` and `2i + 1`, twins of each other.
///
/// The half-edges leaving each vertex are kept in one flat array, vertex by
/// vertex (`out[out_start[v]..out_start[v + 1]]`), rather than one vector
/// per vertex: the same lists in the same order, without an allocation per
/// vertex.
pub(super) struct HalfEdges<'a> {
    coords: &'a [Coord2d],
    from: Vec<usize>,
    /// Where each vertex's list starts in `out`; one entry past the last
    /// vertex holds `out.len()`.
    out_start: Vec<usize>,
    /// The half-edges leaving each vertex, sorted counterclockwise from `+x`.
    out: Vec<usize>,
    /// Each half-edge's position in its origin's list.
    pos: Vec<usize>,
}

impl<'a> HalfEdges<'a> {
    /// The graph of `edges`, given as distinct undirected vertex-index
    /// pairs over `coords`. The edges must not overlap collinearly at a
    /// shared vertex (no two leave a vertex in the same direction), which
    /// every caller has established exactly.
    ///
    /// Each vertex's list is filled in increasing half-edge order and then
    /// sorted by a stable sort, so it is the list a stable sort of the
    /// half-edges pushed in that order gives.
    pub(super) fn new(coords: &'a [Coord2d], edges: &[(usize, usize)]) -> Self {
        let mut from = Vec::with_capacity(2 * edges.len());
        for &(a, b) in edges {
            from.push(a);
            from.push(b);
        }
        let mut out_start = vec![0_usize; coords.len() + 1];
        for &v in &from {
            out_start[v + 1] += 1;
        }
        for v in 0..coords.len() {
            out_start[v + 1] += out_start[v];
        }
        let mut fill = out_start.clone();
        let mut out = vec![0_usize; from.len()];
        for (h, &v) in from.iter().enumerate() {
            out[fill[v]] = h;
            fill[v] += 1;
        }
        let to = |h: usize| from[h ^ 1];
        for v in 0..coords.len() {
            out[out_start[v]..out_start[v + 1]]
                .sort_by(|&h, &k| ccw_cmp(coords[v], coords[to(h)], coords[to(k)]));
        }
        let mut pos = vec![0; from.len()];
        for v in 0..coords.len() {
            for (i, &h) in out[out_start[v]..out_start[v + 1]].iter().enumerate() {
                pos[h] = i;
            }
        }
        Self {
            coords,
            from,
            out_start,
            out,
            pos,
        }
    }

    pub(super) fn len(&self) -> usize {
        self.from.len()
    }

    pub(super) fn coord(&self, v: usize) -> Coord2d {
        self.coords[v]
    }

    pub(super) fn origin(&self, h: usize) -> usize {
        self.from[h]
    }

    /// The number of vertices.
    pub(super) fn vertex_count(&self) -> usize {
        self.coords.len()
    }

    /// The successor of `h = u → v` with the face on the left: the half-edge
    /// leaving `v` first clockwise from the reverse direction `v → u`. It is
    /// the twin's predecessor in `v`'s counterclockwise list.
    fn next(&self, h: usize) -> usize {
        let twin = h ^ 1;
        let v = self.from[twin];
        let list = &self.out[self.out_start[v]..self.out_start[v + 1]];
        list[(self.pos[twin] + list.len() - 1) % list.len()]
    }

    /// Calls `f` with every boundary walk, as the half-edges in order, each
    /// half-edge in exactly one walk, starting from `starts` in order (a
    /// half-edge already walked is skipped). With `starts` covering every
    /// half-edge, this traces every face. One buffer serves every walk.
    pub(super) fn for_each_walk(
        &self,
        starts: impl IntoIterator<Item = usize>,
        mut f: impl FnMut(&[usize]),
    ) {
        let mut seen = vec![false; self.len()];
        let mut walk = Vec::new();
        for s in starts {
            if seen[s] {
                continue;
            }
            walk.clear();
            let mut h = s;
            while !seen[h] {
                seen[h] = true;
                walk.push(h);
                h = self.next(h);
            }
            f(&walk);
        }
    }
}

// ===========================================================================
// Loops and rings
// ===========================================================================

/// Splits boundary walks into their simple closed loops, at every repeated
/// vertex, with scratch space reused from walk to walk.
///
/// Scanning a walk's origin vertices with a stack, a vertex seen earlier
/// closes the loop formed by the stack above its earlier position; that
/// loop is cut out and the scan continues from the vertex. What is left at
/// the end is the last loop. Each loop visits each vertex once. A loop of
/// fewer than three vertices is an edge walked out and back, which no
/// caller's graph has.
///
/// A vertex's stack position is kept in an array indexed by vertex
/// (`NONE` when it is not on the stack); every entry set during a walk is
/// cleared before the next, so the array is all `NONE` between walks.
pub(super) struct LoopSplitter {
    at: Vec<usize>,
    stack: Vec<usize>,
}

impl LoopSplitter {
    const NONE: usize = usize::MAX;

    /// Scratch for a graph of `vertices` vertices.
    pub(super) fn new(vertices: usize) -> Self {
        Self {
            at: vec![Self::NONE; vertices],
            stack: Vec::new(),
        }
    }

    /// Calls `f` with each loop of the walk whose origin vertices are
    /// `walk`, in the order the loops close, the last loop last.
    pub(super) fn split(
        &mut self,
        walk: impl IntoIterator<Item = usize>,
        mut f: impl FnMut(&[usize]),
    ) {
        for v in walk {
            let p = self.at[v];
            if p != Self::NONE {
                for &u in &self.stack[p..] {
                    self.at[u] = Self::NONE;
                }
                f(&self.stack[p..]);
                self.stack.truncate(p);
            }
            self.at[v] = self.stack.len();
            self.stack.push(v);
        }
        for &u in &self.stack {
            self.at[u] = Self::NONE;
        }
        if !self.stack.is_empty() {
            f(&self.stack);
        }
        self.stack.clear();
    }
}

/// A loop's ring, rotated to start at its smallest vertex
/// ([`coord_exact_cmp`], the first of equal ones), orientation kept. `None` if `LinearRing::new` refuses it, which no loop of a noded
/// graph can be (at least three distinct, not all collinear, vertices).
pub(super) fn ring_of(g: &HalfEdges<'_>, lp: &[usize]) -> Option<LinearRing> {
    let start = (0..lp.len())
        .min_by(|&i, &j| coord_exact_cmp(g.coord(lp[i]), g.coord(lp[j])).then(i.cmp(&j)))
        .unwrap_or(0);
    let coords: Vec<Coord2d> = lp[start..]
        .iter()
        .chain(&lp[..start])
        .map(|&v| g.coord(v))
        .collect();
    LinearRing::new(CoordSeq::new(coords)).ok()
}

/// `true` if a loop's ring, as [`ring_of`] returns it, runs
/// counterclockwise, decided exactly without its area.
///
/// The ring starts at its smallest vertex `v` ([`coord_exact_cmp`]) and is
/// simple. At that vertex the ring turns the way it winds, so the sign of
/// `orient(prev, v, next)` is its orientation. The sign is never zero: the
/// two neighbors cannot lie on opposite sides of `v` on one line, since one
/// of them would then be smaller than `v`, and they cannot lie on the same
/// side, since the two edges would then overlap. Used only where the signed
/// area is not finite, so its sign is unavailable.
pub(super) fn runs_counterclockwise(ring: &LinearRing) -> bool {
    let c = ring.coords().as_slice();
    crate::exact::kernel::orient(c[c.len() - 1], c[0], c[1]) == arithmetica::Sign::Positive
}

/// Assigns each hole to the innermost exterior containing it, by exact
/// placement, and builds the polygons. Each exterior comes with its signed
/// area, computed once by the caller to classify it. Returns the polygons (each
/// `Ok(polygon)` or `Err(shell)` if `Polygon::new` refused it) and the holes
/// no exterior contains.
///
/// A hole's candidates are the exteriors whose placement says the hole
/// lies `Inside` (touching at points allowed). Candidates are nested
/// regions, since each contains the hole, and the innermost is the
/// candidate no other candidate lies inside, decided by `place_ring` too.
pub(super) fn assemble(
    shells: Vec<(LinearRing, f64)>,
    holes: Vec<LinearRing>,
) -> (Vec<Result<Polygon, LinearRing>>, Vec<LinearRing>) {
    let rings: Vec<&LinearRing> = shells.iter().map(|(r, _)| r).collect();
    let mut owned: Vec<Vec<LinearRing>> = vec![Vec::new(); shells.len()];
    let mut orphans = Vec::new();
    for hole in holes {
        match innermost(&rings, &hole) {
            Some(i) => owned[i].push(hole),
            None => orphans.push(hole),
        }
    }
    let polys = shells
        .into_iter()
        .zip(owned)
        .map(|((shell, area), mut hs)| {
            if hs.is_empty() {
                // The area is in hand from the classification.
                return Polygon::from_exterior_with_area(shell, area);
            }
            hs.sort_by(|a, b| seq_cmp(a.coords().as_slice(), b.coords().as_slice()));
            Polygon::new_or_exterior(shell, hs)
        })
        .collect();
    (polys, orphans)
}

fn innermost(shells: &[&LinearRing], hole: &LinearRing) -> Option<usize> {
    let he = hole.envelope().ok()?;
    let candidates: Vec<usize> = (0..shells.len())
        .filter(|&i| {
            shells[i].envelope().is_ok_and(|e| e.contains_envelope(&he))
                && place_ring(hole, shells[i]) == RingPlacement::Inside
        })
        .collect();
    candidates.iter().copied().find(|&c| {
        candidates
            .iter()
            .all(|&d| d == c || place_ring(shells[d], shells[c]) != RingPlacement::Inside)
    })
}
