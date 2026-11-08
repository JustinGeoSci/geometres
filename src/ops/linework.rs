// src/ops/linework.rs

//! Operations on line work: [`line_merge`], [`line_merge_directed`],
//! [`line_sequence`] and [`split_by_point`].
//!
//! # Every output coordinate is an input coordinate
//!
//! Merging and sequencing join and reorient whole input line strings at
//! shared endpoints; they never cut one. Splitting cuts a line only at a
//! point the caller supplied, and only where that point lies exactly on the
//! line. The decisions are exact: endpoints are matched by `==` on stored
//! coordinates, and "on the line" is `orient2d` on stored coordinates.
//!
//! # The node graph
//!
//! Each input line string is an edge between its two endpoints; the nodes
//! are the distinct endpoints under `==` (so `-0.0` and `0.0` are one
//! node). A node's degree counts edge ends, so a closed line string adds two
//! to its own node.

use core::cmp::Ordering;
use std::collections::HashMap;

use crate::GeometryError;
use crate::exact::on_segment;
use crate::ops::basic::{key, lex, lex_seq};
use crate::predicates::{Borrowed, Relatable, borrowed};
use crate::types::{Coord2d, CoordSeq, Geometry, LineString, MultiLineString};

// ===========================================================================
// The graph
// ===========================================================================

/// The lineal leaves of a geometry, in order: a line string, a
/// multi-line-string's components, a collection's lineal leaves depth
/// first. Points and polygons contribute nothing.
pub(crate) fn lineal_leaves(g: Borrowed<'_>) -> Vec<&LineString> {
    match g {
        Borrowed::LineString(l) => vec![l],
        Borrowed::MultiLineString(m) => m.line_strings().iter().collect(),
        Borrowed::GeometryCollection(gc) => gc
            .leaves()
            .flat_map(|leaf| match leaf {
                Geometry::LineString(l) => vec![l],
                Geometry::MultiLineString(m) => m.line_strings().iter().collect(),
                _ => Vec::new(),
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// The undirected node graph of a set of line strings.
struct Graph<'a> {
    lines: Vec<&'a LineString>,
    /// Per node: its point (the first occurrence's bits).
    node_pt: Vec<Coord2d>,
    /// Per line: its start node and end node.
    ends: Vec<(usize, usize)>,
    /// Per node: the incident (line, the line leaves this node forward).
    /// A closed line appears twice, once each way.
    star: Vec<Vec<(usize, bool)>>,
}

impl<'a> Graph<'a> {
    fn new(lines: Vec<&'a LineString>) -> Self {
        let mut index: HashMap<(u64, u64), usize> = HashMap::new();
        let mut node_pt: Vec<Coord2d> = Vec::new();
        let mut node = |c: Coord2d, node_pt: &mut Vec<Coord2d>| {
            *index.entry(key(c)).or_insert_with(|| {
                node_pt.push(c);
                node_pt.len() - 1
            })
        };
        let mut ends = Vec::with_capacity(lines.len());
        for l in &lines {
            let s = node(l.start(), &mut node_pt);
            let e = node(l.end(), &mut node_pt);
            ends.push((s, e));
        }
        let mut star = vec![Vec::new(); node_pt.len()];
        for (k, &(s, e)) in ends.iter().enumerate() {
            star[s].push((k, true));
            star[e].push((k, false));
        }
        Self {
            lines,
            node_pt,
            ends,
            star,
        }
    }

    fn degree(&self, v: usize) -> usize {
        self.star[v].len()
    }

    /// The node at the far end of line `k` traversed `forward` (from its
    /// start) or backward.
    fn far(&self, k: usize, forward: bool) -> usize {
        if forward {
            self.ends[k].1
        } else {
            self.ends[k].0
        }
    }
}

/// The coordinates of a chain of (line, forward) steps, joined at their
/// shared nodes, each junction once.
fn chain_coords(lines: &[&LineString], chain: &[(usize, bool)]) -> Vec<Coord2d> {
    let mut out: Vec<Coord2d> = Vec::new();
    for &(k, fwd) in chain {
        let v = lines[k].coords().as_slice();
        let it: Box<dyn Iterator<Item = &Coord2d>> = if fwd {
            Box::new(v.iter())
        } else {
            Box::new(v.iter().rev())
        };
        let skip = usize::from(!out.is_empty());
        out.extend(it.skip(skip).copied());
    }
    out
}

/// Joined line strings, sorted by vertex sequence (stable).
fn finish(mut lines: Vec<Vec<Coord2d>>) -> MultiLineString {
    lines.sort_by(|a, b| lex_seq(a, b));
    // Each joined sequence is input line strings concatenated at shared
    // endpoints with the shared vertex written once: at least two vertices,
    // finite, and no two consecutive equal (each input line has none, and a
    // junction joins two lines that are equal only at the junction vertex,
    // which appears once).
    MultiLineString::new(
        lines
            .into_iter()
            .map(|c| LineString::from_valid(CoordSeq::new(c)))
            .collect(),
    )
}

// ===========================================================================
// line_merge
// ===========================================================================

/// Orients a merged chain by majority: if more of its input lines run
/// backward than forward along it, it is reversed. On a tie, an open chain
/// starts at its lexicographically smaller end; a closed one keeps its start
/// and leaves it toward the smaller of its two neighbors.
fn orient_by_majority(chain: &mut [(usize, bool)], coords: &mut [Coord2d]) {
    let fwd = chain.iter().filter(|s| s.1).count();
    let back = chain.len() - fwd;
    let reverse = match back.cmp(&fwd) {
        Ordering::Greater => true,
        Ordering::Less => false,
        Ordering::Equal => {
            let n = coords.len();
            if coords[0] == coords[n - 1] {
                n >= 3 && lex(coords[n - 2], coords[1]) == Ordering::Less
            } else {
                lex(coords[n - 1], coords[0]) == Ordering::Less
            }
        }
    };
    if reverse {
        chain.reverse();
        for s in chain.iter_mut() {
            s.1 = !s.1;
        }
        coords.reverse();
    }
}

/// Merges line strings that meet end to end at nodes of degree 2 into
/// maximal line strings.
///
/// The input is `g`'s lineal leaves: a `LineString`, a `MultiLineString`'s
/// components, or a collection's lineal leaves. Points and polygons are
/// ignored. The lines are taken as given (merging does not node them):
/// lines are joined only at shared **endpoints**, matched exactly.
///
/// - **Nodes.** A node of degree 2 joins its two lines; any other degree
///   (1, or 3 and more) ends them.
/// - **Orientation**. A merged line runs in the direction of the majority of
///   its input lines. On a tie, an open line starts at its lexicographically
///   smaller end.
/// - **Loops.** A component whose nodes all have degree 2 is a loop. It
///   becomes one closed line starting at its smallest node, oriented by
///   majority; on a tie it leaves its start toward the smaller neighbor.
/// - **Output order:** by vertex sequence, lexicographically.
///
/// Every output vertex is an input vertex, and every input segment appears
/// in exactly one output line, once.
#[must_use]
pub fn line_merge(g: &impl Relatable) -> MultiLineString {
    let graph = Graph::new(lineal_leaves(borrowed(g)));
    let n = graph.lines.len();
    let mut used = vec![false; n];
    let mut out: Vec<Vec<Coord2d>> = Vec::new();
    // Walk from `v` along (k, fwd) through degree-2 nodes.
    let walk = |v: usize, first: (usize, bool), used: &mut [bool]| -> Vec<(usize, bool)> {
        let mut chain = vec![first];
        used[first.0] = true;
        let mut at = graph.far(first.0, first.1);
        while at != v && graph.degree(at) == 2 {
            let Some(&(k, fwd)) = graph.star[at].iter().find(|&&(k, _)| !used[k]) else {
                break;
            };
            used[k] = true;
            chain.push((k, fwd));
            at = graph.far(k, fwd);
        }
        chain
    };
    // Open chains: from every node whose degree is not 2.
    for v in 0..graph.node_pt.len() {
        if graph.degree(v) == 2 {
            continue;
        }
        for i in 0..graph.star[v].len() {
            let first = graph.star[v][i];
            if !used[first.0] {
                let mut chain = walk(v, first, &mut used);
                let mut coords = chain_coords(&graph.lines, &chain);
                orient_by_majority(&mut chain, &mut coords);
                out.push(coords);
            }
        }
    }
    // Loops: every remaining line lies on a cycle of degree-2 nodes. Start at
    // the loop's smallest node.
    for k0 in 0..n {
        if used[k0] {
            continue;
        }
        let probe = walk(graph.ends[k0].0, (k0, true), &mut used.clone());
        let start = probe
            .iter()
            .map(|&(k, fwd)| {
                if fwd {
                    graph.ends[k].0
                } else {
                    graph.ends[k].1
                }
            })
            .min_by(|&a, &b| lex(graph.node_pt[a], graph.node_pt[b]))
            .unwrap_or(graph.ends[k0].0);
        let first = graph.star[start]
            .iter()
            .copied()
            .find(|&(k, _)| !used[k])
            .unwrap_or((k0, true));
        let mut chain = walk(start, first, &mut used);
        let mut coords = chain_coords(&graph.lines, &chain);
        orient_by_majority(&mut chain, &mut coords);
        out.push(coords);
    }
    finish(out)
}

/// As [`line_merge`], but joining only lines that agree in direction: two
/// lines meet at a node only when one ends there and the other starts there
/// and no other line touches it. Every output line keeps its input lines'
/// direction; a directed loop starts at its smallest node.
#[must_use]
pub fn line_merge_directed(g: &impl Relatable) -> MultiLineString {
    let graph = Graph::new(lineal_leaves(borrowed(g)));
    let n = graph.lines.len();
    let nv = graph.node_pt.len();
    let (mut ins, mut outs) = (vec![0usize; nv], vec![0usize; nv]);
    for &(s, e) in &graph.ends {
        outs[s] += 1;
        ins[e] += 1;
    }
    let joins = |v: usize| ins[v] == 1 && outs[v] == 1;
    // The line leaving `v` forward, when `v` joins.
    let next_from = |v: usize| graph.star[v].iter().find(|s| s.1).map(|s| s.0);
    let mut used = vec![false; n];
    let mut out = Vec::new();
    let follow = |k0: usize, used: &mut [bool]| -> Vec<(usize, bool)> {
        let mut chain = vec![(k0, true)];
        used[k0] = true;
        let mut at = graph.ends[k0].1;
        while joins(at) {
            match next_from(at) {
                Some(k) if !used[k] => {
                    used[k] = true;
                    chain.push((k, true));
                    at = graph.ends[k].1;
                }
                _ => break,
            }
        }
        chain
    };
    // Chains start at lines whose start node does not join.
    for k in 0..n {
        if !used[k] && !joins(graph.ends[k].0) {
            let chain = follow(k, &mut used);
            out.push(chain_coords(&graph.lines, &chain));
        }
    }
    // Directed loops: start at the smallest node, on the line leaving it.
    for k0 in 0..n {
        if used[k0] {
            continue;
        }
        let mut probe_used = used.clone();
        let probe = follow(k0, &mut probe_used);
        let start_line = probe
            .iter()
            .map(|s| s.0)
            .min_by(|&a, &b| {
                lex(
                    graph.node_pt[graph.ends[a].0],
                    graph.node_pt[graph.ends[b].0],
                )
            })
            .unwrap_or(k0);
        let chain = follow(start_line, &mut used);
        out.push(chain_coords(&graph.lines, &chain));
    }
    finish(out)
}

// ===========================================================================
// line_sequence
// ===========================================================================

/// Orders and orients `g`'s lineal leaves so that each starts where the
/// previous one ended, one sequence per connected component.
///
/// The lines are not merged: the result holds each input line once,
/// possibly reversed. A component can be sequenced exactly when it has at
/// most two nodes of odd degree (an Euler path exists). The rules that make
/// the result deterministic:
///
/// - **Start node:** the smallest odd-degree node, preferring one where an
///   input line already starts; with no odd-degree node, the smallest node.
/// - **Next line:** at each node, the unused incident line whose far node
///   is smallest, then the one earliest in the input (Hierholzer's
///   algorithm, so the path covers every line once).
/// - **Components** in order of their start nodes.
///
/// # Errors
///
/// [`GeometryError::NotSequenceable`] if some component has more than two
/// nodes of odd degree.
pub fn line_sequence(g: &impl Relatable) -> Result<MultiLineString, GeometryError> {
    let graph = Graph::new(lineal_leaves(borrowed(g)));
    let nv = graph.node_pt.len();
    // Components of nodes.
    let mut comp = vec![usize::MAX; nv];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for v0 in 0..nv {
        if comp[v0] != usize::MAX {
            continue;
        }
        let c = comps.len();
        let mut stack = vec![v0];
        comp[v0] = c;
        let mut members = Vec::new();
        while let Some(v) = stack.pop() {
            members.push(v);
            for &(k, fwd) in &graph.star[v] {
                let w = graph.far(k, fwd);
                if comp[w] == usize::MAX {
                    comp[w] = c;
                    stack.push(w);
                }
            }
        }
        comps.push(members);
    }
    let starts_a_line: Vec<bool> = {
        let mut s = vec![false; nv];
        for &(a, _) in &graph.ends {
            s[a] = true;
        }
        s
    };
    let mut sequences: Vec<(usize, Vec<(usize, bool)>)> = Vec::new();
    for members in &comps {
        let odd: Vec<usize> = members
            .iter()
            .copied()
            .filter(|&v| graph.degree(v) % 2 == 1)
            .collect();
        if odd.len() > 2 {
            return Err(GeometryError::NotSequenceable);
        }
        let smallest = |vs: &mut dyn Iterator<Item = usize>| {
            vs.min_by(|&a, &b| lex(graph.node_pt[a], graph.node_pt[b]))
        };
        let start = if odd.is_empty() {
            smallest(&mut members.iter().copied())
        } else {
            smallest(&mut odd.iter().copied().filter(|&v| starts_a_line[v]))
                .or_else(|| smallest(&mut odd.iter().copied()))
        };
        let Some(start) = start else { continue };
        sequences.push((start, euler_path(&graph, start)));
    }
    sequences.sort_by(|a, b| lex(graph.node_pt[a.0], graph.node_pt[b.0]));
    let mut out = Vec::new();
    for (_, path) in sequences {
        for (k, fwd) in path {
            let l = graph.lines[k];
            out.push(if fwd {
                l.clone()
            } else {
                crate::ops::basic::Reversible::reversed(l)
            });
        }
    }
    Ok(MultiLineString::new(out))
}

/// Hierholzer's algorithm from `start` over the lines of its component,
/// choosing at each node the unused line whose far node is smallest, then
/// the earliest. Returns the path as (line, forward) steps in order.
fn euler_path(graph: &Graph<'_>, start: usize) -> Vec<(usize, bool)> {
    let mut used = vec![false; graph.lines.len()];
    // Each node's incident lines, sorted once by the choice rule.
    let order = |v: usize| {
        let mut s = graph.star[v].clone();
        s.sort_by(|a, b| {
            lex(
                graph.node_pt[graph.far(a.0, a.1)],
                graph.node_pt[graph.far(b.0, b.1)],
            )
            .then(a.0.cmp(&b.0))
            .then(b.1.cmp(&a.1))
        });
        s
    };
    let stars: Vec<Vec<(usize, bool)>> = (0..graph.node_pt.len()).map(order).collect();
    let mut cursor = vec![0usize; graph.node_pt.len()];
    // The stack holds (node, the step that reached it).
    let mut stack: Vec<(usize, Option<(usize, bool)>)> = vec![(start, None)];
    let mut path: Vec<(usize, bool)> = Vec::new();
    while let Some(&(v, via)) = stack.last() {
        let mut advanced = false;
        while cursor[v] < stars[v].len() {
            let (k, fwd) = stars[v][cursor[v]];
            cursor[v] += 1;
            if !used[k] {
                used[k] = true;
                stack.push((graph.far(k, fwd), Some((k, fwd))));
                advanced = true;
                break;
            }
        }
        if !advanced {
            stack.pop();
            if let Some(step) = via {
                path.push(step);
            }
        }
    }
    path.reverse();
    path
}

// ===========================================================================
// split_by_point
// ===========================================================================

/// Splits `line` at `p`: at every passage of the line through `p`, the line
/// is cut there, so a line through `p` once gives two line strings meeting
/// at `p`, and one through it `k` times gives `k + 1`.
///
/// A passage is a vertex of the line equal to `p` other than its first and
/// last, or a segment with `p` exactly on its interior (`orient2d` on stored
/// coordinates, and `p` within the segment's box). `p` is written into the
/// result as given, so every output coordinate is a vertex of `line` or `p`.
///
/// A point *near* the line is not snapped to it: it is not on the line. A
/// caller who needs a point brought onto a line uses a constructive split
/// or snapping, which construct coordinates and so do not belong here.
///
/// # Errors
///
/// - [`GeometryError::PointNotOnLineInterior`] if `p` is the line's first
///   or last vertex (for a closed line, its closing vertex), or not on the
///   line at all.
/// - [`GeometryError::NonFiniteCoordinate`] (position 0) if `p` is not
///   finite: the exact test gives no answer for a point that is not a point.
pub fn split_by_point(line: &LineString, p: Coord2d) -> Result<MultiLineString, GeometryError> {
    let v = line.coords().as_slice();
    let n = v.len();
    if v[0] == p || v[n - 1] == p {
        return Err(GeometryError::PointNotOnLineInterior);
    }
    let mut parts: Vec<Vec<Coord2d>> = Vec::new();
    let mut cur: Vec<Coord2d> = vec![v[0]];
    for i in 0..n - 1 {
        let (a, b) = (v[i], v[i + 1]);
        if a != p && b != p && on_segment(p, a, b)? {
            cur.push(p);
            parts.push(core::mem::replace(&mut cur, vec![p]));
        }
        cur.push(b);
        if b == p && i + 1 < n - 1 {
            parts.push(core::mem::replace(&mut cur, vec![b]));
        }
    }
    if parts.is_empty() {
        return Err(GeometryError::PointNotOnLineInterior);
    }
    parts.push(cur);
    // Each part runs from a vertex or `p` to the next cut: at least two
    // coordinates, no two consecutive equal (`p` is strictly inside the
    // segment it splits, so it differs from both ends).
    Ok(MultiLineString::new(
        parts
            .into_iter()
            .map(|c| LineString::from_valid(CoordSeq::new(c)))
            .collect(),
    ))
}
