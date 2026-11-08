// src/ops/coverage/polygonize.rs

//! The faces of noded linework (the argument is on `polygonize_noded`).

use crate::GeometryError;
use crate::exact::kernel::{on_segment, segments_properly_cross};
use crate::index::BoxIndex;
use crate::ops::distance::{Element, elements};
use crate::ops::line_merge;
use crate::predicates::Borrowed;
use crate::types::{Coord2d, CoordSeq, LineString, LinearRing, MultiLineString, Polygon};

use super::Polygonized;
use super::graph::{
    HalfEdges, LoopSplitter, assemble, coord_exact_cmp, key, ring_of, runs_counterclockwise,
    seq_cmp,
};

/// Polygonizes the linework of `g`.
pub(super) fn run(g: Borrowed<'_>) -> Result<Polygonized, GeometryError> {
    let segments = distinct_segments(g);
    check_noded(&segments)?;

    let (coords, edges) = number_vertices(&segments);

    let (kept, dangles) = remove_dangles(coords.len(), &edges);
    let (graph, index) = subgraph(&coords, &edges, &kept);
    let (kept, cuts) = remove_cut_edges(&graph, &index, kept);
    // With no cut edge removed, the faces' graph is the graph just built
    // (the same kept edges, so the same half-edges); it is rebuilt only
    // when cut edges went.
    let (polygons, invalid) = if cuts.is_empty() {
        faces(&graph)
    } else {
        faces(&subgraph(&coords, &edges, &kept).0)
    };

    let as_lines = |ids: &[usize]| {
        let lines: Vec<LineString> = ids
            .iter()
            .filter_map(|&e| {
                let (a, b) = edges[e];
                LineString::new(CoordSeq::new(vec![coords[a], coords[b]])).ok()
            })
            .collect();
        line_merge(&MultiLineString::new(lines))
    };
    Ok(Polygonized {
        polygons,
        dangles: as_lines(&dangles),
        cut_edges: as_lines(&cuts),
        invalid_rings: MultiLineString::new(invalid),
    })
}

/// The vertices of `segments`, numbered in first-seen order over the
/// endpoint sequence `a₀, b₀, a₁, b₁, …` (each kept as first met, with its
/// stored bits), and each segment as its two vertex numbers.
///
/// Endpoints are sorted by key with a stable sort: equal keys are exactly
/// equal coordinates, so each run of equal keys is one vertex, its entries
/// in position order, and its first entry is its first occurrence.
/// Numbering the runs in order of first occurrence gives the numbers a
/// first-seen map assigns, without hashing. (A key is packed into one
/// `u128` for the comparisons; the packing is injective.)
pub(super) fn number_vertices(
    segments: &[(Coord2d, Coord2d)],
) -> (Vec<Coord2d>, Vec<(usize, usize)>) {
    let packed = |c: Coord2d| {
        let (x, y) = key(c);
        (u128::from(x) << 64) | u128::from(y)
    };
    let mut ends: Vec<(u128, usize)> = Vec::with_capacity(2 * segments.len());
    for (i, &(a, b)) in segments.iter().enumerate() {
        ends.push((packed(a), 2 * i));
        ends.push((packed(b), 2 * i + 1));
    }
    ends.sort_by_key(|e| e.0);
    // Each endpoint's run, and each run's first occurrence.
    let mut run_of = vec![0_usize; ends.len()];
    let mut firsts: Vec<(usize, usize)> = Vec::new();
    for (i, &(k, pos)) in ends.iter().enumerate() {
        if i == 0 || ends[i - 1].0 != k {
            firsts.push((pos, firsts.len()));
        }
        run_of[pos] = firsts.len() - 1;
    }
    firsts.sort_unstable();
    let mut number = vec![0_usize; firsts.len()];
    let mut coords = Vec::with_capacity(firsts.len());
    for (n, &(pos, run)) in firsts.iter().enumerate() {
        number[run] = n;
        let (a, b) = segments[pos / 2];
        coords.push(if pos % 2 == 0 { a } else { b });
    }
    let edges = (0..segments.len())
        .map(|i| (number[run_of[2 * i]], number[run_of[2 * i + 1]]))
        .collect();
    (coords, edges)
}

/// The distinct segments of the linework (points ignored), each with its
/// endpoints in exact lexicographic order, sorted. Identical segments in
/// either direction are one, and the one kept is the first met, with its
/// stored bits (a `-0.0` included).
///
/// # Why one sort does the work of a set and a sort
///
/// Each segment, endpoints ordered by [`coord_exact_cmp`], gets the integer
/// key [`order_key`] of its four values. For finite values, `x + 0.0` maps
/// `-0.0` to `0.0` and changes nothing else, and [`total_key`] is strictly
/// increasing in the value; so two keys compare as the segments compare
/// under `coord_exact_cmp`, and they are equal exactly when the segments
/// have the same two endpoints, which is the identity
/// [`edge_key`](super::graph) gives. Sorting by `(key, position met)` and
/// keeping the first of each run of equal keys therefore keeps, for each
/// distinct segment, the copy a first-seen set kept, in the order a sort of
/// those copies by `coord_exact_cmp` gives (their keys are distinct, so that
/// order is total).
pub(super) fn distinct_segments(g: Borrowed<'_>) -> Vec<(Coord2d, Coord2d)> {
    let mut keyed: Vec<Keyed> = Vec::new();
    for e in elements(g) {
        if let Element::Segment(a, b) = e {
            let s = if coord_exact_cmp(a, b).is_le() {
                (a, b)
            } else {
                (b, a)
            };
            keyed.push((order_key(s), keyed.len(), s));
        }
    }
    keyed.sort_unstable_by_key(|k| (k.0, k.1));
    keyed.dedup_by(|later, kept| later.0 == kept.0);
    keyed.into_iter().map(|(_, _, s)| s).collect()
}

/// A segment with its [`order_key`] and the position it was met at.
type Keyed = ((u128, u128), usize, (Coord2d, Coord2d));

/// The four values of an ordered segment as [`total_key`]s, `-0.0` first
/// mapped to `0.0`, packed two to a `u128` (high half first, so the packing
/// keeps the lexicographic order of the four).
fn order_key((a, b): (Coord2d, Coord2d)) -> (u128, u128) {
    let pair =
        |c: Coord2d| (u128::from(total_key(c.x + 0.0)) << 64) | u128::from(total_key(c.y + 0.0));
    (pair(a), pair(b))
}

/// `x`'s bits mapped so that unsigned order is `f64::total_cmp`'s order,
/// which on finite values other than `-0.0` is the numeric order:
/// negative values (sign bit set) have every bit flipped, the others only
/// the sign bit.
fn total_key(x: f64) -> u64 {
    let b = x.to_bits();
    if b >> 63 == 1 { !b } else { b | (1 << 63) }
}

// ===========================================================================
// The noding check
// ===========================================================================

/// `Ok` when every pair of distinct segments meets only at an endpoint of
/// both; otherwise `NotNoded` for the first offending pair in segment
/// order.
///
/// Two distinct segments that share a point are noded there exactly when
/// the point is an endpoint of both. So a pair is offending when it
/// properly crosses, or an endpoint of one lies in the other's interior.
/// Every collinear overlap of positive length between distinct segments
/// implies the second. All exact: `orient2d` and stored-coordinate order.
fn check_noded(segments: &[(Coord2d, Coord2d)]) -> Result<(), GeometryError> {
    let entries: Vec<((Coord2d, Coord2d), u32)> = segments
        .iter()
        .enumerate()
        .map(|(i, &s)| (s, u32::try_from(i).unwrap_or(u32::MAX)))
        .collect();
    let mut first: Option<(usize, usize)> = None;
    BoxIndex::from_edges(&entries).self_join(|i, j| {
        let (i, j) = (i.min(j) as usize, i.max(j) as usize);
        if first.is_some_and(|f| f <= (i, j)) {
            return;
        }
        if !noded_pair(segments[i], segments[j]) {
            first = Some((i, j));
        }
    });
    match first {
        None => Ok(()),
        Some((i, j)) => Err(GeometryError::NotNoded {
            first_pair: (segments[i], segments[j]),
        }),
    }
}

/// Two distinct segments meet only at an endpoint shared by both: they do
/// not properly cross, and no endpoint of one lies strictly inside the
/// other (which also excludes every collinear overlap of positive length:
/// two distinct overlapping collinear segments always have an endpoint
/// strictly inside the other).
pub(super) fn noded_pair(s: (Coord2d, Coord2d), t: (Coord2d, Coord2d)) -> bool {
    if segments_properly_cross(s.0, s.1, t.0, t.1) {
        return false;
    }
    let strictly_inside =
        |p: Coord2d, (a, b): (Coord2d, Coord2d)| p != a && p != b && on_segment(p, a, b);
    let endpoint_inside = strictly_inside(t.0, s)
        || strictly_inside(t.1, s)
        || strictly_inside(s.0, t)
        || strictly_inside(s.1, t);
    !endpoint_inside
}

// ===========================================================================
// Dangles and cut edges
// ===========================================================================

/// Removes edges with a free end until none is left. Returns the kept edge
/// flags and the removed edges, in removal order.
///
/// Each vertex's incident edges are in one flat array, in increasing edge
/// order (`incident[start[v]..start[v + 1]]`).
fn remove_dangles(nv: usize, edges: &[(usize, usize)]) -> (Vec<bool>, Vec<usize>) {
    let mut start = vec![0_usize; nv + 1];
    for &(a, b) in edges {
        start[a + 1] += 1;
        start[b + 1] += 1;
    }
    for v in 0..nv {
        start[v + 1] += start[v];
    }
    let mut degree: Vec<usize> = (0..nv).map(|v| start[v + 1] - start[v]).collect();
    let mut fill = start.clone();
    let mut incident = vec![0_usize; start[nv]];
    for (e, &(a, b)) in edges.iter().enumerate() {
        for v in [a, b] {
            incident[fill[v]] = e;
            fill[v] += 1;
        }
    }
    let mut kept = vec![true; edges.len()];
    let mut removed = Vec::new();
    let mut queue: Vec<usize> = (0..nv).filter(|&v| degree[v] == 1).collect();
    while let Some(v) = queue.pop() {
        if degree[v] != 1 {
            continue;
        }
        let Some(&e) = incident[start[v]..start[v + 1]].iter().find(|&&e| kept[e]) else {
            continue;
        };
        kept[e] = false;
        removed.push(e);
        let (a, b) = edges[e];
        for u in [a, b] {
            degree[u] -= 1;
            if degree[u] == 1 {
                queue.push(u);
            }
        }
    }
    (kept, removed)
}

/// The graph of the kept edges, with each graph edge's index into `edges`.
fn subgraph<'a>(
    coords: &'a [Coord2d],
    edges: &[(usize, usize)],
    kept: &[bool],
) -> (HalfEdges<'a>, Vec<usize>) {
    let index: Vec<usize> = (0..edges.len()).filter(|&e| kept[e]).collect();
    let pairs: Vec<(usize, usize)> = index.iter().map(|&e| edges[e]).collect();
    (HalfEdges::new(coords, &pairs), index)
}

/// Removes the edges whose two half-edges lie on the same boundary walk:
/// the same face on both sides (a bridge). A non-bridge edge lies on a
/// cycle, which separates the plane, so its sides are different faces.
/// Removing every bridge at once leaves no dangle: a vertex keeps either no
/// edge or the edges of cycles through it.
///
/// `graph` and `index` are [`subgraph`] of the kept edges.
fn remove_cut_edges(
    graph: &HalfEdges<'_>,
    index: &[usize],
    mut kept: Vec<bool>,
) -> (Vec<bool>, Vec<usize>) {
    let mut walk_of = vec![0_usize; graph.len()];
    let mut w = 0;
    graph.for_each_walk(0..graph.len(), |walk| {
        for &h in walk {
            walk_of[h] = w;
        }
        w += 1;
    });
    let mut cuts = Vec::new();
    for (i, &e) in index.iter().enumerate() {
        if walk_of[2 * i] == walk_of[2 * i + 1] {
            kept[e] = false;
            cuts.push(e);
        }
    }
    (kept, cuts)
}

// ===========================================================================
// Faces
// ===========================================================================

/// The bounded faces of the kept edges, as polygons sorted by exterior.
///
/// Every half-edge is walked once. Each walk is split into simple loops;
/// counterclockwise loops are face exteriors, clockwise loops are holes of
/// the innermost exterior containing them, and clockwise loops no exterior
/// contains bound the unbounded face and are dropped.
/// Loops that cannot become polygons (a loop `LinearRing::new` refuses, or
/// an exterior `Polygon::new` refuses) are returned as closed lines. For
/// noded input there are none (see `Polygonized::invalid_rings`); they are
/// reported rather than dropped all the same.
fn faces(graph: &HalfEdges<'_>) -> (Vec<Polygon>, Vec<LineString>) {
    let mut shells: Vec<(LinearRing, f64)> = Vec::new();
    let mut holes: Vec<LinearRing> = Vec::new();
    let mut invalid: Vec<LineString> = Vec::new();
    let mut loops = LoopSplitter::new(graph.vertex_count());
    graph.for_each_walk(0..graph.len(), |walk| {
        loops.split(walk.iter().map(|&h| graph.origin(h)), |lp| {
            let Some(ring) = ring_of(graph, lp) else {
                let mut c: Vec<Coord2d> = lp.iter().map(|&v| graph.coord(v)).collect();
                c.push(c[0]);
                invalid.extend(LineString::new(CoordSeq::new(c)).ok());
                return;
            };
            let area = ring.signed_area();
            if !area.is_finite() {
                // An area beyond the range of `f64`: no polygon can carry the
                // ring. A counterclockwise loop bounds a face, which is
                // reported rather than read as a hole. A clockwise one is the
                // outline of the unbounded face, or a hole of a face larger
                // still, which is itself reported, so it is dropped as any
                // hole without an exterior is.
                if runs_counterclockwise(&ring) {
                    invalid.extend(closed_line(&ring));
                }
            } else if area > 0.0 {
                shells.push((ring, area));
            } else {
                holes.push(ring);
            }
        });
    });
    let (polys, _unbounded) = assemble(shells, holes);
    let mut out: Vec<Polygon> = Vec::new();
    for p in polys {
        match p {
            Ok(p) => out.push(p),
            Err(shell) => invalid.extend(closed_line(&shell)),
        }
    }
    out.sort_by(|a, b| {
        seq_cmp(
            a.exterior().coords().as_slice(),
            b.exterior().coords().as_slice(),
        )
    });
    (out, invalid)
}

fn closed_line(r: &LinearRing) -> Option<LineString> {
    let mut c = r.coords().as_slice().to_vec();
    c.push(c[0]);
    LineString::new(CoordSeq::new(c)).ok()
}
