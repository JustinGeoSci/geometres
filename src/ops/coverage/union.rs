// src/ops/coverage/union.rs

//! The union of a valid coverage by removing shared edges (the argument is
//! on `coverage_union`).

use std::collections::BTreeMap;

use crate::GeometryError;
use crate::types::{Coord2d, LinearRing, MultiPolygon, Polygon};

use super::graph::{HalfEdges, Key, LoopSplitter, assemble, key, ring_of, seq_cmp};
use super::validate::edges_of;

/// The dissolve of `polys`, which must be a valid coverage.
///
/// # Errors
///
/// [`GeometryError::NonFiniteResult`] when a ring of the union has a signed
/// area that is not finite. Each input polygon has a finite area, but their
/// union can exceed the range of `f64` (coordinates near 1e154 and beyond),
/// and `Polygon::new` refuses such a ring. The refusal is reported rather
/// than letting the ring be read as a hole.
///
/// # Panics
///
/// Only on a broken invariant, never on valid input: an edge appearing twice
/// in the same direction (overlapping interiors), a traced loop that is not
/// a ring, a hole no exterior contains, or a polygon or multipolygon its
/// constructor refuses. Each contradicts the validation that ran first, so
/// each is a defect to report: a panic rather than a substitute answer.
pub(super) fn dissolve(polys: &[Polygon]) -> Result<MultiPolygon, GeometryError> {
    let (coords, boundary) = boundary_edges(polys);
    if boundary.is_empty() {
        return Ok(MultiPolygon::empty());
    }
    let graph = HalfEdges::new(&coords, &boundary);
    // Original directions are the even half-edges: the union's interior is
    // on their left, so their walks are exactly the interior faces' walks.
    let starts = (0..graph.len()).step_by(2);
    let mut shells: Vec<(LinearRing, f64)> = Vec::new();
    let mut holes: Vec<LinearRing> = Vec::new();
    let mut overflow = false;
    let mut loops = LoopSplitter::new(graph.vertex_count());
    graph.for_each_walk(starts, |walk| {
        loops.split(walk.iter().map(|&h| graph.origin(h)), |lp| {
            let ring = ring_of(&graph, lp).expect(
                "a dissolved boundary loop is a ring: noded, simple, at least three vertices",
            );
            let area = ring.signed_area();
            if !area.is_finite() {
                overflow = true;
            } else if area > 0.0 {
                shells.push((ring, area));
            } else {
                holes.push(ring);
            }
        });
    });
    if overflow {
        return Err(GeometryError::NonFiniteResult);
    }
    let (polys, orphans) = assemble(shells, holes);
    assert!(
        orphans.is_empty(),
        "every hole of a valid coverage's union lies inside an exterior"
    );
    let mut parts: Vec<Polygon> = polys
        .into_iter()
        .map(|p| p.expect("a dissolved face with exactly placed holes is a valid polygon"))
        .collect();
    parts.sort_by(|a, b| {
        seq_cmp(
            a.exterior().coords().as_slice(),
            b.exterior().coords().as_slice(),
        )
    });
    Ok(MultiPolygon::new(parts).expect("the faces of a union meet at most at points"))
}

/// The union's boundary: every directed ring edge whose reverse does not
/// also appear, as vertex-index pairs over the distinct coordinates (in
/// first-seen order).
///
/// In a valid coverage a shared edge appears exactly twice, once in each
/// direction; an edge appearing twice in the same direction would mean two
/// interiors on the same side of it, which validation excluded.
fn boundary_edges(polys: &[Polygon]) -> (Vec<Coord2d>, Vec<(usize, usize)>) {
    let mut ids: BTreeMap<Key, usize> = BTreeMap::new();
    let mut coords: Vec<Coord2d> = Vec::new();
    let mut id = |c: Coord2d, coords: &mut Vec<Coord2d>| {
        *ids.entry(key(c)).or_insert_with(|| {
            coords.push(c);
            coords.len() - 1
        })
    };
    let mut directed: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for e in edges_of(polys) {
        let (a, b) = (id(e.a, &mut coords), id(e.b, &mut coords));
        let n = directed.entry((a, b)).or_insert(0);
        *n += 1;
        assert!(
            *n == 1,
            "an edge twice in the same direction means overlapping interiors"
        );
    }
    let boundary = directed
        .keys()
        .filter(|&&(a, b)| !directed.contains_key(&(b, a)))
        .copied()
        .collect();
    (coords, boundary)
}
