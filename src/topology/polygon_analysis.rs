// src/topology/polygon_analysis.rs

//! Exact DE-9IM computation for area × area pairs: polygons and
//! multipolygons, in any combination.
//!
//! Every cell is derived from one question asked of every ring edge of each
//! operand: *which parts of this edge lie in the interior, on the boundary,
//! or in the exterior of the other operand?* That question is answered by
//! [`summarize_segment`], the exact sub-interval classifier already used for
//! line × polygon (see [`line_analysis`](super::line_analysis)): every sign
//! decision is `arithmetica::orient2d` on the original `f64` coordinates, and
//! no intersection point is ever constructed or rounded.
//!
//! A multipolygon enters as the ring set of all its parts at once: one
//! classification of all parts of A against all parts of B, rather than
//! pairwise matrices combined. [`Areal`](super::areal) explains why every
//! property used below holds of that ring set. The boundary is topological:
//! a point where parts touch is on the boundary.
//!
//! # Why edges rather than faces
//!
//! A sweep that builds a DCEL and flood-fills it labels only the faces the
//! sweep creates. Where the boundaries touch without crossing (a shared
//! vertex, a vertex on an edge, a shared edge segment, nesting with contact)
//! the sweep produces no split vertices, so the faces and the boundary cells
//! would have to be patched in afterwards, case by case. Classifying ring
//! edges needs no such regime split: crossing, touching, sharing, nesting,
//! and holes are all different classifications of the same edges.
//!
//! # Line cells
//!
//! Let `pa` summarize the ring edges of A against B and `pb` the ring edges
//! of B against A.
//!
//! - `BI` (∂A ∩ int B) is `Line` iff `pa` has an interior piece; `BE` iff `pa`
//!   has an exterior piece. `IB` (int A ∩ ∂B) and `EB` likewise from `pb`.
//!   None of these can be a `Point`: the interior and exterior of an areal
//!   operand are open, so a curve that meets one meets it in a piece of
//!   positive length.
//! - `BB` is `Line` iff either summary has a piece along the other boundary.
//!
//! # Point cells
//!
//! Otherwise `BB` is `Point` iff the boundaries meet at isolated points. The
//! candidates are finite: every vertex of A on ∂B, every vertex of B on ∂A,
//! and every proper crossing of an edge of A with an edge of B at no vertex
//! of either. Any other isolated contact would be two edges touching at an
//! interior point of both, which is collinear overlap (a stretch, already
//! counted) or impossible. The point where two parts of one operand touch is
//! a vertex on that operand's own boundary, so it is found like any other.
//!
//! # Area cells
//!
//! The rings of A and B cut the plane into faces. Each face lies wholly in one
//! of `int A | ext A` times `int B | ext B`, and every face is bounded by at
//! least one open piece of ring edge of positive length. So the set of face
//! kinds present is the union, over all such pieces, of the kinds of the two
//! faces beside it, and each piece's two sides follow from its classification
//! (polygon interiors lie to the left of every directed ring edge, because
//! exteriors are counterclockwise and holes clockwise):
//!
//! | Piece | Faces beside it |
//! |---|---|
//! | ring of B inside int A | (int, int) and (int, ext): `II`, `IE` |
//! | ring of B inside ext A | (ext, int) and (ext, ext): `EI` |
//! | ring of A inside int B | (int, int) and (ext, int): `II`, `EI` |
//! | ring of A inside ext B | (int, ext) and (ext, ext): `IE` |
//! | shared, same direction | both interiors on the left: (int, int), (ext, ext): `II` |
//! | shared, opposite direction | (int, ext) and (ext, int): `IE`, `EI` |
//!
//! `EE` is always `Area`. Nesting needs no special case: a polygon strictly
//! inside another has all its ring edges inside the other's interior, which
//! gives `II` and `IE` (and not `EI`), and a polygon inside a hole of the other
//! has its ring edges in the other's exterior, which gives `EI` and, through
//! the hole ring, `IE`.
//!
//! # Complexity
//!
//! `O(n·m)` for rings of `n` and `m` edges in total over all parts, as for
//! the other exact analyses in this crate, plus `O(v·n)` to locate the `v`
//! contact candidates. A disjoint bounding-box pair returns immediately.

use crate::exact::kernel::{on_segment, segments_properly_cross};
use crate::topology::areal::{AreaOperand, segment_window};
use crate::topology::de9im::IntersectionMatrix;
use crate::topology::label::{Dimension, Location};
use crate::topology::line_analysis::{PieceSummary, summarize_segment};
use crate::types::Coord2d;

/// Classification of every ring edge of `a` against `b`.
///
/// A ring has no boundary points of its own, so every point of it counts as
/// an interior point of the curve.
///
/// Only the edges of `a` whose box meets `b`'s envelope are classified
/// against `b`; any other edge of `a` lies wholly outside `b`'s box, so in
/// `b`'s exterior, and contributes exactly one thing: an `exterior` piece
/// (the far-edge rule). An unindexed `a` finds the same
/// edges by a scan of box comparisons.
fn ring_summary<A: AreaOperand + ?Sized, B: AreaOperand + ?Sized>(a: &A, b: &B) -> PieceSummary {
    let mut pieces = PieceSummary::default();
    let w = b.all_window();
    a.for_edges_near(&w, &mut |p, q| {
        crate::count!(ring_edge_classifications);
        summarize_segment(p, q, b, &mut pieces);
        true
    });
    if a.has_edge_outside(&w) {
        pieces.exterior = true;
    }
    pieces
}

/// Row-major index of the cell for locations `(in A, in B)`.
fn cell_index(la: Location, lb: Location) -> usize {
    let k = |l: Location| match l {
        Location::Interior => 0,
        Location::Boundary => 1,
        Location::Exterior => 2,
    };
    3 * k(la) + k(lb)
}

/// Every isolated point at which the boundaries meet, with its location in
/// A and in B:
///
/// - a vertex of A on B's rings, and a vertex of B on A's rings;
/// - a proper crossing of a ring edge of A with one of B, which lies on both
///   boundaries, unless a vertex of either lies there, which is then the
///   candidate instead.
///
/// Every candidate lies on both boundaries, so this is the `BB = Point`
/// test, taken over all parts' rings.
fn contact_cells<A: AreaOperand + ?Sized, B: AreaOperand + ?Sized>(a: &A, b: &B) -> [bool; 9] {
    let mut out = [false; 9];
    let mut put = |x: Coord2d| out[cell_index(a.locate(x), b.locate(x))] = true;
    // A vertex of one on the other's boundary lies in the other's box, so
    // only the vertices near that box are candidates.
    let (wa, wb) = (a.all_window(), b.all_window());
    let mut on_b: Vec<Coord2d> = Vec::new();
    a.for_vertices_near(&wb, &mut |v| {
        if b.locate(v) == Location::Boundary {
            on_b.push(v);
        }
        true
    });
    let mut on_a: Vec<Coord2d> = Vec::new();
    b.for_vertices_near(&wa, &mut |v| {
        if a.locate(v) == Location::Boundary {
            on_a.push(v);
        }
        true
    });
    for v in on_b.into_iter().chain(on_a) {
        put(v);
    }
    // A proper crossing of an edge of `a` with an edge of `b` lies in both
    // edges' boxes, as does any vertex on both.
    let mut crossing = false;
    a.for_edges_near(&wb, &mut |p, q| {
        let wpq = segment_window(p, q);
        b.for_edges_near(&wpq, &mut |c, d| {
            if segments_properly_cross(p, q, c, d) {
                let mut at_vertex = false;
                a.for_vertices_near(&wpq, &mut |v| {
                    at_vertex = on_segment(v, p, q) && on_segment(v, c, d);
                    !at_vertex
                });
                if !at_vertex {
                    b.for_vertices_near(&wpq, &mut |v| {
                        at_vertex = on_segment(v, p, q) && on_segment(v, c, d);
                        !at_vertex
                    });
                }
                if !at_vertex {
                    crossing = true;
                }
            }
            !crossing
        });
        !crossing
    });
    if crossing {
        out[cell_index(Location::Boundary, Location::Boundary)] = true;
    }
    out
}

/// Computes the DE-9IM matrix of areal `a` against areal `b` exactly. Two
/// polygons, or the parts of valid multipolygons, are classified in one pass
/// over all parts' rings.
///
/// See the module documentation for the derivation of the line and area
/// cells, and [`contact_cells`] for the points.
pub(crate) fn polygons_matrix<A: AreaOperand + ?Sized, B: AreaOperand + ?Sized>(
    a: &A,
    b: &B,
) -> IntersectionMatrix {
    // Strictly separated envelopes: no point is shared, so the matrix is the
    // disjoint one. Strict inequality keeps touching envelopes on the full
    // path. An empty operand never reaches here (it is a point set).
    if let (Some(ea), Some(eb)) = (a.envelope(), b.envelope())
        && (ea.2 < eb.0 || eb.2 < ea.0 || ea.3 < eb.1 || eb.3 < ea.1)
    {
        return IntersectionMatrix {
            interior_interior: Dimension::Empty,
            interior_boundary: Dimension::Empty,
            interior_exterior: Dimension::Area,
            boundary_interior: Dimension::Empty,
            boundary_boundary: Dimension::Empty,
            boundary_exterior: Dimension::Line,
            exterior_interior: Dimension::Area,
            exterior_boundary: Dimension::Line,
            exterior_exterior: Dimension::Area,
        };
    }

    let pa = ring_summary(a, b);
    let pb = ring_summary(b, a);
    let pts = contact_cells(a, b);
    let same = pa.shared_same || pb.shared_same;
    let opposite = pa.shared_opp || pb.shared_opp;

    // A cell takes the greatest dimension any of its sources gives.
    let cell = |d: Dimension, flag: bool, point: bool| {
        if flag {
            d
        } else if point {
            Dimension::Point
        } else {
            Dimension::Empty
        }
    };
    let ii = cell_index(Location::Interior, Location::Interior);
    let ib = cell_index(Location::Interior, Location::Boundary);
    let ie = cell_index(Location::Interior, Location::Exterior);
    let bi = cell_index(Location::Boundary, Location::Interior);
    let bb = cell_index(Location::Boundary, Location::Boundary);
    let be = cell_index(Location::Boundary, Location::Exterior);
    let ei = cell_index(Location::Exterior, Location::Interior);
    let eb = cell_index(Location::Exterior, Location::Boundary);

    IntersectionMatrix {
        interior_interior: cell(Dimension::Area, pa.interior || pb.interior || same, pts[ii]),
        interior_boundary: cell(Dimension::Line, pb.interior, pts[ib]),
        interior_exterior: cell(
            Dimension::Area,
            pb.interior || pa.exterior || opposite,
            pts[ie],
        ),
        boundary_interior: cell(Dimension::Line, pa.interior, pts[bi]),
        boundary_boundary: cell(Dimension::Line, pa.boundary || pb.boundary, pts[bb]),
        boundary_exterior: cell(Dimension::Line, pa.exterior, pts[be]),
        exterior_interior: cell(
            Dimension::Area,
            pa.interior || pb.exterior || opposite,
            pts[ei],
        ),
        exterior_boundary: cell(Dimension::Line, pb.exterior, pts[eb]),
        exterior_exterior: Dimension::Area,
    }
}

// ===========================================================================
// Whether two areas meet, with an early exit
// ===========================================================================

/// Whether the closed regions `a` and `b` share a point: `relate(a, b)`'s
/// `intersects`, decided without building the matrix and stopping at the
/// first witness.
///
/// Three exact searches, in order, each ending at its first hit:
///
/// 1. a vertex of `a` not in `b`'s exterior;
/// 2. a vertex of `b` not in `a`'s exterior;
/// 3. a proper crossing of an edge of `a` with an edge of `b`.
///
/// **Why these three are complete:** let `x` be a common point.
///
/// - If the boundaries meet at some point `y`, then either `y` is a vertex
///   of one operand (found by 1 or 2, since it lies on the other's
///   boundary), or it is interior to an edge of each. Two edges meeting
///   only at interior points of both either cross properly (found by 3) or
///   overlap collinearly, and a collinear overlap of two segments ends at an
///   endpoint of one of them lying on the other (found by 1 or 2).
/// - If the boundaries do not meet, every ring lies wholly in the other
///   operand's interior or wholly in its exterior (it is connected and does
///   not touch the other's boundary). Let `P` be the part of `b` and `Q` the
///   part of `a` containing `x`. If either outer ring lies in the other
///   operand's interior, its vertices are found by 1 or 2. Otherwise each
///   outer ring lies in the other's exterior. The filled outer rings both
///   contain `x`, so one encloses the other; say `Q`'s encloses `P`'s. Then
///   `P`'s outer ring, inside `Q`'s filled outer ring but outside `a`, lies in
///   a hole of `Q`, and so does all of `P` with `x`, which is then not in
///   `Q`: a contradiction. The case with the roles exchanged is the same.
///
/// So no common point escapes the three searches, and each search succeeds
/// only on a common point. Every test is `locate` or a proper-crossing test
/// on stored coordinates, both exact. A strictly separated pair of envelopes
/// returns `false` at once, as [`polygons_matrix`] does.
pub(crate) fn areas_intersect<A: AreaOperand + ?Sized, B: AreaOperand + ?Sized>(
    a: &A,
    b: &B,
) -> bool {
    if let (Some(ea), Some(eb)) = (a.envelope(), b.envelope())
        && (ea.2 < eb.0 || eb.2 < ea.0 || ea.3 < eb.1 || eb.3 < ea.1)
    {
        return false;
    }
    let (wa, wb) = (a.all_window(), b.all_window());
    let mut hit = false;
    // A vertex of `a` outside `b`'s box is in `b`'s exterior.
    a.for_vertices_near(&wb, &mut |v| {
        hit = b.locate(v) != Location::Exterior;
        !hit
    });
    if hit {
        return true;
    }
    b.for_vertices_near(&wa, &mut |v| {
        hit = a.locate(v) != Location::Exterior;
        !hit
    });
    if hit {
        return true;
    }
    // A crossing lies in both edges' boxes.
    a.for_edges_near(&wb, &mut |p, q| {
        b.for_edges_near(&segment_window(p, q), &mut |c, d| {
            hit = segments_properly_cross(p, q, c, d);
            !hit
        });
        !hit
    });
    hit
}
