// src/ops/basic.rs

//! The basic operations that rearrange a geometry's own coordinates:
//! [`boundary`], [`normalize`], [`reverse`] and [`unique_points`].
//!
//! # Every output coordinate is an input coordinate
//!
//! None of these computes a coordinate. `boundary` returns ring vertices and
//! line endpoints; `normalize` and `reverse` reorder the input; and
//! `unique_points` drops repeats. Every decision is a comparison of stored
//! `f64` values (`==`, and `<` in the lexicographic order below) or the exact
//! location the predicates use, so nothing is approximate.
//!
//! # The lexicographic order
//!
//! Points are ordered by `x`, then `y`, comparing the stored values exactly.
//! `-0.0` and `0.0` compare equal; where that leaves two keys equal, the
//! next key decides, and then the original position (every sort here is
//! stable), so every result is deterministic.

use core::cmp::Ordering;
use std::collections::HashSet;

use crate::GeometryError;
use crate::predicates::{Borrowed, Relatable, borrowed};
use crate::topology::line_analysis::LineSet;
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString,
    MultiPoint, MultiPolygon, Point, Polygon,
};

// ===========================================================================
// The order
// ===========================================================================

/// The exact lexicographic order of two points, `x` then `y`.
pub(crate) fn lex(a: Coord2d, b: Coord2d) -> Ordering {
    a.x.partial_cmp(&b.x)
        .unwrap_or(Ordering::Equal)
        .then(a.y.partial_cmp(&b.y).unwrap_or(Ordering::Equal))
}

/// The lexicographic order of two point sequences under [`lex`].
pub(crate) fn lex_seq(a: &[Coord2d], b: &[Coord2d]) -> Ordering {
    for (p, q) in a.iter().zip(b) {
        let o = lex(*p, *q);
        if o != Ordering::Equal {
            return o;
        }
    }
    a.len().cmp(&b.len())
}

/// The key `-0.0` and `0.0` share, for hashing points by `==`.
pub(crate) fn key(c: Coord2d) -> (u64, u64) {
    ((c.x + 0.0).to_bits(), (c.y + 0.0).to_bits())
}

/// Every coordinate of a geometry, in its order: rings exterior then holes
/// (without a closing vertex), parts in order, a collection's leaves depth
/// first.
pub(crate) fn coords_of(g: Borrowed<'_>) -> Vec<Coord2d> {
    let mut out = Vec::new();
    push_coords(g, &mut out);
    out
}

fn push_polygon(p: &Polygon, out: &mut Vec<Coord2d>) {
    out.extend_from_slice(p.exterior().coords().as_slice());
    for h in p.holes() {
        out.extend_from_slice(h.coords().as_slice());
    }
}

fn push_coords(g: Borrowed<'_>, out: &mut Vec<Coord2d>) {
    match g {
        Borrowed::Point(p) => out.push(p.coord()),
        Borrowed::MultiPoint(m) => out.extend(m.points().iter().map(Point::coord)),
        Borrowed::LineString(l) => out.extend_from_slice(l.coords().as_slice()),
        Borrowed::MultiLineString(m) => {
            for l in m.line_strings() {
                out.extend_from_slice(l.coords().as_slice());
            }
        }
        Borrowed::Polygon(p) => push_polygon(p, out),
        Borrowed::MultiPolygon(m) => {
            for p in m.polygons() {
                push_polygon(p, out);
            }
        }
        Borrowed::GeometryCollection(gc) => {
            for leaf in gc.leaves() {
                push_coords(Borrowed::from_geometry(leaf), out);
            }
        }
    }
}

// ===========================================================================
// boundary
// ===========================================================================

/// A ring as a closed line string: its vertices, then the first again.
fn closed(ring: &LinearRing) -> LineString {
    let v = ring.coords().as_slice();
    let mut c = Vec::with_capacity(v.len() + 1);
    c.extend_from_slice(v);
    c.push(v[0]);
    // A ring has at least three distinct vertices and no two consecutive
    // equal, including across the closing edge, so the closed sequence has
    // no two consecutive equal either.
    LineString::from_valid(CoordSeq::new(c))
}

fn points(pts: &[Coord2d]) -> MultiPoint {
    MultiPoint::new(pts.iter().map(|&c| Point::from_finite(c)).collect())
}

fn rings_of(polys: &[&Polygon]) -> Vec<LineString> {
    polys
        .iter()
        .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes().iter()))
        .map(closed)
        .collect()
}

/// The topological boundary of `g`, as `relate` defines it.
///
/// - **`Point`, `MultiPoint`:** the empty `MultiPoint`; a point has no
///   boundary.
/// - **`LineString`, `MultiLineString`:** the OGC mod-2 boundary (the points
///   that are an endpoint of an odd number of line strings) as a
///   `MultiPoint` in exact lexicographic order. A closed line string's
///   endpoints cancel, so its boundary is empty. This is the boundary
///   `relate` uses, computed by the same function.
/// - **`Polygon`, `MultiPolygon`:** every ring as a closed `LineString`
///   (its first vertex repeated last), exterior then holes, part by part, as
///   a `MultiLineString`. Rings keep their orientation.
/// - **`GeometryCollection`:** the boundary `relate` locates for it, made of
///   every ring of its polygonal leaves, and the mod-2 endpoints of its line
///   leaves that lie outside the polygonal leaves (an endpoint inside them
///   is interior to the collection; one on a ring is already on the ring).
///   A `MultiLineString` of the rings when there are only rings, a
///   `MultiPoint` when there are only points, a `GeometryCollection` of the
///   two (rings first) when there are both, and the empty `MultiPoint` when
///   the boundary is empty.
///
/// Every output coordinate is a ring vertex or a line endpoint of the input.
///
/// # Errors
///
/// [`GeometryError::CollectionElementsOverlap`] for a collection whose
/// polygonal leaves overlap, which `relate` refuses for the same reason: its
/// boundary would be the boundary of their union, which this function does
/// not construct.
pub fn boundary(g: &Geometry) -> Result<Geometry, GeometryError> {
    Ok(match g {
        Geometry::Point(_) | Geometry::MultiPoint(_) => Geometry::MultiPoint(MultiPoint::empty()),
        Geometry::LineString(l) => Geometry::MultiPoint(points(LineSet::single(l).boundary())),
        Geometry::MultiLineString(m) => {
            let set = LineSet::new(m.line_strings().iter().collect());
            Geometry::MultiPoint(points(set.boundary()))
        }
        Geometry::Polygon(p) => Geometry::MultiLineString(MultiLineString::new(rings_of(&[p]))),
        Geometry::MultiPolygon(m) => {
            let parts: Vec<&Polygon> = m.polygons().iter().collect();
            Geometry::MultiLineString(MultiLineString::new(rings_of(&parts)))
        }
        Geometry::GeometryCollection(gc) => {
            let (rings, pts) = crate::topology::collection::collection_boundary(gc)?;
            let lines: Vec<LineString> = rings.into_iter().map(closed).collect();
            match (lines.is_empty(), pts.is_empty()) {
                (true, true) => Geometry::MultiPoint(MultiPoint::empty()),
                (false, true) => Geometry::MultiLineString(MultiLineString::new(lines)),
                (true, false) => Geometry::MultiPoint(points(&pts)),
                (false, false) => Geometry::GeometryCollection(GeometryCollection::from_valid(
                    vec![
                        Geometry::MultiLineString(MultiLineString::new(lines)),
                        Geometry::MultiPoint(points(&pts)),
                    ],
                    1,
                )),
            }
        }
    })
}

// ===========================================================================
// normalize
// ===========================================================================

/// The rotation of a ring that starts at its smallest vertex, keeping its
/// orientation. Where the smallest vertex occurs more than once (a ring that
/// touches itself there), the rotation whose whole sequence is smallest is
/// taken, and among equal ones the first.
fn normalize_ring(r: &LinearRing) -> LinearRing {
    let v = r.coords().as_slice();
    let n = v.len();
    let rotation = |k: usize| v[k..].iter().chain(&v[..k]).copied().collect::<Vec<_>>();
    let mut best = 0;
    let mut best_seq = rotation(0);
    for k in 1..n {
        if lex(v[k], v[best]) == Ordering::Greater {
            continue;
        }
        let seq = rotation(k);
        if lex_seq(&seq, &best_seq) == Ordering::Less {
            best = k;
            best_seq = seq;
        }
    }
    LinearRing::from_valid(CoordSeq::new(best_seq))
}

fn normalize_polygon(p: &Polygon) -> Polygon {
    let exterior = normalize_ring(p.exterior());
    let mut holes: Vec<LinearRing> = p.holes().iter().map(normalize_ring).collect();
    holes.sort_by(|a, b| lex_seq(a.coords().as_slice(), b.coords().as_slice()));
    Polygon::from_valid(exterior, holes)
}

fn normalize_line(l: &LineString) -> LineString {
    let v = l.coords().as_slice();
    if v.len() >= 4 && v[0] == v[v.len() - 1] {
        return normalize_closed_line(v);
    }
    let rev: Vec<Coord2d> = v.iter().rev().copied().collect();
    if lex_seq(v, &rev) == Ordering::Greater {
        LineString::from_valid(CoordSeq::new(rev))
    } else {
        l.clone()
    }
}

/// A closed line string's canonical form: its start is not part of the
/// geometry (a closed line has no boundary), so, as for a ring, it starts at
/// its smallest vertex; and, having no orientation to keep, it runs in
/// whichever direction makes its sequence smaller.
fn normalize_closed_line(v: &[Coord2d]) -> LineString {
    let body = &v[..v.len() - 1];
    let n = body.len();
    let close = |mut seq: Vec<Coord2d>| {
        seq.push(seq[0]);
        seq
    };
    let min = body
        .iter()
        .copied()
        .min_by(|a, b| lex(*a, *b))
        .unwrap_or(body[0]);
    let mut best: Option<Vec<Coord2d>> = None;
    for k in 0..n {
        if lex(body[k], min) != Ordering::Equal {
            continue;
        }
        let fwd: Vec<Coord2d> = body[k..].iter().chain(&body[..k]).copied().collect();
        let mut bwd: Vec<Coord2d> = vec![fwd[0]];
        bwd.extend(fwd[1..].iter().rev());
        for cand in [fwd, bwd] {
            let cand = close(cand);
            if best
                .as_ref()
                .is_none_or(|b| lex_seq(&cand, b) == Ordering::Less)
            {
                best = Some(cand);
            }
        }
    }
    // A rotation and possibly a reversal of a valid closed line: the same
    // vertices with the same cyclic neighbors, so still no two consecutive
    // equal, closed, at least four coordinates.
    LineString::from_valid(CoordSeq::new(best.unwrap_or_else(|| v.to_vec())))
}

fn normalize_multi_point(m: &MultiPoint) -> MultiPoint {
    let mut pts = m.points().to_vec();
    pts.sort_by(|a, b| lex(a.coord(), b.coord()));
    MultiPoint::new(pts)
}

fn normalize_multi_line(m: &MultiLineString) -> MultiLineString {
    let mut ls: Vec<LineString> = m.line_strings().iter().map(normalize_line).collect();
    ls.sort_by(|a, b| lex_seq(a.coords().as_slice(), b.coords().as_slice()));
    MultiLineString::new(ls)
}

fn normalize_multi_polygon(m: &MultiPolygon) -> MultiPolygon {
    let mut ps: Vec<Polygon> = m.polygons().iter().map(normalize_polygon).collect();
    ps.sort_by(|a, b| {
        lex_seq(
            a.exterior().coords().as_slice(),
            b.exterior().coords().as_slice(),
        )
    });
    MultiPolygon::from_valid(ps)
}

/// The position of a geometry type in the codec's tag order.
fn type_rank(g: &Geometry) -> u8 {
    use crate::codec::tag;
    match g {
        Geometry::Point(_) => tag::POINT,
        Geometry::LineString(_) => tag::LINE_STRING,
        Geometry::Polygon(_) => tag::POLYGON,
        Geometry::MultiPoint(_) => tag::MULTI_POINT,
        Geometry::MultiLineString(_) => tag::MULTI_LINE_STRING,
        Geometry::MultiPolygon(_) => tag::MULTI_POLYGON,
        Geometry::GeometryCollection(_) => tag::GEOMETRY_COLLECTION,
    }
}

fn normalize_collection(gc: &GeometryCollection) -> GeometryCollection {
    let gs: Vec<Geometry> = gc.geometries().iter().map(normalize).collect();
    // Sort keys computed once: the type's rank, then the whole coordinate
    // sequence (so the first coordinate decides first, as for parts).
    let mut keyed: Vec<(u8, Vec<Coord2d>, Geometry)> = gs
        .into_iter()
        .map(|g| (type_rank(&g), coords_of(Borrowed::from_geometry(&g)), g))
        .collect();
    keyed.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| lex_seq(&a.1, &b.1)));
    GeometryCollection::from_valid(keyed.into_iter().map(|k| k.2).collect(), gc.depth())
}

/// The canonical form of `g`. Equal geometries written differently (rings
/// starting at different vertices, parts or holes in different orders, line
/// strings in either direction) normalize to the same value, byte for byte
/// in the codec.
///
/// - **Rings** start at their lexicographically smallest vertex and keep
///   their orientation: exteriors stay counterclockwise and holes clockwise,
///   the type invariant.
/// - **Holes** are sorted by their rotated vertex sequence.
/// - **Parts** of a multi-geometry are sorted by their (normalized) vertex
///   sequence, a polygon by its exterior.
/// - **A `LineString`** is reversed if its coordinate sequence compares
///   greater than its reverse. A **closed** one, which has no boundary and
///   so no meaningful start, starts at its smallest vertex, as a ring does,
///   and runs in whichever direction makes its sequence smaller.
/// - **A `MultiPoint`'s** points are sorted.
/// - **A `GeometryCollection`'s** elements are normalized, then sorted by
///   type (codec tag order), then by their whole coordinate sequence, so by
///   first coordinate first (an empty element before every other of its
///   type).
///
/// The order is exact (the module documentation): ties under `==`, which
/// include `-0.0` against `0.0`, fall to the next key and then to the
/// original position. `normalize` is idempotent. Every coordinate of the
/// result is a coordinate of `g`, bit for bit.
#[must_use]
pub fn normalize(g: &Geometry) -> Geometry {
    match g {
        Geometry::Point(p) => Geometry::Point(*p),
        Geometry::MultiPoint(m) => Geometry::MultiPoint(normalize_multi_point(m)),
        Geometry::LineString(l) => Geometry::LineString(normalize_line(l)),
        Geometry::MultiLineString(m) => Geometry::MultiLineString(normalize_multi_line(m)),
        Geometry::Polygon(p) => Geometry::Polygon(normalize_polygon(p)),
        Geometry::MultiPolygon(m) => Geometry::MultiPolygon(normalize_multi_polygon(m)),
        Geometry::GeometryCollection(gc) => Geometry::GeometryCollection(normalize_collection(gc)),
    }
}

// ===========================================================================
// reverse
// ===========================================================================

/// A geometry whose direction can be reversed without leaving its type.
///
/// Implemented for the lineal types, whose direction is part of the value,
/// and for the point types, which have none (reversal is the identity). Not
/// for polygons: a reversed ring would invert its winding, which the type's
/// invariant forbids.
pub trait Reversible: Sized {
    /// `self` with its direction reversed.
    #[must_use]
    fn reversed(&self) -> Self;
}

impl Reversible for LineString {
    fn reversed(&self) -> Self {
        // The same coordinates in reverse order: still finite, at least two,
        // no two consecutive equal.
        LineString::from_valid(CoordSeq::new(
            self.coords().as_slice().iter().rev().copied().collect(),
        ))
    }
}

impl Reversible for MultiLineString {
    /// Each line string reversed, in the same order.
    fn reversed(&self) -> Self {
        MultiLineString::new(
            self.line_strings()
                .iter()
                .map(Reversible::reversed)
                .collect(),
        )
    }
}

impl Reversible for Point {
    fn reversed(&self) -> Self {
        *self
    }
}

impl Reversible for MultiPoint {
    fn reversed(&self) -> Self {
        self.clone()
    }
}

/// `g` with its direction reversed: a line string's coordinates in reverse
/// order; a multi-line-string's components each reversed, in the same order;
/// a point unchanged. See [`Reversible`] for why polygons are not offered.
#[must_use]
pub fn reverse<T: Reversible>(g: &T) -> T {
    g.reversed()
}

// ===========================================================================
// unique_points
// ===========================================================================

/// Every distinct coordinate of `g`, in order of first occurrence (rings
/// exterior then holes, without a closing vertex; parts in order; a
/// collection's leaves depth first).
///
/// Distinct means unequal under `==`, so `-0.0` and `0.0` count as one, and
/// the first occurrence's bits are kept.
#[must_use]
pub fn unique_points(g: &impl Relatable) -> MultiPoint {
    let mut seen: HashSet<(u64, u64)> = HashSet::new();
    let pts = coords_of(borrowed(g))
        .into_iter()
        .filter(|&c| seen.insert(key(c)))
        .map(Point::from_finite)
        .collect();
    MultiPoint::new(pts)
}
