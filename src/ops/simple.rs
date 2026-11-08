// src/ops/simple.rs

//! OGC simplicity for every geometry type: [`is_simple`].
//!
//! # The semantics, per type
//!
//! | Type | Simple when |
//! |---|---|
//! | `Point` | always |
//! | `MultiPoint` | no point repeats (coordinates compared by value) |
//! | `LineString` | it passes through no point twice, except that a closed line's first and last vertex coincide |
//! | `MultiLineString` | each element is simple, and two elements meet only at points on the boundary of both: an endpoint of each, neither element closed (a closed line has no boundary under mod-2) |
//! | `Polygon`, `MultiPolygon` | every ring is simple (OGC reads polygon simplicity as ring simplicity) |
//! | `GeometryCollection` | every element is simple; elements are not compared with each other |
//!
//! An empty geometry is simple.
//!
//! # Exactness
//!
//! Every decision is `arithmetica::orient2d` on stored coordinates or an
//! exact comparison of stored coordinates: proper crossings by four
//! orientation signs, every touching configuration by the closed
//! point-on-segment test, repeats and endpoints by value. Two segments that
//! touch without crossing properly meet at an endpoint of one of them, so
//! every contact is located at an input coordinate; none is constructed.

use crate::exact::kernel::{first_self_contact, on_segment, segments_properly_cross};
use crate::predicates::{Borrowed, Relatable, borrowed};
use crate::types::{Coord2d, Geometry, LineString, MultiLineString, MultiPoint, Polygon};

/// `true` if `g` is simple in the OGC sense; the module documentation gives
/// the rule for each type.
///
/// A collection is walked without recursion; each element is judged on its
/// own.
#[must_use]
pub fn is_simple(g: &impl Relatable) -> bool {
    match borrowed(g) {
        Borrowed::Point(_) => true,
        Borrowed::MultiPoint(m) => multi_point_is_simple(m),
        Borrowed::LineString(l) => line_is_simple(l),
        Borrowed::MultiLineString(m) => multi_line_is_simple(m),
        Borrowed::Polygon(p) => polygon_is_simple(p),
        Borrowed::MultiPolygon(m) => m.polygons().iter().all(polygon_is_simple),
        Borrowed::GeometryCollection(c) => c.leaves().all(leaf_is_simple),
    }
}

fn leaf_is_simple(g: &Geometry) -> bool {
    match g {
        Geometry::Point(_) => true,
        Geometry::MultiPoint(m) => multi_point_is_simple(m),
        Geometry::LineString(l) => line_is_simple(l),
        Geometry::MultiLineString(m) => multi_line_is_simple(m),
        Geometry::Polygon(p) => polygon_is_simple(p),
        Geometry::MultiPolygon(m) => m.polygons().iter().all(polygon_is_simple),
        // `leaves` never yields a collection.
        Geometry::GeometryCollection(c) => c.leaves().all(leaf_is_simple),
    }
}

/// No coordinate repeats. Sorting by total order would split `-0.0` from
/// `0.0`, so the comparison is pairwise by value.
fn multi_point_is_simple(m: &MultiPoint) -> bool {
    let pts: Vec<Coord2d> = m.points().iter().map(|p| p.coord()).collect();
    pts.iter()
        .enumerate()
        .all(|(i, p)| !pts[i + 1..].contains(p))
}

/// Every ring of the polygon is simple: the exact self-contact test the
/// vouches use (a crossing, a vertex on a non-adjacent edge, a vertex
/// visited twice, or a fold-back all make a ring non-simple).
fn polygon_is_simple(p: &Polygon) -> bool {
    core::iter::once(p.exterior())
        .chain(p.holes())
        .all(|r| first_self_contact(r.coords().as_slice()).is_none())
}

/// `true` if the line's first and last vertices coincide.
fn is_closed(v: &[Coord2d]) -> bool {
    v.len() > 2 && v.first() == v.last()
}

/// No self-contact except the closing one: the module rule for a line.
fn line_is_simple(l: &LineString) -> bool {
    let v = l.coords().as_slice();
    let m = v.len() - 1; // edges: i runs v[i] → v[i + 1]
    let closed = is_closed(v);
    for i in 0..m {
        let (a, b) = (v[i], v[i + 1]);
        // The next edge shares `b`; any other common point is a fold-back.
        if i + 1 < m && (on_segment(v[i + 2], a, b) || on_segment(a, b, v[i + 2])) {
            return false;
        }
        for j in (i + 2)..m {
            let (c, d) = (v[j], v[j + 1]);
            if closed && i == 0 && j == m - 1 {
                // These share the closing vertex `a == d`; only a fold-back
                // across it is a further contact.
                if on_segment(b, c, d) || on_segment(c, a, b) {
                    return false;
                }
                continue;
            }
            if segments_touch(a, b, c, d) {
                return false;
            }
        }
    }
    true
}

/// Any common point of the closed segments `ab` and `cd`.
fn segments_touch(a: Coord2d, b: Coord2d, c: Coord2d, d: Coord2d) -> bool {
    segments_properly_cross(a, b, c, d)
        || on_segment(a, c, d)
        || on_segment(b, c, d)
        || on_segment(c, a, b)
        || on_segment(d, a, b)
}

/// Each element simple, and elements meeting only at endpoints of both, with
/// neither closed.
fn multi_line_is_simple(m: &MultiLineString) -> bool {
    let lines = m.line_strings();
    if !lines.iter().all(line_is_simple) {
        return false;
    }
    let boxes: Vec<(f64, f64, f64, f64)> = lines.iter().map(line_box).collect();
    for i in 0..lines.len() {
        for j in (i + 1)..lines.len() {
            let (p, q) = (boxes[i], boxes[j]);
            if p.0 <= q.2
                && q.0 <= p.2
                && p.1 <= q.3
                && q.1 <= p.3
                && !meet_only_at_boundaries(&lines[i], &lines[j])
            {
                return false;
            }
        }
    }
    true
}

/// `true` if every common point of two simple lines is an endpoint of each
/// and neither is closed.
fn meet_only_at_boundaries(x: &LineString, y: &LineString) -> bool {
    let (u, w) = (x.coords().as_slice(), y.coords().as_slice());
    let boundary = |v: &[Coord2d], p: Coord2d| -> bool {
        !is_closed(v) && (v.first() == Some(&p) || v.last() == Some(&p))
    };
    for e in u.windows(2) {
        let (a, b) = (e[0], e[1]);
        for f in w.windows(2) {
            let (c, d) = (f[0], f[1]);
            if segments_properly_cross(a, b, c, d) {
                return false;
            }
            // Without a proper crossing, every common point is one of the
            // four endpoints lying on the other segment.
            let mut common: Vec<Coord2d> = Vec::new();
            for (p, s, t) in [(a, c, d), (b, c, d), (c, a, b), (d, a, b)] {
                if on_segment(p, s, t) && !common.contains(&p) {
                    common.push(p);
                }
            }
            // Two distinct common points: the segments overlap along a
            // stretch.
            if common.len() >= 2 {
                return false;
            }
            if let Some(&p) = common.first()
                && !(boundary(u, p) && boundary(w, p))
            {
                return false;
            }
        }
    }
    true
}

fn line_box(l: &LineString) -> (f64, f64, f64, f64) {
    l.bbox()
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used)]
// test-lint-posture: see the workspace lint comment in Cargo.toml
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::types::{CoordSeq, GeometryCollection, LinearRing, MultiPolygon, Point};

    fn line(pts: &[(f64, f64)]) -> LineString {
        LineString::new(CoordSeq::new(
            pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
        ))
        .unwrap()
    }

    fn ring(pts: &[(f64, f64)]) -> LinearRing {
        LinearRing::new(CoordSeq::new(
            pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
        ))
        .unwrap()
    }

    fn pt(x: f64, y: f64) -> Point {
        Point::new(Coord2d::new(x, y)).unwrap()
    }

    #[test]
    fn points_and_multipoints() {
        assert!(is_simple(&pt(1.0, 2.0)));
        assert!(is_simple(&MultiPoint::empty()));
        assert!(is_simple(&MultiPoint::new(vec![
            pt(0.0, 0.0),
            pt(1.0, 0.0)
        ])));
        assert!(!is_simple(&MultiPoint::new(vec![
            pt(0.0, 0.0),
            pt(1.0, 0.0),
            pt(0.0, 0.0)
        ])));
        // -0.0 and 0.0 are the same point.
        assert!(!is_simple(&MultiPoint::new(vec![
            pt(0.0, 1.0),
            pt(-0.0, 1.0)
        ])));
    }

    #[test]
    fn lines_open_and_closed() {
        assert!(is_simple(&line(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])));
        // A closed line: start and end coincide, nothing else touches.
        assert!(is_simple(&line(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0)
        ])));
        // A crossing.
        assert!(!is_simple(&line(&[
            (0.0, 0.0),
            (2.0, 2.0),
            (2.0, 0.0),
            (0.0, 2.0)
        ])));
        // Touching itself at an interior vertex (passes through (1,1) twice).
        assert!(!is_simple(&line(&[
            (0.0, 0.0),
            (1.0, 1.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (1.0, 1.0),
            (0.0, 2.0)
        ])));
        // An endpoint on its own interior (a "6").
        assert!(!is_simple(&line(&[
            (0.0, 0.0),
            (0.0, 4.0),
            (2.0, 4.0),
            (2.0, 2.0),
            (0.0, 2.0)
        ])));
        // A fold-back.
        assert!(!is_simple(&line(&[(0.0, 0.0), (2.0, 0.0), (1.0, 0.0)])));
        // A closed line that also touches itself elsewhere.
        assert!(!is_simple(&line(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (1.0, 1.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (1.0, 1.0),
            (0.0, 0.0)
        ])));
        // Closed, folding back across the closing vertex.
        assert!(!is_simple(&line(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (1.0, 0.0),
            (0.0, 0.0)
        ])));
        // A collinear continuation is not a contact.
        assert!(is_simple(&line(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])));
    }

    #[test]
    fn multilines_under_the_mod_2_rule() {
        let a = line(&[(0.0, 0.0), (2.0, 0.0)]);
        let b = line(&[(2.0, 0.0), (4.0, 1.0)]);
        // Endpoint to endpoint: on both boundaries, simple.
        assert!(is_simple(&MultiLineString::new(vec![a.clone(), b.clone()])));
        // Three lines meeting at one shared endpoint: still simple.
        let c = line(&[(2.0, 0.0), (2.0, -3.0)]);
        assert!(is_simple(&MultiLineString::new(vec![
            a.clone(),
            b.clone(),
            c
        ])));
        // An endpoint on another element's interior: not simple.
        let t = line(&[(1.0, 0.0), (1.0, 2.0)]);
        assert!(!is_simple(&MultiLineString::new(vec![a.clone(), t])));
        // Crossing elements.
        let x = line(&[(1.0, -1.0), (1.0, 1.0)]);
        assert!(!is_simple(&MultiLineString::new(vec![a.clone(), x])));
        // Overlapping along a stretch.
        let o = line(&[(1.0, 0.0), (3.0, 0.0)]);
        assert!(!is_simple(&MultiLineString::new(vec![a.clone(), o])));
        // Meeting a closed element at its closing point: a closed line has no
        // boundary under mod-2, so the contact is interior to it.
        let closed = line(&[(2.0, 0.0), (3.0, 1.0), (3.0, -1.0), (2.0, 0.0)]);
        assert!(!is_simple(&MultiLineString::new(vec![a.clone(), closed])));
        // A non-simple element makes the whole non-simple.
        let bow = line(&[(10.0, 0.0), (12.0, 2.0), (12.0, 0.0), (10.0, 2.0)]);
        assert!(!is_simple(&MultiLineString::new(vec![a, bow])));
        assert!(is_simple(&MultiLineString::empty()));
    }

    #[test]
    fn polygons_by_their_rings() {
        let square =
            Polygon::new_simple(ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])).unwrap();
        assert!(is_simple(&square));
        // A ring with a pinch at (2,2): constructible, not simple.
        let pinch = Polygon::new_simple(ring(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (2.0, 2.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (2.0, 2.0),
        ]))
        .unwrap();
        assert!(!is_simple(&pinch));
        // A hole touching the exterior: each ring is simple, so the polygon
        // is (touches between rings are a validity matter, not simplicity).
        let touching = Polygon::new(
            ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            vec![ring(&[(0.0, 2.0), (2.0, 3.0), (2.0, 1.0)])],
        )
        .unwrap();
        assert!(is_simple(&touching));
        assert!(!is_simple(&MultiPolygon::new(vec![pinch.clone()]).unwrap()));
        assert!(is_simple(&MultiPolygon::new(vec![square.clone()]).unwrap()));

        let gc = GeometryCollection::new(vec![
            Geometry::Polygon(square),
            Geometry::LineString(line(&[(0.0, 0.0), (9.0, 9.0)])),
            // Elements are not compared: this point lies on the line.
            Geometry::Point(pt(1.0, 1.0)),
        ])
        .unwrap();
        assert!(is_simple(&gc));
        let nested = GeometryCollection::new(vec![
            Geometry::GeometryCollection(gc),
            Geometry::Polygon(pinch),
        ])
        .unwrap();
        assert!(!is_simple(&nested));
        assert!(is_simple(&GeometryCollection::empty()));
    }
}
