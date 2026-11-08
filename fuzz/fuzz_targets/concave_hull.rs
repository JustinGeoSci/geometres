// fuzz/fuzz_targets/concave_hull.rs

//! Fuzz target: the concave hull against an exact integer brute force.
//!
//! The input is decoded into up to 64 points with integer coordinates below
//! 2²⁰ (duplicates and axis-aligned runs common), and a threshold, which is
//! either an integer edge length or a ratio in `[0, 1]`.
//!
//! Checked on every input:
//!
//! - the call succeeds (the operation's own certification never fires);
//! - a degenerate result (no Delaunay triangles) is not a polygon, and a
//!   polygon result's ring vertices are input points, bit for bit;
//! - containment, checked independently: every input point is inside or on the
//!   ring, by an `i128` crossing-number walk that shares no code with the
//!   library;
//! - every ring edge is an edge of the Delaunay triangulation.

#![no_main]

use geometres::ops::{ConcaveHullParameter, concave_hull, delaunay};
use geometres::types::{Coord2d, Geometry, MultiPoint, Point};
use libfuzzer_sys::fuzz_target;

struct Bytes<'a> {
    data: &'a [u8],
    at: usize,
}

impl Bytes<'_> {
    fn u8(&mut self) -> u8 {
        let b = self.data.get(self.at).copied().unwrap_or(0);
        self.at += 1;
        b
    }

    fn u32(&mut self) -> u32 {
        u32::from_le_bytes([self.u8(), self.u8(), self.u8(), self.u8()])
    }
}

fn decode(b: &mut Bytes<'_>) -> (Vec<(i64, i64)>, ConcaveHullParameter) {
    let parameter = if b.u8() % 2 == 0 {
        ConcaveHullParameter::MaxEdgeLength(f64::from(b.u32() % (1 << 21)))
    } else {
        ConcaveHullParameter::MaxEdgeLengthRatio(f64::from(b.u8()) / 255.0)
    };
    let n = 3 + usize::from(b.u8()) % 62;
    let mut pts: Vec<(i64, i64)> = Vec::with_capacity(n);
    for _ in 0..n {
        let kind = b.u8() % 4;
        let p = match (kind, pts.last().copied()) {
            (0, Some(_)) => pts[usize::from(b.u8()) % pts.len()],
            (1, Some((x, _))) => (x, i64::from(b.u32() % (1 << 20))),
            (2, Some((_, y))) => (i64::from(b.u32() % (1 << 20)), y),
            _ => (
                i64::from(b.u32() % (1 << 20)),
                i64::from(b.u32() % (1 << 20)),
            ),
        };
        pts.push(p);
    }
    (pts, parameter)
}

fn cross(a: (i64, i64), b: (i64, i64), p: (i64, i64)) -> i128 {
    i128::from(b.0 - a.0) * i128::from(p.1 - a.1) - i128::from(b.1 - a.1) * i128::from(p.0 - a.0)
}

fn on_segment(a: (i64, i64), b: (i64, i64), p: (i64, i64)) -> bool {
    cross(a, b, p) == 0
        && p.0 >= a.0.min(b.0)
        && p.0 <= a.0.max(b.0)
        && p.1 >= a.1.min(b.1)
        && p.1 <= a.1.max(b.1)
}

/// Inside or on the ring, by the exact crossing number: a horizontal ray
/// from `p`, counting edges that straddle its line (half-open in `y`) and
/// cross strictly to its right. The side is read from the exact turn sign.
fn contains(ring: &[(i64, i64)], p: (i64, i64)) -> bool {
    let mut inside = false;
    for k in 0..ring.len() {
        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
        if on_segment(a, b, p) {
            return true;
        }
        if (a.1 > p.1) != (b.1 > p.1) {
            let c = cross(a, b, p);
            if (b.1 > a.1 && c > 0) || (b.1 < a.1 && c < 0) {
                inside = !inside;
            }
        }
    }
    inside
}

#[expect(clippy::cast_possible_truncation, reason = "coordinates are small integers")]
fn as_int(c: Coord2d) -> (i64, i64) {
    (c.x as i64, c.y as i64)
}

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let (pts, parameter) = decode(&mut b);
    let g = MultiPoint::new(
        pts.iter()
            .map(|&(x, y)| Point::new(Coord2d::new(x as f64, y as f64)).unwrap())
            .collect(),
    );
    let hull = concave_hull(&g, parameter).unwrap();
    let t = delaunay(&g).unwrap();
    let Geometry::Polygon(poly) = &hull else {
        assert!(t.triangles().is_empty(), "degenerate result for a triangulated set");
        return;
    };
    assert!(!t.triangles().is_empty());

    let ring: Vec<(i64, i64)> = poly
        .exterior()
        .coords()
        .as_slice()
        .iter()
        .map(|&c| as_int(c))
        .collect();
    let input: std::collections::HashSet<(i64, i64)> = pts.iter().copied().collect();
    for v in &ring {
        assert!(input.contains(v), "a ring vertex is not an input point");
    }
    let mut dedges = std::collections::HashSet::new();
    for tr in t.triangles() {
        for k in 0..3 {
            let (u, w) = (
                as_int(t.points()[tr[k] as usize]),
                as_int(t.points()[tr[(k + 1) % 3] as usize]),
            );
            dedges.insert((u.min(w), u.max(w)));
        }
    }
    for k in 0..ring.len() {
        let (u, w) = (ring[k], ring[(k + 1) % ring.len()]);
        assert!(
            dedges.contains(&(u.min(w), u.max(w))),
            "a ring edge is not a Delaunay edge"
        );
    }
    for &p in &pts {
        assert!(contains(&ring, p), "{p:?} outside the hull");
    }
});
