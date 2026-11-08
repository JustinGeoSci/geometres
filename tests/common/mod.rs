// tests/common/mod.rs

//! Generated geometry for the integration tests.
//!
//! Every generator is deterministic: a fixed seed gives the same geometry
//! on every run and every platform. Coordinates are small integers, so
//! generated geometries share vertices, run along each other's edges and
//! touch at corners far more often than random floating-point input would,
//! and every coordinate, every translation by a power of two and every
//! scaling by a power of two is exact. The near-degenerate cases (slivers,
//! points one unit in the last place off an edge) are built where a test
//! needs them.

#![allow(
    dead_code,
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]

use geometres::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString,
    MultiPoint, MultiPolygon, Point, Polygon,
};

// ===========================================================================
// The generator
// ===========================================================================

/// xorshift64*: small, fast and the same everywhere.
pub struct Rng(u64);

impl Rng {
    /// A generator from a seed; a zero seed is replaced, since xorshift
    /// stays at zero.
    pub fn new(seed: u64) -> Self {
        Self(if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        })
    }

    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// A value in `0..n`.
    pub fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    /// An integer in `lo..=hi`, as `f64`.
    pub fn int(&mut self, lo: i64, hi: i64) -> f64 {
        (lo + self.below((hi - lo + 1) as u64) as i64) as f64
    }

    /// `true` with probability one in `n`.
    pub fn one_in(&mut self, n: u64) -> bool {
        self.below(n) == 0
    }

    /// The elements of `v` in a random order (Fisher–Yates).
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            v.swap(i, self.below(i as u64 + 1) as usize);
        }
    }
}

// ===========================================================================
// Builders
// ===========================================================================

pub fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

pub fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

pub fn line(pts: &[(f64, f64)]) -> LineString {
    LineString::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
}

pub fn point(x: f64, y: f64) -> Point {
    Point::new(c(x, y)).unwrap()
}

/// The axis-aligned rectangle `[x0, x1] × [y0, y1]`, counterclockwise.
pub fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new_simple(ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])).unwrap()
}

/// The hole `[x0, x1] × [y0, y1]`, clockwise.
pub fn hole(x0: f64, y0: f64, x1: f64, y1: f64) -> LinearRing {
    ring(&[(x0, y0), (x0, y1), (x1, y1), (x1, y0)])
}

// ===========================================================================
// Random geometry
// ===========================================================================

/// A jagged polygon with its lower left corner at `(ox, oy)`: a straight
/// bottom edge and a top made of `n` random heights in `3..=h`, so the
/// ring is simple (each top vertex is above the bottom edge and the top is
/// monotone in `x`). One time in two it has a rectangular hole near the
/// bottom, which every height clears.
pub fn jagged(rng: &mut Rng, ox: f64, oy: f64, n: usize, h: i64) -> Polygon {
    let w = n as f64;
    let mut pts = vec![(ox, oy), (ox + w, oy)];
    for i in (0..=n).rev() {
        pts.push((ox + i as f64, oy + rng.int(3, h)));
    }
    let ext = ring(&pts);
    let holes = if n >= 3 && rng.one_in(2) {
        let a = ox + rng.int(1, (n - 2) as i64);
        vec![hole(a, oy + 1.0, a + 1.0, oy + 2.0)]
    } else {
        vec![]
    };
    Polygon::new(ext, holes).unwrap()
}

/// A random walk of `n` vertices on the integer lattice from `(ox, oy)`,
/// with steps in `-3..=3` and no repeated consecutive vertex. It may cross
/// itself, which a `LineString` allows.
pub fn walk(rng: &mut Rng, ox: f64, oy: f64, n: usize) -> LineString {
    let mut pts = vec![(ox, oy)];
    while pts.len() < n {
        let (x, y) = pts[pts.len() - 1];
        let next = (x + rng.int(-3, 3), y + rng.int(-3, 3));
        if next != (x, y) {
            pts.push(next);
        }
    }
    line(&pts)
}

/// `n` lattice points in `[lo, hi]²`, repeats allowed.
pub fn points(rng: &mut Rng, n: usize, lo: i64, hi: i64) -> MultiPoint {
    MultiPoint::new(
        (0..n)
            .map(|_| point(rng.int(lo, hi), rng.int(lo, hi)))
            .collect(),
    )
}

/// One random geometry of the given kind (0 to 6), placed in `[0, 24]²`.
pub fn geometry(rng: &mut Rng, kind: u64) -> Geometry {
    let o = |rng: &mut Rng| (rng.int(0, 14), rng.int(0, 14));
    match kind {
        0 => Geometry::Point(point(rng.int(0, 24), rng.int(0, 24))),
        1 => {
            let n = 1 + rng.below(6) as usize;
            Geometry::MultiPoint(points(rng, n, 0, 24))
        }
        2 => {
            let (x, y) = o(rng);
            let n = 2 + rng.below(6) as usize;
            Geometry::LineString(walk(rng, x, y, n))
        }
        3 => {
            let k = 1 + rng.below(3) as usize;
            Geometry::MultiLineString(MultiLineString::new(
                (0..k)
                    .map(|_| {
                        let (x, y) = o(rng);
                        let n = 2 + rng.below(4) as usize;
                        walk(rng, x, y, n)
                    })
                    .collect(),
            ))
        }
        4 => {
            let (x, y) = o(rng);
            let n = 2 + rng.below(7) as usize;
            Geometry::Polygon(jagged(rng, x, y, n, 9))
        }
        5 => {
            // Two jagged parts side by side with a gap of at least one, so
            // they are disjoint and the multipolygon is valid.
            let (x, y) = o(rng);
            let n1 = 2 + rng.below(4) as usize;
            let a = jagged(rng, x, y, n1, 8);
            let (bx, by) = (x + n1 as f64 + 1.0 + rng.int(0, 2), y + rng.int(-2, 2));
            let n2 = 2 + rng.below(4) as usize;
            let b = jagged(rng, bx, by, n2, 8);
            Geometry::MultiPolygon(MultiPolygon::new(vec![a, b]).unwrap())
        }
        _ => {
            // A collection of one leaf of each dimension; one polygonal
            // leaf, so it is always an accepted operand.
            let (x, y) = o(rng);
            GeometryCollection::new(vec![
                Geometry::Point(point(rng.int(0, 24), rng.int(0, 24))),
                Geometry::LineString(walk(rng, x, y, 3)),
                Geometry::Polygon(jagged(rng, y, x, 3, 6)),
            ])
            .map(Geometry::GeometryCollection)
            .unwrap()
        }
    }
}

/// `n` random geometries, every kind in turn.
pub fn geometries(seed: u64, n: usize) -> Vec<Geometry> {
    let mut rng = Rng::new(seed);
    (0..n).map(|i| geometry(&mut rng, i as u64 % 7)).collect()
}

// ===========================================================================
// Coordinate maps
// ===========================================================================

/// `g` with every coordinate replaced by `f` of it. The map must preserve
/// orientation and validity (a translation, or a scaling by a positive
/// factor), so every constructor accepts the result.
pub fn map(g: &Geometry, f: &impl Fn(Coord2d) -> Coord2d) -> Geometry {
    let seq = |s: &CoordSeq| CoordSeq::new(s.iter().map(|&p| f(p)).collect());
    let ls = |l: &LineString| LineString::new(seq(l.coords())).unwrap();
    let lr = |r: &LinearRing| LinearRing::new(seq(r.coords())).unwrap();
    let pg =
        |p: &Polygon| Polygon::new(lr(p.exterior()), p.holes().iter().map(lr).collect()).unwrap();
    match g {
        Geometry::Point(p) => Geometry::Point(Point::new(f(p.coord())).unwrap()),
        Geometry::MultiPoint(m) => Geometry::MultiPoint(MultiPoint::new(
            m.points()
                .iter()
                .map(|p| Point::new(f(p.coord())).unwrap())
                .collect(),
        )),
        Geometry::LineString(l) => Geometry::LineString(ls(l)),
        Geometry::MultiLineString(m) => Geometry::MultiLineString(MultiLineString::new(
            m.line_strings().iter().map(ls).collect(),
        )),
        Geometry::Polygon(p) => Geometry::Polygon(pg(p)),
        Geometry::MultiPolygon(m) => Geometry::MultiPolygon(
            MultiPolygon::new(m.polygons().iter().map(pg).collect()).unwrap(),
        ),
        Geometry::GeometryCollection(gc) => Geometry::GeometryCollection(
            GeometryCollection::new(gc.geometries().iter().map(|e| map(e, f)).collect()).unwrap(),
        ),
        other => panic!("no map for {other:?}"),
    }
}

/// Every coordinate of `g`, in storage order.
pub fn coords(g: &Geometry) -> Vec<Coord2d> {
    let mut out = Vec::new();
    collect(g, &mut out);
    out
}

fn collect(g: &Geometry, out: &mut Vec<Coord2d>) {
    let rings = |p: &Polygon, out: &mut Vec<Coord2d>| {
        out.extend(p.exterior().coords().iter());
        for h in p.holes() {
            out.extend(h.coords().iter());
        }
    };
    match g {
        Geometry::Point(p) => out.push(p.coord()),
        Geometry::MultiPoint(m) => out.extend(m.points().iter().map(Point::coord)),
        Geometry::LineString(l) => out.extend(l.coords().iter()),
        Geometry::MultiLineString(m) => {
            for l in m.line_strings() {
                out.extend(l.coords().iter());
            }
        }
        Geometry::Polygon(p) => rings(p, out),
        Geometry::MultiPolygon(m) => {
            for p in m.polygons() {
                rings(p, out);
            }
        }
        Geometry::GeometryCollection(gc) => {
            for e in gc.geometries() {
                collect(e, out);
            }
        }
        other => panic!("no coordinates for {other:?}"),
    }
}
