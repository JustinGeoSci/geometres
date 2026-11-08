// fuzz/fuzz_targets/collection_relate.rs

//! Fuzz target: `relate` on `GeometryCollection` operands.
//!
//! The input decodes into two small geometries on an integer grid, each
//! wrapped in collection forms. The contract checked:
//!
//! - `relate` never panics, collection operands included;
//! - a refusal happens **exactly** when the collection's polygonal leaves
//!   fail the multipolygon part rule (checked independently through
//!   `MultiPolygon::new` on the cloned leaves), and only
//!   `CollectionElementsOverlap` is ever the error;
//! - the reductions: a singleton collection, a doubly nested
//!   singleton, and the per-leaf split of a multi-geometry relate exactly
//!   as the unwrapped geometry, on either side;
//! - the transpose: `relate(b, a)` is `relate(a, b)` transposed for every
//!   accepted pair involving a collection.

#![no_main]

use geometres::GeometryError;
use geometres::topology::de9im::IntersectionMatrix;
use geometres::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiPolygon, Point,
    Polygon,
};
use libfuzzer_sys::fuzz_target;

/// A cursor over the fuzz bytes; runs out gracefully (zeros).
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

    fn grid(&mut self) -> f64 {
        f64::from(self.u8() % 12)
    }
}

fn sq(x: f64, y: f64, s: f64) -> Option<Polygon> {
    let ring = LinearRing::new(CoordSeq::new(vec![
        Coord2d::new(x, y),
        Coord2d::new(x + s, y),
        Coord2d::new(x + s, y + s),
        Coord2d::new(x, y + s),
    ]))
    .ok()?;
    Polygon::new_simple(ring).ok()
}

fn gen_line(b: &mut Bytes<'_>) -> Option<LineString> {
    let n = 2 + usize::from(b.u8() % 3);
    let mut pts: Vec<Coord2d> = Vec::new();
    while pts.len() < n {
        let p = Coord2d::new(b.grid(), b.grid());
        if pts.last() != Some(&p) {
            pts.push(p);
        } else if pts.len() >= 2 {
            break;
        } else {
            return None;
        }
    }
    LineString::new(CoordSeq::new(pts)).ok()
}

/// A small geometry: a grid square (possibly overlapping others), a line, a
/// point, or a two-part multipolygon on separate cells.
fn gen_leaf(b: &mut Bytes<'_>) -> Option<Geometry> {
    Some(match b.u8() % 5 {
        0 | 1 => Geometry::Polygon(sq(b.grid(), b.grid(), 1.0 + f64::from(b.u8() % 3))?),
        2 => Geometry::LineString(gen_line(b)?),
        3 => Geometry::Point(Point::new(Coord2d::new(b.grid(), b.grid())).ok()?),
        _ => {
            let (i, j) = (f64::from(b.u8() % 3), f64::from(b.u8() % 3));
            Geometry::MultiPolygon(
                MultiPolygon::new(vec![sq(i * 4.0, 0.0, 3.0)?, sq(j * 4.0, 5.0, 3.0)?]).ok()?,
            )
        }
    })
}

/// A collection of 1–3 leaves, nested one level deeper on demand.
fn gen_collection(b: &mut Bytes<'_>) -> Option<GeometryCollection> {
    let n = 1 + usize::from(b.u8() % 3);
    let mut gs: Vec<Geometry> = Vec::new();
    for _ in 0..n {
        gs.push(gen_leaf(b)?);
    }
    let inner = GeometryCollection::new(gs).ok()?;
    if b.u8().is_multiple_of(4) {
        GeometryCollection::new(vec![Geometry::GeometryCollection(inner)]).ok()
    } else {
        Some(inner)
    }
}

/// The polygonal leaves, for the independent restriction check.
fn polygon_leaves(gc: &GeometryCollection) -> Vec<Polygon> {
    let mut out = Vec::new();
    for leaf in gc.leaves() {
        match leaf {
            Geometry::Polygon(p) => out.push(p.clone()),
            Geometry::MultiPolygon(m) => out.extend(m.polygons().iter().cloned()),
            _ => {}
        }
    }
    out
}

fn relate_gc(a: &Geometry, b: &Geometry) -> Result<IntersectionMatrix, GeometryError> {
    IntersectionMatrix::compute_geometry(a, b)
}

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let Some(gc) = gen_collection(&mut b) else {
        return;
    };
    let Some(other) = gen_leaf(&mut b) else {
        return;
    };
    let g = Geometry::GeometryCollection(gc.clone());

    // The restriction, checked independently: the collection is accepted
    // exactly when its polygonal leaves form a valid multipolygon.
    let leaves_ok = MultiPolygon::new(polygon_leaves(&gc)).is_ok();
    match relate_gc(&g, &other) {
        Ok(m) => {
            assert!(leaves_ok, "accepted a collection whose leaves overlap");
            // Transpose.
            let t = relate_gc(&other, &g).expect("the transpose must be accepted too");
            assert_eq!(t, m.swap_ab(), "transpose");
        }
        Err(GeometryError::CollectionElementsOverlap { .. }) => {
            assert!(!leaves_ok, "refused a collection whose leaves are valid");
        }
        Err(e) => panic!("unexpected error {e:?}"),
    }

    // Reductions on the other side's leaf: a singleton and a nested
    // singleton relate exactly as the leaf, against the plain leaf.
    let w1 = Geometry::GeometryCollection(
        GeometryCollection::new(vec![other.clone()])
            .unwrap_or_else(|_| GeometryCollection::empty()),
    );
    let w2 = Geometry::GeometryCollection(
        GeometryCollection::new(vec![w1.clone()]).unwrap_or_else(|_| GeometryCollection::empty()),
    );
    if let Some(third) = gen_leaf(&mut b) {
        let reference = relate_gc(&other, &third).expect("simple pairs never fail");
        for wrapped in [&w1, &w2] {
            assert_eq!(
                relate_gc(wrapped, &third).expect("singletons never refused"),
                reference,
                "reduction"
            );
        }
        // The per-leaf split of a multipoint: repeated points included.
        if let Geometry::MultiPoint(mp) = &other {
            let split = Geometry::GeometryCollection(
                GeometryCollection::new(mp.points().iter().map(|p| Geometry::Point(*p)).collect())
                    .unwrap_or_else(|_| GeometryCollection::empty()),
            );
            assert_eq!(
                relate_gc(&split, &third).expect("point collections never refused"),
                reference,
                "multipoint split"
            );
        }
    }
});
