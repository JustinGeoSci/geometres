// src/codec/encode.rs

//! Encoding half of the native format.
//!
//! Every implementation here appends to a `Vec<u8>` and never allocates per
//! element beyond that vector's own growth: no `std::io::Write`, no formatter,
//! no intermediate string.

use super::{FORMAT_VERSION, tag};
use crate::topology::de9im::IntersectionMatrix;
use crate::topology::label::{Dimension, Location};
use crate::types::{
    Coord2d, CoordSeq, Envelope, Geometry, GeometryCollection, LineString, LinearRing,
    MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};

#[cfg(feature = "3d")]
use crate::types3d::{
    Coord3d, CoordSeq3d, LineString3d, LinearRing3d, MultiLineString3d, MultiPoint3d,
    MultiPolygon3d, Point3d, Polygon3d,
};

/// A type with a native binary encoding.
///
/// `TAG` identifies the type in a document header and, for the variants of
/// [`Geometry`], inside a `Geometry` body. `encode_body` writes the value
/// without any header, so a type can be nested inside another.
pub trait Encode {
    /// This type's tag byte.
    const TAG: u8;

    /// Append this value's body, without version or tag, to `out`.
    fn encode_body(&self, out: &mut Vec<u8>);
}

/// Encode `value` as a complete document: version, tag, body.
///
/// Deterministic: the same value yields the same bytes on every run and
/// platform.
#[must_use]
pub fn encode<T: Encode>(value: &T) -> Vec<u8> {
    let mut out = Vec::new();
    encode_into(value, &mut out);
    out
}

/// Append a complete document for `value` to an existing buffer.
///
/// Useful when several documents go into one allocation; each remains
/// independently decodable because each carries its own version and tag.
pub fn encode_into<T: Encode>(value: &T, out: &mut Vec<u8>) {
    out.push(FORMAT_VERSION);
    out.push(T::TAG);
    value.encode_body(out);
}

// ======================================================================
// Primitives
// ======================================================================

/// Write a `u32` count, little-endian.
///
/// `usize` is narrowed deliberately: the format fixes counts at four bytes so
/// that a document written on a 64-bit machine is readable on a 32-bit one. A
/// collection longer than `u32::MAX` is not representable and saturates, which
/// is unreachable for any geometry that fits in memory.
///
/// Public so that a crate built on this one can encode its own types in the
/// same format, with the same count encoding
/// [`Reader::count`](super::Reader::count) reads back. Such a type takes a tag
/// in
/// [`DOWNSTREAM_FIRST`](super::tag::DOWNSTREAM_FIRST)`..=`[`DOWNSTREAM_LAST`](super::tag::DOWNSTREAM_LAST).
pub fn put_count(n: usize, out: &mut Vec<u8>) {
    let n = u32::try_from(n).unwrap_or(u32::MAX);
    out.extend_from_slice(&n.to_le_bytes());
}

/// Write an `f64` by its bit pattern, little-endian.
///
/// `to_le_bytes` rather than any decimal rendering: a text form would not be
/// deterministic across platforms and would not round-trip every bit pattern.
/// Every bit pattern is written as is, NaN payloads and signed zero included;
/// whether the value is acceptable is the decoder's constructor's decision.
///
/// Public for the same reason as [`put_count`]; [`Reader::f64`](super::Reader::f64)
/// is its inverse.
pub fn put_f64(v: f64, out: &mut Vec<u8>) {
    out.extend_from_slice(&v.to_le_bytes());
}

// ======================================================================
// 2D coordinates
// ======================================================================

impl Encode for Coord2d {
    const TAG: u8 = tag::COORD2D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_f64(self.x, out);
        put_f64(self.y, out);
    }
}

impl Encode for CoordSeq {
    const TAG: u8 = tag::COORD_SEQ;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.len(), out);
        for c in self.as_slice() {
            c.encode_body(out);
        }
    }
}

// ======================================================================
// 2D geometry
// ======================================================================

impl Encode for Point {
    const TAG: u8 = tag::POINT;

    fn encode_body(&self, out: &mut Vec<u8>) {
        self.coord().encode_body(out);
    }
}

impl Encode for LineString {
    const TAG: u8 = tag::LINE_STRING;

    fn encode_body(&self, out: &mut Vec<u8>) {
        self.coords().encode_body(out);
    }
}

impl Encode for LinearRing {
    const TAG: u8 = tag::LINEAR_RING;

    /// The ring's `n` **distinct** vertices, with closure implied.
    ///
    /// No closing vertex is written. `LinearRing::new` strips one if given, so
    /// writing it here would make a round trip asymmetric and would mix two
    /// conventions for the same ring.
    fn encode_body(&self, out: &mut Vec<u8>) {
        self.coords().encode_body(out);
    }
}

impl Encode for Polygon {
    const TAG: u8 = tag::POLYGON;

    fn encode_body(&self, out: &mut Vec<u8>) {
        self.exterior().encode_body(out);
        put_count(self.holes().len(), out);
        for h in self.holes() {
            h.encode_body(out);
        }
    }
}

impl Encode for MultiPoint {
    const TAG: u8 = tag::MULTI_POINT;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.points().len(), out);
        for p in self.points() {
            p.encode_body(out);
        }
    }
}

impl Encode for MultiLineString {
    const TAG: u8 = tag::MULTI_LINE_STRING;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.line_strings().len(), out);
        for l in self.line_strings() {
            l.encode_body(out);
        }
    }
}

impl Encode for MultiPolygon {
    const TAG: u8 = tag::MULTI_POLYGON;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.polygons().len(), out);
        for p in self.polygons() {
            p.encode_body(out);
        }
    }
}

impl Encode for Geometry {
    const TAG: u8 = tag::GEOMETRY;

    /// A variant byte, then that variant's body.
    ///
    /// The variant byte is the inner type's own `TAG`, so the tag space is
    /// shared between document headers and `Geometry` discriminants.
    ///
    /// The match is deliberately exhaustive with **no wildcard arm**.
    /// `#[non_exhaustive]` does not apply inside the defining crate, so the
    /// compiler requires every variant to be listed here. Adding a variant to
    /// `Geometry` without giving it an encoding is therefore a build failure
    /// rather than a value that silently encodes as a placeholder. That is
    /// the stronger guarantee, and the reason a catch-all was removed.
    fn encode_body(&self, out: &mut Vec<u8>) {
        match self {
            Geometry::Point(p) => {
                out.push(Point::TAG);
                p.encode_body(out);
            }
            Geometry::LineString(l) => {
                out.push(LineString::TAG);
                l.encode_body(out);
            }
            Geometry::Polygon(p) => {
                out.push(Polygon::TAG);
                p.encode_body(out);
            }
            Geometry::MultiPoint(m) => {
                out.push(MultiPoint::TAG);
                m.encode_body(out);
            }
            Geometry::MultiLineString(m) => {
                out.push(MultiLineString::TAG);
                m.encode_body(out);
            }
            Geometry::MultiPolygon(m) => {
                out.push(MultiPolygon::TAG);
                m.encode_body(out);
            }
            Geometry::GeometryCollection(c) => {
                out.push(GeometryCollection::TAG);
                c.encode_body(out);
            }
        }
    }
}

impl Encode for GeometryCollection {
    const TAG: u8 = tag::GEOMETRY_COLLECTION;

    /// `count:u32`, then each element as a `Geometry` body: its variant byte,
    /// then that variant's body. A nested collection is variant
    /// `GEOMETRY_COLLECTION` with this same body.
    ///
    /// Written with an explicit stack, not by recursion: the walk is the same
    /// depth-first order as [`GeometryCollection::leaves`], and a nested
    /// collection's header is written when it is entered. Every element that
    /// is not a collection is encoded by `Geometry`'s own (non-recursive)
    /// arm.
    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.len(), out);
        let mut stack = vec![self.geometries().iter()];
        while let Some(top) = stack.last_mut() {
            match top.next() {
                None => {
                    stack.pop();
                }
                Some(Geometry::GeometryCollection(c)) => {
                    out.push(Self::TAG);
                    put_count(c.len(), out);
                    stack.push(c.geometries().iter());
                }
                Some(g) => g.encode_body(out),
            }
        }
    }
}

impl Encode for Envelope {
    const TAG: u8 = tag::ENVELOPE;

    /// Four bounds in `min_x, min_y, max_x, max_y` order.
    ///
    /// The corners are written rather than a corner plus extents: a width is a
    /// subtraction of two `f64` values and need not be exactly representable,
    /// so reconstructing `max` from `min + width` could land on a different
    /// coordinate than the one encoded.
    fn encode_body(&self, out: &mut Vec<u8>) {
        let (min, max) = (self.min(), self.max());
        put_f64(min.x, out);
        put_f64(min.y, out);
        put_f64(max.x, out);
        put_f64(max.y, out);
    }
}

// ======================================================================
// Topology results
// ======================================================================

impl Encode for Location {
    const TAG: u8 = tag::LOCATION;

    fn encode_body(&self, out: &mut Vec<u8>) {
        out.push(match self {
            Location::Interior => 0,
            Location::Boundary => 1,
            Location::Exterior => 2,
        });
    }
}

impl Encode for Dimension {
    const TAG: u8 = tag::DIMENSION;

    /// The OGC numeric value as a signed byte, so `Empty` is −1 rather than a
    /// renumbered ordinal. Keeping the specification's values means the wire
    /// form is readable against the DE-9IM tables directly.
    fn encode_body(&self, out: &mut Vec<u8>) {
        let v: i8 = match self {
            Dimension::Empty => -1,
            Dimension::Point => 0,
            Dimension::Line => 1,
            Dimension::Area => 2,
        };
        // `cast_unsigned` rather than `as`: the reinterpretation is intended
        // (−1 is written as 0xFF and read back as −1), and naming it keeps
        // `cast_sign_loss` enabled everywhere else. That lint once caught a
        // sign inversion in this encoding.
        out.push(v.cast_unsigned());
    }
}

impl Encode for IntersectionMatrix {
    const TAG: u8 = tag::INTERSECTION_MATRIX;

    /// Nine cells in row-major DE-9IM order, which is the order they are
    /// written in the OGC pattern string.
    fn encode_body(&self, out: &mut Vec<u8>) {
        for d in [
            self.interior_interior,
            self.interior_boundary,
            self.interior_exterior,
            self.boundary_interior,
            self.boundary_boundary,
            self.boundary_exterior,
            self.exterior_interior,
            self.exterior_boundary,
            self.exterior_exterior,
        ] {
            d.encode_body(out);
        }
    }
}

// ======================================================================
// 3D geometry
// ======================================================================

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for Coord3d {
    const TAG: u8 = tag::COORD3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_f64(self.x, out);
        put_f64(self.y, out);
        put_f64(self.z, out);
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for CoordSeq3d {
    const TAG: u8 = tag::COORD_SEQ3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.len(), out);
        for c in self.as_slice() {
            c.encode_body(out);
        }
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for Point3d {
    const TAG: u8 = tag::POINT3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        self.coord().encode_body(out);
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for LineString3d {
    const TAG: u8 = tag::LINE_STRING3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        self.coords().encode_body(out);
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for LinearRing3d {
    const TAG: u8 = tag::LINEAR_RING3D;

    /// Distinct vertices, closure implied, as for the 2D ring.
    fn encode_body(&self, out: &mut Vec<u8>) {
        self.coords().encode_body(out);
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for Polygon3d {
    const TAG: u8 = tag::POLYGON3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        self.exterior().encode_body(out);
        put_count(self.holes().len(), out);
        for h in self.holes() {
            h.encode_body(out);
        }
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for MultiPoint3d {
    const TAG: u8 = tag::MULTI_POINT3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.points().len(), out);
        for p in self.points() {
            p.encode_body(out);
        }
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for MultiLineString3d {
    const TAG: u8 = tag::MULTI_LINE_STRING3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.line_strings().len(), out);
        for l in self.line_strings() {
            l.encode_body(out);
        }
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Encode for MultiPolygon3d {
    const TAG: u8 = tag::MULTI_POLYGON3D;

    fn encode_body(&self, out: &mut Vec<u8>) {
        put_count(self.polygons().len(), out);
        for p in self.polygons() {
            p.encode_body(out);
        }
    }
}
