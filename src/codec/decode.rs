// src/codec/decode.rs

//! Decoding half of the native format.
//!
//! # The rule this module exists to enforce
//!
//! **Every geometry is rebuilt through its public constructor.** Nothing here
//! assembles a type from fields. `Point::new`, `LineString::new`,
//! `LinearRing::new` and `Polygon::new` each validate, and their errors are
//! propagated as [`CodecError::Geometry`]. A byte stream describing a NaN
//! coordinate, a two-vertex ring or a clockwise exterior therefore yields an
//! error and not a geometry. `Point`'s field is private for the same reason:
//! no path, this decoder included, can bypass the constructor.
//!
//! # And the rule that keeps it from panicking
//!
//! Counts come from the stream and are therefore hostile. Nothing here
//! pre-allocates from a count, indexes a slice, or slices by a length it read:
//! every read goes through [`Reader`], which is bounds-checked and returns
//! `Result`. `forbid(unsafe_code)` and `deny(unwrap_used)` rule out memory
//! unsafety and unwraps, but neither prevents an out-of-bounds index or an
//! out-of-memory abort from `Vec::with_capacity(huge)`.

use super::{CodecError, Encode, FORMAT_VERSION, Reader, tag};
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

/// A type that can be rebuilt from the native format.
pub trait Decode: Sized + Encode {
    /// Read this value's body; the header has already been consumed.
    ///
    /// # Errors
    ///
    /// Any [`CodecError`]. Implementations must not panic on any input.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError>;
}

/// Decode a complete document into `T`.
///
/// Checks the version, checks that the tag is `T`'s, decodes the body, and
/// requires the stream to be fully consumed.
///
/// # Errors
///
/// - [`CodecError::UnexpectedEnd`]: the stream ends mid-value.
/// - [`CodecError::UnsupportedVersion`]: the version byte is not
///   [`FORMAT_VERSION`].
/// - [`CodecError::UnexpectedTag`]: the document is not a `T`.
/// - [`CodecError::Geometry`]: the bytes describe an invalid geometry.
/// - [`CodecError::TrailingBytes`]: bytes remained after a complete value.
pub fn decode<T: Decode>(bytes: &[u8]) -> Result<T, CodecError> {
    let mut r = Reader::new(bytes);

    let version = r.u8()?;
    if version != FORMAT_VERSION {
        return Err(CodecError::UnsupportedVersion { found: version });
    }

    let found = r.u8()?;
    if found != T::TAG {
        return Err(CodecError::UnexpectedTag {
            expected: T::TAG,
            found,
        });
    }

    let value = T::decode_body(&mut r)?;

    // A document holds exactly one value. Leftover bytes mean the stream is
    // not what it claims, so it is refused rather than partly believed.
    let remaining = r.remaining();
    if remaining > 0 {
        return Err(CodecError::TrailingBytes { remaining });
    }
    Ok(value)
}

/// Read `n` elements, where `n` came from the stream.
///
/// The capacity hint is clamped: `n` is attacker-controlled, and
/// `Vec::with_capacity(n)` on a four-billion count aborts the process rather
/// than returning an error, which would turn a malformed stream into a crash.
/// Reserving a bounded amount and letting the vector grow costs a few
/// reallocations on genuinely large input and cannot be provoked into an
/// out-of-memory abort by a short hostile document.
fn read_n<T, F>(r: &mut Reader<'_>, n: usize, mut each: F) -> Result<Vec<T>, CodecError>
where
    F: FnMut(&mut Reader<'_>) -> Result<T, CodecError>,
{
    const CAPACITY_CLAMP: usize = 1024;
    let mut out = Vec::with_capacity(n.min(CAPACITY_CLAMP));
    for _ in 0..n {
        out.push(each(r)?);
    }
    Ok(out)
}

// ======================================================================
// 2D coordinates
// ======================================================================

impl Decode for Coord2d {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let x = r.f64()?;
        let y = r.f64()?;
        Ok(Coord2d::new(x, y))
    }
}

impl Decode for CoordSeq {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        Ok(CoordSeq::new(read_n(r, n, Coord2d::decode_body)?))
    }
}

// ======================================================================
// 2D geometry, each through its constructor
// ======================================================================

impl Decode for Point {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        Ok(Point::new(Coord2d::decode_body(r)?)?)
    }
}

impl Decode for LineString {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        Ok(LineString::new(CoordSeq::decode_body(r)?)?)
    }
}

impl Decode for LinearRing {
    /// The sequence read is the ring's distinct vertices; no closing vertex is
    /// expected, and one is refused. `LinearRing::from_distinct` enforces the
    /// minimum length, finiteness and the no-consecutive-duplicates rule,
    /// including the wrap-around pair, without `LinearRing::new`'s stripping
    /// of a closing vertex: a stream that repeats the first vertex at the end
    /// is not what the encoder writes, and accepting it would give the ring
    /// a second encoding that does not re-encode to the bytes read.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        Ok(LinearRing::from_distinct(CoordSeq::decode_body(r)?)?)
    }
}

impl Decode for Polygon {
    /// `Polygon::new` checks exterior and hole winding, so a stream describing
    /// a clockwise exterior is refused here rather than producing a polygon
    /// whose orientation every downstream algorithm would misread.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let exterior = LinearRing::decode_body(r)?;
        let n = r.count()?;
        let holes = read_n(r, n, LinearRing::decode_body)?;
        Ok(Polygon::new(exterior, holes)?)
    }
}

impl Decode for MultiPoint {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        Ok(MultiPoint::new(read_n(r, n, Point::decode_body)?))
    }
}

impl Decode for MultiLineString {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        Ok(MultiLineString::new(read_n(r, n, LineString::decode_body)?))
    }
}

impl Decode for MultiPolygon {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        MultiPolygon::new(read_n(r, n, Polygon::decode_body)?).map_err(CodecError::Geometry)
    }
}

impl Decode for Geometry {
    /// An unknown variant byte is an error.
    ///
    /// `Geometry` is `#[non_exhaustive]`, so a future variant will reach a
    /// decoder that predates it. Refusing is the only safe reading: skipping
    /// the field would require knowing its length, and guessing would produce
    /// a geometry the writer did not describe.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        match r.u8()? {
            tag::GEOMETRY_COLLECTION => Ok(Geometry::GeometryCollection(
                GeometryCollection::decode_body(r)?,
            )),
            variant => decode_leaf(variant, r),
        }
    }
}

/// The body of a `Geometry` of variant `variant`, which is not a collection.
fn decode_leaf(variant: u8, r: &mut Reader<'_>) -> Result<Geometry, CodecError> {
    Ok(match variant {
        tag::POINT => Geometry::Point(Point::decode_body(r)?),
        tag::LINE_STRING => Geometry::LineString(LineString::decode_body(r)?),
        tag::POLYGON => Geometry::Polygon(Polygon::decode_body(r)?),
        tag::MULTI_POINT => Geometry::MultiPoint(MultiPoint::decode_body(r)?),
        tag::MULTI_LINE_STRING => Geometry::MultiLineString(MultiLineString::decode_body(r)?),
        tag::MULTI_POLYGON => Geometry::MultiPolygon(MultiPolygon::decode_body(r)?),
        found => return Err(CodecError::UnknownVariant { found }),
    })
}

/// A collection being decoded: the elements still to read, and those read.
struct Frame {
    remaining: usize,
    elements: Vec<Geometry>,
}

impl Frame {
    fn open(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        const CAPACITY_CLAMP: usize = 1024;
        let remaining = r.count()?;
        Ok(Self {
            remaining,
            elements: Vec::with_capacity(remaining.min(CAPACITY_CLAMP)),
        })
    }
}

impl Decode for GeometryCollection {
    /// Decoded with an explicit stack of open collections, not by recursion,
    /// and **the depth is checked before a level is entered**: a nested
    /// collection's header at depth [`GeometryCollection::MAX_NESTING`] is
    /// refused with [`CodecError::NestingTooDeep`] before anything inside it
    /// is read. A hostile stream therefore cannot nest the decoder deeper
    /// than the limit, and every collection that decodes was built by
    /// [`GeometryCollection::new`], which checks the same limit.
    ///
    /// Counts are hostile: capacity is clamped, and a count larger than the
    /// stream ends in [`CodecError::UnexpectedEnd`].
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        // `stack.len()` is the depth of the innermost open collection.
        let mut stack = vec![Frame::open(r)?];
        loop {
            // The stack is never empty here: the outermost frame is the last
            // popped, and that returns.
            let Some(top) = stack.last_mut() else {
                return Err(CodecError::UnexpectedEnd {
                    needed: 1,
                    available: r.remaining(),
                });
            };
            if top.remaining == 0 {
                let done = stack.pop().map_or_else(Vec::new, |f| f.elements);
                let c = GeometryCollection::new(done)?;
                match stack.last_mut() {
                    None => return Ok(c),
                    Some(parent) => parent.elements.push(Geometry::GeometryCollection(c)),
                }
                continue;
            }
            top.remaining -= 1;
            let variant = r.u8()?;
            if variant == tag::GEOMETRY_COLLECTION {
                if stack.len() >= GeometryCollection::MAX_NESTING {
                    return Err(CodecError::NestingTooDeep {
                        limit: GeometryCollection::MAX_NESTING,
                    });
                }
                stack.push(Frame::open(r)?);
            } else {
                let g = decode_leaf(variant, r)?;
                top.elements.push(g);
            }
        }
    }
}

impl Decode for Envelope {
    /// Rebuilt through `Envelope::new`, so a stream describing inverted or
    /// non-finite bounds is refused rather than producing a rectangle whose
    /// every containment test would be silently false.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let min_x = r.f64()?;
        let min_y = r.f64()?;
        let max_x = r.f64()?;
        let max_y = r.f64()?;
        Ok(Envelope::new(
            Coord2d::new(min_x, min_y),
            Coord2d::new(max_x, max_y),
        )?)
    }
}

// ======================================================================
// Topology results
// ======================================================================

impl Decode for Location {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        match r.u8()? {
            0 => Ok(Location::Interior),
            1 => Ok(Location::Boundary),
            2 => Ok(Location::Exterior),
            found => Err(CodecError::InvalidEnumValue {
                type_name: "Location",
                found: i16::from(found),
            }),
        }
    }
}

impl Decode for Dimension {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        match r.i8()? {
            -1 => Ok(Dimension::Empty),
            0 => Ok(Dimension::Point),
            1 => Ok(Dimension::Line),
            2 => Ok(Dimension::Area),
            found => Err(CodecError::InvalidEnumValue {
                type_name: "Dimension",
                found: i16::from(found),
            }),
        }
    }
}

impl Decode for IntersectionMatrix {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        Ok(IntersectionMatrix {
            interior_interior: Dimension::decode_body(r)?,
            interior_boundary: Dimension::decode_body(r)?,
            interior_exterior: Dimension::decode_body(r)?,
            boundary_interior: Dimension::decode_body(r)?,
            boundary_boundary: Dimension::decode_body(r)?,
            boundary_exterior: Dimension::decode_body(r)?,
            exterior_interior: Dimension::decode_body(r)?,
            exterior_boundary: Dimension::decode_body(r)?,
            exterior_exterior: Dimension::decode_body(r)?,
        })
    }
}

// ======================================================================
// 3D geometry
// ======================================================================

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for Coord3d {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let x = r.f64()?;
        let y = r.f64()?;
        let z = r.f64()?;
        Ok(Coord3d::new(x, y, z))
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for CoordSeq3d {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        Ok(CoordSeq3d::new(read_n(r, n, Coord3d::decode_body)?))
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for Point3d {
    /// Rebuilt through `Point3d::new`, which validates finiteness on all three
    /// axes. A stream describing a NaN or infinite coordinate is therefore
    /// refused rather than decoded, so the validate-on-decode guarantee holds
    /// in 3D as it does in 2D.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        Ok(Point3d::new(Coord3d::decode_body(r)?)?)
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for LineString3d {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        Ok(LineString3d::new(CoordSeq3d::decode_body(r)?)?)
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for LinearRing3d {
    /// Distinct vertices, a closing vertex refused, as for `LinearRing`.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        Ok(LinearRing3d::from_distinct(CoordSeq3d::decode_body(r)?)?)
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for Polygon3d {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let exterior = LinearRing3d::decode_body(r)?;
        let n = r.count()?;
        let holes = read_n(r, n, LinearRing3d::decode_body)?;
        Ok(Polygon3d::new(exterior, holes)?)
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for MultiPoint3d {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        Ok(MultiPoint3d::new(read_n(r, n, Point3d::decode_body)?))
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for MultiLineString3d {
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        Ok(MultiLineString3d::new(read_n(
            r,
            n,
            LineString3d::decode_body,
        )?))
    }
}

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
impl Decode for MultiPolygon3d {
    /// Rebuilt through `MultiPolygon3d::new`, so parts that overlap in plan
    /// on the wire are refused, as they are at construction.
    fn decode_body(r: &mut Reader<'_>) -> Result<Self, CodecError> {
        let n = r.count()?;
        Ok(MultiPolygon3d::new(read_n(r, n, Polygon3d::decode_body)?)?)
    }
}
