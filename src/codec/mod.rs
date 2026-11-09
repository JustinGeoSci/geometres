// src/codec/mod.rs

//! Native binary serialization for `geometres` types.
//!
//! This module is the authoritative definition of a serialized `geometres`
//! value. It replaces the crate's former optional derive-based
//! serialization, and the format described here is part of the crate's
//! public contract.
//!
//! # Why a native format
//!
//! Not dependency count: the previous mechanism was already optional, so the
//! shipped graph was clean. The properties wanted are **deterministic encoding,
//! low allocation overhead, a byte-slice interface with no I/O, and a path to
//! compact or zero-copy representations**, and no single off-the-shelf format
//! crate provides all four.
//!
//! # The four guarantees
//!
//! **Deterministic.** The same value produces byte-identical output on every
//! run, platform and build. No map iteration, no hash seeding, no float
//! formatting: every `f64` is written with [`f64::to_le_bytes`] and read with
//! [`f64::from_le_bytes`], never through a decimal representation.
//!
//! **Round-trip exact.** Decoding an encoding yields a `PartialEq`-equal value,
//! `f64` bit patterns included. A coordinate that comes back with different
//! bits is a defect, not a tolerance.
//!
//! **Validating on decode.** Every geometry is rebuilt through its public
//! constructor (`Point::new`, `LineString::new`, `LinearRing::new`,
//! `Polygon::new`, `MultiPolygon::new`, `GeometryCollection::new` and the 3D
//! counterparts), and the constructor's error is propagated. No type is ever
//! assembled field-by-field from bytes. `Point`'s field is private so that
//! deserialization cannot bypass `Point::new`, and this decoder honors that.
//! **A malformed or hostile byte stream produces a [`CodecError`], never an
//! invalid geometry and never a panic.** `Coord2d` and `CoordSeq` are
//! coordinate values rather than geometries: like their constructors, their
//! decoders accept any `f64`, NaN and infinity included.
//!
//! **No I/O.** Encoding appends to a `Vec<u8>` and decoding reads a `&[u8]`;
//! nothing here uses `std::io`. The codec allocates, and like the rest of the
//! crate it is built on `std`: it is not `no_std`.
//!
//! # Format
//!
//! All integers are little-endian. All lengths are explicit `u32` counts, so a
//! decoder never scans for a boundary.
//!
//! ```text
//! document := version:u8 tag:u8 body
//! ```
//!
//! `version` is [`FORMAT_VERSION`]. A document whose version is not recognized
//! is rejected rather than guessed. The field exists from byte zero so the
//! format can change without orphaning stored artifacts.
//!
//! `tag` identifies the top-level type ([`tag`] module). An unrecognized tag is
//! an error, and so is a tag that does not match the type being decoded.
//!
//! ## Primitives
//!
//! | element | encoding |
//! |---|---|
//! | `f64` | 8 bytes, `to_le_bytes` |
//! | count | `u32`, 4 bytes little-endian |
//! | `Coord2d` | `x:f64 y:f64` |
//! | `Coord3d` | `x:f64 y:f64 z:f64` |
//! | coordinate sequence | `count:u32` then `count` coordinates |
//!
//! ## Bodies
//!
//! | type | body |
//! |---|---|
//! | `Point`, `Point3d` | one coordinate |
//! | `LineString`, `LineString3d` | coordinate sequence |
//! | `LinearRing`, `LinearRing3d` | coordinate sequence; see the closure note |
//! | `Polygon`, `Polygon3d` | exterior sequence, `holes:u32`, then each hole's sequence |
//! | `Envelope` | `min_x:f64 min_y:f64 max_x:f64 max_y:f64` |
//! | `MultiPoint3d` | `count:u32` then each point's body |
//! | `MultiLineString3d` | `count:u32` then each line string's body |
//! | `MultiPolygon3d` | `count:u32` then each polygon's body (tag `0x28`) |
//! | `MultiPoint` | `count:u32` then each point's body |
//! | `MultiLineString` | `count:u32` then each line string's body |
//! | `MultiPolygon` | `count:u32` then each polygon's body |
//! | `Geometry` | `variant:u8` then that variant's body |
//! | `GeometryCollection` | `count:u32` then each element as a `Geometry` body (tag `0x11`) |
//! | `Location` | `u8`: 0 interior, 1 boundary, 2 exterior |
//! | `Dimension` | `i8`: −1 empty, 0 point, 1 line, 2 area (the OGC values) |
//! | `IntersectionMatrix` | nine `Dimension` values in row-major DE-9IM order |
//!
//! The `Geometry` variant byte reuses the top-level [`tag`] values, so one tag
//! space serves both.
//!
//! ## Collections and their depth
//!
//! A `GeometryCollection`'s elements are `Geometry` bodies, so a collection
//! nests by writing variant `0x11` with a collection body. The nesting is at
//! most `GeometryCollection::MAX_NESTING` (32) levels, a permanent part of
//! the format: the decoder refuses a deeper header with
//! [`CodecError::NestingTooDeep`] before reading inside it, and the type
//! refuses a deeper value at construction, so every collection that exists
//! round-trips and no stream can nest the decoder past the limit.
//! `Geometry` is `#[non_exhaustive]`; an unknown variant byte is
//! [`CodecError::UnknownVariant`], never a silently skipped field.
//!
//! ## Ring closure, stated explicitly
//!
//! **A `LinearRing` is encoded as its `n` distinct vertices with closure
//! implied. The encoder does not write a closing vertex and the decoder
//! refuses one** (as the duplicate pair across the closing edge,
//! `GeometryError::ConsecutiveDuplicateVertices`). Refusing it, rather than
//! stripping it as `LinearRing::new` strips wire-format input, keeps every
//! value's encoding unique: whatever decodes re-encodes to the bytes read.
//!
//! This is spelled out because the opposite convention produces conversion
//! bugs at format boundaries: shapefiles and OGC WKT repeat the first vertex
//! at the end, `LinearRing::new` strips it, and a round trip that re-appends
//! it grows the ring on every pass. An independent implementation of this
//! format must not write the closing vertex.
//!
//! # What is not encoded
//!
//! The result wrappers [`Unvouched<T>`](crate::types::Unvouched) and
//! [`MaybeVouched<T, U>`](crate::types::MaybeVouched) have **no encoding**,
//! and that is a design decision rather than an omission. They describe the
//! outcome of a computation, an invalid result, and a stored artifact has
//! no outcome to describe. Encoding one would also create a route by which an
//! invalid geometry could be written to a file and read back as though it
//! were a value worth keeping, which is the opposite of what the wrappers are
//! for. The inner geometry encodes on its own once a caller has unwrapped it,
//! which is the explicit act the wrapper exists to require.
//!
//! `GeometryError` has no encoding. It never carried a serialization derive
//! either, so there is nothing here it replaces, and an error is a returned value
//! rather than a stored artifact. Giving it a wire format would create a
//! compatibility surface with no consumer.

mod decode;
mod encode;

#[cfg(test)]
mod tests;

pub use decode::{Decode, decode};
pub use encode::{Encode, encode, encode_into, put_count, put_f64};

use core::fmt;

/// Version byte written at the start of every document.
///
/// Present from the first release of the format. A decoder rejects any other
/// value rather than attempting to interpret the remainder.
pub const FORMAT_VERSION: u8 = 1;

/// Type tags.
///
/// These identify the top-level type of a document and double as the variant
/// discriminants inside a [`Geometry`](crate::types::Geometry) body, so a
/// value's tag means the same thing wherever it appears.
///
/// Values are fixed by the format. Changing one is a format change and
/// requires a [`FORMAT_VERSION`] bump.
pub mod tag {
    /// `Coord2d`
    pub const COORD2D: u8 = 0x01;
    /// `CoordSeq`
    pub const COORD_SEQ: u8 = 0x02;
    /// `Point`
    pub const POINT: u8 = 0x03;
    /// `LineString`
    pub const LINE_STRING: u8 = 0x04;
    /// `LinearRing`
    pub const LINEAR_RING: u8 = 0x05;
    /// `Polygon`
    pub const POLYGON: u8 = 0x06;
    /// `MultiPoint`
    pub const MULTI_POINT: u8 = 0x07;
    /// `MultiLineString`
    pub const MULTI_LINE_STRING: u8 = 0x08;
    /// `MultiPolygon`
    pub const MULTI_POLYGON: u8 = 0x09;
    /// `Geometry`
    pub const GEOMETRY: u8 = 0x0A;
    /// `IntersectionMatrix`
    pub const INTERSECTION_MATRIX: u8 = 0x0B;
    /// `Dimension`
    pub const DIMENSION: u8 = 0x0C;
    /// `Location`
    pub const LOCATION: u8 = 0x0D;
    // 0x0E is unassigned and stays so: a document carrying it
    // decodes as an unknown tag rather than as another type.
    /// `Envelope`
    pub const ENVELOPE: u8 = 0x0F;
    /// Reserved: no type in this crate takes this value.
    pub const TESSELLATE_METHOD: u8 = 0x10;
    /// `GeometryCollection`.
    pub const GEOMETRY_COLLECTION: u8 = 0x11;

    /// `Coord3d`
    pub const COORD3D: u8 = 0x20;
    /// `CoordSeq3d`
    pub const COORD_SEQ3D: u8 = 0x21;
    /// `Point3d`
    pub const POINT3D: u8 = 0x22;
    /// `LineString3d`
    pub const LINE_STRING3D: u8 = 0x23;
    /// `LinearRing3d`
    pub const LINEAR_RING3D: u8 = 0x24;
    /// `Polygon3d`
    pub const POLYGON3D: u8 = 0x25;
    /// `MultiPoint3d`
    pub const MULTI_POINT3D: u8 = 0x26;
    /// `MultiLineString3d`
    pub const MULTI_LINE_STRING3D: u8 = 0x27;
    /// `MultiPolygon3d`
    pub const MULTI_POLYGON3D: u8 = 0x28;

    /// First tag value reserved for types defined by crates built on this
    /// one. No type in this crate takes a value in
    /// `DOWNSTREAM_FIRST..=DOWNSTREAM_LAST`, so a downstream encoding can
    /// never be read back as a type of this crate. Allocating values inside
    /// the range is the downstream crate's own business.
    pub const DOWNSTREAM_FIRST: u8 = 0x40;
    /// Last tag value reserved for crates built on this one (inclusive).
    pub const DOWNSTREAM_LAST: u8 = 0x7F;
}

/// Why a byte stream could not be decoded.
///
/// Every variant is a refusal. None of them is recoverable into a partial
/// value: a decode either yields a fully valid geometry or yields one of
/// these.
#[non_exhaustive]
#[derive(Clone, Debug, PartialEq)]
pub enum CodecError {
    /// The stream ended before the value did.
    ///
    /// `needed` is how many more bytes the current read wanted, `available`
    /// how many remained.
    UnexpectedEnd {
        /// Bytes the read required.
        needed: usize,
        /// Bytes actually left in the stream.
        available: usize,
    },
    /// The version byte is not [`FORMAT_VERSION`].
    UnsupportedVersion {
        /// The version byte that was read.
        found: u8,
    },
    /// The tag byte does not identify the type being decoded.
    ///
    /// Produced both by an unrecognized tag and by a recognized tag for a
    /// different type, because decoding a `Polygon` document as a `Point` is
    /// as much an error as decoding noise.
    UnexpectedTag {
        /// The tag the decoder required.
        expected: u8,
        /// The tag the stream carried.
        found: u8,
    },
    /// A `Geometry` variant byte that this version does not define.
    ///
    /// `Geometry` is `#[non_exhaustive]`, so an unknown variant is treated as
    /// a refusal rather than skipped.
    UnknownVariant {
        /// The variant byte that was read.
        found: u8,
    },
    /// An enum payload outside its defined range: a `Location` byte above 2,
    /// or a `Dimension` byte outside −1..=2.
    InvalidEnumValue {
        /// Which enum rejected the value.
        type_name: &'static str,
        /// The byte that was read.
        found: i16,
    },
    /// A count that cannot be represented on this platform.
    ///
    /// Only reachable where `usize` is narrower than `u32`.
    CountOverflow {
        /// The count that was read.
        count: u32,
    },
    /// The bytes decoded, but the value they describe is not a valid geometry.
    ///
    /// This is the constructor's own error, propagated unchanged. It is the
    /// variant that carries the doctrine: a byte stream describing a NaN
    /// coordinate, a two-vertex ring or a clockwise exterior produces this,
    /// never a geometry.
    Geometry(crate::error::GeometryError),
    /// A `GeometryCollection` header was met at the deepest nesting the
    /// format allows; the collection it opens would exceed the limit.
    ///
    /// Checked before the nested collection is read, so a hostile stream
    /// cannot nest the decoder past the limit.
    NestingTooDeep {
        /// The limit, `GeometryCollection::MAX_NESTING`.
        limit: usize,
    },
    /// Bytes remained after a complete value was decoded.
    ///
    /// A document describes exactly one value. Trailing bytes mean the stream
    /// is not what it claims to be, so it is refused rather than partly used.
    TrailingBytes {
        /// How many bytes were left over.
        remaining: usize,
    },
}

impl From<crate::error::GeometryError> for CodecError {
    fn from(e: crate::error::GeometryError) -> Self {
        CodecError::Geometry(e)
    }
}

impl fmt::Display for CodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEnd { needed, available } => write!(
                f,
                "unexpected end of stream: needed {needed} byte(s), {available} available"
            ),
            Self::UnsupportedVersion { found } => write!(
                f,
                "unsupported format version {found} (this build writes and reads version {FORMAT_VERSION})"
            ),
            Self::UnexpectedTag { expected, found } => {
                write!(f, "expected type tag {expected:#04x}, found {found:#04x}")
            }
            Self::UnknownVariant { found } => {
                write!(f, "unknown Geometry variant {found:#04x}")
            }
            Self::InvalidEnumValue { type_name, found } => {
                write!(f, "invalid {type_name} value {found}")
            }
            Self::CountOverflow { count } => {
                write!(f, "count {count} exceeds this platform's usize")
            }
            Self::Geometry(e) => write!(f, "decoded bytes are not a valid geometry: {e}"),
            Self::NestingTooDeep { limit } => {
                write!(f, "geometry collection nests more than {limit} levels deep")
            }
            Self::TrailingBytes { remaining } => {
                write!(f, "{remaining} byte(s) remained after the value")
            }
        }
    }
}

impl core::error::Error for CodecError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Geometry(e) => Some(e),
            _ => None,
        }
    }
}

/// A bounds-checked cursor over a byte slice.
///
/// Every read returns `Result`. There is no indexing and no slicing by a
/// length read from the stream, which is the usual way a hand-written decoder
/// acquires a panic on hostile input. `forbid(unsafe_code)` and
/// `deny(unwrap_used)` prevent memory unsafety but not an out-of-bounds index.
pub struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    /// Start reading at the beginning of `bytes`.
    #[must_use]
    pub fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    /// Bytes not yet consumed.
    #[must_use]
    pub fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.pos)
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], CodecError> {
        let end = self.pos.checked_add(n).ok_or(CodecError::UnexpectedEnd {
            needed: n,
            available: self.remaining(),
        })?;
        let slice = self
            .bytes
            .get(self.pos..end)
            .ok_or(CodecError::UnexpectedEnd {
                needed: n,
                available: self.remaining(),
            })?;
        self.pos = end;
        Ok(slice)
    }

    /// Read one byte.
    ///
    /// # Errors
    ///
    /// [`CodecError::UnexpectedEnd`] if the stream is exhausted.
    pub fn u8(&mut self) -> Result<u8, CodecError> {
        Ok(self.take(1)?[0])
    }

    /// Read one signed byte.
    ///
    /// # Errors
    ///
    /// [`CodecError::UnexpectedEnd`] if the stream is exhausted.
    pub fn i8(&mut self) -> Result<i8, CodecError> {
        // The inverse reinterpretation of the `Dimension` encoding: 0xFF is
        // −1. Named rather than cast with `as`, for the same reason.
        self.u8().map(u8::cast_signed)
    }

    /// Read a little-endian `u32` count and widen it to `usize`.
    ///
    /// # Errors
    ///
    /// [`CodecError::UnexpectedEnd`] if fewer than four bytes remain, or
    /// [`CodecError::CountOverflow`] where `usize` is narrower than `u32`.
    pub fn count(&mut self) -> Result<usize, CodecError> {
        let b = self.take(4)?;
        let raw = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        usize::try_from(raw).map_err(|_| CodecError::CountOverflow { count: raw })
    }

    /// Read a little-endian `f64`.
    ///
    /// The bit pattern is preserved exactly, including NaN payloads; validity
    /// is the constructor's decision, not the reader's.
    ///
    /// # Errors
    ///
    /// [`CodecError::UnexpectedEnd`] if fewer than eight bytes remain.
    pub fn f64(&mut self) -> Result<f64, CodecError> {
        let b = self.take(8)?;
        Ok(f64::from_le_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }
}
