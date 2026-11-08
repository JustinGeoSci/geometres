// src/codec/tests.rs

//! Tests for the native binary format.
//!
//! Organized by the property each group establishes, because the properties
//! are the contract: determinism, exact round-trip, validation on decode, and
//! no panic on any input.
//!
//! The last of those is the one a hand-written decoder most often fails, so
//! truncation and corruption are swept exhaustively rather than sampled.

#![allow(clippy::unwrap_used)]
// test-lint-posture: see the workspace lint comment in Cargo.toml
#![allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap
)]

use super::{CodecError, Decode, Encode, FORMAT_VERSION, decode, encode, tag};
use crate::error::GeometryError;
use crate::topology::de9im::IntersectionMatrix;
use crate::topology::label::{Dimension, Location};
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString,
    MultiPoint, MultiPolygon, Point, Polygon,
};

// ======================================================================
// Fixtures
// ======================================================================

fn c(x: f64, y: f64) -> Coord2d {
    Coord2d::new(x, y)
}

fn seq(pts: &[(f64, f64)]) -> CoordSeq {
    CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())
}

fn ring(pts: &[(f64, f64)]) -> LinearRing {
    LinearRing::new(seq(pts)).unwrap()
}

fn square() -> LinearRing {
    ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])
}

fn hole() -> LinearRing {
    // Clockwise, as a hole must be.
    ring(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])
}

/// A square well away from `holed_polygon`, so the two form a valid
/// multipolygon.
fn far_square() -> Polygon {
    Polygon::new_simple(ring(&[
        (100.0, 100.0),
        (101.0, 100.0),
        (101.0, 101.0),
        (100.0, 101.0),
    ]))
    .unwrap()
}

fn holed_polygon() -> Polygon {
    Polygon::new(square(), vec![hole()]).unwrap()
}

/// Coordinates chosen so that a decimal round trip would lose bits: values
/// needing all 17 significant digits, a subnormal, and signed zero.
fn awkward_values() -> [f64; 6] {
    [
        1_450_000.123_456_789_1,
        -520_000.987_654_321_9,
        f64::MIN_POSITIVE / 2.0,
        -0.0,
        f64::MAX,
        f64::from_bits(1),
    ]
}

// ======================================================================
// Round trip, bit for bit
// ======================================================================

/// Assert a value survives encode/decode as an equal value.
fn round_trip<T: Encode + Decode + PartialEq + core::fmt::Debug>(v: &T) {
    let bytes = encode(v);
    let back: T = decode(&bytes).unwrap();
    assert_eq!(&back, v, "round trip changed the value");
}

#[test]
fn every_2d_type_round_trips() {
    round_trip(&c(1.5, -2.25));
    round_trip(&seq(&[(0.0, 0.0), (1.0, 2.0), (3.0, 4.0)]));
    round_trip(&Point::new(c(1.5, -2.25)).unwrap());
    round_trip(&LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap());
    round_trip(&square());
    round_trip(&holed_polygon());
    round_trip(&MultiPoint::new(vec![
        Point::new(c(1.0, 2.0)).unwrap(),
        Point::new(c(3.0, 4.0)).unwrap(),
    ]));
    round_trip(&MultiLineString::new(vec![
        LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap(),
    ]));
    round_trip(&MultiPolygon::new(vec![holed_polygon()]).unwrap());
}

#[test]
fn every_geometry_variant_round_trips() {
    for g in [
        Geometry::Point(Point::new(c(1.0, 2.0)).unwrap()),
        Geometry::LineString(LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap()),
        Geometry::Polygon(holed_polygon()),
        Geometry::MultiPoint(MultiPoint::new(vec![Point::new(c(1.0, 2.0)).unwrap()])),
        Geometry::MultiLineString(MultiLineString::new(vec![
            LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap(),
        ])),
        Geometry::MultiPolygon(MultiPolygon::new(vec![holed_polygon(), far_square()]).unwrap()),
    ] {
        round_trip(&g);
    }
}

#[test]
fn every_topology_result_round_trips() {
    for l in [Location::Interior, Location::Boundary, Location::Exterior] {
        round_trip(&l);
    }
    for d in [
        Dimension::Empty,
        Dimension::Point,
        Dimension::Line,
        Dimension::Area,
    ] {
        round_trip(&d);
    }
    round_trip(&IntersectionMatrix::from_pattern("212101212").unwrap());
}

/// Bit patterns, not values.
///
/// `assert_eq!` on `f64` would accept a value that differs by a signed zero
/// and would reject a NaN that survived perfectly. Comparing `to_bits` is the
/// property actually wanted: the decoder returns the bits the encoder wrote.
#[test]
fn coordinates_round_trip_bit_for_bit() {
    for x in awkward_values() {
        for y in awkward_values() {
            let original = c(x, y);
            let back: Coord2d = decode(&encode(&original)).unwrap();
            assert_eq!(
                back.x.to_bits(),
                x.to_bits(),
                "x lost bits: {x:?} -> {:?}",
                back.x
            );
            assert_eq!(
                back.y.to_bits(),
                y.to_bits(),
                "y lost bits: {y:?} -> {:?}",
                back.y
            );
        }
    }
}

/// Signed zero is the case a naive comparison misses: `-0.0 == 0.0` is true,
/// so only the bit test can tell whether the sign survived.
#[test]
fn signed_zero_survives() {
    let back: Coord2d = decode(&encode(&c(-0.0, 0.0))).unwrap();
    assert_eq!(back.x.to_bits(), (-0.0_f64).to_bits());
    assert_eq!(back.y.to_bits(), 0.0_f64.to_bits());
}

// ======================================================================
// Determinism
// ======================================================================

#[test]
fn encoding_is_byte_identical_across_repetitions() {
    let p = holed_polygon();
    let first = encode(&p);
    for _ in 0..64 {
        assert_eq!(encode(&p), first, "encoding is not deterministic");
    }
}

/// Two structurally equal values built by different routes must encode
/// identically. The output carries no identity, address or insertion order.
#[test]
fn equal_values_built_differently_encode_identically() {
    let a = Polygon::new(square(), vec![hole()]).unwrap();
    let mut coords = vec![c(0.0, 0.0), c(4.0, 0.0), c(4.0, 4.0), c(0.0, 4.0)];
    coords.push(c(0.0, 0.0)); // closing vertex, which LinearRing::new strips
    let b = Polygon::new(
        LinearRing::new(CoordSeq::new(coords)).unwrap(),
        vec![hole()],
    )
    .unwrap();
    assert_eq!(a, b);
    assert_eq!(encode(&a), encode(&b));
}

/// A fixed value encodes to a fixed, hard-coded byte string.
///
/// Determinism across *runs* and *platforms* cannot be shown by encoding twice
/// in one process, which rules out only nondeterminism within a process, such
/// as hash seeding. Pinning the exact bytes makes the property checkable on
/// any machine and in any future build: if endianness, field order, the
/// version byte or a tag ever changes, this fails.
#[test]
fn a_known_value_encodes_to_known_bytes() {
    // Coord2d { x: 1.0, y: -2.0 }
    // 1.0f64  little-endian = 00 00 00 00 00 00 F0 3F
    // -2.0f64 little-endian = 00 00 00 00 00 00 00 C0
    let expected: Vec<u8> = vec![
        FORMAT_VERSION,
        tag::COORD2D,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0xF0,
        0x3F,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0xC0,
    ];
    assert_eq!(encode(&c(1.0, -2.0)), expected);

    // A triangle ring: count 3, then three coordinates.
    let tri = ring(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);
    let mut want: Vec<u8> = vec![FORMAT_VERSION, tag::LINEAR_RING];
    want.extend_from_slice(&3u32.to_le_bytes());
    for (x, y) in [(0.0_f64, 0.0_f64), (1.0, 0.0), (0.0, 1.0)] {
        want.extend_from_slice(&x.to_le_bytes());
        want.extend_from_slice(&y.to_le_bytes());
    }
    assert_eq!(encode(&tri), want);

    // And the golden bytes decode back to the value they came from.
    assert_eq!(decode::<Coord2d>(&expected).unwrap(), c(1.0, -2.0));
    assert_eq!(decode::<LinearRing>(&want).unwrap(), tri);
}

/// The encoded length is exactly what the format description predicts, which
/// catches an accidental extra field or padding byte.
#[test]
fn encoded_length_matches_the_documented_format() {
    // version + tag + count + 4 coords * 16 bytes
    assert_eq!(encode(&square()).len(), 1 + 1 + 4 + 4 * 16);
    // version + tag + exterior(count + 4*16) + holes count + hole(count + 4*16)
    assert_eq!(
        encode(&holed_polygon()).len(),
        1 + 1 + (4 + 64) + 4 + (4 + 64)
    );
    // version + tag + x + y
    assert_eq!(encode(&Point::new(c(0.0, 0.0)).unwrap()).len(), 1 + 1 + 16);
}

// ======================================================================
// Ring closure convention
// ======================================================================

/// The encoder writes `n` distinct vertices and no closing vertex.
///
/// Stated as a test because the opposite convention produces a known class of
/// boundary bugs: a format that writes the closing vertex and a constructor
/// that strips it grow the ring on every round trip.
#[test]
fn a_ring_encodes_its_distinct_vertices_with_no_closing_vertex() {
    let r = square();
    assert_eq!(r.coords().len(), 4);
    let bytes = encode(&r);
    let count = u32::from_le_bytes([bytes[2], bytes[3], bytes[4], bytes[5]]);
    assert_eq!(count, 4, "a closing vertex was written");

    // And it survives repeated round trips without growing.
    let mut current = r;
    for pass in 0..5 {
        current = decode(&encode(&current)).unwrap();
        assert_eq!(current.coords().len(), 4, "ring grew on pass {pass}");
    }
}

// ======================================================================
// Validation on decode
// ======================================================================

/// Hand-build a document, so values the constructors would reject can be put
/// on the wire. Encoding a real value cannot produce them.
fn doc(tag_byte: u8, body: &[u8]) -> Vec<u8> {
    let mut v = vec![FORMAT_VERSION, tag_byte];
    v.extend_from_slice(body);
    v
}

fn coords_body(pts: &[(f64, f64)]) -> Vec<u8> {
    let mut b = (pts.len() as u32).to_le_bytes().to_vec();
    for &(x, y) in pts {
        b.extend_from_slice(&x.to_le_bytes());
        b.extend_from_slice(&y.to_le_bytes());
    }
    b
}

#[test]
fn a_nan_point_is_refused() {
    let body = {
        let mut b = f64::NAN.to_le_bytes().to_vec();
        b.extend_from_slice(&0.0_f64.to_le_bytes());
        b
    };
    assert_eq!(
        decode::<Point>(&doc(tag::POINT, &body)).unwrap_err(),
        CodecError::Geometry(GeometryError::NonFiniteCoordinate {
            coord_index: 0,
            axis: 0
        })
    );
}

#[test]
fn a_two_vertex_ring_is_refused() {
    let body = coords_body(&[(0.0, 0.0), (1.0, 1.0)]);
    assert_eq!(
        decode::<LinearRing>(&doc(tag::LINEAR_RING, &body)).unwrap_err(),
        CodecError::Geometry(GeometryError::InsufficientVertices {
            required: 3,
            found: 2
        })
    );
}

#[test]
fn a_clockwise_exterior_is_refused() {
    // The square's vertices in reverse: a valid ring, wrong winding.
    let body = {
        let mut b = coords_body(&[(0.0, 4.0), (4.0, 4.0), (4.0, 0.0), (0.0, 0.0)]);
        b.extend_from_slice(&0u32.to_le_bytes()); // no holes
        b
    };
    assert_eq!(
        decode::<Polygon>(&doc(tag::POLYGON, &body)).unwrap_err(),
        CodecError::Geometry(GeometryError::IncorrectWindingOrder)
    );
}

#[test]
fn a_counterclockwise_hole_is_refused() {
    let mut body = coords_body(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]);
    body.extend_from_slice(&1u32.to_le_bytes());
    // Counterclockwise, which is wrong for a hole.
    body.extend_from_slice(&coords_body(&[
        (1.0, 1.0),
        (2.0, 1.0),
        (2.0, 2.0),
        (1.0, 2.0),
    ]));
    assert!(matches!(
        decode::<Polygon>(&doc(tag::POLYGON, &body)).unwrap_err(),
        CodecError::Geometry(GeometryError::IncorrectHoleWindingOrder { .. })
    ));
}

#[test]
fn a_duplicate_consecutive_vertex_is_refused() {
    let body = coords_body(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);
    assert!(matches!(
        decode::<LinearRing>(&doc(tag::LINEAR_RING, &body)).unwrap_err(),
        CodecError::Geometry(GeometryError::ConsecutiveDuplicateVertices { .. })
    ));
}

// ======================================================================
// Header handling
// ======================================================================

/// The version field exists and is checked.
///
/// This test fails if the version byte is removed from the format, which is
/// the point: a format without a version cannot be revised without orphaning
/// every stored artifact.
#[test]
fn the_version_field_is_present_and_enforced() {
    let good = encode(&Point::new(c(1.0, 2.0)).unwrap());
    assert_eq!(good[0], FORMAT_VERSION, "byte 0 must be the version");

    for bad_version in [0u8, FORMAT_VERSION.wrapping_add(1), 0xFF] {
        if bad_version == FORMAT_VERSION {
            continue;
        }
        let mut bytes = good.clone();
        bytes[0] = bad_version;
        assert_eq!(
            decode::<Point>(&bytes).unwrap_err(),
            CodecError::UnsupportedVersion { found: bad_version },
            "version {bad_version} must be refused"
        );
    }
}

#[test]
fn a_document_of_the_wrong_type_is_refused() {
    let bytes = encode(&Point::new(c(1.0, 2.0)).unwrap());
    assert_eq!(
        decode::<LineString>(&bytes).unwrap_err(),
        CodecError::UnexpectedTag {
            expected: tag::LINE_STRING,
            found: tag::POINT
        }
    );
}

/// A stream describing two parts that share an edge is refused, not decoded
/// into a multipolygon whose matrix would be meaningless.
#[test]
fn a_multipolygon_stream_with_overlapping_parts_is_refused() {
    let part = encode(&holed_polygon());
    let mut body = 2u32.to_le_bytes().to_vec();
    body.extend_from_slice(&part[2..]);
    body.extend_from_slice(&part[2..]);
    for (t, b) in [
        (tag::MULTI_POLYGON, body.clone()),
        (tag::GEOMETRY, [&[tag::MULTI_POLYGON][..], &body].concat()),
    ] {
        let err = if t == tag::GEOMETRY {
            decode::<Geometry>(&doc(t, &b)).unwrap_err()
        } else {
            decode::<MultiPolygon>(&doc(t, &b)).unwrap_err()
        };
        assert_eq!(
            err,
            CodecError::Geometry(GeometryError::OverlappingParts {
                first: 0,
                second: 1
            })
        );
    }
}

#[test]
fn an_unknown_geometry_variant_is_refused() {
    for unknown in [0x00u8, 0x7F, 0xFE, 0xFF] {
        let bytes = doc(tag::GEOMETRY, &[unknown]);
        assert_eq!(
            decode::<Geometry>(&bytes).unwrap_err(),
            CodecError::UnknownVariant { found: unknown },
            "variant {unknown:#04x} must be refused"
        );
    }
}

#[test]
fn an_out_of_range_enum_value_is_refused() {
    assert!(matches!(
        decode::<Location>(&doc(tag::LOCATION, &[3])).unwrap_err(),
        CodecError::InvalidEnumValue { .. }
    ));
    assert!(matches!(
        decode::<Dimension>(&doc(tag::DIMENSION, &[7])).unwrap_err(),
        CodecError::InvalidEnumValue { .. }
    ));
    // -2 as a byte is 0xFE.
    assert!(matches!(
        decode::<Dimension>(&doc(tag::DIMENSION, &[0xFE])).unwrap_err(),
        CodecError::InvalidEnumValue { .. }
    ));
}

#[test]
fn trailing_bytes_are_refused() {
    let mut bytes = encode(&Point::new(c(1.0, 2.0)).unwrap());
    bytes.push(0);
    assert_eq!(
        decode::<Point>(&bytes).unwrap_err(),
        CodecError::TrailingBytes { remaining: 1 }
    );
}

// ======================================================================
// Hostile input
// ======================================================================

/// Every prefix of a valid document must error, never panic, never yield a
/// value.
///
/// Swept exhaustively rather than sampled: a decoder's panic usually lives at
/// one specific offset, and sampling is how it survives a test suite.
#[test]
fn every_truncation_errors_without_panicking() {
    let documents: Vec<Vec<u8>> = vec![
        encode(&c(1.0, 2.0)),
        encode(&Point::new(c(1.0, 2.0)).unwrap()),
        encode(&square()),
        encode(&holed_polygon()),
        encode(&MultiPolygon::new(vec![holed_polygon(), far_square()]).unwrap()),
        encode(&Geometry::Polygon(holed_polygon())),
        encode(&IntersectionMatrix::from_pattern("212101212").unwrap()),
    ];

    for full in &documents {
        for cut in 0..full.len() {
            let prefix = &full[..cut];
            // The type does not matter for the panic property; Polygon is the
            // deepest nesting and so the most likely to index out of bounds.
            let _ = decode::<Polygon>(prefix);
            let _ = decode::<MultiPolygon>(prefix);
            let _ = decode::<Geometry>(prefix);
            let _ = decode::<Point>(prefix);
        }
        // A prefix of a document can never decode as that same complete value.
        for cut in 0..full.len() {
            assert!(
                decode::<Polygon>(&full[..cut]).is_err() || cut == full.len(),
                "a truncated document decoded at cut {cut}"
            );
        }
    }
}

/// Flipping any single bit must produce either a valid geometry or an error,
/// never a panic.
///
/// Some flips land in a coordinate mantissa and produce a perfectly valid
/// geometry, which is expected: this is a panic-freedom test, not an integrity
/// one. The format carries no checksum and does not claim to detect
/// corruption.
#[test]
fn every_single_bit_flip_errors_or_decodes_without_panicking() {
    let documents: Vec<Vec<u8>> = vec![
        encode(&Point::new(c(1.0, 2.0)).unwrap()),
        encode(&square()),
        encode(&holed_polygon()),
        encode(&Geometry::MultiPoint(MultiPoint::new(vec![
            Point::new(c(1.0, 2.0)).unwrap(),
        ]))),
    ];

    for full in &documents {
        for byte_index in 0..full.len() {
            for bit in 0..8u32 {
                let mut corrupted = full.clone();
                corrupted[byte_index] ^= 1 << bit;
                // Any outcome except a panic is acceptable.
                let _ = decode::<Point>(&corrupted);
                let _ = decode::<LinearRing>(&corrupted);
                let _ = decode::<Polygon>(&corrupted);
                let _ = decode::<Geometry>(&corrupted);
            }
        }
    }
}

/// A huge count must not be believed.
///
/// This is the classic hand-written-decoder crash: `Vec::with_capacity` on an
/// attacker-supplied count aborts the process rather than returning an error.
/// The decoder clamps its capacity hint, so a four-byte lie costs a bounded
/// allocation and then hits end-of-stream.
#[test]
fn an_enormous_count_errors_rather_than_exhausting_memory() {
    for count in [u32::MAX, u32::MAX / 2, 1 << 30, 1 << 20] {
        let body = count.to_le_bytes().to_vec();
        assert!(matches!(
            decode::<CoordSeq>(&doc(tag::COORD_SEQ, &body)).unwrap_err(),
            CodecError::UnexpectedEnd { .. }
        ));
        assert!(matches!(
            decode::<MultiPolygon>(&doc(tag::MULTI_POLYGON, &body)).unwrap_err(),
            CodecError::UnexpectedEnd { .. }
        ));
    }
}

#[test]
fn an_empty_stream_errors() {
    assert!(matches!(
        decode::<Point>(&[]).unwrap_err(),
        CodecError::UnexpectedEnd { .. }
    ));
}

// ======================================================================
// GeometryCollection
// ======================================================================

/// One collection holding an element of every type, empty ones included.
fn every_kind_collection() -> GeometryCollection {
    GeometryCollection::new(vec![
        Geometry::Point(Point::new(c(1.0, 2.0)).unwrap()),
        Geometry::LineString(LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap()),
        Geometry::Polygon(holed_polygon()),
        Geometry::MultiPoint(MultiPoint::new(vec![Point::new(c(-0.0, 5.0)).unwrap()])),
        Geometry::MultiLineString(MultiLineString::new(vec![
            LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap(),
        ])),
        Geometry::MultiPolygon(MultiPolygon::new(vec![holed_polygon(), far_square()]).unwrap()),
        Geometry::MultiPoint(MultiPoint::empty()),
        Geometry::MultiLineString(MultiLineString::empty()),
        Geometry::MultiPolygon(MultiPolygon::empty()),
        Geometry::GeometryCollection(GeometryCollection::empty()),
    ])
    .unwrap()
}

/// A collection `depth` levels deep, `leaf` at the bottom.
fn nested(depth: usize, leaf: Geometry) -> GeometryCollection {
    let mut g = GeometryCollection::new(vec![leaf]).unwrap();
    for _ in 1..depth {
        g = GeometryCollection::new(vec![Geometry::GeometryCollection(g)]).unwrap();
    }
    g
}

/// The bytes of a chain of `levels` collection headers, each holding the
/// next, the last holding one point: a hand-built document, so depths the
/// type would refuse can be put on the wire.
fn collection_chain_body(levels: usize) -> Vec<u8> {
    let mut b = Vec::new();
    for level in 0..levels {
        if level > 0 {
            b.push(tag::GEOMETRY_COLLECTION);
        }
        b.extend_from_slice(&1u32.to_le_bytes());
    }
    b.push(tag::POINT);
    b.extend_from_slice(&1.0f64.to_le_bytes());
    b.extend_from_slice(&2.0f64.to_le_bytes());
    b
}

#[test]
fn the_collection_tag_is_0x11_and_the_others_are_unchanged() {
    assert_eq!(tag::GEOMETRY_COLLECTION, 0x11);
    assert_eq!(GeometryCollection::TAG, 0x11);
    assert_eq!(FORMAT_VERSION, 1);
    for (t, v) in [
        (tag::POINT, 0x03),
        (tag::POLYGON, 0x06),
        (tag::MULTI_POLYGON, 0x09),
        (tag::GEOMETRY, 0x0A),
        (tag::ENVELOPE, 0x0F),
        (tag::TESSELLATE_METHOD, 0x10),
        (tag::COORD3D, 0x20),
        (tag::DOWNSTREAM_FIRST, 0x40),
    ] {
        assert_eq!(t, v);
    }
}

#[test]
fn collections_of_every_element_type_round_trip_byte_for_byte() {
    for gc in [
        GeometryCollection::empty(),
        every_kind_collection(),
        GeometryCollection::new(vec![
            Geometry::GeometryCollection(every_kind_collection()),
            Geometry::Point(Point::new(c(awkward_values()[2], awkward_values()[5])).unwrap()),
        ])
        .unwrap(),
    ] {
        let bytes = encode(&gc);
        let back: GeometryCollection = decode(&bytes).unwrap();
        assert_eq!(back, gc);
        assert_eq!(encode(&back), bytes, "re-encoding changed the bytes");
        // The same value as a `Geometry` variant.
        let g = Geometry::GeometryCollection(gc);
        let gbytes = encode(&g);
        let gback: Geometry = decode(&gbytes).unwrap();
        assert_eq!(gback, g);
        assert_eq!(encode(&gback), gbytes);
    }
}

#[test]
fn a_known_collection_encodes_to_known_bytes() {
    assert_eq!(
        encode(&GeometryCollection::empty()),
        vec![FORMAT_VERSION, 0x11, 0, 0, 0, 0]
    );
    // GEOMETRYCOLLECTION(GEOMETRYCOLLECTION(POINT(1 2))): the outer count,
    // the inner variant byte and count, the point's variant byte and body.
    let gc = nested(2, Geometry::Point(Point::new(c(1.0, 2.0)).unwrap()));
    assert_eq!(encode(&gc), doc(0x11, &collection_chain_body(2)));
}

#[test]
fn a_collection_nested_to_the_limit_round_trips() {
    let gc = nested(
        GeometryCollection::MAX_NESTING,
        Geometry::Polygon(holed_polygon()),
    );
    assert_eq!(gc.depth(), 32);
    round_trip(&gc);
    round_trip(&Geometry::GeometryCollection(gc));
}

#[test]
fn a_document_one_level_deeper_than_the_limit_is_refused() {
    let at_limit = doc(tag::GEOMETRY_COLLECTION, &collection_chain_body(32));
    let gc: GeometryCollection = decode(&at_limit).unwrap();
    assert_eq!(gc.depth(), 32);

    let beyond = doc(tag::GEOMETRY_COLLECTION, &collection_chain_body(33));
    assert_eq!(
        decode::<GeometryCollection>(&beyond).unwrap_err(),
        CodecError::NestingTooDeep { limit: 32 }
    );
    // Through the `Geometry` variant path too.
    let mut body = vec![tag::GEOMETRY_COLLECTION];
    body.extend_from_slice(&collection_chain_body(33));
    assert_eq!(
        decode::<Geometry>(&doc(tag::GEOMETRY, &body)).unwrap_err(),
        CodecError::NestingTooDeep { limit: 32 }
    );
}

#[test]
fn the_depth_is_refused_before_the_deeper_content_is_read() {
    // A chain of headers far past the limit, with nothing after the last:
    // the refusal comes at level 33, not as an end-of-stream error at the
    // bottom.
    let mut body = Vec::new();
    for level in 0..10_000 {
        if level > 0 {
            body.push(tag::GEOMETRY_COLLECTION);
        }
        body.extend_from_slice(&1u32.to_le_bytes());
    }
    assert_eq!(
        decode::<GeometryCollection>(&doc(tag::GEOMETRY_COLLECTION, &body)).unwrap_err(),
        CodecError::NestingTooDeep { limit: 32 }
    );
}

#[test]
fn an_invalid_element_inside_a_collection_is_refused_by_its_constructor() {
    // A collection of one clockwise "polygon".
    let mut body = 1u32.to_le_bytes().to_vec();
    body.push(tag::POLYGON);
    body.extend_from_slice(&coords_body(&[
        (0.0, 0.0),
        (0.0, 1.0),
        (1.0, 1.0),
        (1.0, 0.0),
    ]));
    body.extend_from_slice(&0u32.to_le_bytes());
    assert_eq!(
        decode::<GeometryCollection>(&doc(tag::GEOMETRY_COLLECTION, &body)).unwrap_err(),
        CodecError::Geometry(GeometryError::IncorrectWindingOrder)
    );
}

#[test]
fn an_unknown_variant_inside_a_collection_is_refused() {
    for found in [0x00, 0x0A, 0x0E, 0x10, 0x12, 0x20, 0x40, 0xFF] {
        let mut body = 1u32.to_le_bytes().to_vec();
        body.push(found);
        assert_eq!(
            decode::<GeometryCollection>(&doc(tag::GEOMETRY_COLLECTION, &body)).unwrap_err(),
            CodecError::UnknownVariant { found },
            "variant {found:#04x}"
        );
    }
}

#[test]
fn an_enormous_collection_count_errors_rather_than_exhausting_memory() {
    for count in [u32::MAX, 1 << 30, 1 << 20] {
        let body = count.to_le_bytes().to_vec();
        assert!(matches!(
            decode::<GeometryCollection>(&doc(tag::GEOMETRY_COLLECTION, &body)).unwrap_err(),
            CodecError::UnexpectedEnd { .. }
        ));
        // And inside a nested header.
        let mut nested_body = 1u32.to_le_bytes().to_vec();
        nested_body.push(tag::GEOMETRY_COLLECTION);
        nested_body.extend_from_slice(&count.to_le_bytes());
        assert!(matches!(
            decode::<GeometryCollection>(&doc(tag::GEOMETRY_COLLECTION, &nested_body)).unwrap_err(),
            CodecError::UnexpectedEnd { .. }
        ));
    }
}

#[test]
fn trailing_bytes_after_a_collection_are_refused() {
    let mut bytes = encode(&every_kind_collection());
    bytes.push(0);
    assert_eq!(
        decode::<GeometryCollection>(&bytes).unwrap_err(),
        CodecError::TrailingBytes { remaining: 1 }
    );
}

#[test]
fn collection_documents_survive_every_truncation_and_bit_flip() {
    let documents = [
        encode(&every_kind_collection()),
        encode(&Geometry::GeometryCollection(nested(
            4,
            Geometry::Polygon(holed_polygon()),
        ))),
    ];
    for full in &documents {
        for cut in 0..full.len() {
            assert!(decode::<GeometryCollection>(&full[..cut]).is_err());
            assert!(decode::<Geometry>(&full[..cut]).is_err());
        }
        for byte_index in 0..full.len() {
            for bit in 0..8u32 {
                let mut corrupted = full.clone();
                corrupted[byte_index] ^= 1 << bit;
                for value in [
                    decode::<GeometryCollection>(&corrupted)
                        .ok()
                        .map(Geometry::GeometryCollection),
                    decode::<Geometry>(&corrupted).ok(),
                ]
                .into_iter()
                .flatten()
                {
                    // Whatever decodes re-encodes to a document that decodes
                    // to the same value.
                    let again: Geometry = decode(&encode(&value)).unwrap();
                    assert_eq!(again, value);
                }
            }
        }
    }
}

#[test]
fn a_ring_stream_with_a_closing_vertex_is_refused() {
    // The encoder never writes one; accepting it would give a ring two
    // encodings.
    let closed = coords_body(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 0.0)]);
    assert_eq!(
        decode::<LinearRing>(&doc(tag::LINEAR_RING, &closed)).unwrap_err(),
        CodecError::Geometry(GeometryError::ConsecutiveDuplicateVertices { index: 3 })
    );
    // As a polygon's exterior too.
    let mut body = coords_body(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]);
    body.extend_from_slice(&0u32.to_le_bytes());
    assert_eq!(
        decode::<Polygon>(&doc(tag::POLYGON, &body)).unwrap_err(),
        CodecError::Geometry(GeometryError::ConsecutiveDuplicateVertices { index: 4 })
    );
    // The same closing vertex with the other zero is the same point.
    let signed = coords_body(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (-0.0, 0.0)]);
    assert!(decode::<LinearRing>(&doc(tag::LINEAR_RING, &signed)).is_err());
}

/// Decodes `bytes` as `T`; if it decodes, the value re-encodes to exactly
/// `bytes`.
fn canonical<T: Encode + Decode + PartialEq + core::fmt::Debug>(bytes: &[u8]) {
    if let Ok(v) = decode::<T>(bytes) {
        assert_eq!(encode(&v), bytes, "{v:?} re-encoded differently");
    }
}

/// The fuzz target's property, swept here over every truncation and single
/// bit flip of a document of each type: whatever decodes re-encodes to the
/// bytes it was decoded from, so every value has exactly one encoding.
#[test]
fn whatever_decodes_re_encodes_to_the_bytes_read() {
    let documents: Vec<Vec<u8>> = vec![
        encode(&c(1.0, -2.0)),
        encode(&seq(&[(0.0, 0.0), (1.0, 2.0)])),
        encode(&Point::new(c(1.0, 2.0)).unwrap()),
        encode(&LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap()),
        encode(&square()),
        encode(&holed_polygon()),
        encode(&MultiPoint::new(vec![Point::new(c(1.0, 2.0)).unwrap()])),
        encode(&MultiLineString::new(vec![
            LineString::new(seq(&[(0.0, 0.0), (1.0, 1.0)])).unwrap(),
        ])),
        encode(&MultiPolygon::new(vec![holed_polygon(), far_square()]).unwrap()),
        encode(&Geometry::Polygon(holed_polygon())),
        encode(&every_kind_collection()),
        encode(&crate::types::Envelope::new(c(0.0, -1.0), c(2.0, 3.0)).unwrap()),
        encode(&IntersectionMatrix::from_pattern("212101212").unwrap()),
        encode(&Dimension::Line),
        encode(&Location::Boundary),
    ];
    let check = |bytes: &[u8]| {
        canonical::<Coord2d>(bytes);
        canonical::<CoordSeq>(bytes);
        canonical::<Point>(bytes);
        canonical::<LineString>(bytes);
        canonical::<LinearRing>(bytes);
        canonical::<Polygon>(bytes);
        canonical::<MultiPoint>(bytes);
        canonical::<MultiLineString>(bytes);
        canonical::<MultiPolygon>(bytes);
        canonical::<Geometry>(bytes);
        canonical::<GeometryCollection>(bytes);
        canonical::<crate::types::Envelope>(bytes);
        canonical::<IntersectionMatrix>(bytes);
        canonical::<Dimension>(bytes);
        canonical::<Location>(bytes);
    };
    for full in &documents {
        check(full);
        for cut in 0..full.len() {
            check(&full[..cut]);
        }
        for byte_index in 0..full.len() {
            for bit in 0..8u32 {
                let mut corrupted = full.clone();
                corrupted[byte_index] ^= 1 << bit;
                check(&corrupted);
            }
        }
    }
}

// ======================================================================
// 3D
// ======================================================================

#[cfg(feature = "3d")]
mod three_d {
    use super::*;
    use crate::types3d::{
        Coord3d, CoordSeq3d, LineString3d, LinearRing3d, MultiPoint3d, MultiPolygon3d, Point3d,
        Polygon3d,
    };

    fn c3(x: f64, y: f64, z: f64) -> Coord3d {
        Coord3d::new(x, y, z)
    }

    fn seq3(pts: &[(f64, f64, f64)]) -> CoordSeq3d {
        CoordSeq3d::new(pts.iter().map(|&(x, y, z)| c3(x, y, z)).collect())
    }

    fn ring3() -> LinearRing3d {
        LinearRing3d::new(seq3(&[
            (0.0, 0.0, 0.0),
            (4.0, 0.0, 0.0),
            (4.0, 4.0, 0.0),
            (0.0, 4.0, 0.0),
        ]))
        .unwrap()
    }

    /// A square of side 4 at `(ox, oy)`, height `z`.
    fn square3(ox: f64, oy: f64, z: f64) -> Polygon3d {
        Polygon3d::new_simple(
            LinearRing3d::new(seq3(&[
                (ox, oy, z),
                (ox + 4.0, oy, z),
                (ox + 4.0, oy + 4.0, z),
                (ox, oy + 4.0, z),
            ]))
            .unwrap(),
        )
        .unwrap()
    }

    /// Two parts, one holed, at different heights, with awkward `z` values
    /// (`-0.0`, a subnormal, the smallest normal). The `z` extent was kept
    /// small for the degeneracy threshold of the time, which used the 3D
    /// perimeter; the exact test no longer depends on it.
    fn multi3() -> MultiPolygon3d {
        let hole = LinearRing3d::new(seq3(&[
            (1.0, 1.0, -0.0),
            (1.0, 3.0, 5e-324),
            (3.0, 3.0, 1.5),
            (3.0, 1.0, f64::MIN_POSITIVE),
        ]))
        .unwrap();
        let holed = Polygon3d::new(ring3(), vec![hole]).unwrap();
        MultiPolygon3d::new(vec![holed, square3(10.0, 0.0, 7.5)]).unwrap()
    }

    #[test]
    fn every_3d_type_round_trips() {
        round_trip(&c3(1.0, 2.0, 3.0));
        round_trip(&seq3(&[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]));
        round_trip(&Point3d::new(c3(1.0, 2.0, 3.0)).unwrap());
        round_trip(&LineString3d::new(seq3(&[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)])).unwrap());
        round_trip(&ring3());
        round_trip(&Polygon3d::new_simple(ring3()).unwrap());
        round_trip(&multi3());
        round_trip(&MultiPolygon3d::empty());
    }

    /// `MultiPolygon3d` round-trips bit for bit,
    /// including `-0.0`, the smallest subnormal and the smallest normal in `z`.
    #[test]
    fn multi_polygon3d_round_trips_bit_for_bit() {
        let mp = multi3();
        let back: MultiPolygon3d = decode(&encode(&mp)).unwrap();
        for (p, q) in mp.polygons().iter().zip(back.polygons()) {
            let rings = |p: &Polygon3d| {
                std::iter::once(p.exterior())
                    .chain(p.holes())
                    .flat_map(|r| r.coords().iter().copied())
                    .collect::<Vec<_>>()
            };
            for (a, b) in rings(p).into_iter().zip(rings(q)) {
                assert_eq!(a.x.to_bits(), b.x.to_bits());
                assert_eq!(a.y.to_bits(), b.y.to_bits());
                assert_eq!(a.z.to_bits(), b.z.to_bits());
            }
        }
        // Tag and layout: version, tag 0x28, count 2, then each part's body.
        let bytes = encode(&mp);
        assert_eq!(
            &bytes[..6],
            &[FORMAT_VERSION, tag::MULTI_POLYGON3D, 2, 0, 0, 0]
        );
        let mut expected = Vec::new();
        for p in mp.polygons() {
            p.encode_body(&mut expected);
        }
        assert_eq!(&bytes[6..], &expected[..]);
    }

    /// Parts that overlap in plan are refused on the wire exactly as at
    /// construction: the decoder goes through `MultiPolygon3d::new`.
    #[test]
    fn overlapping_parts_on_the_wire_are_refused() {
        let mut bytes = vec![FORMAT_VERSION, tag::MULTI_POLYGON3D];
        bytes.extend_from_slice(&2u32.to_le_bytes());
        // The same footprint at two heights: disjoint in space, overlapping
        // in plan.
        square3(0.0, 0.0, 0.0).encode_body(&mut bytes);
        square3(0.0, 0.0, 9.0).encode_body(&mut bytes);
        assert_eq!(
            decode::<MultiPolygon3d>(&bytes).unwrap_err(),
            CodecError::Geometry(GeometryError::OverlappingParts {
                first: 0,
                second: 1
            })
        );
        // And an edge-sharing pair.
        let mut bytes = vec![FORMAT_VERSION, tag::MULTI_POLYGON3D];
        bytes.extend_from_slice(&2u32.to_le_bytes());
        square3(0.0, 0.0, 0.0).encode_body(&mut bytes);
        square3(4.0, 0.0, 0.0).encode_body(&mut bytes);
        assert!(matches!(
            decode::<MultiPolygon3d>(&bytes).unwrap_err(),
            CodecError::Geometry(GeometryError::OverlappingParts { .. })
        ));
    }

    /// The truncation and single-bit-flip sweeps of the 2D suite, over a
    /// `MultiPolygon3d` document.
    #[test]
    fn multi_polygon3d_survives_every_truncation_and_bit_flip() {
        let full = encode(&multi3());
        for cut in 0..full.len() {
            assert!(
                decode::<MultiPolygon3d>(&full[..cut]).is_err(),
                "a truncated document decoded at cut {cut}"
            );
            let _ = decode::<Polygon3d>(&full[..cut]);
        }
        for byte_index in 0..full.len() {
            for bit in 0..8u32 {
                let mut corrupted = full.clone();
                corrupted[byte_index] ^= 1 << bit;
                // Any outcome except a panic is acceptable.
                let _ = decode::<MultiPolygon3d>(&corrupted);
                let _ = decode::<Polygon3d>(&corrupted);
            }
        }
    }

    #[test]
    fn an_enormous_3d_count_errors_rather_than_exhausting_memory() {
        for count in [u32::MAX, 1 << 30] {
            let body = count.to_le_bytes().to_vec();
            assert!(matches!(
                decode::<MultiPolygon3d>(&doc(tag::MULTI_POLYGON3D, &body)).unwrap_err(),
                CodecError::UnexpectedEnd { .. }
            ));
        }
    }

    #[test]
    fn three_d_coordinates_round_trip_bit_for_bit() {
        for z in awkward_values() {
            let back: Coord3d = decode(&encode(&c3(1.0, 2.0, z))).unwrap();
            assert_eq!(back.z.to_bits(), z.to_bits(), "z lost bits");
        }
    }

    /// A NaN `Point3d` is refused on decode, as its 2D counterpart is.
    ///
    /// The decoder refuses through `Point3d::new`, which checks every axis, so
    /// the validate-on-decode guarantee holds in 3D as it does in 2D. The
    /// document is hand-built, because encoding a real `Point3d` cannot
    /// produce a non-finite one.
    #[test]
    fn a_non_finite_point3d_is_refused_on_decode() {
        for (x, y, z, axis) in [
            (f64::NAN, 0.0, 0.0, 0u8),
            (0.0, f64::INFINITY, 0.0, 1),
            (0.0, 0.0, f64::NEG_INFINITY, 2),
        ] {
            let mut bytes = vec![FORMAT_VERSION, tag::POINT3D];
            for v in [x, y, z] {
                bytes.extend_from_slice(&v.to_le_bytes());
            }
            assert_eq!(
                decode::<Point3d>(&bytes).unwrap_err(),
                CodecError::Geometry(GeometryError::NonFiniteCoordinate {
                    coord_index: 0,
                    axis
                }),
                "a non-finite z on axis {axis} must be refused, not decoded"
            );
        }
    }

    /// The same holds through a collection: an invalid element cannot be
    /// smuggled in via `MultiPoint3d`, because each element decodes through
    /// `Point3d::new`.
    #[test]
    fn a_non_finite_element_is_refused_through_multi_point3d() {
        let mut bytes = vec![FORMAT_VERSION, tag::MULTI_POINT3D];
        bytes.extend_from_slice(&2u32.to_le_bytes());
        for (x, y, z) in [(1.0_f64, 2.0_f64, 3.0_f64), (f64::NAN, 0.0, 0.0)] {
            for v in [x, y, z] {
                bytes.extend_from_slice(&v.to_le_bytes());
            }
        }
        assert_eq!(
            decode::<MultiPoint3d>(&bytes).unwrap_err(),
            CodecError::Geometry(GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 0
            })
        );
    }

    #[test]
    fn three_d_truncation_errors_without_panicking() {
        let full = encode(&Polygon3d::new_simple(ring3()).unwrap());
        for cut in 0..full.len() {
            assert!(decode::<Polygon3d>(&full[..cut]).is_err());
        }
    }
}
