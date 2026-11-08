// fuzz/fuzz_targets/codec_decode.rs

//! Fuzz target: the native codec's decoder.
//!
//! The input's first byte picks a type, and the rest is that type's body
//! under a correct version and tag, so the fuzzer reaches every decoder
//! rather than spending its time on headers. The whole input is also decoded
//! as every type, header included.
//!
//! The contract checked, for every type:
//!
//! - decoding never panics, whatever the bytes;
//! - whatever decodes re-encodes to exactly the bytes it was decoded from
//!   (every value has one encoding), and those bytes decode again to the
//!   same value, bit for bit.
//!
//! "The same value" is judged by the encoding, not by `PartialEq`. The raw
//! coordinate types (`Coord2d`, `CoordSeq` and their 3D forms) are not
//! validated, so they can hold NaN, and NaN is not `==` to itself. The
//! encoding writes every `f64` by its bit pattern, so equal encodings are
//! bitwise-equal values; for every validated type the two notions agree.

#![no_main]

use geometres::codec::{Decode, Encode, FORMAT_VERSION, decode, encode};
use geometres::topology::de9im::IntersectionMatrix;
use geometres::topology::label::{Dimension, Location};
use geometres::types::{
    Coord2d, CoordSeq, Envelope, Geometry, GeometryCollection, LineString, LinearRing,
    MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};
use geometres::types3d::{
    Coord3d, CoordSeq3d, LineString3d, LinearRing3d, MultiLineString3d, MultiPoint3d,
    MultiPolygon3d, Point3d, Polygon3d,
};
use libfuzzer_sys::fuzz_target;

/// Decodes `bytes` as `T`; a value that decodes must re-encode to `bytes`
/// and decode again to itself, bit for bit (compared by encoding; see the
/// module documentation for why not by `PartialEq`).
fn check<T: Encode + Decode + PartialEq + core::fmt::Debug>(bytes: &[u8]) {
    if let Ok(v) = decode::<T>(bytes) {
        let again = encode(&v);
        assert_eq!(again, bytes, "{v:?} re-encoded to different bytes");
        let w: T = decode(&again).expect("a re-encoding must decode");
        assert_eq!(encode(&w), again, "{v:?} decoded again to {w:?}");
        // Where `PartialEq` is reflexive (no NaN inside), it agrees.
        // `v == v` is exactly the NaN test, so the lint is waived here.
        #[allow(clippy::eq_op)]
        if v == v {
            assert_eq!(w, v);
        }
    }
}

/// `body` under a correct header for `T`, then checked.
fn check_body<T: Encode + Decode + PartialEq + core::fmt::Debug>(body: &[u8]) {
    let mut doc = Vec::with_capacity(body.len() + 2);
    doc.push(FORMAT_VERSION);
    doc.push(T::TAG);
    doc.extend_from_slice(body);
    check::<T>(&doc);
}

macro_rules! every_type {
    ($f:ident, $bytes:expr) => {{
        $f::<Coord2d>($bytes);
        $f::<CoordSeq>($bytes);
        $f::<Point>($bytes);
        $f::<LineString>($bytes);
        $f::<LinearRing>($bytes);
        $f::<Polygon>($bytes);
        $f::<MultiPoint>($bytes);
        $f::<MultiLineString>($bytes);
        $f::<MultiPolygon>($bytes);
        $f::<Geometry>($bytes);
        $f::<GeometryCollection>($bytes);
        $f::<Envelope>($bytes);
        $f::<IntersectionMatrix>($bytes);
        $f::<Dimension>($bytes);
        $f::<Location>($bytes);
        $f::<Coord3d>($bytes);
        $f::<CoordSeq3d>($bytes);
        $f::<Point3d>($bytes);
        $f::<LineString3d>($bytes);
        $f::<LinearRing3d>($bytes);
        $f::<Polygon3d>($bytes);
        $f::<MultiPoint3d>($bytes);
        $f::<MultiLineString3d>($bytes);
        $f::<MultiPolygon3d>($bytes);
    }};
}

fuzz_target!(|data: &[u8]| {
    every_type!(check, data);
    if let Some((&which, body)) = data.split_first() {
        // One type per input for the header-correct path, chosen by the
        // first byte, so each body reaches its decoder at full depth.
        match which % 24 {
            0 => check_body::<Coord2d>(body),
            1 => check_body::<CoordSeq>(body),
            2 => check_body::<Point>(body),
            3 => check_body::<LineString>(body),
            4 => check_body::<LinearRing>(body),
            5 => check_body::<Polygon>(body),
            6 => check_body::<MultiPoint>(body),
            7 => check_body::<MultiLineString>(body),
            8 => check_body::<MultiPolygon>(body),
            9 => check_body::<Geometry>(body),
            10 => check_body::<GeometryCollection>(body),
            11 => check_body::<Envelope>(body),
            12 => check_body::<IntersectionMatrix>(body),
            13 => check_body::<Dimension>(body),
            14 => check_body::<Location>(body),
            15 => check_body::<Coord3d>(body),
            16 => check_body::<CoordSeq3d>(body),
            17 => check_body::<Point3d>(body),
            18 => check_body::<LineString3d>(body),
            19 => check_body::<LinearRing3d>(body),
            20 => check_body::<Polygon3d>(body),
            21 => check_body::<MultiPoint3d>(body),
            22 => check_body::<MultiLineString3d>(body),
            _ => check_body::<MultiPolygon3d>(body),
        }
    }
});
