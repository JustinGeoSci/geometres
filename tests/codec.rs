// tests/codec.rs

//! The binary codec: exact round trips, determinism, and decoding of
//! hostile input.
//!
//! Every generated geometry decodes to a value equal to the original and
//! re-encodes to the same bytes. Every proper prefix of an encoding is
//! refused, and corrupted encodings decode or are refused but never panic.
//! Whatever decodes has passed the constructors' checks and round-trips
//! itself.

#![allow(clippy::unwrap_used, clippy::cast_possible_truncation)]

mod common;

use common::{Rng, c, geometries};
use geometres::codec::{decode, encode};
use geometres::topology::de9im::IntersectionMatrix;
use geometres::types::{Envelope, Geometry};

#[test]
fn every_geometry_round_trips_exactly() {
    for g in geometries(0x5EED_0501, 280) {
        let bytes = encode(&g);
        let back: Geometry = decode(&bytes).unwrap();
        assert_eq!(back, g);
        assert_eq!(encode(&back), bytes, "re-encoding is byte-identical");
    }
}

#[test]
fn coordinates_round_trip_bit_for_bit() {
    // Signed zeros, subnormals, extremes and values with no short decimal
    // form keep every bit.
    let values = [
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        -f64::from_bits(3),
        1.0 / 3.0,
        f64::MAX,
        f64::MIN,
        1e-300,
        123_456_789.123_456_79,
    ];
    for &x in &values {
        for &y in &values {
            let p = c(x, y);
            let back: geometres::types::Coord2d = decode(&encode(&p)).unwrap();
            assert_eq!(
                (back.x.to_bits(), back.y.to_bits()),
                (x.to_bits(), y.to_bits())
            );
        }
    }
}

#[test]
fn envelopes_and_matrices_round_trip() {
    let e = Envelope::new(c(-1.5, 2.0), c(3.0, 2.0)).unwrap();
    assert_eq!(decode::<Envelope>(&encode(&e)).unwrap(), e);
    for pattern in ["FF2FF1212", "0F1FF0102", "212101212", "FFFFFFFFF"] {
        let m = IntersectionMatrix::from_pattern(pattern).unwrap();
        assert_eq!(decode::<IntersectionMatrix>(&encode(&m)).unwrap(), m);
    }
}

#[test]
fn every_proper_prefix_is_refused() {
    for g in geometries(0x5EED_0502, 70) {
        let bytes = encode(&g);
        for len in 0..bytes.len() {
            assert!(
                decode::<Geometry>(&bytes[..len]).is_err(),
                "prefix of {len} bytes decoded: {g:?}"
            );
        }
    }
}

#[test]
fn trailing_bytes_are_refused() {
    for g in geometries(0x5EED_0503, 35) {
        let mut bytes = encode(&g);
        bytes.push(0);
        assert!(decode::<Geometry>(&bytes).is_err(), "{g:?}");
    }
}

#[test]
fn corrupted_input_never_panics_and_what_decodes_round_trips() {
    let mut rng = Rng::new(0x5EED_0504);
    let mut decoded = 0;
    for g in geometries(0x5EED_0505, 140) {
        let bytes = encode(&g);
        for _ in 0..50 {
            let mut b = bytes.clone();
            for _ in 0..=rng.below(3) {
                let i = rng.below(b.len() as u64) as usize;
                b[i] ^= 1 << rng.below(8);
            }
            if let Ok(h) = decode::<Geometry>(&b) {
                // A flipped coordinate bit can give another value that the
                // constructors accept; it must round-trip like any other.
                assert_eq!(decode::<Geometry>(&encode(&h)).unwrap(), h);
                decoded += 1;
            }
        }
    }
    // Some corruptions land in coordinate bits and still decode.
    assert!(decoded > 0);
}

#[test]
fn random_bytes_never_panic() {
    let mut rng = Rng::new(0x5EED_0506);
    for _ in 0..20_000 {
        let n = rng.below(64) as usize;
        let b: Vec<u8> = (0..n).map(|_| rng.next() as u8).collect();
        let _ = decode::<Geometry>(&b);
    }
}
