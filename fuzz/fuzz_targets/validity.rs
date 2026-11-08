// fuzz/fuzz_targets/validity.rs

//! Fuzz target: validity of raw rings against the constructors.
//!
//! The input is decoded into an exterior and up to three holes of up to
//! nine vertices each. Most coordinates lie on a small integer grid, where
//! touches, shared stretches, collinear vertices, repeats and closing
//! vertices are common; a few are drawn from the edges of the `f64` range
//! (signed zero, a subnormal, values whose shoelace products overflow, NaN
//! and infinity).
//!
//! The contract checked:
//!
//! - `validate` never panics;
//! - `validate` reports an issue a constructor refuses if and only if
//!   building the rings and then the polygon fails;
//! - when the polygon builds, the issues `validate` reports are exactly the
//!   issues `validity_issues` reports on the polygon (self-intersection and
//!   a disconnected interior, which construction does not check). So
//!   `validate` reports nothing if and only if the polygon builds and
//!   `validity_issues` is empty.

#![no_main]

use geometres::ops::{ValidityIssue, is_valid, validate, validity_issues};
use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon};
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

    /// One coordinate value: usually a grid value in 0..=16, sometimes an
    /// edge of the `f64` range.
    fn value(&mut self) -> f64 {
        let b = self.u8();
        match b {
            0..=239 => f64::from(b % 17),
            240 => -0.0,
            241 => f64::from_bits(1),
            242 => 1e200,
            243 => -1e200,
            244 => 0.5,
            245 => 16.000_000_000_000_004,
            246 => f64::NAN,
            247 => f64::INFINITY,
            _ => f64::from(b % 17) + 0.25,
        }
    }
}

fn rings_from(data: &[u8]) -> Vec<CoordSeq> {
    let mut b = Bytes { data, at: 0 };
    let n_rings = 1 + usize::from(b.u8() % 4);
    (0..n_rings)
        .map(|_| {
            let n = usize::from(b.u8() % 10);
            CoordSeq::new(
                (0..n)
                    .map(|_| {
                        let x = b.value();
                        let y = b.value();
                        Coord2d::new(x, y)
                    })
                    .collect(),
            )
        })
        .collect()
}

fuzz_target!(|data: &[u8]| {
    let rings = rings_from(data);
    let issues = validate(&rings);
    let refused = issues.iter().any(ValidityIssue::constructor_refuses);

    let built: Option<Polygon> = (|| {
        let mut lrs = rings
            .iter()
            .map(|c| LinearRing::new(c.clone()).ok())
            .collect::<Option<Vec<_>>>()?;
        if lrs.is_empty() {
            return None;
        }
        let exterior = lrs.remove(0);
        Polygon::new(exterior, lrs).ok()
    })();

    match built {
        None => assert!(
            refused,
            "the constructors refuse, validate does not: {issues:?} {rings:?}"
        ),
        Some(p) => {
            assert!(
                !refused,
                "the constructors build, validate refuses: {issues:?} {rings:?}"
            );
            let located: Vec<ValidityIssue> =
                validity_issues(&p).into_iter().map(|l| l.issue).collect();
            assert_eq!(located, issues, "{rings:?}");
            assert_eq!(is_valid(&p), issues.is_empty());
        }
    }
});
