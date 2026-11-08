// fuzz/fuzz_targets/simplify_topology.rs

//! Fuzz target: topology-preserving simplification.
//!
//! The input is decoded into a polygon or a multipolygon: one to three
//! parts, each an exterior and up to two holes of up to twelve vertices.
//! Coordinates lie on a half-unit grid in `[0, 24]`, where touches, shared
//! vertices, vertices on edges, collinear runs and narrow gaps are common.
//! The last byte picks a tolerance in `[0, 8]`, in steps of an eighth.
//! Inputs the constructors refuse are skipped.
//!
//! The contract checked, for every constructed input:
//!
//! - neither entry point panics;
//! - promise 1: each output ring is a subsequence of its input ring, bit
//!   for bit, in the input's order;
//! - promise 3: for every pair of rings (across all parts), the DE-9IM
//!   matrix of the rings as closed lines is the same before and after, and
//!   each ring is as simple as it was;
//! - promise 5: a valid input is simplified and vouched, and the output has
//!   no validity issue; an invalid one is returned unvouched, never
//!   vouched.

#![no_main]

use geometres::ops::{
    is_simple, simplify_preserve_topology, simplify_preserve_topology_unvouched, validity_issues,
};
use geometres::predicates::relate;
use geometres::types::{
    Coord2d, CoordSeq, LineString, LinearRing, MaybeVouched, MultiPolygon, Polygon,
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

    fn value(&mut self) -> f64 {
        f64::from(self.u8() % 49) * 0.5
    }

    fn ring(&mut self) -> Option<LinearRing> {
        let n = 3 + usize::from(self.u8() % 10);
        let v: Vec<Coord2d> = (0..n)
            .map(|_| {
                let x = self.value();
                let y = self.value();
                Coord2d::new(x, y)
            })
            .collect();
        LinearRing::new(CoordSeq::new(v)).ok()
    }

    fn polygon(&mut self) -> Option<Polygon> {
        let mut ext = self.ring()?;
        if ext.signed_area() < 0.0 {
            ext = reversed(&ext)?;
        }
        let holes = (0..self.u8() % 3)
            .filter_map(|_| {
                let h = self.ring()?;
                if h.signed_area() > 0.0 {
                    reversed(&h)
                } else {
                    Some(h)
                }
            })
            .collect();
        Polygon::new(ext, holes).ok()
    }
}

fn reversed(r: &LinearRing) -> Option<LinearRing> {
    let mut v = r.coords().as_slice().to_vec();
    v.reverse();
    LinearRing::new(CoordSeq::new(v)).ok()
}

fn rings_of(m: &MultiPolygon) -> Vec<Vec<Coord2d>> {
    m.polygons()
        .iter()
        .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes()))
        .map(|r| r.coords().as_slice().to_vec())
        .collect()
}

fn closed(v: &[Coord2d]) -> LineString {
    let mut w = v.to_vec();
    w.push(v[0]);
    LineString::new(CoordSeq::new(w)).expect("a ring closed is a line")
}

fn bits(p: Coord2d) -> (u64, u64) {
    (p.x.to_bits(), p.y.to_bits())
}

/// Promises 1 and 3 between an input and an output multipolygon.
fn check(input: &MultiPolygon, output: &MultiPolygon) {
    let (ri, ro) = (rings_of(input), rings_of(output));
    assert_eq!(ri.len(), ro.len(), "ring count changed");
    for (a, b) in ri.iter().zip(&ro) {
        let mut k = 0;
        for &p in b {
            while k < a.len() && bits(a[k]) != bits(p) {
                k += 1;
            }
            assert!(
                k < a.len(),
                "output vertex {p:?} is not an input vertex in order"
            );
            k += 1;
        }
    }
    let (li, lo): (Vec<LineString>, Vec<LineString>) = (
        ri.iter().map(|r| closed(r)).collect(),
        ro.iter().map(|r| closed(r)).collect(),
    );
    for i in 0..li.len() {
        assert_eq!(
            is_simple(&li[i]),
            is_simple(&lo[i]),
            "ring {i} simplicity changed"
        );
        for j in i + 1..li.len() {
            assert_eq!(
                relate(&li[i], &li[j]).expect("lines relate"),
                relate(&lo[i], &lo[j]).expect("lines relate"),
                "rings {i} and {j} changed relation"
            );
        }
    }
}

fuzz_target!(|data: &[u8]| {
    let mut b = Bytes { data, at: 0 };
    let parts: Vec<Polygon> = (0..1 + b.u8() % 3).filter_map(|_| b.polygon()).collect();
    let tolerance = f64::from(data.last().copied().unwrap_or(0) % 65) / 8.0;
    let Ok(input) = MultiPolygon::new(parts) else {
        return;
    };
    if input.polygons().is_empty() {
        return;
    }
    let valid = validity_issues(&input).is_empty();
    match simplify_preserve_topology_unvouched(&input, tolerance) {
        Ok(MaybeVouched::Vouched(out)) => {
            assert!(valid, "an invalid input was vouched for");
            assert!(
                validity_issues(&out).is_empty(),
                "a vouched output is invalid"
            );
            check(&input, &out);
            assert_eq!(
                simplify_preserve_topology(&input, tolerance).ok(),
                Some(out)
            );
        }
        Ok(MaybeVouched::Unvouched(u)) => {
            assert!(!valid, "a valid input was refused: {:?}", u.reason());
            check(&input, u.inner());
        }
        // A constructor refusing the output: ruled out for a valid input.
        Err(e) => assert!(!valid, "a valid input failed: {e:?}"),
    }
    // A single part goes through the `Polygon` entry point the same way.
    if let [p] = input.polygons() {
        let single = simplify_preserve_topology_unvouched(p, tolerance);
        if valid {
            assert!(
                matches!(single, Ok(MaybeVouched::Vouched(_))),
                "polygon entry point refused"
            );
        }
    }
});
