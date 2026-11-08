// src/ops/pinch.rs

//! The exact *pinch* split: a closed coordinate sequence cut at every vertex it
//! passes through twice (ADR-26, ADR-32).
//!
//! A *pinch* is a ring that visits one of its own vertices twice, such as
//! ESRI's "inverted hole" (an exterior that loops in to draw a hole touching
//! it at a point) or a figure eight. It is not an OGC-simple ring, but it
//! bounds a well-defined region. Cutting it at the repeated vertex into
//! loops that each visit no vertex twice gives that region in OGC form
//! **without moving a point**: every loop is a run of the ring's own edges, and
//! together the loops use each edge once.
//!
//! **The split does not decide** which loop is a hole and which an exterior,
//! whether a loop is a valid ring, or whether the fill rules agree on the
//! region. That remains the caller's job. A caller can settle it by
//! certified signed area and exact placement, and refuse what it cannot
//! decide.
//!
//! # Exactness
//!
//! Vertices are compared by their stored `f64` bits, after `−0.0` is mapped
//! to `0.0` (the two are the same point); no tolerance, no arithmetic. A
//! sequence with a non-finite coordinate is returned whole: it is not a
//! ring, and the constructors refuse it.

use crate::types::{Coord2d, CoordSeq};

/// Splits a closed coordinate sequence at every vertex it passes through
/// twice, exactly.
///
/// The sequence is first normalized: consecutive duplicate vertices are
/// collapsed (a zero-length edge, which removes no point of the ring), and
/// the OGC closing vertex, if present, is dropped. Every returned sequence
/// is a run of the normalized ring's edges that visits no vertex twice, in
/// the ring's own direction, starting at a pinch vertex (except, possibly,
/// the last one, which starts where the input does).
///
/// # The order, and how to put the ring back together
///
/// One pass with a stack: each vertex is pushed; when a vertex is already on
/// the stack, everything from its earlier position up is a closed loop,
/// popped off whole. The ring's closing edge closes the last loop. Loops
/// come out in that pop order, the last being whatever remained on the
/// stack. To reproduce the input's (normalized) vertex sequence, start from
/// the last loop and, taking the others in **reverse** order, insert each
/// one immediately before the first occurrence of its first vertex. The
/// test `rejoining_the_loops_reproduces_the_ring` does exactly that on
/// 10,000 generated pinched rings.
///
/// An unpinched ring comes back as one sequence; an empty one as none.
///
/// ```
/// use geometres::ops::split_pinches;
/// use geometres::types::{Coord2d, CoordSeq};
///
/// // An inverted hole: the unit square loops in through (1, 1) to draw a
/// // diamond, and comes back out through (1, 1).
/// let ring = CoordSeq::new(
///     [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (1.0, 1.0), (1.5, 1.5), (1.0, 1.9), (1.0, 1.1), (1.0, 1.0), (0.0, 2.0)]
///         .map(|(x, y)| Coord2d::new(x, y))
///         .to_vec(),
/// );
/// let loops = split_pinches(&ring);
/// assert_eq!(loops.len(), 2);
/// assert_eq!(loops[0].len(), 4); // the diamond, starting at the pinch
/// assert_eq!(loops[1].len(), 5); // the square, with the pinch vertex once
/// ```
#[must_use]
pub fn split_pinches(ring: &CoordSeq) -> Vec<CoordSeq> {
    let v = normalized(ring.as_slice());
    if v.is_empty() {
        return Vec::new();
    }
    if v.iter().any(|c| !c.x.is_finite() || !c.y.is_finite()) {
        return vec![CoordSeq::new(v)];
    }
    split_at_repeats(&v)
        .into_iter()
        .map(CoordSeq::new)
        .collect()
}

/// `ring` with consecutive duplicates collapsed and the closing vertex
/// dropped, so `[a, b, c, a, a]` and `[a, a, b, c, a]` both become
/// `[a, b, c]`. Duplicates are exact: the same bits, `−0.0` and `0.0` alike.
fn normalized(seq: &[Coord2d]) -> Vec<Coord2d> {
    let mut out: Vec<Coord2d> = Vec::with_capacity(seq.len());
    for &c in seq {
        if out.last().is_none_or(|last| key(last) != key(&c)) {
            out.push(c);
        }
    }
    if out.len() > 1 && key(&out[0]) == key(&out[out.len() - 1]) {
        out.pop();
    }
    out
}

/// The identity of a vertex: its bits, with `−0.0` mapped to `0.0`. Adding
/// `0.0` maps `−0.0` to `0.0` and changes no other value.
fn key(c: &Coord2d) -> (u64, u64) {
    ((c.x + 0.0).to_bits(), (c.y + 0.0).to_bits())
}

/// The stack pass (see [`split_pinches`]).
fn split_at_repeats(v: &[Coord2d]) -> Vec<Vec<Coord2d>> {
    let mut loops = Vec::new();
    let mut stack: Vec<Coord2d> = Vec::with_capacity(v.len());
    let mut at: std::collections::HashMap<(u64, u64), usize> = std::collections::HashMap::new();
    for c in v {
        if let Some(&p) = at.get(&key(c)) {
            let lp: Vec<Coord2d> = stack.drain(p..).collect();
            for u in &lp {
                at.remove(&key(u));
            }
            loops.push(lp);
        }
        at.insert(key(c), stack.len());
        stack.push(*c);
    }
    if !stack.is_empty() {
        loops.push(stack);
    }
    loops
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests {
    use super::*;

    fn seq(pts: &[(f64, f64)]) -> CoordSeq {
        CoordSeq::new(pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect())
    }

    fn pts(s: &CoordSeq) -> Vec<(f64, f64)> {
        s.iter().map(|c| (c.x, c.y)).collect()
    }

    /// The documented reconstruction: start from the last loop, and insert
    /// the others in reverse order before the first occurrence of their
    /// first vertex.
    fn rejoin(loops: &[CoordSeq]) -> Vec<Coord2d> {
        let Some((last, rest)) = loops.split_last() else {
            return Vec::new();
        };
        let mut out: Vec<Coord2d> = last.as_slice().to_vec();
        for lp in rest.iter().rev() {
            let first = lp.as_slice()[0];
            let at = out
                .iter()
                .position(|c| key(c) == key(&first))
                .expect("a loop's first vertex is in the ring");
            out.splice(at..at, lp.as_slice().iter().copied());
        }
        out
    }

    fn visits_no_vertex_twice(s: &CoordSeq) -> bool {
        let mut seen = std::collections::HashSet::new();
        s.iter().all(|c| seen.insert(key(c)))
    }

    // --- the extract tests' cases (arrangement/extract.rs::split_tests) ---

    #[test]
    fn a_ring_without_repeats_is_one_loop() {
        let r = seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let loops = split_pinches(&r);
        assert_eq!(loops.len(), 1);
        assert_eq!(pts(&loops[0]), pts(&r));
    }

    #[test]
    fn one_pinch_splits_into_two_loops_sharing_the_pinch() {
        // Vertex ids 0,1,6,7,8,6,9 as coordinates: the inner run 6,7,8
        // becomes its own loop starting at the pinch; the rest keeps it.
        let r = seq(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (6.0, 0.0),
            (7.0, 1.0),
            (8.0, 2.0),
            (6.0, 0.0),
            (9.0, 5.0),
        ]);
        let loops = split_pinches(&r);
        assert_eq!(pts(&loops[0]), vec![(6.0, 0.0), (7.0, 1.0), (8.0, 2.0)]);
        assert_eq!(
            pts(&loops[1]),
            vec![(0.0, 0.0), (1.0, 0.0), (6.0, 0.0), (9.0, 5.0)]
        );
    }

    #[test]
    fn two_pinches_give_three_loops_in_pop_order() {
        // 0,6,1,6,2,6,3: the same vertex pinched twice.
        let r = seq(&[
            (0.0, 0.0),
            (6.0, 6.0),
            (1.0, 0.0),
            (6.0, 6.0),
            (2.0, 0.0),
            (6.0, 6.0),
            (3.0, 0.0),
        ]);
        let loops = split_pinches(&r);
        assert_eq!(loops.len(), 3);
        assert_eq!(pts(&loops[0]), vec![(6.0, 6.0), (1.0, 0.0)]);
        assert_eq!(pts(&loops[1]), vec![(6.0, 6.0), (2.0, 0.0)]);
        assert_eq!(pts(&loops[2]), vec![(0.0, 0.0), (6.0, 6.0), (3.0, 0.0)]);
        assert_eq!(rejoin(&loops), r.as_slice());
    }

    #[test]
    fn nested_pinches_rejoin_in_the_documented_order() {
        // 0,6,7,8,7,9,6,3: a pinch at 7 inside a pinch at 6.
        let r = seq(&[
            (0.0, 0.0),
            (6.0, 0.0),
            (7.0, 0.0),
            (8.0, 0.0),
            (7.0, 0.0),
            (9.0, 0.0),
            (6.0, 0.0),
            (3.0, 0.0),
        ]);
        let loops = split_pinches(&r);
        assert_eq!(loops.len(), 3);
        assert_eq!(rejoin(&loops), r.as_slice());
    }

    // --- Cases from a consumer's materialization, at the split's own level ---

    #[test]
    fn the_inverted_hole_splits_at_its_pinch() {
        let r = seq(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (3.0, 1.0),
            (2.0, 2.0),
            (1.0, 1.0),
            (2.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
        ]);
        let loops = split_pinches(&r);
        assert_eq!(loops.len(), 2);
        assert_eq!(
            pts(&loops[0]),
            vec![(2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (1.0, 1.0)]
        );
        assert_eq!(
            pts(&loops[1]),
            vec![(0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
        );
    }

    #[test]
    fn a_pinch_at_the_first_vertex_is_split_the_same() {
        // The figure eight started at its pinch vertex, with the closing
        // vertex present.
        let r = seq(&[
            (1.0, 1.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (1.0, 2.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
        ]);
        let loops = split_pinches(&r);
        assert_eq!(loops.len(), 2);
        assert_eq!(
            pts(&loops[0]),
            vec![(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]
        );
        assert_eq!(
            pts(&loops[1]),
            vec![(1.0, 1.0), (0.0, 1.0), (0.0, 0.0), (1.0, 0.0)]
        );
    }

    #[test]
    fn a_pinched_hole_becomes_two_loops() {
        let r = seq(&[
            (2.0, 2.0),
            (2.0, 3.0),
            (3.0, 3.0),
            (3.0, 4.0),
            (4.0, 4.0),
            (4.0, 3.0),
            (3.0, 3.0),
            (3.0, 2.0),
        ]);
        let loops = split_pinches(&r);
        assert_eq!(loops.len(), 2);
        assert!(loops.iter().all(visits_no_vertex_twice));
        assert_eq!(rejoin(&loops), r.as_slice());
    }

    /// The split does not judge the loops: a bowtie loop (which a caller
    /// would refuse afterwards) and a two-vertex spike are returned as cut.
    #[test]
    fn the_split_cuts_without_judging_the_loops() {
        let bowtie = seq(&[
            (0.0, 0.0),
            (5.0, 0.0),
            (7.0, 2.0),
            (7.0, 1.0),
            (5.0, 1.0),
            (5.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
        ]);
        let loops = split_pinches(&bowtie);
        assert_eq!(loops.len(), 2);
        assert_eq!(loops[0].len(), 4);
        let spike = seq(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (3.0, 1.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (0.0, 2.0),
        ]);
        let loops = split_pinches(&spike);
        assert_eq!(loops.len(), 2);
        assert_eq!(pts(&loops[0]), vec![(2.0, 1.0), (3.0, 1.0)]);
    }

    // --- normalization and exactness ---

    #[test]
    fn consecutive_duplicates_and_the_closing_vertex_are_dropped() {
        let a = seq(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0), (0.0, 0.0)]);
        let b = seq(&[(0.0, 0.0), (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]);
        assert_eq!(
            pts(&split_pinches(&a)[0]),
            vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]
        );
        assert_eq!(
            pts(&split_pinches(&b)[0]),
            vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]
        );
    }

    #[test]
    fn negative_zero_is_the_same_vertex_and_an_ulp_is_not() {
        let r = seq(&[(0.0, 0.0), (1.0, 0.0), (-0.0, 0.0), (0.0, 1.0)]);
        assert_eq!(split_pinches(&r).len(), 2, "−0.0 repeats 0.0");
        let ulp = f64::from_bits(0.0f64.to_bits() + 1);
        let r = seq(&[(0.0, 0.0), (1.0, 0.0), (ulp, 0.0), (0.0, 1.0)]);
        assert_eq!(
            split_pinches(&r).len(),
            1,
            "one ulp apart is another vertex"
        );
    }

    #[test]
    fn empty_and_non_finite_input() {
        assert!(split_pinches(&CoordSeq::new(Vec::new())).is_empty());
        let r = seq(&[(0.0, 0.0), (f64::NAN, 0.0), (0.0, 0.0), (1.0, 1.0)]);
        assert_eq!(split_pinches(&r).len(), 1);
    }

    // --- generated: rejoin reproduces the ring; no loop visits a vertex twice ---

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> usize {
            usize::try_from(self.next() % n).unwrap_or(0)
        }
    }

    /// A ring of up to 40 vertices drawn from a pool of 12 distinct points,
    /// so repeats are frequent, nested and at the same vertex several
    /// times; consecutive repeats are allowed in, to exercise the
    /// normalization.
    fn generated(g: &mut Lcg) -> CoordSeq {
        let pool: Vec<Coord2d> = (0..12)
            .map(|i| Coord2d::new(f64::from(i % 4) * 1.5 - 0.0, f64::from(i / 4) * 2.25))
            .collect();
        let n = 3 + g.below(38);
        let v: Vec<Coord2d> = (0..n).map(|_| pool[g.below(12)]).collect();
        CoordSeq::new(v)
    }

    #[test]
    fn rejoining_the_loops_reproduces_the_ring() {
        let mut g = Lcg(0x9E37_79B9_7F4A_7C15);
        let mut pinched = 0usize;
        for _ in 0..10_000 {
            let r = generated(&mut g);
            let loops = split_pinches(&r);
            let norm = normalized(r.as_slice());
            assert_eq!(rejoin(&loops), norm, "{:?}", pts(&r));
            assert!(loops.iter().all(visits_no_vertex_twice), "{:?}", pts(&r));
            let total: usize = loops.iter().map(CoordSeq::len).sum();
            assert_eq!(total, norm.len(), "every vertex in exactly one loop");
            pinched += usize::from(loops.len() > 1);
        }
        assert!(pinched > 5_000, "{pinched} pinched rings");
    }
}
