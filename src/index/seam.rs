// src/index/seam.rs

//! The crate's seam with `situs`.
//!
//! Everything the index module needs from `situs` passes through here:
//! the conversion of an [`Envelope`] to a `situs::Rect` and back, the box of
//! an edge, and a thin wrapper around `situs::Index` keyed by `Envelope`.
//! `situs` types appear in no public signature of this crate (ADR-31); the
//! public module is typed in `Envelope`, the geometry types and `usize`. A
//! search of the crate for `situs::` finds it in this file only.
//!
//! # Why the conversions cannot fail
//!
//! `Envelope` and `situs::Rect` have the same invariants: four finite
//! bounds, `min <= max` on each axis, zero extent allowed. So every
//! `Envelope` converts to a `Rect`, and every `Rect` to an `Envelope`,
//! without a refusal path. An edge box is the minimum and maximum of two
//! coordinates of a validated ring, which are finite, so `min <= max` holds
//! by construction. Each conversion returns the value, with the invariant
//! stated in its `expect`: a broken invariant would panic, never substitute
//! a box, because a wrong box could lose a candidate silently. The seam's
//! property tests hold the equivalence in place. They check that
//! `Rect::new` and `Envelope::new` accept and refuse the same inputs, and
//! that the two types' closed overlap and containment tests agree.
//!
//! The conversions copy stored values; no arithmetic.

use crate::types::{Coord2d, Envelope};

/// The `situs::Rect` of an envelope. Cannot fail: see the module doc.
///
/// The `expect` is the invariant's statement, not a guess: both types
/// refuse exactly non-finite and inverted bounds, and the envelope passed
/// its constructor. Were the invariant ever broken, a panic is the right
/// outcome; a substituted box could lose a candidate silently.
fn to_rect(e: &Envelope) -> situs::Rect {
    let (min, max) = (e.min(), e.max());
    situs::Rect::new(min.x, min.y, max.x, max.y)
        .expect("an Envelope satisfies Rect's invariants: finite, min <= max")
}

/// The box of the segment `a → b`: the minimum and maximum of the two
/// endpoints on each axis, exact. Cannot fail for the coordinates of a
/// validated ring, which are finite, so the chosen minimum is at most the
/// chosen maximum.
fn edge_rect(a: Coord2d, b: Coord2d) -> situs::Rect {
    let (min_x, max_x) = if b.x < a.x { (b.x, a.x) } else { (a.x, b.x) };
    let (min_y, max_y) = if b.y < a.y { (b.y, a.y) } else { (a.y, b.y) };
    situs::Rect::new(min_x, min_y, max_x, max_y)
        .expect("an edge of a validated ring has finite endpoints, so its box is a Rect")
}

/// The envelope of a `situs::Rect`. Cannot fail: see the module doc.
fn from_rect(r: &situs::Rect) -> Envelope {
    Envelope::new(
        Coord2d::new(r.min_x(), r.min_y()),
        Coord2d::new(r.max_x(), r.max_y()),
    )
    .expect("a Rect satisfies Envelope's invariants: finite, min <= max")
}

/// A `situs::Index` keyed by [`Envelope`].
///
/// The public index types hold one of these. Its entries are
/// `(Envelope, u32)`; the `u32` is the caller's id, and the public types
/// map it back to a `usize` part index or an `EdgeRef`.
#[derive(Clone, Debug)]
pub(crate) struct BoxIndex {
    inner: situs::Index,
}

impl BoxIndex {
    /// Bulk-loads the boxes (STR packing; see `situs`).
    pub(super) fn from_envelopes(entries: &[(Envelope, u32)]) -> Self {
        let entries: Vec<(situs::Rect, u32)> =
            entries.iter().map(|(e, id)| (to_rect(e), *id)).collect();
        Self {
            inner: situs::Index::bulk_load(entries),
        }
    }

    /// Bulk-loads edge boxes given as endpoint pairs.
    pub(crate) fn from_edges(entries: &[((Coord2d, Coord2d), u32)]) -> Self {
        let entries: Vec<(situs::Rect, u32)> = entries
            .iter()
            .map(|((a, b), id)| (edge_rect(*a, *b), *id))
            .collect();
        Self {
            inner: situs::Index::bulk_load(entries),
        }
    }

    /// Every id whose box meets `window` (closed), through `visit`.
    pub(crate) fn candidates(&self, window: &Envelope, visit: impl FnMut(u32)) {
        self.inner.query(&to_rect(window), visit);
    }

    /// [`Self::candidates`] with an early exit: the search stops at the
    /// first `Break` and returns it. Same candidates, same no-false-negative
    /// obligation; built for `ops::distance`'s `dwithin`, which needs only
    /// one passing pair.
    pub(crate) fn try_candidates<B>(
        &self,
        window: &Envelope,
        visit: impl FnMut(u32) -> core::ops::ControlFlow<B>,
    ) -> core::ops::ControlFlow<B> {
        self.inner.try_query(&to_rect(window), visit)
    }

    /// Every unordered pair of distinct entries whose boxes meet (closed),
    /// through `visit`, each pair once (`situs`' dual-tree self-join). Built
    /// for `ops::coverage`, whose candidate pairs are pairs of edges or of
    /// polygons; the obligation is the same: no false negatives.
    pub(crate) fn self_join(&self, visit: impl FnMut(u32, u32)) {
        self.inner.self_join(visit);
    }

    /// Bulk-loads boxes given as envelopes; [`Self::from_envelopes`] under
    /// the crate-visible name `ops::coverage` uses.
    pub(crate) fn from_boxes(entries: &[(Envelope, u32)]) -> Self {
        Self::from_envelopes(entries)
    }

    pub(super) fn len(&self) -> usize {
        self.inner.len()
    }

    /// The envelope of every entry, or `None` if there are none.
    pub(super) fn bounds(&self) -> Option<Envelope> {
        self.inner.bounds().map(|r| from_rect(&r))
    }
}

// ===========================================================================
// SegmentIndex: the dynamic tree, for a changing set of segments
// ===========================================================================

/// A `situs::DynamicIndex` (the R\*-tree) over segments and points, keyed by
/// a caller's `u32` id.
///
/// Built for topology-preserving simplification (`ops::simplify_topology`),
/// whose set of current segments changes as runs are replaced by chords: a
/// run's segments are removed and the chord inserted. Each entry's box is
/// [`edge_rect`] of its two endpoints (a point is the segment from itself to
/// itself), so the box contains the segment exactly, and the tree's one
/// obligation, no false negatives, carries over: every entry whose box meets
/// a window is reported.
///
/// Removal names the entry by the same endpoints and id it was inserted
/// with. `situs` removes exactly the entry with that box and id, so a caller
/// that removes what it inserted leaves the tree holding exactly the current
/// set; the seam's property test below holds that in place against brute
/// force.
#[derive(Clone, Debug, Default)]
pub(crate) struct SegmentIndex {
    inner: situs::DynamicIndex,
}

impl SegmentIndex {
    /// An empty index.
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Adds the segment `a → b` under `id`.
    pub(crate) fn insert(&mut self, a: Coord2d, b: Coord2d, id: u32) {
        self.inner.insert(edge_rect(a, b), id);
    }

    /// Removes the segment `a → b` inserted under `id`; `true` if it was
    /// there.
    pub(crate) fn remove(&mut self, a: Coord2d, b: Coord2d, id: u32) -> bool {
        self.inner.remove(&edge_rect(a, b), id)
    }

    /// Every id whose box meets the closed box spanned by `lo` and `hi`
    /// (any two corners, in either order), through `visit`.
    pub(crate) fn candidates(&self, lo: Coord2d, hi: Coord2d, visit: impl FnMut(u32)) {
        self.inner.query(&edge_rect(lo, hi), visit);
    }

    /// The number of entries.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.inner.len()
    }
}

// ===========================================================================
// The seam
// ===========================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::error::GeometryError;
    use proptest::prelude::*;

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    /// `Rect::new` and `Envelope::new` accept and refuse the same bounds,
    /// including non-finite and inverted ones.
    fn same_verdict(a: f64, b: f64, x: f64, y: f64) {
        let rect = situs::Rect::new(a, b, x, y);
        let env = Envelope::new(c(a, b), c(x, y));
        match (rect, env) {
            (Ok(r), Ok(e)) => {
                assert_eq!(to_rect(&e), r);
                assert_eq!(from_rect(&r), e);
            }
            (Err(situs::RectError::NonFinite), Err(GeometryError::NonFiniteCoordinate { .. }))
            | (Err(situs::RectError::Inverted), Err(GeometryError::InvertedEnvelope)) => {}
            (r, e) => {
                panic!("Rect::new gave {r:?} but Envelope::new gave {e:?} for ({a}, {b}, {x}, {y})")
            }
        }
    }

    #[test]
    fn fixed_bounds_get_the_same_verdict() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            same_verdict(bad, 0.0, 1.0, 1.0);
            same_verdict(0.0, bad, 1.0, 1.0);
            same_verdict(0.0, 0.0, bad, 1.0);
            same_verdict(0.0, 0.0, 1.0, bad);
            same_verdict(bad, 5.0, 0.0, 0.0);
        }
        same_verdict(1.0, 0.0, 0.0, 1.0);
        same_verdict(0.0, 1.0, 1.0, 0.0);
        same_verdict(2.0, 3.0, 2.0, 3.0);
        same_verdict(-0.0, 0.0, 0.0, -0.0);
    }

    #[test]
    fn edge_rect_is_the_endpoints_box_in_either_order() {
        let r = edge_rect(c(3.0, -1.0), c(-2.0, 4.0));
        assert_eq!(
            (r.min_x(), r.min_y(), r.max_x(), r.max_y()),
            (-2.0, -1.0, 3.0, 4.0)
        );
        assert_eq!(edge_rect(c(-2.0, 4.0), c(3.0, -1.0)), r);
        let p = edge_rect(c(1.0, 1.0), c(1.0, 1.0));
        assert_eq!((p.min_x(), p.max_x()), (1.0, 1.0));
    }

    /// Mostly finite, often degenerate, sometimes not finite.
    fn bound() -> impl Strategy<Value = f64> {
        prop_oneof![
            6 => (-20i32..=20).prop_map(|v| f64::from(v) * 0.5),
            2 => prop::num::f64::ANY,
            1 => prop_oneof![Just(f64::NAN), Just(f64::INFINITY), Just(f64::NEG_INFINITY), Just(-0.0)],
        ]
    }

    /// Finite, on a small grid, so boxes touch and coincide often.
    fn bound_finite() -> impl Strategy<Value = f64> {
        (-12i32..=12).prop_map(|v| f64::from(v) * 0.5)
    }

    fn envelope() -> impl Strategy<Value = Envelope> {
        (bound(), bound(), bound(), bound()).prop_filter_map("finite", |(a, b, x, y)| {
            let (lo_x, hi_x) = if a <= x { (a, x) } else { (x, a) };
            let (lo_y, hi_y) = if b <= y { (b, y) } else { (y, b) };
            Envelope::new(c(lo_x, lo_y), c(hi_x, hi_y)).ok()
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(20_000))]

        /// `Rect::new` and `Envelope::new` accept and refuse the same
        /// bounds.
        #[test]
        fn new_agrees(a in bound(), b in bound(), x in bound(), y in bound()) {
            same_verdict(a, b, x, y);
        }

        /// `intersects` agrees on accepted pairs, and so does
        /// containment of one box in another.
        #[test]
        fn intersects_agrees(e in envelope(), f in envelope()) {
            let (r, s) = (to_rect(&e), to_rect(&f));
            prop_assert_eq!(r.intersects(&s), e.intersects(&f));
            prop_assert_eq!(r.contains_rect(&s), e.contains_envelope(&f));
            prop_assert_eq!(from_rect(&r.union(&s)).contains_envelope(&e), true);
        }

        /// `contains_point` agrees with `contains_coord`,
        /// including for NaN and infinite coordinates.
        #[test]
        fn contains_agrees(e in envelope(), x in bound(), y in bound()) {
            prop_assert_eq!(to_rect(&e).contains_point(x, y), e.contains_coord(c(x, y)));
        }

        /// The dynamic index, after any sequence of inserts and removals,
        /// reports exactly the current entries whose boxes meet a window,
        /// as brute force over the current set does.
        #[test]
        fn segment_index_tracks_inserts_and_removals(
            segs in prop::collection::vec((bound_finite(), bound_finite(), bound_finite(), bound_finite()), 1..60),
            removals in prop::collection::vec(any::<prop::sample::Index>(), 0..60),
            windows in prop::collection::vec((bound_finite(), bound_finite(), bound_finite(), bound_finite()), 1..8),
        ) {
            let mut index = SegmentIndex::new();
            let mut live: Vec<(Coord2d, Coord2d, u32)> = Vec::new();
            for (k, &(a, b, x, y)) in segs.iter().enumerate() {
                let id = u32::try_from(k).unwrap();
                index.insert(c(a, b), c(x, y), id);
                live.push((c(a, b), c(x, y), id));
            }
            for r in &removals {
                if live.is_empty() {
                    break;
                }
                let (a, b, id) = live.swap_remove(r.index(live.len()));
                prop_assert!(index.remove(a, b, id));
                // A second removal of the same entry finds nothing.
                prop_assert!(!index.remove(a, b, id));
            }
            prop_assert_eq!(index.len(), live.len());
            for &(a, b, x, y) in &windows {
                let window = edge_rect(c(a, b), c(x, y));
                let mut got = Vec::new();
                index.candidates(c(a, b), c(x, y), |id| got.push(id));
                got.sort_unstable();
                let mut want: Vec<u32> = live
                    .iter()
                    .filter(|(p, q, _)| edge_rect(*p, *q).intersects(&window))
                    .map(|&(_, _, id)| id)
                    .collect();
                want.sort_unstable();
                prop_assert_eq!(got, want);
            }
        }

        /// The round trip is the identity, bit for bit.
        #[test]
        fn round_trip_is_exact(e in envelope()) {
            let back = from_rect(&to_rect(&e));
            prop_assert_eq!(back.min().x.to_bits(), e.min().x.to_bits());
            prop_assert_eq!(back.min().y.to_bits(), e.min().y.to_bits());
            prop_assert_eq!(back.max().x.to_bits(), e.max().x.to_bits());
            prop_assert_eq!(back.max().y.to_bits(), e.max().y.to_bits());
        }
    }
}
