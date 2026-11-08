// src/index/mod.rs

//! Candidate filters over the parts and the edges of a polygonal value.
//!
//! # What an index is, here
//!
//! An index **filters candidates**. [`PartIndex::candidates`] returns every
//! part whose box meets a window; [`EdgeIndex::candidates`] returns every
//! ring edge whose box meets it. Neither decides anything geometric: the
//! caller tests each candidate exactly, with this crate's predicates or
//! its own arrangement. The index's one obligation is therefore **no false
//! negatives**: nothing whose box meets the window is left out. It meets that
//! obligation exactly, because:
//!
//! 1. every entry's box contains its geometry: a part's box is its exterior
//!    ring's envelope, which contains the holes too (`Polygon::new` checks
//!    that each hole lies inside the exterior), and an edge's box is the
//!    minimum and maximum of its two endpoints, both exact selections of
//!    stored coordinates;
//! 2. every node's box contains its children's, by exact union; and
//! 3. the overlap test is the closed-interval comparison of stored values,
//!    the same one [`Envelope::intersects`] makes.
//!
//! The argument in full, and the tests that attack each condition, are in
//! `situs`, the crate that holds the tree. This module never names it in a
//! public signature: the API is typed in [`Envelope`], [`Polygon`] and
//! `usize`, and the one private module `seam` converts (ADR-31). An
//! `Envelope` and a `situs` box have the same validity rules, so the
//! conversion has no failure path; the seam's property tests hold that in
//! place.
//!
//! # Candidates are exactly the entries whose boxes meet the window
//!
//! Not a superset: a leaf entry is reported only after its own box is
//! tested. So a caller may rely on `candidates` returning precisely the
//! parts, or edges, whose boxes intersect the window (closed, so touching
//! counts), in an order that is deterministic but not meaningful.
//!
//! # Ids
//!
//! Entries are identified by `u32` inside the tree. More than `u32::MAX`
//! parts or edges is [`GeometryError::TooManyEntries`], never a truncation.
//! Edges are numbered through a side table of [`EdgeRef`]s, so the public
//! API speaks in `(part, ring, edge)` triples and `usize`.
//!
//! # Not here
//!
//! Prepared predicates build on these, and a part filter is a consumer
//! of [`PartIndex`].

mod seam;

pub(crate) use seam::{BoxIndex, SegmentIndex};

use crate::error::GeometryError;
use crate::types::{Coord2d, Envelope, Polygon};

#[cfg(test)]
mod tests;

/// Converts an entry count to the index's `u32` id space.
///
/// # Errors
///
/// [`GeometryError::TooManyEntries`] if `count` exceeds `u32::MAX`.
pub(crate) fn check_count(count: usize) -> Result<u32, GeometryError> {
    u32::try_from(count).map_err(|_| GeometryError::TooManyEntries { count })
}

// ===========================================================================
// PartIndex
// ===========================================================================

/// A candidate filter over the parts of a polygonal value.
///
/// Built once over a slice of parts (a `MultiPolygon`'s polygons, or any
/// `&[Polygon]`); each part's box is its exterior ring's envelope.
/// [`Self::candidates`] then returns the indices of exactly the parts whose
/// boxes meet a window. Exact: see the module doc.
///
/// ```
/// use geometres::index::PartIndex;
/// use geometres::types::{Coord2d, CoordSeq, Envelope, LinearRing, Polygon};
///
/// let square = |x: f64, y: f64| {
///     let ring = LinearRing::new(CoordSeq::new(vec![
///         Coord2d::new(x, y), Coord2d::new(x + 1.0, y),
///         Coord2d::new(x + 1.0, y + 1.0), Coord2d::new(x, y + 1.0),
///     ]))?;
///     Polygon::new_simple(ring)
/// };
/// let parts = vec![square(0.0, 0.0)?, square(5.0, 5.0)?, square(1.0, 0.0)?];
/// let index = PartIndex::new(&parts)?;
///
/// let window = Envelope::new(Coord2d::new(1.0, 0.5), Coord2d::new(1.5, 0.5))?;
/// let mut hits = Vec::new();
/// index.candidates(&window, |part| hits.push(part));
/// hits.sort_unstable();
/// assert_eq!(hits, [0, 2]); // part 0 touches the window at x = 1; part 2 contains it
/// # Ok::<(), geometres::GeometryError>(())
/// ```
#[derive(Clone, Debug)]
pub struct PartIndex {
    index: seam::BoxIndex,
}

impl PartIndex {
    /// Indexes `parts` by their exterior envelopes.
    ///
    /// `O(n log n)`. The part at position `i` of `parts` is reported as
    /// candidate `i`.
    ///
    /// # Errors
    ///
    /// - [`GeometryError::TooManyEntries`] if `parts.len() > u32::MAX`.
    /// - A `LinearRing::envelope` error, which is unreachable for a valid
    ///   ring and is propagated rather than unwrapped.
    pub fn new(parts: &[Polygon]) -> Result<Self, GeometryError> {
        check_count(parts.len())?;
        let mut entries: Vec<(Envelope, u32)> = Vec::with_capacity(parts.len());
        for (i, p) in parts.iter().enumerate() {
            // `check_count` passed, so `i < u32::MAX` and this cannot fail.
            let id = check_count(i)?;
            entries.push((p.exterior().envelope()?, id));
        }
        Ok(Self {
            index: seam::BoxIndex::from_envelopes(&entries),
        })
    }

    /// Calls `visit` with the index of every part whose box meets `window`.
    ///
    /// Exactly those parts, each once, closed (a part touching the window's
    /// edge is a candidate). The order is deterministic but not meaningful.
    pub fn candidates(&self, window: &Envelope, mut visit: impl FnMut(usize)) {
        self.index.candidates(window, |id| visit(id as usize));
    }

    /// The number of parts indexed.
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// `true` if no part is indexed.
    pub fn is_empty(&self) -> bool {
        self.index.len() == 0
    }

    /// The envelope of every indexed part, or `None` if there are none.
    ///
    /// Every bound is one of the parts' coordinates; nothing is computed.
    pub fn bounds(&self) -> Option<Envelope> {
        self.index.bounds()
    }
}

// ===========================================================================
// EdgeIndex
// ===========================================================================

/// One ring edge of a polygonal value, by position.
///
/// `ring` 0 is the part's exterior; `ring` `k ≥ 1` is `holes()[k − 1]`.
/// `edge` is the position in `LinearRing::edges`, so edge `i` runs from
/// vertex `i` to vertex `i + 1`, and the last edge closes the ring.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeRef {
    /// Position of the part in the indexed slice.
    pub part: usize,
    /// 0 for the exterior ring, `k` for hole `k − 1`.
    pub ring: usize,
    /// Position of the edge in the ring.
    pub edge: usize,
}

/// A candidate filter over every ring edge of a polygonal value.
///
/// Each edge's box is the minimum and maximum of its two endpoints.
/// [`Self::candidates`] returns exactly the edges whose boxes meet a
/// window, as [`EdgeRef`]s. Exact: see the module doc.
#[derive(Clone, Debug)]
pub struct EdgeIndex {
    index: seam::BoxIndex,
    /// The edge each `u32` id names, in entry order: part by part, exterior
    /// first and then holes, each ring in `edges` order.
    refs: Vec<EdgeRef>,
}

impl EdgeIndex {
    /// Indexes every edge of every ring of `parts`.
    ///
    /// `O(E log E)` for `E` edges.
    ///
    /// # Errors
    ///
    /// [`GeometryError::TooManyEntries`] if the edge count exceeds
    /// `u32::MAX`.
    pub fn new(parts: &[Polygon]) -> Result<Self, GeometryError> {
        let mut refs: Vec<EdgeRef> = Vec::new();
        let mut entries: Vec<((Coord2d, Coord2d), u32)> = Vec::new();
        for (part, p) in parts.iter().enumerate() {
            let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
            for (ring, r) in rings.enumerate() {
                for (edge, (a, b)) in r.edges().enumerate() {
                    let id = check_count(refs.len())?;
                    refs.push(EdgeRef { part, ring, edge });
                    entries.push(((a, b), id));
                }
            }
        }
        Ok(Self {
            index: seam::BoxIndex::from_edges(&entries),
            refs,
        })
    }

    /// Calls `visit` with every edge whose box meets `window`.
    ///
    /// Exactly those edges, each once, closed. The order is deterministic
    /// but not meaningful.
    pub fn candidates(&self, window: &Envelope, mut visit: impl FnMut(EdgeRef)) {
        self.index
            .candidates(window, |id| visit(self.refs[id as usize]));
    }

    /// The number of edges indexed.
    pub fn len(&self) -> usize {
        self.refs.len()
    }

    /// `true` if no edge is indexed.
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    /// The envelope of every indexed edge, or `None` if there are none.
    pub fn bounds(&self) -> Option<Envelope> {
        self.index.bounds()
    }
}
