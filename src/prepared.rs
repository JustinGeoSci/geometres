// src/prepared.rs

//! Prepared polygonal masks for repeated predicate queries (ADR-29).
//!
//! A mask prepared once answers repeated relate queries in less time than
//! the unprepared operation, and every answer is **identical** to the one the
//! unprepared operation returns. **Preparation changes the cost of a query
//! rather than its result.**
//!
//! # The prepared-geometry contract
//!
//! 1. **Identity.** For every input, a prepared operation returns exactly
//!    what the unprepared operation returns: the same `IntersectionMatrix`
//!    cell for cell, the same output bytes and the same refusal with the
//!    same reason and `inputs`.
//!
//! 2. **Simplicity (ADR-32).** Preparation runs
//!    [`exact::first_self_contact`](crate::exact::first_self_contact) on
//!    every ring of every part, once. If every ring is simple, the
//!    preparation is `Simple`, and the fast structures may assume simple
//!    rings. Otherwise it is `Delegating`: every query calls the unprepared
//!    operation, and is correct and not faster. The verdict is visible
//!    ([`PreparedPolygonal::is_delegating`]), so a caller can see why a mask
//!    is slow. The constructors are unchanged (ADR-27): a non-simple
//!    ring is still a `Polygon` if `Polygon::new` accepted it.
//!
//! 3. Exactness of a rounded result. A predicate evaluated on a rounded
//!    overlay result is exact for that result rather than for the operands
//!    that produced it. `relate(A ∪ B, B)` can report `EI = 2` when a
//!    rounded crossing of the union lies 7·10⁻¹⁴ inside `B`. Preparation
//!    does not change this. A caller must not read such a disagreement as a wrong
//!    answer, or test for exact containment of operands in rounded results.
//!
//! 4. **Immutability.** A prepared value is immutable and `Sync`, so a
//!    consumer can share it across a parallel feature loop.
//!
//! # What preparation builds
//!
//! - the parts, their exterior envelopes, and a `PartIndex` over them;
//! - an `EdgeIndex` over every ring edge of every part;
//! - the simplicity verdict.
//!
//! # How a prepared `relate` is computed, and why it is the same matrix
//!
//! `relate(mask, g)` for a `Simple` mask runs in three steps, each with its
//! own identity argument.
//!
//! ## 1. The envelope early-out
//!
//! If the mask's envelope and `g`'s are strictly separated, the matrix is the
//! disjoint one. This is the early-out `polygons_matrix` has for areal `g`; for
//! a line or point `g` the unprepared code has none, so here it is an addition,
//! exact because strictly separated closed boxes share no point: the mask's
//! interior and boundary lie wholly in `g`'s exterior (`IE = 2`, `BE = 1`),
//! `g`'s interior and boundary wholly in the mask's exterior (`EI`, `EB` of
//! `g`'s own dimensions), and nothing else meets.
//!
//! ## 2. The candidate parts
//!
//! The parts whose closed exterior box meets `g`'s closed box are the
//! candidates `C`; the rest, `F`, are far. The matrix is computed over `C`
//! alone, and `F`'s contribution is added by rule. The rule, and its proof: let
//! the mask be `M = C ∪ F`. The parts are interior-disjoint and meet at most at
//! points, and the boundary is topological (every part's rings), so
//!
//! - `I(M) = I(C) ∪ I(F)`, `B(M) = B(C) ∪ B(F)`, and `E(M) = E(C) ∩ E(F)`.
//!
//! A far part's box misses `g`'s box, so `F ∩ B(g) = ∅` as closed boxes,
//! hence `F` lies in `E(g)` with its closure. Then, cell by cell:
//!
//! - `I(M) ∩ I(g) = I(C) ∩ I(g)`, and likewise `IB`, `BI`, `BB`, because
//!   `I(F)` and `B(F)` miss `g`'s closure;
//! - `I(M) ∩ E(g) = (I(C) ∩ E(g)) ∪ I(F)`: the `C` cell, raised to `Area`
//!   when `F` is not empty; and `B(M) ∩ E(g)` likewise raised to `Line`;
//! - `E(M) ∩ I(g) = E(C) ∩ E(F) ∩ I(g) = E(C) ∩ I(g)`, because
//!   `I(g) ⊆ E(F)`; and `EB` likewise;
//! - `EE` is always `Area`.
//!
//! So the matrix over `M` is the matrix over `C` with `IE` raised to `Area`
//! and `BE` to `Line` when any part is far. When no part is a candidate,
//! one far part stands in for `C`: its own matrix against `g` is the
//! disjoint matrix, which is what the rule yields for an empty `C`, and it
//! keeps the computation on the unprepared code path. A subset of a valid
//! multipolygon's parts is valid (the same pairwise condition over fewer
//! pairs), so `C` is a legitimate areal operand.
//!
//! ## 3. The candidate edges
//!
//! Over `C`, the analyses visit only the edges whose closed box meets each
//! question's window, through the `EdgeIndex` (the crate-private
//! `topology::areal::IndexedAreal`): the edges containing a point, crossing a
//! segment, or crossing a point's ray; and the far-edge rule, which classifies
//! a mask edge wholly outside `g`'s box as an `exterior` piece without visiting
//! `g`. The argument that skipping a far edge changes no answer is in
//! `topology::areal`. The analyses themselves are written once, over a trait,
//! so the prepared and unprepared paths share every line of logic; what differs
//! is only which edges each accessor visits.
//!
//! # The predicates
//!
//! `intersects`, `disjoint`, `covers`, `contains` and `within` answer from
//! the full matrix through its predicate methods, after one envelope early-out
//! each: separated boxes decide `intersects` and `disjoint` (as in step 1),
//! and a `g` whose box is not within the mask's box cannot be covered or
//! contained, because some point of `g` is outside the mask's box. Each
//! early-out returns what the matrix would.
//!
//! # Hole-to-face association
//!
//! None is needed. The ray parity is computed per ring and combined per
//! part, exactly as the unindexed `locate` does, so no containment tree is
//! built.

use std::sync::OnceLock;

use crate::error::GeometryError;
use crate::exact::first_self_contact;
use crate::index::EdgeRef;
use crate::index::{EdgeIndex, PartIndex};
use crate::ops::distance as dist;
use crate::predicates;
use crate::topology::areal::{AreaOperand, Areal, IndexedAreal, WindowedAreal};
use crate::topology::de9im::IntersectionMatrix;
use crate::topology::exact_point;
pub use crate::topology::exact_point::ExactPoint;
use crate::topology::label::Dimension;
use crate::topology::label::Location;
use crate::topology::line_analysis::{SegmentEvent, segment_events};
use crate::topology::polygon_analysis::areas_intersect;
use crate::types::{Coord2d, Envelope, Geometry, LineString, MultiPolygon, Polygon, Polygonal};

/// Why a prepared value delegates: the first non-simple ring.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DelegationReason {
    /// The part, by index.
    pub part: usize,
    /// The ring within the part: `0` is the exterior, `k` the hole `k − 1`.
    pub ring: usize,
    /// The first pair of edges of that ring that touch other than at a
    /// shared vertex, by edge index (`exact::first_self_contact`).
    pub edges: (usize, usize),
}

/// The simplicity verdict.
#[derive(Clone, Debug)]
enum Kind {
    /// Every ring is simple; the fast structures are built.
    Simple {
        part_index: PartIndex,
        edge_index: EdgeIndex,
        /// Each part's exterior envelope.
        boxes: Vec<Envelope>,
    },
    /// The mask has at most [`SMALL_MASK_EDGES`] edges: no index is built.
    /// The simplicity verdict is taken on first need and kept (`None`: every ring
    /// is simple; `Some`: the first ring that is not).
    ///
    /// # Why the verdict can wait
    ///
    /// The verdict is a function of the parts alone, so taking it later
    /// gives the same verdict. And the questions a point asks (`locate`,
    /// through [`PreparedPolygonal::operand`]) run over the unindexed
    /// operand whatever the verdict: an unindexed simple mask and a
    /// delegating mask answer them by the same code. Only the questions
    /// that do depend on it (`relate`, the early outs, `edges_near`,
    /// `is_delegating`) take it. So a parcel prepared for point queries
    /// alone never pays the simplicity check, which was most of its
    /// preparation.
    Small(OnceLock<Option<DelegationReason>>),
    /// Some ring is not simple; every query runs the unprepared operation.
    Delegating(DelegationReason),
}

/// Ring edges up to which a simple mask is kept unindexed. A performance
/// choice only: the indexed and unindexed operands give the same answers,
/// and below this size building the index costs more than it saves.
const SMALL_MASK_EDGES: usize = 64;

/// A polygonal mask prepared for repeated `relate`.
///
/// Built once from a `Polygon` or a valid `MultiPolygon`; immutable; `Sync`.
/// Every method returns what the unprepared `relate(mask, other)` returns.
///
/// ```
/// use geometres::prepared::PreparedPolygonal;
/// use geometres::predicates::covers;
/// use geometres::topology::de9im::IntersectionMatrix;
/// use geometres::types::{Coord2d, CoordSeq, Geometry, LinearRing, Point, Polygon};
///
/// let ring = LinearRing::new(CoordSeq::new(vec![
///     Coord2d::new(0.0, 0.0), Coord2d::new(4.0, 0.0), Coord2d::new(4.0, 4.0), Coord2d::new(0.0, 4.0),
/// ]))?;
/// let mask = Polygon::new_simple(ring)?;
/// let prepared = PreparedPolygonal::new(&mask)?;
/// assert!(!prepared.is_delegating());
///
/// let p = Geometry::Point(Point::new(Coord2d::new(1.0, 1.0))?);
/// let unprepared = IntersectionMatrix::compute_geometry(&Geometry::Polygon(mask.clone()), &p)?;
/// assert_eq!(prepared.relate(&p)?, unprepared);
/// assert_eq!(prepared.covers(&p)?, unprepared.covers());
/// assert_eq!(prepared.covers(&p)?, covers(&mask, &p)?);
/// # Ok::<(), geometres::GeometryError>(())
/// ```
#[derive(Clone, Debug)]
pub struct PreparedPolygonal {
    kind: Kind,
    /// The parts, as given.
    mask: MultiPolygon,
    /// The envelope of every part's exterior, `None` for no parts.
    envelope: Option<Envelope>,
}

impl PreparedPolygonal {
    /// Prepares `geom`: the simplicity verdict, the part boxes and index, and the
    /// edge index.
    ///
    /// # Errors
    ///
    /// [`GeometryError::TooManyEntries`] if the mask has more than
    /// `u32::MAX` parts or edges; the only failure.
    pub fn new(geom: &impl Polygonal) -> Result<Self, GeometryError> {
        let parts = geom.parts();
        let mask = MultiPolygon::from_retained_parts(geom, |_, _| true);
        let envelope = Envelope::from_points(
            parts
                .iter()
                .flat_map(|p| p.exterior().coords().as_slice().iter().copied()),
        )?;
        let edges: usize = parts
            .iter()
            .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes().iter()))
            .map(crate::types::LinearRing::num_vertices)
            .sum();
        let kind = if edges > 0 && edges <= SMALL_MASK_EDGES {
            Kind::Small(OnceLock::new())
        } else if let Some(reason) = first_non_simple_ring(parts) {
            Kind::Delegating(reason)
        } else {
            let boxes = parts
                .iter()
                .map(|p| p.exterior().envelope())
                .collect::<Result<Vec<_>, _>>()?;
            Kind::Simple {
                part_index: PartIndex::new(parts)?,
                edge_index: EdgeIndex::new(parts)?,
                boxes,
            }
        };
        Ok(Self {
            kind,
            mask,
            envelope,
        })
    }

    /// Whether every query runs the unprepared operation, because some
    /// ring is not simple.
    pub fn is_delegating(&self) -> bool {
        self.delegation_reason().is_some()
    }

    /// The first non-simple ring, when delegating.
    pub fn delegation_reason(&self) -> Option<DelegationReason> {
        match &self.kind {
            Kind::Delegating(r) => Some(*r),
            Kind::Simple { .. } => None,
            Kind::Small(verdict) => {
                *verdict.get_or_init(|| first_non_simple_ring(self.mask.polygons()))
            }
        }
    }

    /// The mask's parts.
    pub fn parts(&self) -> &[Polygon] {
        self.mask.polygons()
    }

    /// The envelope of every part, or `None` for an empty mask.
    pub fn envelope(&self) -> Option<Envelope> {
        self.envelope
    }

    /// The DE-9IM matrix of the mask (A) against `other` (B): exactly
    /// `IntersectionMatrix::compute_geometry(mask, other)`.
    ///
    /// # Errors
    ///
    /// As `compute_geometry`: a collection operand whose polygonal leaves
    /// overlap is refused (`CollectionElementsOverlap`);
    /// no other pair fails. Prepared answers over collections are identical
    /// to unprepared `relate`, as for every other operand.
    pub fn relate(&self, other: &Geometry) -> Result<IntersectionMatrix, GeometryError> {
        let (part_index, edge_index, boxes) = match &self.kind {
            Kind::Simple {
                part_index,
                edge_index,
                boxes,
            } => (part_index, edge_index, boxes),
            Kind::Small(_) if self.is_delegating() => return self.unprepared(other),
            // The unindexed operand, by the same code `relate` runs.
            Kind::Small(_) => {
                return IntersectionMatrix::compute_area_against(&Areal::multi(&self.mask), other);
            }
            Kind::Delegating(_) => return self.unprepared(other),
        };
        let (Some(mask_env), Some(other_env)) = (self.envelope, envelope_of(other)?) else {
            // An empty mask, or an empty `other`: a point set, with no
            // early-out and no candidates to choose.
            return self.unprepared(other);
        };
        if separated(&mask_env, &other_env) {
            crate::count!(prepared_early_outs);
            return Ok(disjoint_matrix(other));
        }

        // The candidate parts.
        let parts = self.mask.polygons();
        let mut active = vec![false; parts.len()];
        part_index.candidates(&other_env, |i| active[i] = true);
        let far = active.iter().any(|a| !a);
        if !active.iter().any(|a| *a) {
            // No candidate: one far part stands in (see the module doc).
            active[0] = true;
        }
        let env = active
            .iter()
            .zip(boxes)
            .filter(|(a, _)| **a)
            .map(|(_, b)| (b.min().x, b.min().y, b.max().x, b.max().y))
            .reduce(|e, f| (e.0.min(f.0), e.1.min(f.1), e.2.max(f.2), e.3.max(f.3)));

        // The candidate edges.
        // The candidate edges near `other`, taken from the index once, so
        // every question of the analysis is answered from them (and the
        // `+x` strip for `locate`) rather than by a query of its own.
        let operand = WindowedAreal::of_indexed(
            IndexedAreal::new(parts, edge_index, Some(&active), env),
            (other_env.min(), other_env.max()),
        );
        let mut m = IntersectionMatrix::compute_area_against(&operand, other)?;
        if far {
            m.interior_exterior = Dimension::Area;
            m.boundary_exterior = Dimension::Line;
        }
        Ok(m)
    }

    /// [`IntersectionMatrix::intersects`] of [`relate`](Self::relate).
    ///
    /// # Errors
    ///
    /// As [`relate`](Self::relate).
    pub fn intersects(&self, other: &Geometry) -> Result<bool, GeometryError> {
        Ok(!self.disjoint(other)?)
    }

    /// [`IntersectionMatrix::disjoint`] of [`relate`](Self::relate).
    ///
    /// # Errors
    ///
    /// As [`relate`](Self::relate).
    pub fn disjoint(&self, other: &Geometry) -> Result<bool, GeometryError> {
        if let Geometry::Point(p) = other {
            // A point meets the mask exactly when it is not in its exterior.
            return Ok(self.locate(p.coord()) == Location::Exterior);
        }
        if self.boxes_separated(other)? {
            return Ok(true);
        }
        // Two areas: the early-exit search for one common point
        // (`polygon_analysis::areas_intersect`), over the mask's operand.
        let polygonal = match other {
            Geometry::Polygon(p) => Some(Areal::single(p)),
            Geometry::MultiPolygon(m) if !m.is_empty() => Some(Areal::multi(m)),
            _ => None,
        };
        if let (Some(b), false, false) = (
            polygonal,
            self.is_delegating(),
            self.mask.polygons().is_empty(),
        ) {
            let meet = match self.operand() {
                Operand::Indexed(a) => areas_intersect(&a, &b),
                Operand::Plain(a) => areas_intersect(&a, &b),
            };
            return Ok(!meet);
        }
        Ok(self.relate(other)?.disjoint())
    }

    /// [`IntersectionMatrix::covers`] of [`relate`](Self::relate).
    ///
    /// # Errors
    ///
    /// As [`relate`](Self::relate).
    pub fn covers(&self, other: &Geometry) -> Result<bool, GeometryError> {
        if let Geometry::Point(p) = other {
            // A point is covered exactly when it is not in the exterior.
            return Ok(self.locate(p.coord()) != Location::Exterior);
        }
        if self.box_escapes(other)? || self.vertex_escapes(other) {
            return Ok(false);
        }
        Ok(self.relate(other)?.covers())
    }

    /// [`IntersectionMatrix::contains`] of [`relate`](Self::relate).
    ///
    /// # Errors
    ///
    /// As [`relate`](Self::relate).
    pub fn contains(&self, other: &Geometry) -> Result<bool, GeometryError> {
        if let Geometry::Point(p) = other {
            // A point's interior is the point: contained exactly when it is
            // in the mask's interior.
            return Ok(self.locate(p.coord()) == Location::Interior);
        }
        if self.box_escapes(other)? || self.vertex_escapes(other) {
            return Ok(false);
        }
        Ok(self.relate(other)?.contains())
    }

    /// [`IntersectionMatrix::within`] of [`relate`](Self::relate): the mask within
    /// `other`.
    ///
    /// # Errors
    ///
    /// As [`relate`](Self::relate).
    pub fn within(&self, other: &Geometry) -> Result<bool, GeometryError> {
        Ok(self.relate(other)?.within())
    }

    // =======================================================================
    // The exact questions a prepared overlay asks
    // =======================================================================

    /// Exact location of a stored point in the mask: on any part's ring is
    /// `Boundary`; inside some part is `Interior`; otherwise `Exterior`.
    /// Exactly `exact::locate_in_areal` of the mask; through the edge index
    /// for a `Simple` mask.
    pub fn locate(&self, p: Coord2d) -> Location {
        match self.operand() {
            Operand::Indexed(a) => a.locate(p),
            Operand::Plain(a) => a.locate(p),
        }
    }

    /// The location of the points immediately to the **left** of the ray
    /// from `at` toward `toward`, just after `at`: `Interior` or
    /// `Exterior`, never `Boundary`. `at` may be a constructed rational
    /// point (a crossing); `toward` is a stored coordinate that fixes the
    /// direction exactly. This is the winding side of the face to the left
    /// of an arrangement half-edge leaving `at`; the argument is in
    /// `topology::exact_point`.
    pub fn location_left_of(&self, at: ExactPoint<'_>, toward: Coord2d) -> Location {
        match self.operand() {
            Operand::Indexed(a) => exact_point::location_left_of(&at, toward, &a),
            Operand::Plain(a) => exact_point::location_left_of(&at, toward, &a),
        }
    }

    /// Every event on the segment `p → q` against the mask, exactly
    /// `exact::segment_events` of the mask; through the edge index for a
    /// `Simple` mask.
    pub fn segment_events(&self, p: Coord2d, q: Coord2d) -> Vec<SegmentEvent> {
        match self.operand() {
            Operand::Indexed(a) => segment_events(p, q, &a),
            Operand::Plain(a) => segment_events(p, q, &a),
        }
    }

    /// Every ring edge whose closed box meets `window`, as `(part, ring,
    /// edge)` references in ascending order, each once. These are the edges
    /// a windowed arrangement admits. For a `Delegating` mask every edge is
    /// reported, whatever the window.
    pub fn edges_near(&self, window: &Envelope, mut visit: impl FnMut(EdgeRef)) {
        match &self.kind {
            Kind::Simple { edge_index, .. } => {
                let mut refs: Vec<EdgeRef> = Vec::new();
                edge_index.candidates(window, |r| refs.push(r));
                refs.sort_unstable();
                for r in refs {
                    visit(r);
                }
            }
            // The edges whose box meets the window, found by a scan, in
            // the same ascending order the index's sorted candidates have.
            Kind::Small(_) if !self.is_delegating() => {
                for (part, p) in self.mask.polygons().iter().enumerate() {
                    let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
                    for (ring, r) in rings.enumerate() {
                        for (edge, (a, b)) in r.edges().enumerate() {
                            if a.x.min(b.x) <= window.max().x
                                && a.x.max(b.x) >= window.min().x
                                && a.y.min(b.y) <= window.max().y
                                && a.y.max(b.y) >= window.min().y
                            {
                                visit(EdgeRef { part, ring, edge });
                            }
                        }
                    }
                }
            }
            Kind::Small(_) | Kind::Delegating(_) => {
                for (part, p) in self.mask.polygons().iter().enumerate() {
                    let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
                    for (ring, r) in rings.enumerate() {
                        for edge in 0..r.num_vertices() {
                            visit(EdgeRef { part, ring, edge });
                        }
                    }
                }
            }
        }
    }

    /// The stored endpoints of an edge, in ring direction.
    ///
    /// # Panics
    ///
    /// If `r` does not name an edge of this mask; references come from
    /// [`edges_near`](Self::edges_near), which only yields valid ones.
    #[must_use]
    pub fn edge(&self, r: EdgeRef) -> (Coord2d, Coord2d) {
        let p = &self.mask.polygons()[r.part];
        let ring = if r.ring == 0 {
            p.exterior()
        } else {
            &p.holes()[r.ring - 1]
        };
        let c = ring.coords().as_slice();
        (c[r.edge], c[(r.edge + 1) % c.len()])
    }

    /// [`ops::dwithin`](crate::ops::dwithin) of this mask and `other`,
    /// decided **exactly** and with the same answer: the prepared form reuses the
    /// mask's edge index as the candidate filter, and every decision (the
    /// intersection test, each element-pair comparison) is the same exact
    /// predicate the free function makes. For a `Delegating` mask every edge
    /// is a candidate, which is correct and not faster.
    ///
    /// # Errors
    ///
    /// As [`ops::dwithin`](crate::ops::dwithin): an invalid `d`, or a
    /// collection `other` whose polygonal leaves overlap.
    pub fn dwithin(&self, other: &Geometry, d: f64) -> Result<bool, GeometryError> {
        dist::check_distance(d)?;
        if self.intersects(other)? {
            return Ok(true);
        }
        let probe = dist::elements(predicates::Borrowed::from_geometry(other));
        if probe.is_empty() || self.mask.polygons().is_empty() {
            return Ok(false);
        }
        for &e in &probe {
            let window = dist::grow_outward(e.corners(), d);
            let mut hit = false;
            self.edges_near(&window, |r| {
                if !hit {
                    let (a, b) = self.edge(r);
                    hit = dist::pair_within(e, dist::Element::Segment(a, b), d);
                }
            });
            if hit {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// [`ops::distance`](crate::ops::distance()) of this mask and `other`,
    /// **bit for bit**: exactly `0.0` if and only if they intersect, and
    /// otherwise the same `f64` the free function returns, because pruning
    /// through the edge index never discards a pair whose computed value
    /// could be the minimum (the argument is in
    /// [`ops::distance`](mod@crate::ops::distance)'s module documentation).
    ///
    /// # Errors
    ///
    /// As [`ops::distance`](crate::ops::distance()):
    /// [`GeometryError::EmptyOperand`] when `other` has no points or the
    /// mask no parts, and a collection refusal as for
    /// [`Self::dwithin`].
    pub fn distance(&self, other: &Geometry) -> Result<f64, GeometryError> {
        let probe = dist::elements(predicates::Borrowed::from_geometry(other));
        let touching = self.intersects(other)?;
        if probe.is_empty() || self.mask.polygons().is_empty() {
            return Err(GeometryError::EmptyOperand);
        }
        if touching {
            return Ok(0.0);
        }
        let mask_elements = dist::elements(predicates::Borrowed::MultiPolygon(&self.mask));
        let m = dist::coordinate_magnitude(&probe).max(dist::coordinate_magnitude(&mask_elements));
        let delta = dist::DIST_ERR_C * dist::U * m;
        let mut best = dist::pair_distance(probe[0], mask_elements[0]);
        for &e in &probe {
            let r = (best + delta).next_up();
            if !r.is_finite() {
                for &other_e in &mask_elements {
                    best = best.min(dist::pair_distance(e, other_e));
                }
                continue;
            }
            let window = dist::grow_outward(e.corners(), r);
            self.edges_near(&window, |er| {
                let (a, b) = self.edge(er);
                best = best.min(dist::pair_distance(e, dist::Element::Segment(a, b)));
            });
        }
        // The free function's zero guard, so the two agree bit for bit.
        Ok(dist::certified_positive(best))
    }

    /// The operand the exact questions run over: indexed when `Simple`.
    fn operand(&self) -> Operand<'_> {
        match &self.kind {
            Kind::Simple { edge_index, .. } => Operand::Indexed(IndexedAreal::new(
                self.mask.polygons(),
                edge_index,
                None,
                self.envelope
                    .map(|e| (e.min().x, e.min().y, e.max().x, e.max().y)),
            )),
            Kind::Small(_) | Kind::Delegating(_) => Operand::Plain(Areal::multi(&self.mask)),
        }
    }

    /// Some vertex of `other` lies in the mask's exterior: then some point
    /// of `other` is outside the mask, so the mask neither covers nor
    /// contains it. Exact (`locate` on a stored coordinate); stops at the
    /// first such vertex. `false` for a delegating mask, which answers
    /// through the unprepared matrix alone.
    fn vertex_escapes(&self, other: &Geometry) -> bool {
        if self.is_delegating() || self.mask.polygons().is_empty() {
            return false;
        }
        let mut out = false;
        let mut check = |c: Coord2d| {
            if !out {
                out = self.locate(c) == Location::Exterior;
            }
        };
        match other {
            Geometry::Point(p) => check(p.coord()),
            Geometry::MultiPoint(mp) => mp.points().iter().for_each(|p| check(p.coord())),
            Geometry::LineString(l) => l.coords().as_slice().iter().copied().for_each(check),
            Geometry::MultiLineString(ml) => ml
                .line_strings()
                .iter()
                .flat_map(|l| l.coords().as_slice().iter().copied())
                .for_each(check),
            Geometry::Polygon(p) => p
                .exterior()
                .coords()
                .as_slice()
                .iter()
                .copied()
                .for_each(check),
            Geometry::MultiPolygon(mp) => mp
                .polygons()
                .iter()
                .flat_map(|p| p.exterior().coords().as_slice().iter().copied())
                .for_each(check),
            // A collection answers through the matrix.
            Geometry::GeometryCollection(_) => {}
        }
        out
    }

    fn unprepared(&self, other: &Geometry) -> Result<IntersectionMatrix, GeometryError> {
        crate::count!(prepared_delegations);
        IntersectionMatrix::compute_geometry(&Geometry::MultiPolygon(self.mask.clone()), other)
    }

    /// Strictly separated boxes: no shared point, so disjoint.
    fn boxes_separated(&self, other: &Geometry) -> Result<bool, GeometryError> {
        if self.is_delegating() {
            return Ok(false);
        }
        Ok(match (self.envelope, envelope_of(other)?) {
            (Some(a), Some(b)) => separated(&a, &b),
            _ => false,
        })
    }

    /// `other`'s box not within the mask's: some point of `other` is
    /// outside the mask's box, so outside the mask.
    fn box_escapes(&self, other: &Geometry) -> Result<bool, GeometryError> {
        if self.is_delegating() {
            return Ok(false);
        }
        Ok(match (self.envelope, envelope_of(other)?) {
            (Some(a), Some(b)) => !a.contains_envelope(&b),
            _ => false,
        })
    }
}

/// The mask as an areal operand for the exact questions.
enum Operand<'a> {
    Indexed(IndexedAreal<'a>),
    Plain(Areal<'a>),
}

/// The first ring, in part and ring order, that is not simple.
fn first_non_simple_ring(parts: &[Polygon]) -> Option<DelegationReason> {
    for (part, p) in parts.iter().enumerate() {
        let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
        for (ring, r) in rings.enumerate() {
            if let Some(edges) = first_self_contact(r) {
                return Some(DelegationReason { part, ring, edges });
            }
        }
    }
    None
}

/// The closed box of every coordinate of `g`, or `None` for an empty
/// geometry, from the cached boxes of its lines and rings (holes lie inside
/// exteriors). Its bounds are only compared, so the sign of a zero bound,
/// the one way it can differ from a pass over the coordinates
/// ([`bbox_fold`](crate::types::bbox_fold)), does not matter.
fn envelope_of(g: &Geometry) -> Result<Option<Envelope>, GeometryError> {
    let bbox = match g {
        Geometry::Point(p) => {
            let c = p.coord();
            Some((c.x, c.y, c.x, c.y))
        }
        Geometry::MultiPoint(mp) => {
            return Envelope::from_points(mp.points().iter().map(|p| p.coord()));
        }
        Geometry::LineString(l) => Some(l.bbox()),
        Geometry::MultiLineString(ml) => {
            crate::types::bbox_fold(ml.line_strings().iter().map(LineString::bbox))
        }
        Geometry::Polygon(p) => Some(p.exterior().bbox()),
        Geometry::MultiPolygon(mp) => {
            crate::types::bbox_fold(mp.polygons().iter().map(|p| p.exterior().bbox()))
        }
        // A collection's envelope covers every leaf's bounding coordinates
        // (holes lie inside exteriors), which is what the early-outs need.
        Geometry::GeometryCollection(gc) => return Ok(gc.envelope()),
    };
    bbox.map(|(x0, y0, x1, y1)| Envelope::new(Coord2d::new(x0, y0), Coord2d::new(x1, y1)))
        .transpose()
}

/// Strict separation of two closed boxes: no shared point.
fn separated(a: &Envelope, b: &Envelope) -> bool {
    a.max().x < b.min().x || b.max().x < a.min().x || a.max().y < b.min().y || b.max().y < a.min().y
}

/// The matrix of a non-empty areal A disjoint from a non-empty B of any
/// type: `IE = 2`, `BE = 1`, `EI` and `EB` of B's own dimensions, `EE = 2`.
fn disjoint_matrix(other: &Geometry) -> IntersectionMatrix {
    let (ei, eb) = match predicates::dimension(other) {
        Dimension::Area => (Dimension::Area, Dimension::Line),
        Dimension::Line => (Dimension::Line, boundary_dimension_of_lines(other)),
        _ => (Dimension::Point, Dimension::Empty),
    };
    IntersectionMatrix {
        interior_interior: Dimension::Empty,
        interior_boundary: Dimension::Empty,
        interior_exterior: Dimension::Area,
        boundary_interior: Dimension::Empty,
        boundary_boundary: Dimension::Empty,
        boundary_exterior: Dimension::Line,
        exterior_interior: ei,
        exterior_boundary: eb,
        exterior_exterior: Dimension::Area,
    }
}

/// The dimension of a linear geometry's mod-2 boundary: `Point` if some
/// endpoint is the end of an odd number of components, else `Empty` (a
/// closed line, or a set of lines whose ends all pair up).
///
/// A dimension-`Line` collection has no areal leaves, so its boundary is
/// its line leaves' joint mod-2 boundary, the same rule over the leaves
/// (point leaves have no boundary and toggle nothing).
fn boundary_dimension_of_lines(g: &Geometry) -> Dimension {
    fn gather<'g>(g: &'g Geometry, lines: &mut Vec<&'g crate::types::LineString>) {
        match g {
            Geometry::LineString(l) => lines.push(l),
            Geometry::MultiLineString(ml) => lines.extend(ml.line_strings()),
            _ => {}
        }
    }
    let mut lines: Vec<&crate::types::LineString> = Vec::new();
    match g {
        Geometry::GeometryCollection(gc) => {
            for leaf in gc.leaves() {
                gather(leaf, &mut lines);
            }
        }
        _ => gather(g, &mut lines),
    }
    let mut ends: Vec<crate::types::Coord2d> = Vec::new();
    for l in lines {
        let c = l.coords().as_slice();
        ends.push(c[0]);
        ends.push(c[c.len() - 1]);
    }
    let odd = ends
        .iter()
        .any(|e| ends.iter().filter(|f| *f == e).count() % 2 == 1);
    if odd {
        Dimension::Point
    } else {
        Dimension::Empty
    }
}

mod line;

#[cfg(test)]
mod tests;
