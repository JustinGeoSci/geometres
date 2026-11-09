// src/types/result_wrappers.rs

//! Wrappers for results that are not a plain geometry: [`Unvouched<T>`] for a
//! geometry that is **invalid**, and [`MaybeVouched<T, U>`] for an operation
//! that returns either a vouched value or an unvouched one.
//!
//! # Why a wrapper, and why only for invalid results
//!
//! A caller holding an `Unvouched` value has a broken geometry, and must not
//! mistake it for a valid one. A result that is valid but has components of
//! several dimensions (an intersection that is an area, a line and a point)
//! needs no wrapper. It is a `GeometryCollection`, a valid geometry like any
//! other, and the operations that can produce one return a `Geometry`.
//!
//! # The two invariants
//!
//! Both are properties of the API. They are stated here because they are the
//! reason the wrapper exists rather than incidental facts about its definition.
//!
//! ## 1. The inner value requires an explicit unwrap
//!
//! `Unvouched` implements no `Deref`, `AsRef<T>`, `Borrow<T>` or `Into<T>`, and
//! has no public field. The only routes in are [`Unvouched::into_inner`] and
//! [`Unvouched::inner`]. Reading the value is therefore always a visible
//! decision at the call site.
//!
//! ## 2. It does not satisfy a signature expecting a bare geometry
//!
//! There is no coercion to `T` and no blanket conversion, so
//!
//! ```compile_fail
//! use geometres::types::{Polygon, Unvouched, UnvouchedReason};
//! fn takes_polygon(_: &Polygon) {}
//!
//! let u: Unvouched<Polygon> = // ...
//! #     unreachable!();
//! takes_polygon(&u); // error[E0308]: mismatched types
//! ```
//!
//! does not compile, so an invalid result cannot reenter an exact computation
//! without notice. The entry points of the arrangement take `&Polygon`, and a
//! wrapped value reaches one only through an explicit `into_inner()` call,
//! which is easy to search for and review.
//!
//! # The field is private
//!
//! ```compile_fail
//! use geometres::types::{Polygon, Unvouched};
//! let u: Unvouched<Polygon> = unreachable!();
//! let _ = u.value; // error[E0616]: field `value` is private
//! ```
//!
//! And it does not deref, so a method expecting the inner type is not
//! reachable by autoderef:
//!
//! ```compile_fail
//! use geometres::types::{Polygon, Unvouched};
//! let u: Unvouched<Polygon> = unreachable!();
//! let _ = u.exterior(); // error[E0599]: no method named `exterior`
//! ```

use crate::GeometryError;
use crate::ops::CoverageIssue;
use crate::ops::validate::ValidityIssue;

// ======================================================================
// Unvouched
// ======================================================================

/// A geometry that could be produced but **not vouched for**, because it is
/// invalid.
///
/// Returned by the opt-in unvouched entry points of operations whose default
/// entry point refuses instead, such as `simplify`. A caller that asks for an
/// unvouched result is stating that it would rather have a broken answer it can
/// inspect than no answer; the wrapper makes that choice explicit and keeps the
/// broken value out of every signature that expects a valid one.
///
/// # Invariants
///
/// The inner value is unreachable without [`Unvouched::into_inner`] or
/// [`Unvouched::inner`], and `Unvouched<T>` does not satisfy a signature
/// expecting `T`. The documentation of src/types/result_wrappers.rs
/// gives the reasoning.
///
/// ```
/// use geometres::GeometryError;
/// use geometres::types::{Coord2d, CoordSeq, LinearRing, MaybeVouched, Polygon, UnvouchedReason};
///
/// let ring = LinearRing::new(CoordSeq::new(
///     [(6.0, 5.0), (5.0, 6.0), (2.0, 2.0), (5.0, 5.0), (4.0, 1.0)]
///         .iter()
///         .map(|&(x, y)| Coord2d::new(x, y))
///         .collect(),
/// ))?;
/// let pentagon = Polygon::new_simple(ring)?;
///
/// // Simplifying at tolerance 1.5 drops (5, 6), and the edge (6, 5)-(2, 2)
/// // then crosses the edge (5, 5)-(4, 1). The vouched entry point refuses.
/// assert!(matches!(pentagon.simplify(1.5), Err(GeometryError::Unvouched { .. })));
///
/// // The unvouched entry point returns the invalid polygon, wrapped.
/// let MaybeVouched::Unvouched(u) = pentagon.simplify_unvouched(1.5)? else {
///     unreachable!("the simplified pentagon crosses itself");
/// };
///
/// // The reason is readable without unwrapping the value.
/// assert_eq!(*u.reason(), UnvouchedReason::SelfIntersecting { edges: (0, 2) });
/// assert!(u.reason().is_locally_invalid());
///
/// // Taking the value is an explicit act.
/// let _polygon = u.into_inner();
/// # Ok::<(), GeometryError>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct Unvouched<T> {
    value: T,
    reason: UnvouchedReason,
}

impl<T> Unvouched<T> {
    /// Wrap a value that is known to be invalid, with the reason it is.
    pub fn new(value: T, reason: UnvouchedReason) -> Self {
        Self { value, reason }
    }

    /// Why this value could not be vouched for.
    ///
    /// Readable without touching the value, so a caller can decide whether to
    /// accept the geometry at all before extracting it.
    pub fn reason(&self) -> &UnvouchedReason {
        &self.reason
    }

    /// Borrow the invalid value.
    ///
    /// Named rather than reached through `Deref` so that every use of an
    /// unvouched geometry is visible at its call site.
    pub fn inner(&self) -> &T {
        &self.value
    }

    /// Take the invalid value, discarding the reason.
    pub fn into_inner(self) -> T {
        self.value
    }

    /// Split into the value and its reason.
    pub fn into_parts(self) -> (T, UnvouchedReason) {
        (self.value, self.reason)
    }
}

/// The answer of an entry point that returns an invalid result rather than
/// refusing it, when the result turned out to be valid after all.
///
/// `simplify_unvouched` runs the same vouch as `simplify` and reports what it
/// found. A result the vouch passes is `Vouched`, and one it fails only on
/// simplicity is `Unvouched` with the reason. Wrapping a valid result in
/// [`Unvouched`] would require inventing a reason, and the wrapper means
/// "known to be invalid", so the two outcomes are kept apart.
///
/// Both invariants of [`Unvouched<T>`] hold: the inner value is reached
/// only by matching or by [`MaybeVouched::vouched`], and the type does not
/// satisfy a signature expecting `T`.
///
/// The unvouched value may have a different type `U` from the vouched one. An
/// offset that fails the vouch has no `MultiPolygon` to return, since the
/// refused rings cannot be assembled into one, so it carries them in a
/// structure of its own. `U` defaults to `T` for the cases, like `simplify`,
/// where the refused value is still an instance of the type.
#[derive(Clone, Debug, PartialEq)]
pub enum MaybeVouched<T, U = T> {
    /// The result passed every check.
    Vouched(T),
    /// The result exists but did not pass the vouch; the reason says why.
    Unvouched(Unvouched<U>),
}

impl<T, U> MaybeVouched<T, U> {
    /// The vouched value, or the refusal the vouched entry point would have
    /// returned: [`GeometryError::Unvouched`] with the reason and no inputs.
    ///
    /// # Errors
    ///
    /// When the result was not vouched for.
    pub fn vouched(self) -> Result<T, GeometryError> {
        match self {
            Self::Vouched(v) => Ok(v),
            Self::Unvouched(u) => Err(GeometryError::Unvouched {
                reason: *u.reason(),
                inputs: Vec::new(),
            }),
        }
    }

    /// `true` for a vouched result.
    pub fn is_vouched(&self) -> bool {
        matches!(self, Self::Vouched(_))
    }
}

impl<V> MaybeVouched<V, V> {
    /// The value, whether vouched or not, with the reason if it was not.
    /// Available only where both arms hold the same type.
    pub fn into_parts(self) -> (V, Option<UnvouchedReason>) {
        match self {
            Self::Vouched(v) => (v, None),
            Self::Unvouched(u) => {
                let (v, r) = u.into_parts();
                (v, Some(r))
            }
        }
    }
}

/// Why a geometry could not be vouched for.
///
/// # The two classes, and why the distinction is required
///
/// **Locally invalid**. The geometry is roughly the right answer and locally
/// broken. A consumer may be able to repair it, snap it, or accept it with a
/// warning.
///
/// **Structurally wrong**. The true answer has a different *shape* than the
/// return type can express. A single ring where the truth is two disjoint
/// polygons is not a degraded version of the right answer; it is a different
/// answer, and it is worse than no answer at all.
///
/// A consumer may accept the first class while still refusing the second, so
/// the classes must be distinguishable without string matching. Use
/// [`UnvouchedReason::is_locally_invalid`] and
/// [`UnvouchedReason::is_structurally_wrong`].
///
/// # Edge indices, not constructed points
///
/// `SelfIntersecting` carries the indices of the two offending edges rather
/// than the point where they cross. An index is exact and requires no
/// construction, while an intersection point requires a division and is
/// therefore approximate. A diagnostic that is exact is worth more than one
/// that is convenient, and an approximate point in a correctness-critical
/// diagnostic is the kind of plausible-looking value this project avoids.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnvouchedReason {
    // ---- locally invalid ----
    /// Two edges of the result cross. **Locally invalid.**
    ///
    /// `edges` holds the indices of the two crossing edges, where edge `i`
    /// runs from vertex `i` to vertex `i + 1` (and the last wraps). Exact: no
    /// intersection point is constructed.
    SelfIntersecting {
        /// Indices of the two crossing edges.
        edges: (usize, usize),
    },
    /// A ring of the result has an exact signed area of zero. **Locally
    /// invalid.**
    DegenerateRing {
        /// Index of the offending ring: 0 is the exterior, 1.. are holes.
        ring: usize,
    },

    // ---- structurally wrong ----
    /// The true answer is several disjoint polygons. **Structurally wrong.**
    ///
    /// Returning one of them, or their union as a single ring, would be a
    /// different geometry rather than an approximation of the right one.
    RequiresMultiplePolygons,
    /// The true answer has a hole this result cannot express. **Structurally
    /// wrong.**
    ///
    /// The trapped-void case: a region genuinely outside the result is
    /// enclosed by its boundary.
    RequiresHole,
    /// A hole of the result contains its own exterior. **Structurally wrong.**
    ///
    /// The true answer is the exterior minus the grown hole, a subtraction
    /// this return type cannot represent. `hole` names which hole, because
    /// "some hole is wrong" is a much weaker diagnostic than "hole 2 is".
    HoleExceedsExterior {
        /// Index of the offending hole, counting from 0.
        hole: usize,
    },

    // ---- rounding-boundary refusals ----
    //
    // These name a ring by its index in the extraction's ring list, which
    // mixes outers and holes in boundary order, rather than by the
    // exterior/hole indexing the variants above use. An extraction that is
    // refused has not yet been assembled into polygons, so there is no
    // exterior to count from.
    /// Rounding collapsed a ring: consecutive vertices became equal, or fewer
    /// than three distinct vertices remain. **Locally invalid.**
    RingCollapsed {
        /// Index of the ring in extraction order.
        ring: usize,
    },
    /// Two edges of a ring touch other than at the vertex consecutive edges
    /// share: a crossing, a pinch, a vertex on a non-adjacent edge, or a
    /// fold-back. **Locally invalid.**
    ///
    /// `edges` are the indices of the two edges in contact, numbered as in
    /// [`Self::SelfIntersecting`]. Always a real pair: every contact
    /// configuration is identified by the edges involved, so none is
    /// invented.
    ///
    /// `ring` is the producer's ring index. The overlay numbers rings in
    /// extraction order; `Polygon::simplify` uses it for a hole, numbered 1.. in the
    /// polygon's hole order, with the exterior reported as
    /// [`Self::SelfIntersecting`].
    RingSelfIntersecting {
        /// Index of the ring: extraction order (overlay) or hole order
        /// (simplify); see above.
        ring: usize,
        /// Indices of the two edges in contact within that ring.
        edges: (usize, usize),
    },
    /// Two rings that must be interior-disjoint are not, after rounding.
    /// **Structurally wrong**. The result is not a valid multipolygon, and a
    /// consumer cannot repair two overlapping outers locally.
    RingsOverlap {
        /// Index of the first ring in extraction order.
        first: usize,
        /// Index of the second ring in extraction order.
        second: usize,
    },
    /// A clockwise ring lies inside no counterclockwise ring, or crosses the
    /// one it should lie inside, after rounding. **Structurally wrong.**
    HoleMisplaced {
        /// Index of the hole ring in extraction order.
        ring: usize,
    },

    // ---- simplification refusals ----
    /// A ring of one part of a simplified multipolygon touches or crosses
    /// itself. **Locally invalid.**
    ///
    /// `part` indexes the multipolygon's parts; `ring` the part's rings, 0
    /// the exterior and 1.. the holes; `edges` as in
    /// [`Self::SelfIntersecting`]. The single-polygon entry points report
    /// the exterior as `SelfIntersecting` and a hole as
    /// `RingSelfIntersecting`; this variant adds the part, which those
    /// cannot carry.
    PartSelfIntersecting {
        /// Index of the part.
        part: usize,
        /// Index of the ring within the part: 0 is the exterior, 1.. holes.
        ring: usize,
        /// Indices of the two edges in contact within that ring.
        edges: (usize, usize),
    },

    // ---- lineal refusals (clipping; overlay residue) ----
    /// Rounding a constructed vertex of a line string in the result made two
    /// consecutive vertices equal: the line is shorter, at that point, than
    /// the rounding can resolve. **Locally invalid.**
    ///
    /// Raised for a clipped run, and for a lower-dimensional residue line of
    /// a complete overlay result.
    RunCollapsed {
        /// Index of the line string in the order the result would have
        /// listed it: for a clip, along the input line, except that on a
        /// closed line a run passing through the closing point is listed
        /// last; for overlay residue, in the residue's canonical order.
        run: usize,
    },

    // ---- topology-preserving simplification ----
    /// The input has a validity issue. The operation computed its result
    /// anyway, since operations do not refuse invalid input, but its
    /// argument holds only for a valid input, so the result carries no
    /// vouch. `issue` is the input's first issue, in the order
    /// [`validity_issues`](crate::ops::validity_issues) reports them.
    /// **Locally invalid.**
    InvalidInput {
        /// The input's first validity issue.
        issue: ValidityIssue,
    },
    /// The result has a validity issue, although the operation's argument
    /// says a valid input cannot produce one. The check exists as a defense;
    /// it is expected never to fire. `issue` is the result's first issue.
    /// **Locally invalid.**
    InvalidResult {
        /// The result's first validity issue.
        issue: ValidityIssue,
    },

    // ---- Delaunay triangulation ----
    /// The triangulation failed its own certification: an output triangle
    /// is not strictly counterclockwise, or an input point lies strictly
    /// inside an output triangle's circumcircle, or an internal invariant
    /// of the construction broke. **Structurally wrong.**
    ///
    /// The checks run on every call (the canonicalization pass computes the
    /// same predicates anyway), and the construction's argument says they
    /// cannot fail, so this is a defense that is expected never to fire. It
    /// carries no location: no vouched output exists to index into.
    NotDelaunay,
    /// The concave hull failed its certification: the eroded region's
    /// boundary did not trace as one simple counterclockwise ring, the ring
    /// failed the polygon constructor's checks, or an input point lies
    /// strictly outside the result. **Structurally wrong.**
    ///
    /// The removal rule's combinatorial argument says none of these can
    /// happen; the checks run on every call as a defense and are expected
    /// never to fire. No location is carried, because no vouched output
    /// exists.
    NotAHull,

    // ---- Visvalingam–Whyatt simplification ----
    /// The simplified geometry has a validity issue. Visvalingam–Whyatt
    /// elimination does not preserve topology (removing vertices can make a
    /// ring cross itself, push a hole outside its exterior, or make parts
    /// overlap), so this is the ordinary refusal of such a result, not a
    /// defect. `issue` is the result's first issue, in the order
    /// [`validity_issues`](crate::ops::validity_issues) reports them.
    /// **Locally invalid.**
    SimplifiedInvalid {
        /// The simplified geometry's first validity issue.
        issue: ValidityIssue,
    },

    // ---- coverage simplification ----
    /// The simplified polygons failed the closing coverage check
    /// ([`coverage_validate`](crate::ops::coverage_validate)), although the
    /// operation's argument says a valid coverage cannot produce that. The
    /// check exists as a defense and is expected never to fire. `first` is
    /// the first issue found. **Structurally wrong.**
    NotACoverage {
        /// The first coverage issue of the simplified polygons.
        first: CoverageIssue,
    },
}

impl UnvouchedReason {
    /// Is this the locally invalid class: roughly right, locally broken?
    ///
    /// A consumer that can repair or tolerate local damage may accept these.
    pub fn is_locally_invalid(self) -> bool {
        matches!(
            self,
            Self::SelfIntersecting { .. }
                | Self::DegenerateRing { .. }
                | Self::RingCollapsed { .. }
                | Self::RingSelfIntersecting { .. }
                | Self::PartSelfIntersecting { .. }
                | Self::RunCollapsed { .. }
                | Self::InvalidInput { .. }
                | Self::InvalidResult { .. }
                | Self::SimplifiedInvalid { .. }
        )
    }

    /// Is this the structurally wrong class: a different answer rather than
    /// a degraded one?
    ///
    /// A consumer should refuse these regardless of its tolerance for local
    /// damage.
    ///
    /// Written as the complement of [`Self::is_locally_invalid`] rather than as
    /// its own `matches!` list: the enum is `#[non_exhaustive]`, and a variant
    /// added without updating both lists would otherwise be silently neither
    /// class. Here it defaults to structurally wrong, which is the refusing
    /// answer and so the safe one.
    pub fn is_structurally_wrong(self) -> bool {
        !self.is_locally_invalid()
    }
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used)]
// test-lint-posture: see the workspace lint comment in Cargo.toml
#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::range_plus_one,
    clippy::no_effect_underscore_binding
)]
mod tests {
    use super::*;
    use crate::types::{Coord2d, CoordSeq, LinearRing, Polygon};

    fn triangle() -> Polygon {
        Polygon::new_simple(
            LinearRing::new(CoordSeq::new(vec![
                Coord2d::new(0.0, 0.0),
                Coord2d::new(1.0, 0.0),
                Coord2d::new(1.0, 1.0),
            ]))
            .unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn unvouched_carries_its_value_and_reason() {
        let u = Unvouched::new(
            triangle(),
            UnvouchedReason::SelfIntersecting { edges: (0, 2) },
        );
        assert_eq!(
            *u.reason(),
            UnvouchedReason::SelfIntersecting { edges: (0, 2) }
        );
        assert_eq!(u.inner(), &triangle());
        let (v, r) = u.into_parts();
        assert_eq!(v, triangle());
        assert!(r.is_locally_invalid());
    }

    /// Every variant belongs to exactly one class.
    ///
    /// The two predicates must partition the enum: a variant in neither would
    /// be invisible to a consumer switching on the class, and a variant in
    /// both would make the distinction meaningless.
    #[test]
    fn every_reason_is_in_exactly_one_class() {
        let all = [
            UnvouchedReason::SelfIntersecting { edges: (0, 1) },
            UnvouchedReason::DegenerateRing { ring: 0 },
            UnvouchedReason::RequiresMultiplePolygons,
            UnvouchedReason::RequiresHole,
            UnvouchedReason::HoleExceedsExterior { hole: 0 },
            UnvouchedReason::InvalidInput {
                issue: ValidityIssue::DegenerateRing { ring: 0 },
            },
            UnvouchedReason::InvalidResult {
                issue: ValidityIssue::DegenerateRing { ring: 0 },
            },
            UnvouchedReason::SimplifiedInvalid {
                issue: ValidityIssue::DegenerateRing { ring: 0 },
            },
            UnvouchedReason::NotACoverage {
                first: crate::ops::CoverageIssue::Overlap { a: 0, b: 1 },
            },
        ];
        for r in all {
            assert_ne!(
                r.is_locally_invalid(),
                r.is_structurally_wrong(),
                "{r:?} is in both classes or neither"
            );
        }
    }

    #[test]
    fn the_two_classes_are_populated_as_documented() {
        assert!(UnvouchedReason::SelfIntersecting { edges: (0, 1) }.is_locally_invalid());
        assert!(UnvouchedReason::DegenerateRing { ring: 0 }.is_locally_invalid());

        assert!(UnvouchedReason::RequiresMultiplePolygons.is_structurally_wrong());
        assert!(UnvouchedReason::RequiresHole.is_structurally_wrong());
        assert!(UnvouchedReason::HoleExceedsExterior { hole: 0 }.is_structurally_wrong());
    }

    /// `SelfIntersecting` carries indices, not coordinates.
    ///
    /// Pinned as a test because the tempting "improvement" is to add the
    /// crossing point, which would make an exact diagnostic approximate.
    #[test]
    fn self_intersecting_carries_edge_indices() {
        let r = UnvouchedReason::SelfIntersecting { edges: (3, 7) };
        let UnvouchedReason::SelfIntersecting { edges } = r else {
            panic!("wrong variant")
        };
        assert_eq!(edges, (3usize, 7usize));
    }

    #[test]
    fn hole_exceeds_exterior_names_the_hole() {
        let r = UnvouchedReason::HoleExceedsExterior { hole: 2 };
        let UnvouchedReason::HoleExceedsExterior { hole } = r else {
            panic!("wrong variant")
        };
        assert_eq!(hole, 2);
    }
}
