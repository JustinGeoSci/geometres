// src/error.rs

//! Geometry construction and validation errors.

/// An error produced by a geometry constructor or validator.
///
/// All variants carry enough context to identify the specific coordinate
/// or index that caused the failure without requiring the caller to
/// re-examine the input.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in minor
/// versions as additional geometry types are introduced. Match arms should
/// include a `_` catch-all or match only the variants they intend to handle.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum GeometryError {
    /// The coordinate sequence has fewer vertices than the type requires.
    InsufficientVertices {
        /// Minimum number of vertices required by this geometry type.
        required: usize,
        /// Number of vertices actually present in the input.
        found: usize,
    },
    /// Two consecutive vertices in the sequence are identical.
    ConsecutiveDuplicateVertices {
        /// Index of the first of the two duplicate vertices.
        index: usize,
    },
    /// A coordinate value is non-finite (NaN, positive infinity, or negative infinity).
    NonFiniteCoordinate {
        /// Index of the coordinate containing the non-finite value.
        coord_index: usize,
        /// Axis that is non-finite: `0` for x, `1` for y, `2` for z
        /// (z is reported only by the 3D types under `feature = "3d"`).
        axis: u8,
    },
    /// The exterior ring of a polygon has incorrect winding order.
    /// Exterior rings must be counterclockwise.
    IncorrectWindingOrder,
    /// A hole ring has incorrect winding order.
    /// Hole rings must be clockwise.
    IncorrectHoleWindingOrder {
        /// Index of the hole ring with incorrect winding order.
        hole_index: usize,
    },
    /// A hole ring does not lie inside the exterior ring.
    ///
    /// A polygon is the region its exterior bounds minus the regions its
    /// holes bound, so every hole must be carved from material the exterior
    /// actually encloses. This is returned when a hole is disjoint from the
    /// exterior, when its boundary crosses the exterior's (properly, or by
    /// passing through it at a point), when it shares a stretch of the
    /// exterior's boundary, or when it surrounds the exterior. A hole that
    /// lies inside and *touches* the exterior at isolated points is accepted,
    /// however many: touches that cut the interior in two (a hole touching
    /// the exterior twice, or a triangle whose vertices all lie on the
    /// exterior) make a disconnected interior, which construction does not
    /// check and `ops::validity_issues` reports.
    ///
    /// Decided exactly, by `arithmetica::orient2d` on the stored coordinates;
    /// no tolerance is involved.
    HoleNotInsideExterior {
        /// Index of the offending hole ring, counting from 0.
        hole_index: usize,
    },
    /// Two hole rings have overlapping interiors.
    ///
    /// Holes must be interior-disjoint: one may not lie inside another, nor
    /// cross it, nor coincide with it, or the net area subtracts the shared
    /// region twice. Holes that touch at a single vertex are accepted. Decided
    /// exactly.
    OverlappingHoles {
        /// Index of the first hole of the pair, counting from 0.
        first: usize,
        /// Index of the second hole of the pair; always greater than `first`.
        second: usize,
    },
    /// A ring has an exact signed area of zero (every vertex collinear, or
    /// lobes that cancel exactly) and cannot be oriented.
    DegenerateRing,
    /// Two polygonal leaves of a `GeometryCollection` overlap: their
    /// interiors meet, or their boundaries share a stretch.
    ///
    /// `relate` and the predicates accept a collection whose polygonal
    /// leaves together form a valid `MultiPolygon`, whose parts meet at most
    /// at finitely many points, decided exactly by the same placement test
    /// `MultiPolygon::new` applies. A collection outside that restriction
    /// is refused rather than answered approximately; full union semantics
    /// over overlapping elements is future state. The collection itself is
    /// still *valid* (OGC puts no condition between elements), so this is a
    /// relate-level refusal, not a validity issue. A caller that wants the
    /// union's relation unions the elements first and relates the vouched
    /// result.
    ///
    /// `first` and `second` are the two leaves' paths, in the convention of
    /// [`ops::validate::LocatedIssue`](crate::ops::LocatedIssue): element
    /// indices through nested collections, then, for a part of a
    /// multipolygon leaf, the part index.
    CollectionElementsOverlap {
        /// Path of the earlier offending polygonal leaf.
        first: Vec<usize>,
        /// Path of the later offending polygonal leaf.
        second: Vec<usize>,
    },
    /// An intermediate or final value of the computation was not finite
    /// (overflow from coordinates of extreme magnitude). No finite answer
    /// can be vouched for, so none is returned.
    NonFiniteResult,
    /// A tolerance or distance parameter was negative or non-finite.
    ///
    /// [`Polygon::simplify`](crate::types::Polygon::simplify) and the other
    /// `simplify` methods take a distance; zero is the identity, and a value
    /// that is not a distance is rejected rather than read as "keep every
    /// vertex".
    /// [`dwithin`](crate::ops::dwithin)'s `d` is held to the same rule.
    InvalidTolerance {
        /// The rejected value.
        tolerance: f64,
    },
    /// [`coverage_union`](crate::ops::coverage_union) was given polygons that
    /// are not a valid coverage. Carries the first issue
    /// [`coverage_validate`](crate::ops::coverage_validate) reports; for
    /// general input, use a constructive union.
    InvalidCoverage {
        /// The first coverage issue found.
        first: crate::ops::CoverageIssue,
    },
    /// [`polygonize_noded`](crate::ops::polygonize_noded) was given
    /// linework that is not noded: the two segments named cross, overlap
    /// collinearly, or one has an endpoint in the other's interior. Node the
    /// input first with a constructive operation.
    NotNoded {
        /// The first offending pair of segments, each by its endpoints.
        first_pair: (
            (crate::types::Coord2d, crate::types::Coord2d),
            (crate::types::Coord2d, crate::types::Coord2d),
        ),
    },
    /// A [`distance`](crate::ops::distance()) operand has no points.
    ///
    /// The distance to nothing is undefined. Returning 0 for an empty
    /// operand, which reads as "they touch", is a plausible wrong answer,
    /// and this crate refuses to give it. (`dwithin` with an empty operand is simply
    /// `false`: no point of one geometry is within anything of the other.)
    EmptyOperand,
    /// Two parts of a `MultiPolygon` overlap: their interiors meet, or their
    /// boundaries share a stretch of positive length. OGC allows parts to
    /// meet at finitely many points only.
    OverlappingParts {
        /// Index of the first part of the offending pair.
        first: usize,
        /// Index of the second part.
        second: usize,
    },
    /// A `MultiPolygon` could not be converted to a single `Polygon` because it
    /// does not contain exactly one polygon (it is empty or holds several).
    NotASinglePolygon,
    /// `line_sequence` found a connected component of the line work with
    /// more than two nodes of odd degree: no path runs through each of its
    /// lines once, so the lines cannot be put in sequence.
    NotSequenceable,
    /// `split_by_point` was given a point that is not on the line's
    /// interior: it is the line's first or last vertex, or not on the line
    /// at all. A point near the line is not snapped to it.
    PointNotOnLineInterior,
    /// An [`Envelope`](crate::types::Envelope) was given bounds where
    /// `min > max` on at least one axis.
    ///
    /// A zero-extent envelope (`min == max`) is *valid*. It is a legitimate
    /// point-sized or line-shaped query rectangle. Only an inverted one is an
    /// error, because every containment test against it would be silently
    /// false.
    InvertedEnvelope,
    /// A spatial index was asked to hold more entries than it can address.
    ///
    /// Index entries are identified by `u32`, so a set of more than
    /// `u32::MAX` parts or edges cannot be indexed. This is an error, never a
    /// truncation: an index that silently dropped entries would have false
    /// negatives. The native codec caps every count at `u32`, so no geometry
    /// that round-trips through it reaches this.
    TooManyEntries {
        /// The number of entries requested.
        count: usize,
    },

    /// A `GeometryCollection` would nest more than its limit of collection
    /// levels.
    ///
    /// The limit, `GeometryCollection::MAX_NESTING` (32), is a permanent
    /// contract: it bounds the recursion of the derived traits and of any
    /// caller's code that walks a collection, and every collection that
    /// exists decodes.
    NestingTooDeep {
        /// The limit that was exceeded.
        limit: usize,
    },

    /// An operation refused its result: the result failed an exact vouch.
    ///
    /// Carries the same reason the unvouched entry point would attach, so a
    /// strict caller and a tolerant one see one diagnostic.
    Unvouched {
        /// The check that failed.
        reason: crate::types::UnvouchedReason,
        /// The inputs implicated in the refusal, by index into the input
        /// collection, sorted and deduplicated. For an operation over a
        /// collection, these are the polygons whose boundaries form the
        /// offending edges or ring. Filled by
        /// `overlay::union_all`; empty from the binary operations, whose two
        /// operands are already named by the call.
        inputs: Vec<usize>,
    },
}

impl core::fmt::Display for GeometryError {
    // One arm per variant; the length is the variant count, not complexity.
    #[allow(clippy::too_many_lines)]
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InsufficientVertices { required, found } => {
                write!(f, "insufficient vertices: need {required}, found {found}")
            }
            Self::ConsecutiveDuplicateVertices { index } => {
                write!(f, "consecutive duplicate vertices at index {index}")
            }
            Self::NonFiniteCoordinate { coord_index, axis } => {
                let ax = match axis {
                    0 => "x",
                    1 => "y",
                    _ => "z",
                };
                write!(f, "non-finite {ax} coordinate at index {coord_index}")
            }
            Self::IncorrectWindingOrder => write!(f, "exterior ring must be counterclockwise"),
            Self::IncorrectHoleWindingOrder { hole_index } => {
                write!(f, "hole {hole_index} must be clockwise")
            }
            Self::DegenerateRing => write!(f, "ring has zero area"),
            Self::HoleNotInsideExterior { hole_index } => {
                write!(f, "hole {hole_index} does not lie inside the exterior ring")
            }
            Self::OverlappingHoles { first, second } => {
                write!(f, "holes {first} and {second} have overlapping interiors")
            }
            Self::NonFiniteResult => write!(f, "computation produced a non-finite value"),
            Self::Unvouched {
                reason: crate::types::UnvouchedReason::NotDelaunay,
                ..
            } => {
                write!(
                    f,
                    "triangulation failed its Delaunay certification (a defect in this \
                     library, not in the input)"
                )
            }
            Self::Unvouched {
                reason: crate::types::UnvouchedReason::NotAHull,
                ..
            } => {
                write!(
                    f,
                    "concave hull failed its containment certification (a defect in this \
                     library, not in the input)"
                )
            }
            Self::Unvouched { reason, inputs } => {
                write!(f, "result failed the exact vouch: {reason:?}")?;
                if inputs.is_empty() {
                    Ok(())
                } else {
                    write!(f, " (implicated inputs: {inputs:?})")
                }
            }
            Self::InvalidCoverage { first } => {
                write!(f, "the polygons are not a valid coverage: {first:?}")
            }
            Self::NotNoded {
                first_pair: ((a, b), (c, d)),
            } => write!(
                f,
                "the linework is not noded: segment ({}, {})–({}, {}) meets segment \
                 ({}, {})–({}, {}) away from a shared endpoint",
                a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y
            ),
            Self::EmptyOperand => {
                write!(f, "distance is undefined for an empty operand")
            }
            Self::InvalidTolerance { tolerance } => {
                write!(
                    f,
                    "tolerance must be a finite non-negative distance, got {tolerance}"
                )
            }
            Self::OverlappingParts { first, second } => write!(
                f,
                "multi-polygon parts {first} and {second} overlap (OGC allows parts to meet at points only)"
            ),
            Self::NotASinglePolygon => {
                write!(f, "multi-polygon does not contain exactly one polygon")
            }
            Self::NotSequenceable => write!(
                f,
                "line work has a connected component with more than two nodes of odd degree; it cannot be sequenced"
            ),
            Self::PointNotOnLineInterior => write!(
                f,
                "point is not on the line's interior (it is an endpoint, or not on the line)"
            ),
            Self::InvertedEnvelope => {
                write!(
                    f,
                    "envelope bounds are inverted: min exceeds max on an axis"
                )
            }
            Self::CollectionElementsOverlap { first, second } => {
                write!(
                    f,
                    "collection elements at paths {first:?} and {second:?} overlap \
                     (relate accepts collections whose polygonal leaves meet at most \
                     at finitely many points)"
                )
            }
            Self::NestingTooDeep { limit } => {
                write!(f, "geometry collection nests more than {limit} levels deep")
            }
            Self::TooManyEntries { count } => {
                write!(
                    f,
                    "spatial index cannot hold {count} entries (the limit is {})",
                    u32::MAX
                )
            }
        }
    }
}

impl std::error::Error for GeometryError {}
