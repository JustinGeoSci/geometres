// src/types/mod.rs

//! Geometry types owned by `geometres`.
//!
//! # The construction contract (permissive constructors)
//!
//! Every constructor validates, so **a value that exists has passed its
//! checks**:
//!
//! - every coordinate is finite;
//! - every ring and line has enough distinct vertices, with no consecutive
//!   duplicates;
//! - every ring's winding is correct and its exact signed area is nonzero,
//!   with exteriors counterclockwise and holes clockwise, decided by
//!   certified area signs;
//! - every hole lies inside its exterior, holes are interior-disjoint, and a
//!   multipolygon's parts meet at most at finitely many points, all decided
//!   exactly.
//!
//! The constructors are deliberately **permissive** beyond that. A value
//! that exists may still be invalid in the OGC sense, in exactly two ways:
//!
//! - **a ring may touch or cross itself**, since simplicity is not checked at
//!   construction;
//! - **a polygon's interior may be disconnected**, since holes touching the
//!   exterior and each other at points can chain into a cut.
//!
//! [`ops::is_valid`](crate::ops::is_valid) and
//! [`ops::validity_issues`](crate::ops::validity_issues) decide both,
//! exactly, and report every issue;
//! [`ops::split_pinches`](crate::ops::split_pinches) repairs a ring that
//! passes through one of its vertices twice. Operations never refuse on the
//! grounds that their input might be invalid: a caller who chooses to
//! operate on an unchecked value gets output that is correct for that value
//! as a point set, or an output that is explicitly not vouched for.

mod coord;
mod envelope;
mod geometry;
mod geometry_collection;
mod line_string;
mod linear_ring;
mod multi_line_string;
mod multi_point;
mod multi_polygon;
mod operand;
mod point;
mod polygon;
mod result_wrappers;

pub(crate) use coord::{Bbox, bbox_fold};
pub use coord::{Coord2d, CoordSeq};
pub use envelope::Envelope;
pub use geometry::Geometry;
pub use geometry_collection::{GeometryCollection, Leaves};
pub use line_string::LineString;
pub use linear_ring::LinearRing;
pub use multi_line_string::MultiLineString;
pub use multi_point::MultiPoint;
pub use multi_polygon::MultiPolygon;
pub(crate) use multi_polygon::{first_overlapping_parts_ref, parts_meet_only_at_points};
// `MultiPolygon3d` applies the 2D part rule to its parts' projections.
#[cfg(feature = "3d")]
pub(crate) use multi_polygon::first_overlapping_parts;
pub use operand::{Lineal, Polygonal};
pub use point::Point;
pub use polygon::Polygon;
// `Polygon3d` validates its holes on their XY projection with the same
// check.
pub use polygon::ring_is_degenerate;
#[cfg(feature = "3d")]
pub(crate) use polygon::validate_hole_placement;
pub use result_wrappers::{MaybeVouched, Unvouched, UnvouchedReason};
