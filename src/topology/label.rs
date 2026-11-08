// src/topology/label.rs

//! Topological labeling types: `Location` and `Dimension`.
//!
//! [`Location`] records where a point or region lies relative to one geometry.
//! [`Dimension`] records the maximum dimension of an intersection and is the
//! cell type of [`IntersectionMatrix`](crate::topology::de9im::IntersectionMatrix).
//!
//! Neither type has any connection to a doubly-connected edge list, so they
//! live directly under `topology` rather than in a directory of their own.

/// The topological location of a region relative to one input geometry.
///
/// Returned by point location and used as the cell semantics of the DE-9IM
/// matrix: a cell records which pair of locations the row and column refer to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Location {
    /// The region lies inside the geometry (in its interior).
    Interior,
    /// The region lies on the boundary of the geometry.
    Boundary,
    /// The region lies outside the geometry (in its exterior).
    Exterior,
}

/// The maximum dimension of a geometric intersection.
///
/// Values correspond to the DE-9IM specification:
///
/// | Value | Meaning |
/// |-------|---------|
/// | `Empty` (−1) | The intersection is empty. |
/// | `Point` (0) | At least one point, no line or area. |
/// | `Line` (1) | At least one line segment, no area. |
/// | `Area` (2) | The intersection contains area. |
///
/// The numeric discriminant values match the OGC DE-9IM convention; the
/// `PartialOrd`/`Ord` derive is therefore meaningful and returns the
/// expected comparison (e.g., `Area > Line`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Dimension {
    /// The intersection is empty (OGC notation: `F`).
    Empty = -1,
    /// The intersection contains at least one point but no line or area
    /// (OGC notation: `0`).
    Point = 0,
    /// The intersection contains at least one line segment but no area
    /// (OGC notation: `1`).
    Line = 1,
    /// The intersection contains area (OGC notation: `2`).
    Area = 2,
}
