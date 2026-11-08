// src/topology/mod.rs

//! Exact topology: DE-9IM computation and the labeling types it produces.
//!
//! # Layout
//!
//! - [`de9im`] holds [`IntersectionMatrix`](de9im::IntersectionMatrix), the
//!   primary public output of this layer. It is computed by exact ring-edge
//!   and sub-interval classification in the three `*_analysis` modules,
//!   rather than by a sweep, a DCEL or a flood fill. That pipeline was
//!   correct only when boundaries did not meet.
//! - [`label`] holds [`Location`](label::Location) and
//!   [`Dimension`](label::Dimension).
//!
//! Nothing in this module constructs an intersection coordinate: every cell of
//! the matrix is decided by `orient2d` on stored coordinates and comparisons
//! of stored coordinates, and, for a collection operand (the `collection` module),
//! including exact rational comparisons of named crossings, which construct nothing.

pub(crate) mod areal;
pub(crate) mod collection;
pub mod de9im;
pub(crate) mod exact_point;
pub mod label;
pub(crate) mod line_analysis;
pub(crate) mod point_analysis;
pub(crate) mod polygon_analysis;
