// src/types3d/multi_line_string3d.rs

//! Collection of 3D line strings (`MultiLineString3d`).
//!
//! Follows [`MultiLineString`](crate::types::MultiLineString) exactly, with
//! an infallible constructor over already-validated elements.

use crate::types3d::LineString3d;

/// A collection of zero or more 3D line strings.
///
/// Semantics match OGC `MultiLineString`: no minimum count, no ordering
/// requirement, components may touch or cross.
///
/// # Invariants
///
/// - The collection itself may be empty.
/// - Every element is a `LineString3d`, so every component has at least two
///   vertices and all finite coordinates. `LineString3d::new` validates both.
///
/// # Boundary
///
/// No boundary accessor is provided. The 2D multipart boundary rule is mod-2
/// over component endpoints, and the 3D rule should extend a decided one
/// rather than be invented here. No current operation needs it.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiLineString3d {
    line_strings: Vec<LineString3d>,
}

impl MultiLineString3d {
    /// Construct a `MultiLineString3d` from a `Vec<LineString3d>`.
    ///
    /// Infallible. An empty vector is accepted and represents an empty
    /// geometry. Element order is preserved.
    pub fn new(line_strings: Vec<LineString3d>) -> Self {
        Self { line_strings }
    }

    /// Construct an empty `MultiLineString3d` (zero components).
    pub fn empty() -> Self {
        Self {
            line_strings: vec![],
        }
    }

    /// Return a slice over the constituent line strings, in construction order.
    pub fn line_strings(&self) -> &[LineString3d] {
        &self.line_strings
    }

    /// Return the number of line strings in the collection.
    pub fn num_line_strings(&self) -> usize {
        self.line_strings.len()
    }

    /// Return `true` if the collection contains no line strings.
    pub fn is_empty(&self) -> bool {
        self.line_strings.is_empty()
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
    use crate::types3d::{Coord3d, CoordSeq3d};

    fn line(pts: &[(f64, f64, f64)]) -> LineString3d {
        LineString3d::new(CoordSeq3d::new(
            pts.iter().map(|&(x, y, z)| Coord3d::new(x, y, z)).collect(),
        ))
        .unwrap()
    }

    #[test]
    fn round_trips_through_its_accessors() {
        let a = line(&[(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]);
        let b = line(&[(2.0, 2.0, 2.0), (3.0, 3.0, 3.0)]);
        let m = MultiLineString3d::new(vec![a.clone(), b.clone()]);
        assert_eq!(m.line_strings(), [a, b].as_slice());
        assert_eq!(m.num_line_strings(), 2);
        assert!(!m.is_empty());
    }

    #[test]
    fn empty_is_empty_and_matches_new_with_an_empty_vec() {
        assert!(MultiLineString3d::empty().is_empty());
        assert_eq!(MultiLineString3d::empty().num_line_strings(), 0);
        assert_eq!(MultiLineString3d::empty(), MultiLineString3d::new(vec![]));
    }

    #[test]
    fn element_order_is_preserved() {
        let a = line(&[(9.0, 0.0, 0.0), (8.0, 0.0, 0.0)]);
        let b = line(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]);
        let m = MultiLineString3d::new(vec![a.clone(), b.clone()]);
        assert_eq!(m.line_strings()[0], a);
        assert_eq!(m.line_strings()[1], b);
    }
}
