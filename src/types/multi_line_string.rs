// src/types/multi_line_string.rs

//! Multi-line-string geometry (`MultiLineString`).

use crate::types::line_string::LineString;

/// A collection of zero or more [`LineString`] values.
///
/// `MultiLineString::new` is infallible: each `LineString` is already
/// validated by [`LineString::new`]. An empty `MultiLineString` is valid and
/// represents an empty geometry.
///
/// # Why construction cannot fail
///
/// Unlike `LineString::new` and `Polygon::new`, which validate a raw
/// coordinate sequence, `MultiLineString` is assembled from complete,
/// already-validated geometries. The collection adds no constraint of its
/// own: there is no minimum count, constituent line strings may touch,
/// cross, overlap, or repeat, and the constituents need not be connected.
///
/// # Invariants
///
/// - Every element is a valid `LineString` (at least 2 vertices, all
///   coordinates finite, no consecutive duplicate vertices), enforced at
///   `LineString` construction time.
/// - The collection itself may be empty.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiLineString {
    line_strings: Vec<LineString>,
}

impl MultiLineString {
    /// Construct a `MultiLineString` from a `Vec<LineString>`.
    ///
    /// Infallible. An empty vector is accepted and represents an empty
    /// geometry. Element order is preserved.
    pub fn new(line_strings: Vec<LineString>) -> Self {
        Self { line_strings }
    }

    /// Construct an empty `MultiLineString` (zero line strings).
    ///
    /// Equivalent to `MultiLineString::new(vec![])`.
    pub fn empty() -> Self {
        Self {
            line_strings: vec![],
        }
    }

    /// Return a slice over the constituent line strings, in construction order.
    pub fn line_strings(&self) -> &[LineString] {
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
    use crate::types::coord::{Coord2d, CoordSeq};

    fn ls(pts: &[(f64, f64)]) -> LineString {
        LineString::new(CoordSeq::new(
            pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
        ))
        .unwrap()
    }

    #[test]
    fn single_line_string_round_trips() {
        let a = ls(&[(0.0, 0.0), (1.0, 1.0)]);
        let mls = MultiLineString::new(vec![a.clone()]);
        assert_eq!(mls.line_strings(), &[a]);
        assert_eq!(mls.num_line_strings(), 1);
        assert!(!mls.is_empty());
    }

    #[test]
    fn multiple_line_strings_preserve_order() {
        let a = ls(&[(0.0, 0.0), (1.0, 0.0)]);
        let b = ls(&[(5.0, 5.0), (6.0, 6.0), (7.0, 5.0)]);
        let mls = MultiLineString::new(vec![a.clone(), b.clone()]);
        assert_eq!(mls.num_line_strings(), 2);
        assert_eq!(mls.line_strings()[0], a);
        assert_eq!(mls.line_strings()[1], b);
    }

    #[test]
    fn crossing_and_repeated_members_are_permitted() {
        let a = ls(&[(0.0, 0.0), (2.0, 2.0)]);
        let b = ls(&[(0.0, 2.0), (2.0, 0.0)]);
        let mls = MultiLineString::new(vec![a.clone(), b, a]);
        assert_eq!(mls.num_line_strings(), 3);
    }

    #[test]
    fn empty_is_empty() {
        let mls = MultiLineString::empty();
        assert!(mls.is_empty());
        assert_eq!(mls.num_line_strings(), 0);
        assert_eq!(mls.line_strings(), &[]);
    }

    #[test]
    fn empty_equals_new_with_empty_vec() {
        assert_eq!(MultiLineString::empty(), MultiLineString::new(vec![]));
        assert!(MultiLineString::new(vec![]).is_empty());
    }
}
