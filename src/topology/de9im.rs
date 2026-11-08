// src/topology/de9im.rs

//! The Dimensionally Extended 9-Intersection Model (DE-9IM) matrix.
//!
//! This module provides [`IntersectionMatrix`], which records the maximum
//! dimension of the intersection between each pair of topological components
//! (interior, boundary, exterior) of two geometries, and is the primary
//! public entry point for the topology module.
//!
//! # OGC string notation
//!
//! The `Display` implementation renders the matrix as a 9-character string
//! in row-major order (Interior/Boundary/Exterior on rows for A, same columns
//! for B):
//!
//! ```text
//! II IB IE
//! BI BB BE
//! EI EB EE
//! ```
//!
//! Flattened: `II IB IE BI BB BE EI EB EE`
//!
//! Each character is `F` (Empty = −1), `0` (Point), `1` (Line), or `2` (Area).
//!
//! **Reference:** OGC 06-103r4 §6.1.2.3; JTS `IntersectionMatrix`.
//!
//! # Computation
//!
//! [`IntersectionMatrix::compute`] classifies every ring edge of each polygon
//! against the other polygon (interior, boundary, exterior, and, along a
//! shared boundary, the relative direction) and reads all nine cells from
//! those classifications. The method, its exactness argument, and the face
//! analysis behind the area cells are in
//! `polygon_analysis`.
//!
//! The matrix is not computed from a planar subdivision. A subdivision
//! pipeline is correct only when the boundaries do not meet, and the edge
//! classification has no such limit.

use super::areal::{AreaOperand, Areal, Window, WindowedAreal};
use super::label::Dimension;
use super::line_analysis::{LineSet, lines_lines_matrix, lines_polygon_matrix};
use super::point_analysis::{Target, points_matrix};
use super::polygon_analysis::polygons_matrix;
use crate::GeometryError;
use crate::types::{
    Coord2d, Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};

/// The Dimensionally Extended 9-Intersection Model matrix.
///
/// Records the maximum dimension of the intersection between each pair of
/// topological components (interior, boundary, exterior) of two geometries.
///
/// Named fields correspond directly to OGC terminology. The matrix is computed
/// by exact ring-edge classification (see [`IntersectionMatrix::compute`]).
///
/// # OGC string notation
///
/// The [`Display`](std::fmt::Display) implementation renders the matrix as a
/// 9-character string in the order:
///
/// `II IB IE BI BB BE EI EB EE`
///
/// where each character is `F` (empty / −1), `0`, `1`, or `2`.
///
/// # Pattern parsing
///
/// The [`from_pattern`](IntersectionMatrix::from_pattern) constructor parses
/// the same format, with `*` (wildcard) matching any dimension. Wildcards are
/// only meaningful for pattern matching in tests; a matrix produced by
/// [`compute`](IntersectionMatrix::compute) never contains wildcards.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IntersectionMatrix {
    /// Interior(A) ∩ Interior(B)
    pub interior_interior: Dimension,
    /// Interior(A) ∩ Boundary(B)
    pub interior_boundary: Dimension,
    /// Interior(A) ∩ Exterior(B)
    pub interior_exterior: Dimension,
    /// Boundary(A) ∩ Interior(B)
    pub boundary_interior: Dimension,
    /// Boundary(A) ∩ Boundary(B)
    pub boundary_boundary: Dimension,
    /// Boundary(A) ∩ Exterior(B)
    pub boundary_exterior: Dimension,
    /// Exterior(A) ∩ Interior(B)
    pub exterior_interior: Dimension,
    /// Exterior(A) ∩ Boundary(B)
    pub exterior_boundary: Dimension,
    /// Exterior(A) ∩ Exterior(B)
    pub exterior_exterior: Dimension,
}

impl IntersectionMatrix {
    /// Computes the DE-9IM matrix for two polygons, exactly.
    ///
    /// Every ring edge of each polygon (exterior and holes) is classified
    /// against the other polygon, and all nine cells follow from those
    /// classifications; see `polygon_analysis` for
    /// the derivation. The result is correct for every valid pair: disjoint,
    /// crossing, touching at a vertex or along an edge segment, nested with
    /// or without contact, equal, and with holes.
    ///
    /// All sign decisions use `arithmetica::orient2d` on the input
    /// coordinates, so a vertex exactly on an edge is `Boundary` and a vertex
    /// one ulp off it is not.
    pub fn compute(a: &Polygon, b: &Polygon) -> Self {
        polygons_matrix(&Areal::single(a), &Areal::single(b))
    }

    /// Computes the DE-9IM matrix for a pair of heterogeneous geometries.
    ///
    /// Every ordered pair of the six [`Geometry`] variants is supported. The
    /// pair is reduced to one of three operand classes, and the matrix comes
    /// from the first row that applies; the reversed order uses
    /// [`swap_ab`](Self::swap_ab):
    ///
    /// | Operand classes | Variants | Computation |
    /// |---|---|---|
    /// | area × area | `Polygon`, `MultiPolygon` × the same | exact ring-edge classification over all parts at once |
    /// | lines × area | `LineString`, `MultiLineString` × `Polygon`, `MultiPolygon` | exact line × area analysis |
    /// | lines × lines | `LineString`, `MultiLineString` × the same | exact line × line analysis |
    /// | points × any | `Point`, `MultiPoint`, and any empty geometry × anything | exact point-set analysis |
    ///
    /// Linear geometry uses the OGC mod-2 boundary rule: a point is in the
    /// boundary of a `MultiLineString` when it is the endpoint of an odd
    /// number of its components. Areal geometry uses the topological
    /// boundary, as OGC does: every part's rings, including a point
    /// where parts touch.
    /// Empty `MultiPoint`, `MultiLineString` and `MultiPolygon` values are
    /// the empty point set.
    ///
    /// A `MultiPolygon`'s parts are valid by construction (`MultiPolygon::new`
    /// rejects overlapping parts), so no matrix is ever computed over an
    /// invalid one.
    ///
    /// # Errors
    ///
    /// No current pair fails. Every geometry type validates finiteness when
    /// it is built (`Point::new`, and the line and polygon constructors), so a
    /// non-finite coordinate, which would make every orientation test report
    /// "collinear" and produce a plausible-looking wrong matrix, cannot reach
    /// this function.
    ///
    /// [`GeometryError::CollectionElementsOverlap`] if either operand is a
    /// `GeometryCollection` whose polygonal leaves overlap: `relate` accepts
    /// a collection whose polygonal leaves together form a valid
    /// `MultiPolygon`; the union
    /// of overlapping elements is future state, reachable explicitly by
    /// unioning the elements first.
    pub fn compute_geometry(a: &Geometry, b: &Geometry) -> Result<Self, GeometryError> {
        Ok(Self::compute_operands(
            &Operand::from_geometry(a)?,
            &Operand::from_geometry(b)?,
        ))
    }

    /// The matrix of two reduced operands: the table of
    /// [`compute_geometry`](Self::compute_geometry), which delegates here,
    /// as do the generic predicates over borrowed geometries.
    pub(crate) fn compute_operands(oa: &Operand<'_>, ob: &Operand<'_>) -> Self {
        use crate::topology::collection::{Strata, merged_matrix};
        match (oa, ob) {
            // Every pair with a collection in it goes to the merged
            // assembly; a pure-stratum collection reduced to one of the
            // other kinds at construction and never reaches this arm.
            (Operand::Collection(_), _) | (_, Operand::Collection(_)) => {
                merged_matrix(&Strata::of(oa), &Strata::of(ob))
            }
            (Operand::Area(a), Operand::Area(b)) => {
                let (wa, wb) = windowed_pair(a, b);
                polygons_matrix(
                    wa.as_ref().map_or(a as &dyn AreaOperand, |w| w),
                    wb.as_ref().map_or(b as &dyn AreaOperand, |w| w),
                )
            }
            (Operand::Lines(a), Operand::Area(b)) => {
                let wb = WindowedAreal::of(b, lines_window(a));
                lines_polygon_matrix(a, wb.as_ref().map_or(b as &dyn AreaOperand, |w| w))
            }
            (Operand::Area(a), Operand::Lines(b)) => {
                let wa = WindowedAreal::of(a, lines_window(b));
                lines_polygon_matrix(b, wa.as_ref().map_or(a as &dyn AreaOperand, |w| w)).swap_ab()
            }
            (Operand::Lines(a), Operand::Lines(b)) => lines_lines_matrix(a, b),
            (Operand::Points(p), other) => points_matrix(p, &other.as_target()),
            (other, Operand::Points(p)) => points_matrix(p, &other.as_target()).swap_ab(),
        }
    }

    /// The matrix of an areal operand `a`, indexed or not, against any
    /// geometry `b`: the rows of [`compute_geometry`](Self::compute_geometry)
    /// whose A is an area, with `a` in that place. For an unindexed `a` the
    /// two functions compute the same thing by the same code (`prepared`
    /// relies on it).
    ///
    /// # Errors
    ///
    /// As [`compute_geometry`](Self::compute_geometry): a collection `b`
    /// whose polygonal leaves overlap.
    pub(crate) fn compute_area_against(
        a: &dyn AreaOperand,
        b: &Geometry,
    ) -> Result<IntersectionMatrix, GeometryError> {
        use crate::topology::collection::{Strata, merged_matrix};
        Ok(match Operand::from_geometry(b)? {
            Operand::Area(ob) => polygons_matrix(a, &ob),
            Operand::Lines(ob) => lines_polygon_matrix(&ob, a).swap_ab(),
            Operand::Points(ob) => points_matrix(&ob, &Target::Area(a)).swap_ab(),
            ref ob @ Operand::Collection(_) => merged_matrix(&Strata::of_area(a), &Strata::of(ob)),
        })
    }

    /// Returns the matrix with the roles of A and B exchanged (transpose).
    ///
    /// `relate(B, A)` is the transpose of `relate(A, B)`: the cell that
    /// records `X(A) ∩ Y(B)` for `relate(A, B)` records `Y(A) ∩ X(B)` for the
    /// exchanged call, because the intersection of two sets does not depend
    /// on which is named first. So `II`, `BB`, and `EE` are unchanged, and
    /// `IB ↔ BI`, `IE ↔ EI`, `BE ↔ EB`. Applying it twice is the identity.
    #[must_use]
    pub fn swap_ab(&self) -> Self {
        Self {
            interior_interior: self.interior_interior,
            interior_boundary: self.boundary_interior,
            interior_exterior: self.exterior_interior,
            boundary_interior: self.interior_boundary,
            boundary_boundary: self.boundary_boundary,
            boundary_exterior: self.exterior_boundary,
            exterior_interior: self.interior_exterior,
            exterior_boundary: self.boundary_exterior,
            exterior_exterior: self.exterior_exterior,
        }
    }

    /// Parses a 9-character OGC DE-9IM pattern string.
    ///
    /// # Character mapping
    ///
    /// | Char | `Dimension` |
    /// |------|-------------|
    /// | `F`  | `Empty` (−1) |
    /// | `0`  | `Point` |
    /// | `1`  | `Line` |
    /// | `2`  | `Area` |
    /// | `*`  | wildcard, stored as `Empty`; semantics handled by the caller |
    ///
    /// # Errors
    ///
    /// Returns `Err(PatternError)` if:
    /// - The string is not exactly 9 characters long.
    /// - Any character is not one of `F`, `0`, `1`, `2`, `*`.
    ///
    /// # Example
    ///
    /// ```
    /// use geometres::topology::de9im::IntersectionMatrix;
    ///
    /// let m = IntersectionMatrix::from_pattern("FF2FF1***").unwrap();
    /// assert_eq!(m.interior_exterior, geometres::topology::label::Dimension::Area);
    /// ```
    pub fn from_pattern(pattern: &str) -> Result<Self, PatternError> {
        let chars: Vec<char> = pattern.chars().collect();
        if chars.len() != 9 {
            return Err(PatternError {
                message: "pattern must be exactly 9 characters",
            });
        }

        let mut dims = [Dimension::Empty; 9];
        // The characters `F` and `*` map to the same value but are distinct in
        // the DE-9IM grammar: `F` asserts the cell is empty, and `*` asserts
        // nothing. The arms are kept separate to mirror the specification, and
        // to keep visible that this parser collapses a don't-care into Empty
        // rather than representing it.
        #[allow(clippy::match_same_arms)]
        for (i, &ch) in chars.iter().enumerate() {
            dims[i] = match ch {
                'F' => Dimension::Empty,
                '0' => Dimension::Point,
                '1' => Dimension::Line,
                '2' => Dimension::Area,
                '*' => Dimension::Empty,
                _ => {
                    return Err(PatternError {
                        message: "invalid character in pattern: expected F, 0, 1, 2, or *",
                    });
                }
            };
        }

        Ok(Self {
            interior_interior: dims[0],
            interior_boundary: dims[1],
            interior_exterior: dims[2],
            boundary_interior: dims[3],
            boundary_boundary: dims[4],
            boundary_exterior: dims[5],
            exterior_interior: dims[6],
            exterior_boundary: dims[7],
            exterior_exterior: dims[8],
        })
    }

    /// Returns all nine dimensions as an array in OGC row-major order.
    ///
    /// The order is: `II IB IE BI BB BE EI EB EE`.
    pub fn cells(&self) -> [Dimension; 9] {
        [
            self.interior_interior,
            self.interior_boundary,
            self.interior_exterior,
            self.boundary_interior,
            self.boundary_boundary,
            self.boundary_exterior,
            self.exterior_interior,
            self.exterior_boundary,
            self.exterior_exterior,
        ]
    }
}

/// A geometry reduced to the three classes the matrix computation needs: a
/// borrowed view, built without cloning the geometry.
pub(crate) enum Operand<'a> {
    /// A finite point set, possibly empty (`Point`, `MultiPoint`, empty
    /// `MultiLineString` and `MultiPolygon`, and a collection with only
    /// point leaves or none).
    Points(Vec<Coord2d>),
    /// A non-empty linear geometry with its mod-2 boundary (`LineString`,
    /// `MultiLineString`, or a collection with only line leaves).
    Lines(LineSet<'a>),
    /// A polygon, or a non-empty valid multipolygon.
    Area(Areal<'a>),
    /// A restricted `GeometryCollection` with polygonal leaves or several
    /// strata: its decomposition into interior-disjoint
    /// polygonal leaves, line leaves, and point leaves.
    Collection(crate::topology::collection::GcParts<'a>),
}

impl<'a> Operand<'a> {
    /// Reduces `g`, by the reduction of the type it holds.
    ///
    /// # Errors
    ///
    /// [`GeometryError::CollectionElementsOverlap`] for a collection whose
    /// polygonal leaves overlap.
    pub(crate) fn from_geometry(g: &'a Geometry) -> Result<Self, GeometryError> {
        Ok(match g {
            Geometry::Point(p) => Self::point(p),
            Geometry::MultiPoint(mp) => Self::multi_point(mp),
            Geometry::LineString(l) => Self::line_string(l),
            Geometry::MultiLineString(ml) => Self::multi_line_string(ml),
            Geometry::Polygon(p) => Self::polygon(p),
            Geometry::MultiPolygon(mp) => Self::multi_polygon(mp),
            Geometry::GeometryCollection(gc) => {
                return crate::topology::collection::collection_operand(gc);
            }
        })
    }

    /// A point. Its coordinates are finite by the `Point` invariant.
    pub(crate) fn point(p: &Point) -> Self {
        Operand::Points(vec![p.coord()])
    }

    /// A multipoint's points, possibly none.
    pub(crate) fn multi_point(mp: &MultiPoint) -> Self {
        Operand::Points(mp.points().iter().map(Point::coord).collect())
    }

    /// One line string.
    pub(crate) fn line_string(l: &'a LineString) -> Self {
        Operand::Lines(LineSet::single(l))
    }

    /// A multi line string; the empty point set if it has no components.
    pub(crate) fn multi_line_string(ml: &'a MultiLineString) -> Self {
        if ml.is_empty() {
            Operand::Points(Vec::new())
        } else {
            Operand::Lines(LineSet::new(ml.line_strings().iter().collect()))
        }
    }

    /// One polygon.
    pub(crate) fn polygon(p: &'a Polygon) -> Self {
        Operand::Area(Areal::single(p))
    }

    /// A multipolygon; the empty point set if it has no parts.
    pub(crate) fn multi_polygon(mp: &'a MultiPolygon) -> Self {
        if mp.is_empty() {
            Operand::Points(Vec::new())
        } else {
            Operand::Area(Areal::multi(mp))
        }
    }

    /// The topological dimension of the point set, `Empty` for an empty one.
    pub(crate) fn dimension(&self) -> Dimension {
        match self {
            Operand::Points(p) if p.is_empty() => Dimension::Empty,
            Operand::Points(_) => Dimension::Point,
            Operand::Lines(_) => Dimension::Line,
            Operand::Area(_) => Dimension::Area,
            Operand::Collection(gc) => crate::topology::collection::parts_dimension(gc),
        }
    }

    fn as_target(&self) -> Target<'_> {
        match self {
            Operand::Points(p) => Target::Points(p),
            Operand::Lines(l) => Target::Lines(l),
            Operand::Area(a) => Target::Area(a),
            // Unreachable: every pair with a collection dispatches to the
            // merged assembly before any target view is taken. Kept total,
            // and correct if ever reached, rather than a panic path.
            Operand::Collection(gc) => Target::Points(&gc.pts),
        }
    }
}

impl std::fmt::Display for IntersectionMatrix {
    /// Renders the matrix as a 9-character OGC string.
    ///
    /// Characters: `F` = Empty, `0` = Point, `1` = Line, `2` = Area.
    /// Order: `II IB IE BI BB BE EI EB EE`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for dim in self.cells() {
            let ch = match dim {
                Dimension::Empty => 'F',
                Dimension::Point => '0',
                Dimension::Line => '1',
                Dimension::Area => '2',
            };
            write!(f, "{ch}")?;
        }
        Ok(())
    }
}

/// Error returned by [`IntersectionMatrix::from_pattern`] on invalid input.
#[derive(Debug, PartialEq)]
pub struct PatternError {
    /// A static description of the parse failure.
    pub message: &'static str,
}

impl std::fmt::Display for PatternError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "DE-9IM pattern parse error: {}", self.message)
    }
}

// ===========================================================================
// Windowing a large unindexed operand
// ===========================================================================

/// Each areal operand windowed to the other's envelope, where it is large
/// enough for that to pay (`WindowedAreal`). A windowed operand answers
/// every query exactly as the unindexed one does, so the matrix is the
/// same either way.
pub(crate) fn windowed_pair<'a>(
    a: &Areal<'a>,
    b: &Areal<'a>,
) -> (Option<WindowedAreal<'a>>, Option<WindowedAreal<'a>>) {
    // Each window is the other operand's envelope; it is computed only when
    // the operand it windows is large, so a small operand against a large
    // one costs one pass over the large one's edges.
    let wa = WindowedAreal::worth_it(a)
        .then(|| WindowedAreal::of(a, b.all_window()))
        .flatten();
    let wb = WindowedAreal::worth_it(b)
        .then(|| {
            let w = wa
                .as_ref()
                .map_or_else(|| a.all_window(), |x| x.all_window());
            WindowedAreal::of(b, w)
        })
        .flatten();
    (wa, wb)
}

/// The box of a non-empty linear operand's vertices.
fn lines_window(l: &LineSet<'_>) -> Window {
    l.bbox().map_or(
        (
            crate::types::Coord2d::new(0.0, 0.0),
            crate::types::Coord2d::new(0.0, 0.0),
        ),
        |(x0, y0, x1, y1)| {
            (
                crate::types::Coord2d::new(x0, y0),
                crate::types::Coord2d::new(x1, y1),
            )
        },
    )
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

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    fn seq(pts: &[(f64, f64)]) -> CoordSeq {
        CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())
    }

    fn ccw_square(ox: f64, oy: f64, side: f64) -> Polygon {
        let ring = LinearRing::new(seq(&[
            (ox, oy),
            (ox + side, oy),
            (ox + side, oy + side),
            (ox, oy + side),
        ]))
        .unwrap();
        Polygon::new_simple(ring).unwrap()
    }

    // ---- existing matrix tests ----

    #[test]
    fn disjoint_polygons_interior_interior_is_empty() {
        let a = ccw_square(0.0, 0.0, 1.0);
        let b = ccw_square(2.0, 0.0, 1.0);
        let m = IntersectionMatrix::compute(&a, &b);
        assert_eq!(m.interior_interior, Dimension::Empty);
        assert_eq!(m.exterior_exterior, Dimension::Area);
    }

    #[test]
    fn from_pattern_valid_with_wildcards() {
        let result = IntersectionMatrix::from_pattern("FF2FF1***");
        assert!(result.is_ok());
        let m = result.unwrap();
        assert_eq!(m.interior_interior, Dimension::Empty);
        assert_eq!(m.interior_exterior, Dimension::Area);
        assert_eq!(m.boundary_exterior, Dimension::Line);
    }

    #[test]
    fn from_pattern_wrong_length_is_err() {
        assert!(IntersectionMatrix::from_pattern("FF2FF1**").is_err());
        assert!(IntersectionMatrix::from_pattern("FF2FF1****").is_err());
        assert!(IntersectionMatrix::from_pattern("").is_err());
    }

    #[test]
    fn from_pattern_invalid_char_is_err() {
        let result = IntersectionMatrix::from_pattern("FF2FF1**X");
        assert!(result.is_err());
        assert!(result.unwrap_err().message.contains("invalid character"));
    }

    #[test]
    fn display_renders_ogc_string() {
        let m = IntersectionMatrix {
            interior_interior: Dimension::Area,
            interior_boundary: Dimension::Line,
            interior_exterior: Dimension::Area,
            boundary_interior: Dimension::Line,
            boundary_boundary: Dimension::Point,
            boundary_exterior: Dimension::Line,
            exterior_interior: Dimension::Area,
            exterior_boundary: Dimension::Line,
            exterior_exterior: Dimension::Area,
        };
        assert_eq!(m.to_string(), "212101212");
    }

    #[test]
    fn display_renders_all_dimension_chars() {
        let m = IntersectionMatrix {
            interior_interior: Dimension::Empty,
            interior_boundary: Dimension::Point,
            interior_exterior: Dimension::Line,
            boundary_interior: Dimension::Area,
            boundary_boundary: Dimension::Empty,
            boundary_exterior: Dimension::Point,
            exterior_interior: Dimension::Line,
            exterior_boundary: Dimension::Area,
            exterior_exterior: Dimension::Empty,
        };
        assert_eq!(m.to_string(), "F012F012F");
    }

    // ---- Boundary cell tests ----

    #[test]
    fn compute_d4_edge_touching_produces_bb_line() {
        // Squares sharing a full edge: the shared boundary is a segment, so
        // BB is Line (the earlier test asserted Point, which was wrong).
        let a = ccw_square(0.0, 0.0, 1.0);
        let b = ccw_square(1.0, 0.0, 1.0);
        let m = IntersectionMatrix::compute(&a, &b);
        assert_eq!(m.boundary_boundary, Dimension::Line);
        assert_eq!(m.interior_interior, Dimension::Empty);
        assert_eq!(m.to_string(), "FF2F11212");
    }

    #[test]
    fn compute_d5_vertex_touching_produces_bb_point() {
        let a = ccw_square(0.0, 0.0, 1.0);
        let b = ccw_square(1.0, 1.0, 1.0);
        let m = IntersectionMatrix::compute(&a, &b);
        assert_eq!(m.boundary_boundary, Dimension::Point);
        assert_eq!(m.interior_interior, Dimension::Empty);
    }

    #[test]
    fn compute_crossing_squares_bb_is_point() {
        // Boundaries cross at two points, so BB is Point (the earlier test
        // asserted Empty, which was wrong); each boundary enters the other
        // polygon's interior, so IB and BI are Line.
        let a = ccw_square(0.0, 0.0, 4.0);
        let b = ccw_square(2.0, 2.0, 4.0);
        let m = IntersectionMatrix::compute(&a, &b);
        assert_eq!(m.boundary_boundary, Dimension::Point);
        assert_eq!(m.interior_interior, Dimension::Area);
        assert_eq!(m.to_string(), "212101212");
    }

    // ---- Equal polygon detection ----

    /// Equal squares produce the correct "2FFF1FFF2" matrix.
    #[test]
    fn compute_d3_equal_squares_produces_correct_matrix() {
        let a = ccw_square(0.0, 0.0, 1.0);
        let b = ccw_square(0.0, 0.0, 1.0);
        let m = IntersectionMatrix::compute(&a, &b);

        assert_eq!(m.interior_interior, Dimension::Area, "II");
        assert_eq!(m.interior_boundary, Dimension::Empty, "IB");
        assert_eq!(m.interior_exterior, Dimension::Empty, "IE");
        assert_eq!(m.boundary_interior, Dimension::Empty, "BI");
        assert_eq!(m.boundary_boundary, Dimension::Line, "BB: shared boundary");
        assert_eq!(m.boundary_exterior, Dimension::Empty, "BE");
        assert_eq!(m.exterior_interior, Dimension::Empty, "EI");
        assert_eq!(m.exterior_boundary, Dimension::Empty, "EB");
        assert_eq!(m.exterior_exterior, Dimension::Area, "EE");
        assert_eq!(m.to_string(), "2FFF1FFF2", "OGC string for equal polygons");
    }

    // ---- Nested polygon integration tests ----

    /// B = [1,3]² nested inside A = [0,4]².
    /// compute(A, B) should produce the "B inside A" matrix.
    #[test]
    fn compute_d1_b_inside_a_produces_correct_matrix() {
        let a = ccw_square(0.0, 0.0, 4.0);
        let b = ccw_square(1.0, 1.0, 2.0);
        let m = IntersectionMatrix::compute(&a, &b);

        assert_eq!(m.interior_interior, Dimension::Area, "II");
        assert_eq!(
            m.interior_boundary,
            Dimension::Line,
            "IB: B's ring inside A"
        );
        assert_eq!(m.interior_exterior, Dimension::Area, "IE: annular region");
        assert_eq!(m.boundary_interior, Dimension::Empty, "BI");
        assert_eq!(m.boundary_boundary, Dimension::Empty, "BB");
        assert_eq!(
            m.boundary_exterior,
            Dimension::Line,
            "BE: A's ring outside B"
        );
        assert_eq!(m.exterior_interior, Dimension::Empty, "EI");
        assert_eq!(m.exterior_boundary, Dimension::Empty, "EB");
        assert_eq!(m.exterior_exterior, Dimension::Area, "EE");
    }

    /// A = [1,3]² nested inside B = [0,4]².
    /// compute(A, B) should produce the "A inside B" matrix.
    #[test]
    fn compute_d2_a_inside_b_produces_correct_matrix() {
        let a = ccw_square(1.0, 1.0, 2.0);
        let b = ccw_square(0.0, 0.0, 4.0);
        let m = IntersectionMatrix::compute(&a, &b);

        assert_eq!(m.interior_interior, Dimension::Area, "II");
        assert_eq!(m.interior_boundary, Dimension::Empty, "IB");
        assert_eq!(m.interior_exterior, Dimension::Empty, "IE");
        assert_eq!(
            m.boundary_interior,
            Dimension::Line,
            "BI: A's ring inside B"
        );
        assert_eq!(m.boundary_boundary, Dimension::Empty, "BB");
        assert_eq!(m.boundary_exterior, Dimension::Empty, "BE");
        assert_eq!(m.exterior_interior, Dimension::Area, "EI: annular region");
        assert_eq!(
            m.exterior_boundary,
            Dimension::Line,
            "EB: B's ring in A's exterior"
        );
        assert_eq!(m.exterior_exterior, Dimension::Area, "EE");
    }
}
