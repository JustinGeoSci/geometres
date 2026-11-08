// tests/predicates.rs

//! Integration tests for `geometres::predicates`.
//!
//! Pattern-based tests construct `IntersectionMatrix` directly via
//! `from_pattern`, with no topology computation, to isolate predicate logic.
//! Each test group also includes at least one geometry-level test that calls
//! the generic predicate function with real `Polygon` values.
//!
//! # Position reference (row-major OGC order)
//!
//! ```text
//! index: 0  1  2  3  4  5  6  7  8
//! cell:  II IB IE BI BB BE EI EB EE
//! ```
//!
//! Characters: `F` = Empty, `0` = Point, `1` = Line, `2` = Area, `*` = wildcard (stored as Empty).

// Test-code lint posture, stated rather than left implicit:
//
//   float_cmp:           an exact-value assertion is a float comparison. These
//                        tests assert bit-exact results against a rational
//                        oracle or a hand-computed value; a tolerance would
//                        defeat their purpose.
//   unreadable_literal:  bit patterns and reference coordinates are
//                        transcribed as written in the source they come from,
//                        so that they can be checked against it.
//   cast_*:              fixture construction, where the ranges are known by
//                        construction.
//   unwrap_used:         in a test, a panic is the failure signal.
#![allow(
    clippy::float_cmp,
    clippy::unreadable_literal,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_lossless,
    clippy::unwrap_used,
    clippy::type_complexity,
    clippy::no_effect_underscore_binding,
    clippy::if_not_else,
    clippy::comparison_chain
)]

#[allow(clippy::unwrap_used)]
mod helpers {
    use geometres::types::{Coord2d, CoordSeq, LinearRing, Polygon};

    pub fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    pub fn seq(pts: &[(f64, f64)]) -> CoordSeq {
        CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())
    }

    /// Axis-aligned CCW square with lower-left corner `(ox, oy)` and given `side`.
    pub fn ccw_square(ox: f64, oy: f64, side: f64) -> Polygon {
        let ring = LinearRing::new(seq(&[
            (ox, oy),
            (ox + side, oy),
            (ox + side, oy + side),
            (ox, oy + side),
        ]))
        .unwrap();
        Polygon::new_simple(ring).unwrap()
    }
}

// ===========================================================================
// equals
// ===========================================================================

mod equals_tests {
    use geometres::predicates::equals;
    use geometres::topology::de9im::IntersectionMatrix;

    /// Identical-interior matrix: II=Area, IE=Empty, EI=Empty → equals is true.
    #[test]
    fn equals_true_shared_interior_no_exterior() {
        // II=2, IB=F, IE=F, BI=F, BB=F, BE=F, EI=F, EB=F, EE=2
        let m = IntersectionMatrix::from_pattern("2FFFFFFF2").unwrap();
        assert!(
            m.equals(),
            "equals should be true when II≠Empty, IE=Empty, EI=Empty"
        );
    }

    /// IE=Empty, EI=Empty but II=Point (e.g. touching points); still equals.
    #[test]
    fn equals_true_point_intersection_no_exterior() {
        // II=0, IE=F, EI=F
        let m = IntersectionMatrix::from_pattern("0FFFF*F*2").unwrap();
        assert!(m.equals());
    }

    /// IE≠Empty means A has parts outside B → not equal.
    #[test]
    fn equals_false_a_has_exterior() {
        // II=2, IE=2 → A has interior outside B
        let m = IntersectionMatrix::from_pattern("2F2FF*FF2").unwrap();
        assert!(!m.equals(), "equals false when IE≠Empty");
    }

    /// EI≠Empty means B has parts outside A → not equal.
    #[test]
    fn equals_false_b_has_exterior() {
        // II=2, IB=F, IE=F, BI=F, BB=F, BE=*, EI=2, EB=F, EE=2
        // Positions: 0=II 1=IB 2=IE 3=BI 4=BB 5=BE 6=EI 7=EB 8=EE
        let m = IntersectionMatrix::from_pattern("2FFFF*2F2").unwrap();
        assert!(!m.equals(), "equals false when EI≠Empty");
    }

    /// Two disjoint squares: II=Empty, IE=Area, EI=Area → not equal.
    ///
    /// The identical-polygon case uses coincident-edge detection; it is covered
    /// at the pattern level by the `equals_true_*` tests above.
    #[test]
    fn polygon_equals_disjoint_squares_not_equal() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(5.0, 0.0, 1.0);
        assert!(!equals(&a, &b).unwrap());
    }

    /// Identical squares are geometrically equal.
    #[test]
    fn polygon_equals_d3_identical_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(0.0, 0.0, 1.0);
        assert!(
            equals(&a, &b).unwrap(),
            "equals must return true for geometrically identical polygons"
        );
    }
}

// ===========================================================================
// disjoint
// ===========================================================================

mod disjoint_tests {
    use geometres::predicates::disjoint;
    use geometres::topology::de9im::IntersectionMatrix;

    /// All four checked cells (II, IB, BI, BB) empty → disjoint.
    #[test]
    fn disjoint_true_all_empty() {
        // II=F IB=F IE=2 BI=F BB=F BE=2 EI=2 EB=2 EE=2
        let m = IntersectionMatrix::from_pattern("FF2FF2222").unwrap();
        assert!(m.disjoint());
    }

    /// Classic disjoint pattern: no contact at all.
    #[test]
    fn disjoint_true_ff2ff2t22() {
        let m = IntersectionMatrix::from_pattern("FF2FF2*22").unwrap();
        assert!(m.disjoint());
    }

    /// II≠Empty → not disjoint.
    #[test]
    fn disjoint_false_interiors_intersect() {
        let m = IntersectionMatrix::from_pattern("2F2FF*FF2").unwrap();
        assert!(!m.disjoint());
    }

    /// BB≠Empty (boundary contact) → not disjoint.
    #[test]
    fn disjoint_false_boundary_contact() {
        // II=F IB=F IE=2 BI=F BB=1 BE=1 EI=2 EB=1 EE=2
        let m = IntersectionMatrix::from_pattern("FF2F11212").unwrap();
        assert!(!m.disjoint());
    }

    /// Two unit squares placed far apart.
    #[test]
    fn polygon_disjoint_far_apart() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(5.0, 0.0, 1.0);
        assert!(disjoint(&a, &b).unwrap());
    }

    /// Edge-touching squares: BB = Point → not disjoint.
    #[test]
    fn polygon_disjoint_d4_edge_touching_is_false() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(1.0, 0.0, 1.0);
        assert!(
            !disjoint(&a, &b).unwrap(),
            "disjoint must return false for edge-touching squares \
             (BB=Point from shared vertices)"
        );
    }

    /// Vertex-touching squares: BB = Point → not disjoint.
    #[test]
    fn polygon_disjoint_d5_vertex_touching_is_false() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(1.0, 1.0, 1.0);
        assert!(
            !disjoint(&a, &b).unwrap(),
            "disjoint must return false for vertex-touching squares \
             (BB=Point from shared vertex (1,1))"
        );
    }
}

// ===========================================================================
// intersects
// ===========================================================================

mod intersects_tests {
    use geometres::predicates::intersects;
    use geometres::topology::de9im::IntersectionMatrix;

    /// Interior overlap → intersects.
    #[test]
    fn intersects_true_interiors_overlap() {
        let m = IntersectionMatrix::from_pattern("2F2FF*FF2").unwrap();
        assert!(m.intersects());
    }

    /// Boundary-only contact → intersects (complement of disjoint).
    #[test]
    fn intersects_true_boundary_contact() {
        let m = IntersectionMatrix::from_pattern("FF2F11212").unwrap();
        assert!(m.intersects());
    }

    /// Fully disjoint → not intersects.
    #[test]
    fn intersects_false_fully_disjoint() {
        let m = IntersectionMatrix::from_pattern("FF2FF2222").unwrap();
        assert!(!m.intersects());
    }

    /// Another fully-disjoint matrix.
    #[test]
    fn intersects_false_no_contact() {
        let m = IntersectionMatrix::from_pattern("FF2FF2*22").unwrap();
        assert!(!m.intersects());
    }

    /// Two overlapping squares: they intersect.
    #[test]
    fn polygon_intersects_overlapping() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(2.0, 2.0, 4.0);
        assert!(intersects(&a, &b).unwrap());
    }

    /// Edge-touching squares: BB = Point → intersects (complement of disjoint).
    #[test]
    fn polygon_intersects_d4_edge_touching_is_true() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(1.0, 0.0, 1.0);
        assert!(
            intersects(&a, &b).unwrap(),
            "intersects must return true for edge-touching squares"
        );
    }

    /// Vertex-touching squares: BB = Point → intersects.
    #[test]
    fn polygon_intersects_d5_vertex_touching_is_true() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(1.0, 1.0, 1.0);
        assert!(
            intersects(&a, &b).unwrap(),
            "intersects must return true for vertex-touching squares"
        );
    }
}

// ===========================================================================
// touches
// ===========================================================================

mod touches_tests {
    use geometres::predicates::touches;
    use geometres::topology::de9im::IntersectionMatrix;

    /// Polygon pair sharing a boundary edge: II=F, BB=Line → touches.
    #[test]
    fn touches_true_boundary_edge_contact() {
        // II=F IB=F IE=2 BI=F BB=1 BE=1 EI=2 EB=1 EE=2
        let m = IntersectionMatrix::from_pattern("FF2F11212").unwrap();
        assert!(m.touches(), "boundary edge contact: II=Empty, BB=Line");
    }

    /// Polygon pair sharing only a vertex: II=F, BB=Point → touches.
    #[test]
    fn touches_true_vertex_only_contact() {
        // II=F IB=F IE=2 BI=F BB=0 BE=2 EI=2 EB=2 EE=2
        let m = IntersectionMatrix::from_pattern("FF2F02222").unwrap();
        assert!(m.touches(), "vertex contact: II=Empty, BB=Point");
    }

    /// Interior overlap: II=Area → touches is false.
    #[test]
    fn touches_false_interior_overlap() {
        let m = IntersectionMatrix::from_pattern("212111212").unwrap();
        assert!(!m.touches(), "interior overlap: touches is false");
    }

    /// Fully disjoint: II=Empty but IB/BI/BB all Empty → touches is false.
    #[test]
    fn touches_false_fully_disjoint() {
        let m = IntersectionMatrix::from_pattern("FF2FF2222").unwrap();
        assert!(
            !m.touches(),
            "disjoint: no boundary contact means touches is false"
        );
    }

    /// Squares sharing the edge x = 1.
    ///
    /// Boundary contact analysis detects the two shared vertices (1,0) and
    /// (1,1) as vertex–vertex coincidences → BB = Point.
    /// `touches()` requires II = Empty (no interior overlap) and at least
    /// one of IB, BI, BB to be non-empty. With BB = Point the predicate
    /// returns true.
    #[test]
    fn polygon_touches_d4_edge_sharing() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(1.0, 0.0, 1.0);
        assert!(
            touches(&a, &b).unwrap(),
            "touches must return true for side-by-side squares \
             sharing the edge x=1 (BB=Point from shared vertices)"
        );
    }

    /// Squares touching at exactly one vertex (1,1).
    ///
    /// Boundary contact analysis detects the shared vertex (1,1) → BB = Point.
    #[test]
    fn polygon_touches_d5_vertex_sharing() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(1.0, 1.0, 1.0);
        assert!(
            touches(&a, &b).unwrap(),
            "touches must return true for squares sharing only \
             the vertex (1,1) (BB=Point from shared vertex)"
        );
    }
}

// ===========================================================================
// crosses, line × line
// ===========================================================================

mod crosses_line_line_tests {
    use geometres::topology::de9im::IntersectionMatrix;
    use geometres::topology::label::Dimension;

    /// Classic X-crossing: II=Point, IE≠Empty, EI≠Empty.
    #[test]
    fn crosses_line_line_true_x_crossing() {
        // II=0 IB=F IE=1 BI=F BB=F BE=1 EI=1 EB=F EE=2
        let m = IntersectionMatrix::from_pattern("0F1FF1102").unwrap();
        assert!(m.crosses(Dimension::Line, Dimension::Line));
    }

    /// II=Point with both IE and EI non-empty: another true case.
    #[test]
    fn crosses_line_line_true_general() {
        // II=0 IE=2 EI=2
        let m = IntersectionMatrix::from_pattern("0F2FF*2F2").unwrap();
        assert!(m.crosses(Dimension::Line, Dimension::Line));
    }

    /// II=Line (collinear overlap) → not a crossing (II must be Point).
    #[test]
    fn crosses_line_line_false_collinear_overlap() {
        // II=1 (line segment overlap, not a point crossing)
        let m = IntersectionMatrix::from_pattern("1F2FF*2F2").unwrap();
        assert!(!m.crosses(Dimension::Line, Dimension::Line));
    }

    /// The OGC line × line pattern for crosses is `0********`, so `IE` and
    /// `EI` are not part of it. A reading that also requires both non-empty
    /// agrees with it on every matrix two lines can produce, because
    /// interiors meeting only at a point leave the rest of each line's
    /// interior outside the other. This matrix (`II = 0`, `EI = F`) is one no
    /// two lines produce, and the full definition reads it as a crossing.
    #[test]
    fn crosses_line_line_reads_only_ii() {
        // II=0 IE=2 EI=F
        let m = IntersectionMatrix::from_pattern("0F2FFFFF2").unwrap();
        assert!(m.crosses(Dimension::Line, Dimension::Line));
    }
}

// ===========================================================================
// crosses, line × area
// ===========================================================================

mod crosses_line_polygon_tests {
    use geometres::topology::de9im::IntersectionMatrix;
    use geometres::topology::label::Dimension;

    /// Line enters and exits polygon: II≠Empty, IE≠Empty.
    #[test]
    fn crosses_line_polygon_true_enters_and_exits() {
        // II=1 (line inside polygon), IE=1 (line outside polygon)
        let m = IntersectionMatrix::from_pattern("1F2FF*2F2").unwrap();
        assert!(m.crosses(Dimension::Line, Dimension::Area));
    }

    /// II=Area (treated as non-empty), IE≠Empty → true.
    #[test]
    fn crosses_line_polygon_true_area_interior() {
        let m = IntersectionMatrix::from_pattern("2F2FF*FF2").unwrap();
        assert!(m.crosses(Dimension::Line, Dimension::Area));
    }

    /// Line entirely inside polygon: IE=Empty → not crossing.
    #[test]
    fn crosses_line_polygon_false_line_within_polygon() {
        // II=1, IE=F
        let m = IntersectionMatrix::from_pattern("1FFFF*FF2").unwrap();
        assert!(!m.crosses(Dimension::Line, Dimension::Area));
    }

    /// Line entirely outside polygon: II=Empty → not crossing.
    #[test]
    fn crosses_line_polygon_false_line_disjoint() {
        let m = IntersectionMatrix::from_pattern("FF2FF2*22").unwrap();
        assert!(!m.crosses(Dimension::Line, Dimension::Area));
    }
}

// ===========================================================================
// within
// ===========================================================================

mod within_tests {
    use geometres::predicates::within;
    use geometres::topology::de9im::IntersectionMatrix;

    /// A fully inside B: II≠Empty, IE=Empty, BE=Empty.
    #[test]
    fn within_true_fully_contained() {
        // II=2, IB=F, IE=F, BI=F, BB=F, BE=F, EI=2, EB=2, EE=2
        let m = IntersectionMatrix::from_pattern("2FFFFFFF2").unwrap();
        // IE=pos[2]=F, BE=pos[5]=F: within is true
        assert!(m.within());
    }

    /// II=Point, IE=Empty, BE=Empty: still within (e.g. a point inside a polygon).
    #[test]
    fn within_true_point_inside() {
        // II=0, IE=F, BE=F
        let m = IntersectionMatrix::from_pattern("0FFFF*FF2").unwrap();
        assert!(m.within());
    }

    /// IE≠Empty → A has parts outside B → not within.
    #[test]
    fn within_false_a_has_exterior() {
        let m = IntersectionMatrix::from_pattern("2F2FF*FF2").unwrap();
        assert!(!m.within());
    }

    /// BE≠Empty → A's boundary crosses outside B → not within.
    #[test]
    fn within_false_boundary_outside() {
        // IE=F but BE=1
        let m = IntersectionMatrix::from_pattern("2FFFF1FF2").unwrap();
        assert!(!m.within());
    }

    /// Two partially-overlapping squares: A=[0,4]² has interior outside B=[2,2]×[6,6],
    /// so A is not within B.
    #[test]
    fn polygon_within_false_overlapping_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(2.0, 2.0, 4.0);
        // A has area outside B (IE=Area), so within(A, B) is false.
        assert!(!within(&a, &b).unwrap());
    }

    /// B = [1,3]² is inside A = [0,4]². A is not within B.
    #[test]
    fn polygon_within_d1_outer_a_not_within_inner_b() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(1.0, 1.0, 2.0);
        // A is the outer polygon; it is not within B (the inner polygon).
        assert!(
            !within(&a, &b).unwrap(),
            "the outer polygon A is not within the inner polygon B"
        );
    }

    /// A = [1,3]² is inside B = [0,4]². A is within B.
    #[test]
    fn polygon_within_d2_a_within_b() {
        let a = super::helpers::ccw_square(1.0, 1.0, 2.0);
        let b = super::helpers::ccw_square(0.0, 0.0, 4.0);
        assert!(
            within(&a, &b).unwrap(),
            "within must return true when A is entirely inside B"
        );
    }

    /// A within B when A = B.
    #[test]
    fn polygon_within_d3_equal_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(0.0, 0.0, 1.0);
        assert!(within(&a, &b).unwrap(), "equal polygons satisfy within");
    }
}

// ===========================================================================
// contains
// ===========================================================================

mod contains_tests {
    use geometres::predicates::contains;
    use geometres::topology::de9im::IntersectionMatrix;

    /// A contains B: II≠Empty, EI=Empty, EB=Empty.
    #[test]
    fn contains_true_b_fully_inside() {
        // II=2, EI=F, EB=F: positions [0]=2, [6]=F, [7]=F
        let m = IntersectionMatrix::from_pattern("2FFFFFF*2").unwrap();
        // EI=pos[6]=F, EB=pos[7]=wildcard→Empty: contains is true
        assert!(m.contains());
    }

    /// Another contains case, with explicit F for EI and EB.
    #[test]
    fn contains_true_explicit_cells() {
        let m = IntersectionMatrix::from_pattern("2FFFFFFF2").unwrap();
        assert!(m.contains());
    }

    /// EI≠Empty → part of B's interior lies outside A → not contains.
    #[test]
    fn contains_false_ei_nonempty() {
        // II=2, IB=F, IE=F, BI=F, BB=F, BE=F, EI=2, EB=F, EE=2
        // Positions: 0=II 1=IB 2=IE 3=BI 4=BB 5=BE 6=EI 7=EB 8=EE
        let m = IntersectionMatrix::from_pattern("2FFFFF2F2").unwrap();
        assert!(!m.contains());
    }

    /// EB≠Empty → part of B's boundary lies outside A → not contains.
    #[test]
    fn contains_false_eb_nonempty() {
        // EB=1
        let m = IntersectionMatrix::from_pattern("2FFFFF1F2").unwrap();
        assert!(!m.contains());
    }

    /// Two partially-overlapping squares: each has area outside the other,
    /// so neither contains the other.
    #[test]
    fn polygon_contains_false_overlapping_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(2.0, 2.0, 4.0);
        // B has area outside A (EI=Area), so contains(A, B) is false.
        assert!(!contains(&a, &b).unwrap());
    }

    /// A = [0,4]² contains B = [1,3]².
    #[test]
    fn polygon_contains_d1_a_contains_b() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(1.0, 1.0, 2.0);
        assert!(
            contains(&a, &b).unwrap(),
            "contains must return true when B is entirely inside A"
        );
    }

    /// B = [0,4]² contains A = [1,3]². A does not contain B.
    #[test]
    fn polygon_contains_d2_b_contains_a_not_vice_versa() {
        let a = super::helpers::ccw_square(1.0, 1.0, 2.0);
        let b = super::helpers::ccw_square(0.0, 0.0, 4.0);
        assert!(
            !contains(&a, &b).unwrap(),
            "the inner polygon A does not contain the outer polygon B"
        );
    }

    /// A contains B when A = B.
    #[test]
    fn polygon_contains_d3_equal_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(0.0, 0.0, 1.0);
        assert!(contains(&a, &b).unwrap(), "equal polygons satisfy contains");
    }
}

// ===========================================================================
// overlaps, area × area
// ===========================================================================

mod overlaps_tests {
    use geometres::predicates::overlaps;
    use geometres::topology::de9im::IntersectionMatrix;
    use geometres::topology::label::Dimension;

    /// Partial overlap: II=Area, IE≠Empty, EI≠Empty.
    #[test]
    fn overlaps_true_partial_overlap() {
        // II=2, IE=2, EI=2: each polygon has area inside and outside the other
        let m = IntersectionMatrix::from_pattern("2F2FF*2F2").unwrap();
        assert!(m.overlaps(Dimension::Area, Dimension::Area));
    }

    /// Another partial overlap case with explicit boundary cells set.
    #[test]
    fn overlaps_true_general_overlap() {
        let m = IntersectionMatrix::from_pattern("212111212").unwrap();
        assert!(m.overlaps(Dimension::Area, Dimension::Area));
    }

    /// A contained in B: IE=Empty → neither has area outside the other (from A's side).
    #[test]
    fn overlaps_false_a_within_b() {
        // II=2, IE=F → A is fully inside B
        let m = IntersectionMatrix::from_pattern("2FFFFFFF2").unwrap();
        assert!(!m.overlaps(Dimension::Area, Dimension::Area));
    }

    /// Disjoint: II=Empty → no shared area.
    #[test]
    fn overlaps_false_disjoint() {
        let m = IntersectionMatrix::from_pattern("FF2FF2222").unwrap();
        assert!(!m.overlaps(Dimension::Area, Dimension::Area));
    }

    /// Two partially-overlapping squares: A=[0,4]², B=[2,2]×[6,6].
    #[test]
    fn polygon_overlaps_crossing_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(2.0, 2.0, 4.0);
        assert!(overlaps(&a, &b).unwrap());
    }
}

// ===========================================================================
// covers
// ===========================================================================

mod covers_tests {
    use geometres::predicates::covers;
    use geometres::topology::de9im::IntersectionMatrix;

    /// A covers B: IE=Empty, BE=Empty: B is fully inside A (including boundary).
    #[test]
    fn covers_true_b_fully_inside() {
        // IE=pos[2]=F, BE=pos[5]=F
        let m = IntersectionMatrix::from_pattern("2FFFFFFF2").unwrap();
        assert!(m.covers());
    }

    /// EI=Empty, EB=Empty with BB contact (boundary sharing allowed).
    #[test]
    fn covers_true_with_boundary_sharing() {
        // IE=F, BE=F: covers allows boundary sharing (unlike strict contains)
        let m = IntersectionMatrix::from_pattern("2FF2FF*F2").unwrap();
        // IE=pos[2]=F, BE=pos[5]=F → covers is true
        assert!(m.covers());
    }

    /// EI≠Empty: B's interior lies in A's exterior → A does not cover B.
    ///
    /// This is the geometrically correct "not covers" test: B has area outside
    /// A (EI = Area at position 6), so `covers(A, B)` = false.
    #[test]
    fn covers_false_b_interior_outside_a() {
        // pos[6] = EI = 2 (B's interior partly outside A)
        let m = IntersectionMatrix::from_pattern("2F2FF*2F2").unwrap();
        assert!(!m.covers());
    }

    /// EB≠Empty: B's boundary lies in A's exterior → A does not cover B.
    ///
    /// `covers(A, B)` checks EI = Empty AND EB = Empty.
    /// When EB = Line (position 7), B's boundary has a portion outside A,
    /// so `covers` returns false.
    #[test]
    fn covers_false_b_boundary_outside_a() {
        // pos[7] = EB = 1 (B's boundary partly in A's exterior)
        let m = IntersectionMatrix::from_pattern("2FFFFFF12").unwrap();
        assert!(!m.covers());
    }

    /// Two partially-overlapping squares: each has area outside the other.
    ///
    /// `covers(A, B)` checks EI = Empty AND EB = Empty. For overlapping
    /// squares B = [2,6]² has interior outside A = [0,4]² (EI = Area),
    /// so covers(A, B) = false.
    #[test]
    fn polygon_covers_false_overlapping_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(2.0, 2.0, 4.0);
        // B's interior lies in A's exterior (EI=Area), so covers(A, B) is false.
        assert!(!covers(&a, &b).unwrap());
    }

    /// A = [0,4]² covers B = [1,3]².
    ///
    /// No point of B lies in A's exterior: EI = Empty, EB = Empty.
    #[test]
    fn polygon_covers_d1_a_covers_b() {
        let a = super::helpers::ccw_square(0.0, 0.0, 4.0);
        let b = super::helpers::ccw_square(1.0, 1.0, 2.0);
        assert!(
            covers(&a, &b).unwrap(),
            "covers must return true when B is entirely inside A"
        );
    }

    /// A = [1,3]² is inside B = [0,4]². B does not cover A in this call.
    #[test]
    fn polygon_covers_d2_inner_does_not_cover_outer() {
        let a = super::helpers::ccw_square(1.0, 1.0, 2.0);
        let b = super::helpers::ccw_square(0.0, 0.0, 4.0);
        // covers(a, b).unwrap() = "does A cover B?" A is inside B, so no.
        assert!(
            !covers(&a, &b).unwrap(),
            "the inner polygon A does not cover the outer polygon B"
        );
    }

    /// A covers B when A = B.
    #[test]
    fn polygon_covers_d3_equal_squares() {
        let a = super::helpers::ccw_square(0.0, 0.0, 1.0);
        let b = super::helpers::ccw_square(0.0, 0.0, 1.0);
        assert!(covers(&a, &b).unwrap(), "equal polygons satisfy covers");
    }
}
