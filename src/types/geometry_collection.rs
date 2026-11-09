// src/types/geometry_collection.rs

//! Heterogeneous collection geometry (`GeometryCollection`).
//!
//! # Nesting, and why it is limited
//!
//! A collection may hold other collections, as in OGC Simple Features and
//! WKB. Its *depth* counts collection levels: a collection whose
//! elements are all non-collections has depth 1, and one holding such a
//! collection has depth 2. The empty collection has depth 1.
//!
//! **The depth is at most [`GeometryCollection::MAX_NESTING`], 32. This is a
//! permanent contract of the type and of the native format.**
//!
//! - **What the limit protects:** The derived `Clone`, `PartialEq` and
//!   `Debug`, and the compiler's drop glue for `Vec<Geometry>`, each recurse
//!   once per level, and so does any caller's code that matches its way down.
//!   With the depth bounded, the worst-case recursion anywhere is a few dozen
//!   small frames, on every thread stack and every target. This crate's own
//!   algorithms do not recurse at all: they walk a collection with
//!   [`GeometryCollection::leaves`] or an explicit stack.
//! - **Why it is checked at construction:** Every value that exists then also
//!   encodes and decodes, so the codec still round-trips. The decoder checks
//!   the same limit before it enters a level, so a hostile byte stream is
//!   refused before it can nest any deeper.
//! - **Why it never changes:** Lowering it would make existing values
//!   unconstructible and existing documents undecodable; raising it would
//!   make older decoders refuse the documents a newer encoder writes.
//!
//! # No condition between elements
//!
//! OGC imposes none: elements may overlap, repeat or touch. Every element is
//! valid by its own construction, so the only check here is the depth. The
//! measures over a collection (`ops::collection_area`,
//! `ops::collection_length`) are sums over its elements, so overlapping
//! elements count twice.

use crate::GeometryError;
use crate::topology::label::Dimension;
use crate::types::{Coord2d, Envelope, Geometry};

/// A collection of geometries of any types, including other collections,
/// nested at most [`MAX_NESTING`](Self::MAX_NESTING) levels deep.
///
/// # Invariants
///
/// - Every element is a valid geometry of its own type (its constructor
///   checked it), so every coordinate is finite.
/// - The depth (see the documentation of src/types/geometry_collection.rs)
///   is at most [`MAX_NESTING`](Self::MAX_NESTING).
#[derive(Clone, Debug, PartialEq)]
pub struct GeometryCollection {
    geometries: Vec<Geometry>,
    /// The nesting depth, cached so that building a nest bottom-up costs
    /// `O(n)` rather than a walk of every child on each construction. A
    /// function of `geometries`, so the derived `PartialEq` stays correct.
    depth: u8,
}

impl GeometryCollection {
    /// The deepest nesting a collection may have: a permanent contract (see the
    /// documentation of src/types/geometry_collection.rs)
    pub const MAX_NESTING: usize = 32;

    /// Construct a collection from its elements, in order.
    ///
    /// # Errors
    ///
    /// [`GeometryError::NestingTooDeep`] if the result would be more than
    /// [`MAX_NESTING`](Self::MAX_NESTING) collection levels deep, that is, if
    /// an element is a collection of depth `MAX_NESTING`.
    pub fn new(geometries: Vec<Geometry>) -> Result<Self, GeometryError> {
        let deepest_child = geometries
            .iter()
            .map(|g| match g {
                Geometry::GeometryCollection(c) => usize::from(c.depth),
                _ => 0,
            })
            .max()
            .unwrap_or(0);
        let depth = deepest_child + 1;
        if depth > Self::MAX_NESTING {
            return Err(GeometryError::NestingTooDeep {
                limit: Self::MAX_NESTING,
            });
        }
        // `depth <= MAX_NESTING = 32`, so it fits; the fallback is unreachable.
        let depth = u8::try_from(depth).unwrap_or(u8::MAX);
        Ok(Self { geometries, depth })
    }

    /// A collection whose depth is already known within the limit: a
    /// validated collection's elements, each replaced by one of the same
    /// nesting depth and permuted. Crate-private.
    pub(crate) fn from_valid(geometries: Vec<Geometry>, depth: usize) -> Self {
        debug_assert!((1..=Self::MAX_NESTING).contains(&depth));
        Self {
            geometries,
            depth: u8::try_from(depth).unwrap_or(u8::MAX),
        }
    }

    /// The collection with no elements, of depth 1.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            geometries: Vec::new(),
            depth: 1,
        }
    }

    /// The elements, in construction order.
    #[must_use]
    pub fn geometries(&self) -> &[Geometry] {
        &self.geometries
    }

    /// The elements, by value.
    #[must_use]
    pub fn into_geometries(self) -> Vec<Geometry> {
        self.geometries
    }

    /// The number of direct elements (a nested collection counts as one).
    #[must_use]
    pub fn len(&self) -> usize {
        self.geometries.len()
    }

    /// The nesting depth: 1 for a collection of non-collections (or of
    /// nothing), one more for each level of collections inside it.
    #[must_use]
    pub fn depth(&self) -> usize {
        usize::from(self.depth)
    }

    /// `true` if the collection has no points at all: no elements, or only
    /// empty ones (an empty multipart geometry, or an empty collection), at
    /// any depth. OGC's `IsEmpty`.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.leaves().all(geometry_is_empty)
    }

    /// The highest dimension of any non-empty element at any depth, or
    /// [`Dimension::Empty`] if the collection [is empty](Self::is_empty).
    #[must_use]
    pub fn dimension(&self) -> Dimension {
        self.leaves()
            .filter(|g| !geometry_is_empty(g))
            .map(leaf_dimension)
            .max()
            .unwrap_or(Dimension::Empty)
    }

    /// The closed box of every element's coordinates, or `None` if the
    /// collection is empty.
    ///
    /// Each bound is an input coordinate, selected by exact comparison
    /// ([`Envelope::from_points`]); nothing is computed. A polygon
    /// contributes its exterior's vertices only, since its holes lie inside
    /// the exterior by construction.
    #[must_use]
    pub fn envelope(&self) -> Option<Envelope> {
        // Every coordinate is finite by the element invariants, so
        // `from_points` cannot return its non-finite error.
        Envelope::from_points(self.leaves().flat_map(leaf_coords)).unwrap_or(None)
    }

    /// The non-collection elements, depth first in element order, at every
    /// depth. Nested collections are entered, not yielded.
    ///
    /// Walks an explicit stack of at most [`MAX_NESTING`](Self::MAX_NESTING)
    /// slice iterators; it never recurses.
    #[must_use]
    pub fn leaves(&self) -> Leaves<'_> {
        Leaves {
            stack: vec![self.geometries.iter()],
        }
    }
}

/// The iterator of [`GeometryCollection::leaves`].
#[derive(Clone, Debug)]
pub struct Leaves<'a> {
    stack: Vec<core::slice::Iter<'a, Geometry>>,
}

impl<'a> Iterator for Leaves<'a> {
    type Item = &'a Geometry;

    fn next(&mut self) -> Option<&'a Geometry> {
        loop {
            let top = self.stack.last_mut()?;
            match top.next() {
                None => {
                    self.stack.pop();
                }
                Some(Geometry::GeometryCollection(c)) => self.stack.push(c.geometries.iter()),
                Some(g) => return Some(g),
            }
        }
    }
}

/// `true` for an empty multipart geometry or an empty collection.
fn geometry_is_empty(g: &Geometry) -> bool {
    match g {
        Geometry::Point(_) | Geometry::LineString(_) | Geometry::Polygon(_) => false,
        Geometry::MultiPoint(m) => m.is_empty(),
        Geometry::MultiLineString(m) => m.is_empty(),
        Geometry::MultiPolygon(m) => m.is_empty(),
        Geometry::GeometryCollection(c) => c.is_empty(),
    }
}

/// The dimension of a non-collection, non-empty element.
fn leaf_dimension(g: &Geometry) -> Dimension {
    match g {
        Geometry::Point(_) | Geometry::MultiPoint(_) => Dimension::Point,
        Geometry::LineString(_) | Geometry::MultiLineString(_) => Dimension::Line,
        Geometry::Polygon(_) | Geometry::MultiPolygon(_) => Dimension::Area,
        Geometry::GeometryCollection(c) => c.dimension(),
    }
}

/// Every coordinate of a non-collection element that bounds it: all of a
/// point's or line's, a polygon's exterior.
fn leaf_coords(g: &Geometry) -> Box<dyn Iterator<Item = Coord2d> + '_> {
    match g {
        Geometry::Point(p) => Box::new(core::iter::once(p.coord())),
        Geometry::MultiPoint(m) => Box::new(m.points().iter().map(crate::types::Point::coord)),
        Geometry::LineString(l) => Box::new(l.coords().as_slice().iter().copied()),
        Geometry::MultiLineString(m) => Box::new(
            m.line_strings()
                .iter()
                .flat_map(|l| l.coords().as_slice().iter().copied()),
        ),
        Geometry::Polygon(p) => Box::new(p.exterior().coords().as_slice().iter().copied()),
        Geometry::MultiPolygon(m) => Box::new(
            m.polygons()
                .iter()
                .flat_map(|p| p.exterior().coords().as_slice().iter().copied()),
        ),
        // `leaves` never yields a collection; an empty iterator keeps this
        // total rather than panicking.
        Geometry::GeometryCollection(_) => Box::new(core::iter::empty()),
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
    use crate::types::{
        CoordSeq, LineString, LinearRing, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
    };

    fn pt(x: f64, y: f64) -> Geometry {
        Geometry::Point(Point::new(Coord2d::new(x, y)).unwrap())
    }

    fn line(pts: &[(f64, f64)]) -> Geometry {
        Geometry::LineString(
            LineString::new(CoordSeq::new(
                pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
            ))
            .unwrap(),
        )
    }

    fn square(x0: f64, y0: f64, s: f64) -> Polygon {
        let pts = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s)];
        Polygon::new_simple(
            LinearRing::new(CoordSeq::new(
                pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
            ))
            .unwrap(),
        )
        .unwrap()
    }

    /// A collection nested `depth` levels deep, a point at the bottom.
    fn nest(depth: usize) -> GeometryCollection {
        let mut c = GeometryCollection::new(vec![pt(1.0, 2.0)]).unwrap();
        for _ in 1..depth {
            c = GeometryCollection::new(vec![Geometry::GeometryCollection(c)]).unwrap();
        }
        c
    }

    #[test]
    fn construction_keeps_elements_in_order() {
        let gs = vec![pt(0.0, 0.0), line(&[(0.0, 0.0), (1.0, 1.0)])];
        let c = GeometryCollection::new(gs.clone()).unwrap();
        assert_eq!(c.geometries(), gs.as_slice());
        assert_eq!(c.len(), 2);
        assert_eq!(c.depth(), 1);
        assert_eq!(c.clone().into_geometries(), gs);
    }

    #[test]
    fn the_empty_collection_is_empty_of_depth_one_with_no_envelope() {
        let e = GeometryCollection::empty();
        assert_eq!(e, GeometryCollection::new(vec![]).unwrap());
        assert!(e.is_empty());
        assert_eq!(e.len(), 0);
        assert_eq!(e.depth(), 1);
        assert_eq!(e.dimension(), Dimension::Empty);
        assert_eq!(e.envelope(), None);
        assert_eq!(e.leaves().count(), 0);
    }

    #[test]
    fn a_collection_of_empty_elements_is_empty_at_any_depth() {
        let inner = GeometryCollection::new(vec![
            Geometry::MultiPoint(MultiPoint::empty()),
            Geometry::GeometryCollection(GeometryCollection::empty()),
        ])
        .unwrap();
        let c = GeometryCollection::new(vec![
            Geometry::MultiLineString(MultiLineString::empty()),
            Geometry::MultiPolygon(MultiPolygon::empty()),
            Geometry::GeometryCollection(inner),
        ])
        .unwrap();
        assert!(c.is_empty());
        assert_eq!(c.len(), 3);
        assert_eq!(c.depth(), 3);
        assert_eq!(c.dimension(), Dimension::Empty);
        assert_eq!(c.envelope(), None);
    }

    #[test]
    fn dimension_is_the_highest_non_empty_element() {
        let points = GeometryCollection::new(vec![pt(0.0, 0.0)]).unwrap();
        assert_eq!(points.dimension(), Dimension::Point);
        let mixed = GeometryCollection::new(vec![
            Geometry::GeometryCollection(points),
            line(&[(0.0, 0.0), (1.0, 0.0)]),
            // An empty areal element does not raise the dimension.
            Geometry::MultiPolygon(MultiPolygon::empty()),
        ])
        .unwrap();
        assert!(!mixed.is_empty());
        assert_eq!(mixed.dimension(), Dimension::Line);
        let with_area = GeometryCollection::new(vec![
            Geometry::GeometryCollection(mixed),
            Geometry::Polygon(square(0.0, 0.0, 1.0)),
        ])
        .unwrap();
        assert_eq!(with_area.dimension(), Dimension::Area);
    }

    #[test]
    fn the_envelope_selects_input_coordinates_across_depths() {
        let inner = GeometryCollection::new(vec![
            pt(-3.0, 7.0),
            Geometry::MultiPolygon(MultiPolygon::new(vec![square(10.0, -2.0, 1.0)]).unwrap()),
        ])
        .unwrap();
        let c = GeometryCollection::new(vec![
            line(&[(0.0, 0.0), (1.0, 1.0)]),
            Geometry::GeometryCollection(inner),
        ])
        .unwrap();
        let e = c.envelope().unwrap();
        assert_eq!(e.min(), Coord2d::new(-3.0, -2.0));
        assert_eq!(e.max(), Coord2d::new(11.0, 7.0));
    }

    #[test]
    fn leaves_are_depth_first_in_element_order() {
        let a = pt(1.0, 0.0);
        let b = pt(2.0, 0.0);
        let c = pt(3.0, 0.0);
        let d = pt(4.0, 0.0);
        let e = pt(5.0, 0.0);
        let deep = GeometryCollection::new(vec![c.clone()]).unwrap();
        let mid = GeometryCollection::new(vec![
            b.clone(),
            Geometry::GeometryCollection(deep),
            Geometry::GeometryCollection(GeometryCollection::empty()),
            d.clone(),
        ])
        .unwrap();
        let top = GeometryCollection::new(vec![
            a.clone(),
            Geometry::GeometryCollection(mid),
            e.clone(),
        ])
        .unwrap();
        let got: Vec<Geometry> = top.leaves().cloned().collect();
        assert_eq!(got, vec![a, b, c, d, e]);
        assert_eq!(top.depth(), 3);
    }

    #[test]
    fn nesting_is_accepted_at_the_limit_and_refused_one_beyond() {
        let at_limit = nest(GeometryCollection::MAX_NESTING);
        assert_eq!(at_limit.depth(), 32);
        assert_eq!(at_limit.leaves().count(), 1);
        let err = GeometryCollection::new(vec![Geometry::GeometryCollection(at_limit.clone())])
            .unwrap_err();
        assert_eq!(err, GeometryError::NestingTooDeep { limit: 32 });
        // A shallow sibling does not rescue a deep one, and order does not
        // matter.
        let err =
            GeometryCollection::new(vec![pt(0.0, 0.0), Geometry::GeometryCollection(at_limit)])
                .unwrap_err();
        assert_eq!(err, GeometryError::NestingTooDeep { limit: 32 });
    }

    #[test]
    fn depth_is_one_more_than_the_deepest_child() {
        let c = GeometryCollection::new(vec![
            Geometry::GeometryCollection(nest(2)),
            Geometry::GeometryCollection(nest(5)),
            pt(0.0, 0.0),
        ])
        .unwrap();
        assert_eq!(c.depth(), 6);
    }

    // --- measures (ops::measure) ---

    use crate::accum::Accum;
    use crate::ops::{area, collection_area, collection_length, length, polygon_length};

    fn holed() -> Polygon {
        let ring = |pts: &[(f64, f64)]| {
            LinearRing::new(CoordSeq::new(
                pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect(),
            ))
            .unwrap()
        };
        Polygon::new(
            ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            vec![ring(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])],
        )
        .unwrap()
    }

    #[test]
    fn polygon_length_counts_every_ring() {
        assert_eq!(polygon_length(&holed()), 16.0 + 4.0);
        assert_eq!(polygon_length(&square(0.0, 0.0, 1.0)), 4.0);
    }

    #[test]
    fn collection_measures_are_the_accum_sum_over_components_in_leaf_order() {
        let a = square(0.0, 0.0, 0.1);
        let b = holed();
        let c = square(10.0, 10.0, 0.3);
        let l1 = LineString::new(CoordSeq::new(vec![
            Coord2d::new(0.0, 0.0),
            Coord2d::new(0.1, 0.7),
        ]))
        .unwrap();
        let l2 = LineString::new(CoordSeq::new(vec![
            Coord2d::new(5.0, 5.0),
            Coord2d::new(5.3, 5.9),
            Coord2d::new(6.0, 5.0),
        ]))
        .unwrap();
        let inner = GeometryCollection::new(vec![
            Geometry::MultiPolygon(MultiPolygon::new(vec![b.clone(), c.clone()]).unwrap()),
            Geometry::MultiLineString(MultiLineString::new(vec![l2.clone()])),
            pt(9.0, 9.0),
        ])
        .unwrap();
        let gc = GeometryCollection::new(vec![
            Geometry::Polygon(a.clone()),
            Geometry::LineString(l1.clone()),
            Geometry::GeometryCollection(inner),
        ])
        .unwrap();

        let mut acc = Accum::new();
        for p in [&a, &b, &c] {
            acc.add(area(p));
        }
        assert_eq!(collection_area(&gc).to_bits(), acc.value().to_bits());

        let mut acc = Accum::new();
        acc.add(polygon_length(&a));
        acc.add(length(&l1));
        acc.add(polygon_length(&b));
        acc.add(polygon_length(&c));
        acc.add(length(&l2));
        assert_eq!(collection_length(&gc).to_bits(), acc.value().to_bits());
    }

    #[test]
    fn overlapping_elements_count_twice() {
        let p = Geometry::Polygon(square(0.0, 0.0, 2.0));
        let gc = GeometryCollection::new(vec![p.clone(), p]).unwrap();
        assert_eq!(collection_area(&gc), 8.0);
        assert_eq!(collection_length(&gc), 16.0);
    }

    #[test]
    fn a_one_polygon_collection_measures_exactly_as_the_polygon() {
        let p = holed();
        let gc = GeometryCollection::new(vec![Geometry::Polygon(p.clone())]).unwrap();
        assert_eq!(collection_area(&gc).to_bits(), area(&p).to_bits());
        assert_eq!(
            collection_length(&gc).to_bits(),
            polygon_length(&p).to_bits()
        );
    }

    #[test]
    fn measures_are_zero_exactly_when_there_is_nothing_to_measure() {
        let points = GeometryCollection::new(vec![
            pt(0.0, 0.0),
            Geometry::MultiPolygon(MultiPolygon::empty()),
            Geometry::GeometryCollection(GeometryCollection::empty()),
        ])
        .unwrap();
        assert_eq!(collection_area(&points), 0.0);
        assert_eq!(collection_length(&points), 0.0);
        assert_eq!(collection_area(&GeometryCollection::empty()), 0.0);
        let lines = GeometryCollection::new(vec![line(&[(0.0, 0.0), (3.0, 4.0)])]).unwrap();
        assert_eq!(collection_area(&lines), 0.0);
        assert_eq!(collection_length(&lines), 5.0);
        // A tiny polygon among large ones still makes the area nonzero.
        let tiny = square(1e6, 1e6, 1e-6);
        let only_tiny = GeometryCollection::new(vec![Geometry::Polygon(tiny)]).unwrap();
        assert!(collection_area(&only_tiny) > 0.0);
    }

    // --- the predicates on collections ---

    #[test]
    fn relate_accepts_collections_and_refuses_overlapping_polygon_leaves() {
        use crate::predicates::{dimension, intersects, relate};
        use crate::topology::de9im::IntersectionMatrix;
        // A singleton collection relates exactly as its element, on either
        // side, through every entry point.
        let gc = GeometryCollection::new(vec![Geometry::Polygon(square(0.0, 0.0, 1.0))]).unwrap();
        let p = square(0.0, 0.0, 1.0);
        let reference = IntersectionMatrix::compute(&p, &p);
        assert_eq!(relate(&gc, &p).unwrap().matrix(), reference);
        assert_eq!(relate(&p, &gc).unwrap().matrix(), reference);
        assert!(intersects(&Geometry::GeometryCollection(gc.clone()), &p).unwrap());
        assert_eq!(
            IntersectionMatrix::compute_geometry(
                &Geometry::Polygon(p.clone()),
                &Geometry::GeometryCollection(gc.clone())
            )
            .unwrap(),
            reference
        );
        let mask = crate::prepared::PreparedPolygonal::new(&p).unwrap();
        assert_eq!(
            mask.relate(&Geometry::GeometryCollection(gc.clone()))
                .unwrap(),
            reference
        );
        // Overlapping polygonal leaves are refused with both leaf paths.
        let overlapping = GeometryCollection::new(vec![
            Geometry::Polygon(square(0.0, 0.0, 2.0)),
            Geometry::Polygon(square(1.0, 1.0, 2.0)),
        ])
        .unwrap();
        assert_eq!(
            relate(&overlapping, &p).unwrap_err(),
            GeometryError::CollectionElementsOverlap {
                first: vec![0],
                second: vec![1],
            }
        );
        assert_eq!(dimension(&gc), Dimension::Area);
        assert_eq!(dimension(&GeometryCollection::empty()), Dimension::Empty);
    }
}
