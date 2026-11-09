// src/predicates/mod.rs

//! The OGC spatial predicates, for any pair of geometries.
//!
//! # Entry points
//!
//! - [`relate`] computes one exact DE-9IM matrix for two geometries and
//!   returns it as a [`Relation`], from which every predicate is read. It is
//!   the entry point for a caller that wants "the spatial relation" rather
//!   than one predicate.
//! - One free function per predicate ([`equals`], [`disjoint`],
//!   [`intersects`], [`touches`], [`crosses`], [`within`], [`contains`],
//!   [`overlaps`], [`covers`] and [`covered_by`]), each `relate(a, b)` read
//!   one way.
//!
//! Every one of them takes its operands as `&impl Relatable`: a borrowed
//! [`Point`], [`LineString`], [`Polygon`], [`MultiPoint`],
//! [`MultiLineString`], [`MultiPolygon`], [`GeometryCollection`] or
//! [`Geometry`], in any combination. A collection operand is accepted when
//! its polygonal leaves together form a valid `MultiPolygon` (line and point
//! leaves are unrestricted), and refused with
//! `GeometryError::CollectionElementsOverlap` otherwise. A caller that wants
//! the union of overlapping elements unions them first and relates the vouched
//! result. For `crosses` and `overlaps` a collection's dimension is its
//! highest non-empty leaf's ([`GeometryCollection::dimension`]).
//!
//! # The matrix level
//!
//! Each predicate is also an [`IntersectionMatrix`] method, for a caller
//! that already holds a matrix (a prepared geometry's, or one parsed with
//! [`IntersectionMatrix::from_pattern`]). `crosses` and `overlaps` are
//! defined by the operands' dimensions, so their methods take the two
//! dimensions ([`dimension`] gives them).
//!
//! | Method | OGC DE-9IM pattern |
//! |---|---|
//! | [`equals`](IntersectionMatrix::equals) | `T*F**FFF*` |
//! | [`disjoint`](IntersectionMatrix::disjoint) | `FF*FF****` |
//! | [`intersects`](IntersectionMatrix::intersects) | complement of `disjoint` |
//! | [`touches`](IntersectionMatrix::touches) | `FT*******` / `F**T*****` / `F***T****` |
//! | [`crosses`](IntersectionMatrix::crosses) | by dimension pair; see the method |
//! | [`within`](IntersectionMatrix::within) | `T*F**F***` |
//! | [`contains`](IntersectionMatrix::contains) | `T*****FF*` |
//! | [`overlaps`](IntersectionMatrix::overlaps) | by dimension pair; see the method |
//! | [`covers`](IntersectionMatrix::covers) | `T*****FF*` / `*T****FF*` / `***T**FF*` / `****T*FF*` |
//! | [`covered_by`](IntersectionMatrix::covered_by) | `T*F**F***` / `*TF**F***` / `**FT*F***` / `**F*TF***` |
//!
//! # Exactness
//!
//! The matrix is exact: every sign decision is `arithmetica::orient2d` on
//! the input coordinates (`topology`). Each predicate is a fixed boolean
//! function of the matrix and the two dimensions, so each is exact too.

use crate::GeometryError;
use crate::topology::de9im::{IntersectionMatrix, Operand};
use crate::topology::label::Dimension;
use crate::types::{
    Geometry, GeometryCollection, LineString, MultiLineString, MultiPoint, MultiPolygon, Point,
    Polygon,
};

// ===========================================================================
// Operands
// ===========================================================================

mod sealed {
    /// The crate-private borrowed view behind [`super::Relatable`]: the
    /// geometry, borrowed where it lies, with its type. Opaque outside this
    /// crate.
    pub struct View<'a>(pub(crate) super::Borrowed<'a>);

    /// Seals [`super::Relatable`] and carries its view.
    pub trait Sealed {
        /// The view of `self`.
        #[doc(hidden)]
        fn view(&self) -> View<'_>;
    }
}

/// A relatable geometry, borrowed, with its type: what the predicates reduce
/// to an [`Operand`], and what `is_simple` and `validity_issues` dispatch on.
#[derive(Clone, Copy)]
pub(crate) enum Borrowed<'a> {
    Point(&'a Point),
    MultiPoint(&'a MultiPoint),
    LineString(&'a LineString),
    MultiLineString(&'a MultiLineString),
    Polygon(&'a Polygon),
    MultiPolygon(&'a MultiPolygon),
    GeometryCollection(&'a GeometryCollection),
}

impl<'a> Borrowed<'a> {
    /// The variant `g` holds.
    pub(crate) fn from_geometry(g: &'a Geometry) -> Self {
        match g {
            Geometry::Point(p) => Self::Point(p),
            Geometry::MultiPoint(m) => Self::MultiPoint(m),
            Geometry::LineString(l) => Self::LineString(l),
            Geometry::MultiLineString(m) => Self::MultiLineString(m),
            Geometry::Polygon(p) => Self::Polygon(p),
            Geometry::MultiPolygon(m) => Self::MultiPolygon(m),
            Geometry::GeometryCollection(c) => Self::GeometryCollection(c),
        }
    }

    /// The reduction the matrix computation needs.
    ///
    /// # Errors
    ///
    /// [`GeometryError::CollectionElementsOverlap`] for a collection whose
    /// polygonal leaves overlap: the predicates accept a
    /// collection whose polygonal leaves together form a valid
    /// `MultiPolygon`, and refuse the rest rather than answer for a union
    /// they cannot certify.
    fn operand(self) -> Result<Operand<'a>, GeometryError> {
        Ok(match self {
            Self::Point(p) => Operand::point(p),
            Self::MultiPoint(m) => Operand::multi_point(m),
            Self::LineString(l) => Operand::line_string(l),
            Self::MultiLineString(m) => Operand::multi_line_string(m),
            Self::Polygon(p) => Operand::polygon(p),
            Self::MultiPolygon(m) => Operand::multi_polygon(m),
            Self::GeometryCollection(gc) => {
                return crate::topology::collection::collection_operand(gc);
            }
        })
    }
}

/// A geometry the predicates accept: [`Point`], [`LineString`],
/// [`Polygon`], [`MultiPoint`], [`MultiLineString`], [`MultiPolygon`],
/// [`GeometryCollection`] or [`Geometry`].
///
/// Sealed: implemented only in this crate, because the matrix computation
/// relies on what the types promise (finite coordinates; a multipolygon's
/// parts meeting at most at finitely many points). A predicate call borrows
/// its operands; it never clones them into a [`Geometry`].
pub trait Relatable: sealed::Sealed {}

macro_rules! relatable {
    ($($ty:ty => $variant:ident),* $(,)?) => {$(
        impl sealed::Sealed for $ty {
            fn view(&self) -> sealed::View<'_> {
                sealed::View(Borrowed::$variant(self))
            }
        }
        impl Relatable for $ty {}
    )*};
}

relatable! {
    Point => Point,
    MultiPoint => MultiPoint,
    LineString => LineString,
    MultiLineString => MultiLineString,
    Polygon => Polygon,
    MultiPolygon => MultiPolygon,
    GeometryCollection => GeometryCollection,
}

impl sealed::Sealed for Geometry {
    fn view(&self) -> sealed::View<'_> {
        sealed::View(Borrowed::from_geometry(self))
    }
}
impl Relatable for Geometry {}

/// A reference to a relatable geometry is relatable, so a `&&Geometry` (from
/// iterating a slice of references, say) can be passed as it is.
impl<T: Relatable + ?Sized> sealed::Sealed for &T {
    fn view(&self) -> sealed::View<'_> {
        (**self).view()
    }
}
impl<T: Relatable + ?Sized> Relatable for &T {}

/// The borrowed view of `g`, for the other modules that dispatch on type.
pub(crate) fn borrowed(g: &impl Relatable) -> Borrowed<'_> {
    g.view().0
}

/// The topological dimension of a geometry, `Empty` for an empty one. It
/// governs `crosses` and `overlaps`.
#[must_use]
pub fn dimension(g: &impl Relatable) -> Dimension {
    match g.view().0 {
        // A collection's dimension is defined for every collection, the
        // refused ones included, so it is read from the type directly.
        Borrowed::GeometryCollection(c) => c.dimension(),
        // No other operand's reduction can fail.
        b => b.operand().map_or(Dimension::Empty, |o| o.dimension()),
    }
}

/// Classify a coordinate against a polygon: interior, boundary or exterior.
///
/// Exact. Every sign decision goes through `arithmetica::orient2d` on the
/// original coordinates, so the verdict is correct rather than
/// tolerance-dependent, and a point exactly on an edge is reported as
/// `Boundary` rather than falling to one side of it.
///
/// Public because constructive algorithms built on this crate state their
/// guarantees in terms of it ("the result is verified
/// `locate != Location::Exterior`"); a contract a consumer cannot check is a
/// claim rather than a guarantee.
pub use crate::topology::line_analysis::locate;

// ===========================================================================
// relate and the per-predicate functions
// ===========================================================================

/// The spatial relation of `a` to `b`: one exact DE-9IM matrix, from which
/// every predicate is read.
///
/// Every ordered pair of geometry types is supported, multipart included. A
/// `MultiPolygon`'s boundary is topological: a point where its parts touch
/// is boundary, as in OGC.
///
/// # Errors
///
/// [`GeometryError::CollectionElementsOverlap`] if either operand is a
/// [`GeometryCollection`] whose polygonal leaves overlap; no other pair
/// fails. A collection whose polygonal leaves form a valid `MultiPolygon`
/// is computed exactly.
pub fn relate(a: &impl Relatable, b: &impl Relatable) -> Result<Relation, GeometryError> {
    let (oa, ob) = (a.view().0.operand()?, b.view().0.operand()?);
    let m = IntersectionMatrix::compute_operands(&oa, &ob);
    Ok(Relation::from_matrix(&m, oa.dimension(), ob.dimension()))
}

/// `a` and `b` have at least one point in common: [`Relation::intersects`].
///
/// The same answer as `relate(a, b)?.intersects()` for every pair. When both
/// operands are areal (a polygon or a multipolygon), the matrix is not
/// built: the answer is decided by exact searches for one common point,
/// which stop at the first one found.
///
/// # Errors
///
/// As [`relate`]: a collection operand whose polygonal leaves overlap.
pub fn intersects(a: &impl Relatable, b: &impl Relatable) -> Result<bool, GeometryError> {
    // A point against a line or an area: the point's exact location, with
    // no operand built.
    if let Some(v) = point_meets(a.view().0, b.view().0) {
        return Ok(v);
    }
    let (oa, ob) = (a.view().0.operand()?, b.view().0.operand()?);
    if let (Operand::Area(x), Operand::Area(y)) = (&oa, &ob) {
        use crate::topology::areal::AreaOperand;
        let (wx, wy) = crate::topology::de9im::windowed_pair(x, y);
        return Ok(crate::topology::polygon_analysis::areas_intersect(
            wx.as_ref().map_or(x as &dyn AreaOperand, |w| w),
            wy.as_ref().map_or(y as &dyn AreaOperand, |w| w),
        ));
    }
    let m = IntersectionMatrix::compute_operands(&oa, &ob);
    Ok(Relation::from_matrix(&m, oa.dimension(), ob.dimension()).intersects())
}

/// Whether a single point and a line or areal geometry meet, in either
/// order, decided by the point's exact location: on some segment of a line
/// (`on_segment`, the closed segment, so endpoints count), or not in the
/// exterior of a polygon or a valid multipolygon (`Areal::locate`). `None`
/// for every other pair, which `intersects` answers through the matrix.
fn point_meets(a: Borrowed<'_>, b: Borrowed<'_>) -> Option<bool> {
    use crate::exact::kernel::on_segment;
    use crate::topology::areal::{AreaOperand, Areal};
    use crate::topology::label::Location;
    let (p, other) = match (a, b) {
        (Borrowed::Point(p), o) | (o, Borrowed::Point(p)) => (p.coord(), o),
        _ => return None,
    };
    let on_line = |l: &LineString| {
        l.coords()
            .as_slice()
            .windows(2)
            .any(|w| on_segment(p, w[0], w[1]))
    };
    match other {
        Borrowed::LineString(l) => Some(on_line(l)),
        Borrowed::MultiLineString(m) => Some(m.line_strings().iter().any(on_line)),
        Borrowed::Polygon(q) => Some(Areal::single(q).locate(p) != Location::Exterior),
        Borrowed::MultiPolygon(m) if !m.is_empty() => {
            Some(Areal::multi(m).locate(p) != Location::Exterior)
        }
        _ => None,
    }
}

/// `a` and `b` have no point in common: [`Relation::disjoint`]. The
/// negation of [`intersects`], decided the same way.
///
/// # Errors
///
/// As [`relate`]: a collection operand whose polygonal leaves overlap.
pub fn disjoint(a: &impl Relatable, b: &impl Relatable) -> Result<bool, GeometryError> {
    Ok(!intersects(a, b)?)
}

/// Defines one free function per predicate, each `relate(a, b)` read through
/// the [`Relation`] method of the same name.
macro_rules! predicate_fns {
    ($($(#[doc = $doc:literal])* $name:ident),* $(,)?) => {$(
        $(#[doc = $doc])*
        ///
        /// # Errors
        ///
        /// As [`relate`]: a collection operand whose polygonal leaves
        /// overlap.
        pub fn $name(a: &impl Relatable, b: &impl Relatable) -> Result<bool, GeometryError> {
            Ok(relate(a, b)?.$name())
        }
    )*};
}

predicate_fns! {
    /// `a` and `b` are the same point set: [`Relation::equals`].
    equals,
    /// The boundaries of `a` and `b` meet and their interiors do not:
    /// [`Relation::touches`].
    touches,
    /// `a` crosses `b`, by the definition for their dimensions:
    /// [`Relation::crosses`].
    crosses,
    /// `a` lies in `b`, their interiors meeting: [`Relation::within`].
    within,
    /// `b` lies in `a`, their interiors meeting: [`Relation::contains`].
    contains,
    /// `a` and `b` overlap, by the definition for their dimensions:
    /// [`Relation::overlaps`].
    overlaps,
    /// No point of `b` lies outside `a`, and they meet: [`Relation::covers`].
    covers,
    /// No point of `a` lies outside `b`, and they meet:
    /// [`Relation::covered_by`]. The same as `covers(b, a)`.
    covered_by,
}

// ===========================================================================
// Relation
// ===========================================================================

/// The spatial relation of one ordered pair: its exact DE-9IM matrix and the
/// two operands' dimensions, with every predicate read from them.
///
/// Returned by [`relate`]. Every predicate is authoritative for every
/// ordered pair of geometry types: `crosses` and `overlaps` are defined for
/// all of them, by dimension.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Relation {
    matrix: IntersectionMatrix,
    dim_a: Dimension,
    dim_b: Dimension,
}

impl Relation {
    /// The relation given by a matrix and the two geometries' dimensions.
    #[must_use]
    pub fn from_matrix(m: &IntersectionMatrix, dim_a: Dimension, dim_b: Dimension) -> Self {
        Self {
            matrix: *m,
            dim_a,
            dim_b,
        }
    }

    /// The exact DE-9IM matrix of `a` against `b`.
    #[must_use]
    pub fn matrix(&self) -> IntersectionMatrix {
        self.matrix
    }

    /// `T*F**FFF*`: the same point set. Two empty geometries are also equal.
    /// See [`IntersectionMatrix::equals`].
    #[must_use]
    pub fn equals(&self) -> bool {
        self.matrix.equals()
    }

    /// `FF*FF****`: no point in common.
    #[must_use]
    pub fn disjoint(&self) -> bool {
        self.matrix.disjoint()
    }

    /// The complement of [`disjoint`](Self::disjoint).
    #[must_use]
    pub fn intersects(&self) -> bool {
        self.matrix.intersects()
    }

    /// Boundaries meet, interiors do not.
    #[must_use]
    pub fn touches(&self) -> bool {
        self.matrix.touches()
    }

    /// By the two dimensions; see [`IntersectionMatrix::crosses`].
    #[must_use]
    pub fn crosses(&self) -> bool {
        self.matrix.crosses(self.dim_a, self.dim_b)
    }

    /// `T*F**F***`: `a` lies in `b`, interiors meeting.
    #[must_use]
    pub fn within(&self) -> bool {
        self.matrix.within()
    }

    /// `T*****FF*`: `b` lies in `a`, interiors meeting.
    #[must_use]
    pub fn contains(&self) -> bool {
        self.matrix.contains()
    }

    /// By the two dimensions; see [`IntersectionMatrix::overlaps`].
    #[must_use]
    pub fn overlaps(&self) -> bool {
        self.matrix.overlaps(self.dim_a, self.dim_b)
    }

    /// They meet, and no point of `b` is outside `a`.
    #[must_use]
    pub fn covers(&self) -> bool {
        self.matrix.covers()
    }

    /// They meet, and no point of `a` is outside `b`.
    #[must_use]
    pub fn covered_by(&self) -> bool {
        self.matrix.covered_by()
    }
}

// ===========================================================================
// The matrix-level predicates
// ===========================================================================

impl IntersectionMatrix {
    /// Tests whether two geometries are spatially equal.
    ///
    /// **OGC DE-9IM pattern:** `T*F**FFF*`
    ///
    /// Both interiors must intersect (`II ≠ Empty`) and no part of either,
    /// interior or boundary, may lie outside the other (`IE`, `BE`, `EI`,
    /// `EB` all `Empty`). Two empty geometries are equal as point sets, so
    /// this method returns `true` for them, although their matrix,
    /// `FFFFFFFF2`, does not match the pattern.
    #[must_use]
    pub fn equals(&self) -> bool {
        // Two empty geometries are equal as point sets. Their matrix,
        // `FFFFFFFF2`, fails `T*F**FFF*` on `II`; it is also the only
        // matrix with all eight non-`EE` cells empty, because a non-empty
        // geometry has a non-empty interior.
        let all_empty = self.cells()[..8].iter().all(|&d| d == Dimension::Empty);
        if all_empty {
            return true;
        }
        self.interior_interior != Dimension::Empty
            && self.interior_exterior == Dimension::Empty
            && self.boundary_exterior == Dimension::Empty
            && self.exterior_interior == Dimension::Empty
            && self.exterior_boundary == Dimension::Empty
    }

    /// Tests whether two geometries have no point in common.
    ///
    /// **OGC DE-9IM pattern:** `FF*FF****`: all of `II`, `IB`, `BI`, `BB`
    /// empty. Equivalent to `!self.intersects()`.
    #[must_use]
    pub fn disjoint(&self) -> bool {
        self.interior_interior == Dimension::Empty
            && self.interior_boundary == Dimension::Empty
            && self.boundary_interior == Dimension::Empty
            && self.boundary_boundary == Dimension::Empty
    }

    /// Tests whether two geometries share at least one point: the
    /// complement of [`disjoint`](Self::disjoint). At least one of `II`,
    /// `IB`, `BI`, `BB` is non-empty.
    #[must_use]
    pub fn intersects(&self) -> bool {
        !self.disjoint()
    }

    /// Tests whether two geometries share boundary points but not interior
    /// points.
    ///
    /// **OGC DE-9IM pattern:** `FT*******` / `F**T*****` / `F***T****`
    ///
    /// - `II = Empty`: the interiors do not meet.
    /// - At least one of `IB`, `BI`, `BB` non-empty: some boundary contact.
    #[must_use]
    pub fn touches(&self) -> bool {
        self.interior_interior == Dimension::Empty
            && (self.interior_boundary != Dimension::Empty
                || self.boundary_interior != Dimension::Empty
                || self.boundary_boundary != Dimension::Empty)
    }

    /// Tests whether geometry A is completely within geometry B.
    ///
    /// **OGC DE-9IM pattern:** `T*F**F***`
    ///
    /// - `II ≠ Empty`: A's interior meets B's interior.
    /// - `IE = Empty`: no part of A's interior lies outside B.
    /// - `BE = Empty`: no part of A's boundary lies outside B.
    ///
    /// `within` of a matrix is `contains` of its transpose.
    #[must_use]
    pub fn within(&self) -> bool {
        self.interior_interior != Dimension::Empty
            && self.interior_exterior == Dimension::Empty
            && self.boundary_exterior == Dimension::Empty
    }

    /// Tests whether geometry A completely contains geometry B.
    ///
    /// **OGC DE-9IM pattern:** `T*****FF*`
    ///
    /// - `II ≠ Empty`: B's interior meets A's interior.
    /// - `EI = Empty`: no part of B's interior lies outside A.
    /// - `EB = Empty`: no part of B's boundary lies outside A.
    ///
    /// `contains` of a matrix is `within` of its transpose.
    #[must_use]
    pub fn contains(&self) -> bool {
        self.interior_interior != Dimension::Empty
            && self.exterior_interior == Dimension::Empty
            && self.exterior_boundary == Dimension::Empty
    }

    /// Tests whether geometry A covers geometry B: they meet, and no point
    /// of B lies outside A.
    ///
    /// **OGC DE-9IM pattern:** `T*****FF*` / `*T****FF*` / `***T**FF*` /
    /// `****T*FF*`
    ///
    /// - One of `II`, `IB`, `BI`, `BB` non-empty: each alternative needs a
    ///   `T` in one of them, which is [`intersects`](Self::intersects). An
    ///   empty B is therefore covered by nothing.
    /// - `EI = Empty`: no part of B's interior lies in A's exterior.
    /// - `EB = Empty`: no part of B's boundary lies in A's exterior.
    ///
    /// More permissive than [`contains`](Self::contains), which also needs
    /// `II ≠ Empty`.
    #[must_use]
    pub fn covers(&self) -> bool {
        self.intersects()
            && self.exterior_interior == Dimension::Empty
            && self.exterior_boundary == Dimension::Empty
    }

    /// Tests whether geometry A is covered by geometry B: they meet, and no
    /// point of A lies outside B.
    ///
    /// **OGC DE-9IM pattern:** `T*F**F***` / `*TF**F***` / `**FT*F***` /
    /// `**F*TF***`
    ///
    /// [`covers`](Self::covers) of the transpose: `IE` and `BE` take the
    /// place of `EI` and `EB`, and `intersects` is symmetric.
    #[must_use]
    pub fn covered_by(&self) -> bool {
        self.intersects()
            && self.interior_exterior == Dimension::Empty
            && self.boundary_exterior == Dimension::Empty
    }

    /// `crosses` for any pair, by the OGC definition for its dimensions:
    ///
    /// | dimensions | pattern |
    /// |---|---|
    /// | `dim A < dim B` (point/line, point/area, line/area) | `T*T******` |
    /// | `dim A > dim B` | `T*****T**` |
    /// | line × line | `0********` |
    /// | point × point, area × area, or either empty | `false` |
    #[must_use]
    pub fn crosses(&self, dim_a: Dimension, dim_b: Dimension) -> bool {
        if dim_a == Dimension::Empty || dim_b == Dimension::Empty {
            return false;
        }
        let t = |d: Dimension| d != Dimension::Empty;
        if dim_a < dim_b {
            t(self.interior_interior) && t(self.interior_exterior)
        } else if dim_a > dim_b {
            t(self.interior_interior) && t(self.exterior_interior)
        } else if dim_a == Dimension::Line {
            self.interior_interior == Dimension::Point
        } else {
            false
        }
    }

    /// `overlaps` for any pair, by the OGC definition: the two geometries
    /// have the same dimension, their interiors meet **in that dimension**,
    /// and neither contains the other.
    ///
    /// | dimensions | pattern |
    /// |---|---|
    /// | point × point | `0*T***T**` (= `T*T***T**`) |
    /// | line × line | `1*T***T**` |
    /// | area × area | `2*T***T**` |
    /// | different dimensions, or either empty | `false` |
    ///
    /// Requiring `II` to equal the common dimension is the definition,
    /// `dim(I(a) ∩ I(b)) = dim(I(a)) = dim(I(b))`; for computed matrices it
    /// coincides with the matrix forms, because two areas' interiors are
    /// open and can only meet in an area.
    #[must_use]
    pub fn overlaps(&self, dim_a: Dimension, dim_b: Dimension) -> bool {
        if dim_a != dim_b || dim_a == Dimension::Empty {
            return false;
        }
        let t = |d: Dimension| d != Dimension::Empty;
        self.interior_interior == dim_a && t(self.interior_exterior) && t(self.exterior_interior)
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod fast_path_tests;

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::types::{Coord2d, CoordSeq, LinearRing};

    const DIMS: [Dimension; 4] = [
        Dimension::Empty,
        Dimension::Point,
        Dimension::Line,
        Dimension::Area,
    ];

    /// Every one of the 4⁹ matrices over `{F, 0, 1, 2}`.
    fn every_matrix() -> impl Iterator<Item = IntersectionMatrix> {
        (0..4usize.pow(9)).map(|mut k| {
            let mut c = [Dimension::Empty; 9];
            for cell in &mut c {
                *cell = DIMS[k % 4];
                k /= 4;
            }
            IntersectionMatrix {
                interior_interior: c[0],
                interior_boundary: c[1],
                interior_exterior: c[2],
                boundary_interior: c[3],
                boundary_boundary: c[4],
                boundary_exterior: c[5],
                exterior_interior: c[6],
                exterior_boundary: c[7],
                exterior_exterior: c[8],
            }
        })
    }

    /// The relations between the matrix-level predicates that hold by
    /// definition, on every matrix: `covered_by` is `covers` of the
    /// transpose, `within` is `contains` of the transpose, `disjoint` is the
    /// complement of `intersects`, and `equals`, `intersects`, `disjoint`
    /// and `touches` are symmetric.
    #[test]
    fn predicates_relate_as_defined_on_every_matrix() {
        for m in every_matrix() {
            let t = m.swap_ab();
            assert_eq!(m.covered_by(), t.covers(), "{m}");
            assert_eq!(m.within(), t.contains(), "{m}");
            assert_eq!(m.disjoint(), !m.intersects(), "{m}");
            assert_eq!(m.equals(), t.equals(), "{m}");
            assert_eq!(m.touches(), t.touches(), "{m}");
            if m.contains() {
                assert!(m.covers(), "{m}");
            }
            for da in DIMS {
                for db in DIMS {
                    assert_eq!(m.crosses(da, db), t.crosses(db, da), "{m} {da:?} {db:?}");
                    assert_eq!(m.overlaps(da, db), t.overlaps(db, da), "{m} {da:?} {db:?}");
                }
            }
        }
    }

    fn square(x: f64) -> Polygon {
        let c = |a: f64, b: f64| Coord2d::new(a, b);
        let ring = LinearRing::new(CoordSeq::new(vec![
            c(x, 0.0),
            c(x + 2.0, 0.0),
            c(x + 2.0, 2.0),
            c(x, 2.0),
        ]))
        .unwrap();
        Polygon::new_simple(ring).unwrap()
    }

    /// A concrete type, its `Geometry` wrapper and a reference to either
    /// give the same relation, read by the same methods.
    #[test]
    fn every_operand_form_gives_the_same_relation() {
        let (a, b) = (square(0.0), square(1.0));
        let (ga, gb) = (Geometry::Polygon(a.clone()), Geometry::Polygon(b.clone()));
        let r = relate(&a, &b).unwrap();
        assert_eq!(r, relate(&ga, &gb).unwrap());
        assert_eq!(r, relate(&&a, &gb).unwrap());
        assert_eq!(r, relate(&ga, &&&b).unwrap());
        assert_eq!(r.matrix(), IntersectionMatrix::compute(&a, &b));
        assert!(r.overlaps() && r.intersects() && !r.touches() && !r.covers());
        assert_eq!(overlaps(&a, &b), Ok(true));
        assert_eq!(dimension(&a), Dimension::Area);
        assert_eq!(dimension(&MultiPolygon::empty()), Dimension::Empty);
    }

    /// `Relation::from_matrix` keeps the dimensions that `crosses` and
    /// `overlaps` read.
    #[test]
    fn a_relation_reads_crosses_and_overlaps_by_its_dimensions() {
        let m = IntersectionMatrix::from_pattern("1F1FF0102").unwrap();
        let lines = Relation::from_matrix(&m, Dimension::Line, Dimension::Line);
        let areas = Relation::from_matrix(&m, Dimension::Area, Dimension::Area);
        assert!(lines.overlaps() && !areas.overlaps());
        assert_eq!(lines.matrix(), m);
    }
}
