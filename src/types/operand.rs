// src/types/operand.rs

//! Operand traits for the overlay entry points: a polygonal operand is a
//! `Polygon` or a `MultiPolygon`, a lineal one a `LineString` or a
//! `MultiLineString`.
//!
//! # Why traits, and why sealed
//!
//! Taking `&impl Polygonal` lets one function accept a single part or many
//! without a second family of `*_multi` entry points, and a call with a
//! `&Polygon` reads the same as one with a `&MultiPolygon`.
//!
//! The traits are **sealed**, implementable only here, because the operations
//! rely on what they promise: a `Polygonal`'s parts meet at most at
//! finitely many points. `Polygon` and `MultiPolygon` guarantee that by
//! construction. An outside implementation could hand back overlapping
//! parts, and the exact winding guard would then refuse every call. The
//! seal closes that door rather than leaving it to the guard.

use crate::types::{LineString, MultiLineString, MultiPolygon, Polygon};

mod sealed {
    pub trait Sealed {}
    impl Sealed for crate::types::Polygon {}
    impl Sealed for crate::types::MultiPolygon {}
    impl Sealed for crate::types::LineString {}
    impl Sealed for crate::types::MultiLineString {}
}

/// A polygonal operand: one polygon, or the parts of a valid multipolygon.
pub trait Polygonal: sealed::Sealed {
    /// The parts. Interior-disjoint, meeting at most at finitely many points.
    fn parts(&self) -> &[Polygon];
}

impl Polygonal for Polygon {
    fn parts(&self) -> &[Polygon] {
        core::slice::from_ref(self)
    }
}

impl Polygonal for MultiPolygon {
    fn parts(&self) -> &[Polygon] {
        self.polygons()
    }
}

/// A lineal operand: one line string, or the components of a multi line
/// string.
pub trait Lineal: sealed::Sealed {
    /// The components, in order.
    fn components(&self) -> &[LineString];
}

impl Lineal for LineString {
    fn components(&self) -> &[LineString] {
        core::slice::from_ref(self)
    }
}

impl Lineal for MultiLineString {
    fn components(&self) -> &[LineString] {
        self.line_strings()
    }
}
