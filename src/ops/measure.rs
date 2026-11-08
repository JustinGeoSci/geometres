// src/ops/measure.rs

//! The length of a polygon, and the area and length of a geometry
//! collection.
//!
//! # The conventions
//!
//! - **A collection's measure is the sum over its elements**. A collection
//!   imposes no condition between elements, so overlapping elements count
//!   twice, and an element repeated is counted again. Where the areal
//!   elements are interior-disjoint, the sum *is* the area of their union.
//! - **The length of a polygon is the length of its whole boundary**, every
//!   ring, holes included. [`super::perimeter`] remains the exterior alone.
//!
//! # How the sum *is* formed, and what it guarantees
//!
//! Each **component** of a non-collection element contributes one term: a
//! polygon's [`area`](super::area) or [`polygon_length`], a line's
//! [`length`](super::length). A multipart element contributes one term per
//! part, and a point contributes nothing. The terms are added with
//! [`Accum`], in the order [`GeometryCollection::leaves`] visits them, so
//! the result is deterministic, and a collection of one polygon gives
//! exactly that polygon's own value (`Accum` returns a single term
//! unchanged).
//!
//! **Sign.** Every term is non-negative, and for non-negative terms
//! `Σ|tᵢ|` is the exact sum `T`. `Accum`'s bound then reads
//! `|r − T| ≤ (ε + γ_k²)·T`, which is less than `T` for any `k` with
//! `k·ε < ½`. So `r` is zero if and only if `T` is, that is, **if and only
//! if every term is zero**. (The terms are added, not multiplied, so
//! `Accum`'s product-underflow exception does not arise here.) A valid
//! polygon's area is positive (its exact net area is, and the sign is
//! certified), and every ring or line has a positive length, so a
//! collection's area is zero exactly when it has no non-empty areal element,
//! and its length is zero exactly when it has no linear or areal element.
//! The positivity of a polygon's area rests on `ops/area.rs`'s exactness
//! argument, which assumes its products do not underflow; for coordinates
//! below about `2⁻⁴⁸⁴` in magnitude that assumption can fail.
//!
//! **Magnitude.** For the area, with `bᵢ` each polygon's own bound
//! (`ops/area.rs`) and `aᵢ` its computed area, the error of the sum is at
//! most `Σbᵢ + (ε + γ_k²)·Σaᵢ`. For the length the per-term accuracy is
//! measured, not proven ([`super::edge_length`]: `hypot`'s accuracy is the
//! platform's), so the collection's length is measured too; the composition
//! adds at most `(ε + γ_k²)·Σℓᵢ`.

use crate::accum::Accum;
use crate::types::{Geometry, GeometryCollection, Polygon};

use super::{area, edge_length, length};

/// The length of a polygon's whole boundary: every edge of every ring,
/// holes included, in one compensated sum.
///
/// [`perimeter`](super::perimeter) is the exterior ring alone. The accuracy
/// is that of [`length`]: a compensated sum of `hypot` edge lengths, measured
/// rather than proven. Positive for every polygon (each ring has at least
/// three distinct vertices).
#[must_use]
pub fn polygon_length(polygon: &Polygon) -> f64 {
    let mut acc = Accum::new();
    for ring in core::iter::once(polygon.exterior()).chain(polygon.holes()) {
        for (a, b) in ring.edges() {
            acc.add(edge_length(a, b));
        }
    }
    acc.value()
}

/// The sum of the areas of a collection's polygons, at every depth: each
/// polygon, alone or a part of a multipolygon, is one term of an [`Accum`]
/// sum. Points and lines contribute no area.
///
/// **Overlapping elements count twice**: a collection imposes no condition
/// between its elements. The module documentation gives the sign argument
/// (zero exactly when there is no polygon) and the bound.
#[must_use]
pub fn collection_area(collection: &GeometryCollection) -> f64 {
    let mut acc = Accum::new();
    for g in collection.leaves() {
        match g {
            Geometry::Polygon(p) => acc.add(area(p)),
            Geometry::MultiPolygon(m) => {
                for p in m.polygons() {
                    acc.add(area(p));
                }
            }
            Geometry::Point(_)
            | Geometry::MultiPoint(_)
            | Geometry::LineString(_)
            | Geometry::MultiLineString(_)
            | Geometry::GeometryCollection(_) => {}
        }
    }
    acc.value()
}

/// The sum of the lengths of a collection's lines and polygon boundaries, at
/// every depth: each line is one term ([`length`]) and each polygon one term
/// ([`polygon_length`], every ring). Points contribute nothing.
///
/// **Overlapping elements count twice.** The module documentation gives the
/// sign argument (zero exactly when there is no line or polygon) and the
/// accuracy, which is measured, not proven.
#[must_use]
pub fn collection_length(collection: &GeometryCollection) -> f64 {
    let mut acc = Accum::new();
    for g in collection.leaves() {
        match g {
            Geometry::LineString(l) => acc.add(length(l)),
            Geometry::MultiLineString(m) => {
                for l in m.line_strings() {
                    acc.add(length(l));
                }
            }
            Geometry::Polygon(p) => acc.add(polygon_length(p)),
            Geometry::MultiPolygon(m) => {
                for p in m.polygons() {
                    acc.add(polygon_length(p));
                }
            }
            Geometry::Point(_) | Geometry::MultiPoint(_) | Geometry::GeometryCollection(_) => {}
        }
    }
    acc.value()
}
