// src/lib.rs

//! Exact, non-constructive computational geometry over `f64` coordinates.

#![forbid(unsafe_code)]
#![cfg_attr(docsrs, feature(doc_cfg))]
#![deny(clippy::unwrap_used)]
#![warn(missing_docs)]

//! Provides the validated geometry types, exact DE-9IM topological
//! predicates, validation, and certified measurements, backed by adaptive
//! floating-point arithmetic from the `arithmetica` crate.
//!
//! # What this crate does, and what it does not
//!
//! Every public function does one of three things: it validates, it decides
//! a predicate exactly, or it returns a quantity whose sign is certified and
//! whose magnitude has a documented bound ([`ops::area`], [`ops::length`],
//! [`accum::Accum`]). No public function returns a coordinate that is not an
//! input coordinate.
//!
//! Algorithms that construct new coordinates (the planar arrangement,
//! rounding to `f64`, overlay, offset, centroid, tessellation) are built on
//! top of this crate and vouch their results with its primitives. The
//! [`exact`], [`accum`] and [`codec`] modules are the public surface they use.
//!
//! # Format neutrality
//!
//! This crate takes coordinates as `f64` values and returns geometric
//! results. It has no knowledge of GIS file formats, coordinate reference
//! systems, or geographic features.

/// Increments a counter of the `counters` module by one, or by `$n`, on the
/// calling thread.
///
/// Without the `unstable-counters` feature this expands to an empty block:
/// the argument is not evaluated and nothing is compiled. Crate-internal; the
/// counters are read through `counters::snapshot`.
macro_rules! count {
    ($name:ident) => {
        count!($name, 1u64)
    };
    ($name:ident, $n:expr) => {{
        #[cfg(feature = "unstable-counters")]
        {
            $crate::counters::add(|c| &c.$name, u64::try_from($n).unwrap_or(u64::MAX));
        }
    }};
}
pub(crate) use count;

pub mod accum;
pub mod codec;
#[cfg(feature = "unstable-counters")]
pub mod counters;

/// Not API: outside semver, and may change or disappear in any release.
///
/// It exists only to hold this `compile_fail` doctest, which proves that
/// without the `unstable-counters` feature the counters do not exist: the
/// item is compiled only then, and the doctest must fail to compile because
/// `geometres::counters` cannot be found.
///
/// ```compile_fail,E0433
/// let _ = geometres::counters::snapshot();
/// ```
#[cfg(not(feature = "unstable-counters"))]
#[doc(hidden)]
pub mod __counters_absent_check {}
pub mod error;
pub mod exact;
pub mod index;
pub mod ops;
pub mod predicates;
pub mod prepared;
pub mod topology;
pub mod types;

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
pub mod types3d;

pub use error::GeometryError;
pub use types::{
    Coord2d, CoordSeq, LineString, LinearRing, MultiLineString, MultiPoint, MultiPolygon, Point,
    Polygon,
};

#[cfg(feature = "3d")]
#[cfg_attr(docsrs, doc(cfg(feature = "3d")))]
pub use types3d::{
    Coord3d, CoordSeq3d, LineString3d, LinearRing3d, MultiPolygon3d, Point3d, Polygon3d,
};
