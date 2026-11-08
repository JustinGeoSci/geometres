// src/types3d/mod.rs

//! 3D geometry types, available only under `feature = "3d"`.
//!
//! These types mirror the 2D types in [`crate::types`] with an added `z`
//! component. They are separate types rather than 2D types with an
//! optional Z.
//! All topology computation in `geometres` uses only the XY projection;
//! Z is stored and round-tripped but takes no part in any predicate.

mod coord3d;
mod line_string3d;
mod linear_ring3d;
mod multi_line_string3d;
mod multi_point3d;
mod multi_polygon3d;
mod point3d;
mod polygon3d;

pub use coord3d::{Coord3d, CoordSeq3d};
pub use line_string3d::LineString3d;
pub use linear_ring3d::LinearRing3d;
pub use multi_line_string3d::MultiLineString3d;
pub use multi_point3d::MultiPoint3d;
pub use multi_polygon3d::MultiPolygon3d;
pub use point3d::Point3d;
pub use polygon3d::Polygon3d;
