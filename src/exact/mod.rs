// src/exact/mod.rs

//! Exact ring- and segment-level primitives, public for crates built on this
//! one.
//!
//! # What every item here guarantees
//!
//! Every decision is `arithmetica::orient2d` on stored input coordinates, or
//! a comparison of stored input coordinates. Nothing is approximate, no
//! tolerance is involved, and **no coordinate is constructed**. A crossing is
//! named by the two edges that make it ([`EventAt::Crossing`]), never
//! computed. These are the same primitives `Polygon::new`,
//! `MultiPolygon::new`, `LinearRing::is_simple` and the DE-9IM matrix use, so
//! a caller that vouches its own output with them applies exactly the checks
//! the types apply.
//!
//! # Finiteness
//!
//! `orient2d` is exact for finite inputs and refuses non-finite ones. Inside
//! this crate every coordinate comes from a validated type, so the internal
//! forms of these primitives (in the private `kernel` submodule) treat that
//! refusal as unreachable. A caller outside the crate has no such guarantee:
//! [`Coord2d::new`] does not validate. So each public item either takes only
//! validated types ([`place_ring`], [`first_self_contact`], [`rings`],
//! [`Areal::new`]) or checks every coordinate it is given and returns
//! [`GeometryError::NonFiniteCoordinate`] rather than a sign or a location
//! derived from a NaN.
//!
//! For the checked items, `coord_index` in that error is the position of the
//! offending coordinate among the function's coordinate arguments, counting
//! from 0 (`p` is 0, `a` is 1, `b` is 2 for [`on_segment`]), and `axis` is
//! 0 for x and 1 for y.
//!
//! # Where the items live
//!
//! [`rings`], [`Areal`], [`EventAt`] and [`SegmentEvent`] are defined in the
//! topology layer and re-exported here, which keeps their bodies beside the
//! DE-9IM code that uses them.
//! `ring_is_degenerate` is not here: it is the constructor's area-threshold
//! guard, not an exact predicate, and it is public as
//! [`types::ring_is_degenerate`](crate::types::ring_is_degenerate).

pub(crate) mod kernel;

use crate::GeometryError;
use crate::topology::areal::AreaOperand;
use crate::topology::label::Location;
use crate::topology::line_analysis;
use crate::types::{Coord2d, LinearRing, Polygon};

pub use crate::topology::areal::Areal;
pub use crate::topology::line_analysis::{EventAt, SegmentEvent, rings};
pub use kernel::{RingPlacement, place_ring};

/// `Ok(())` if every coordinate in `coords` is finite; otherwise the first
/// non-finite one, by position and axis.
fn require_finite(coords: &[Coord2d]) -> Result<(), GeometryError> {
    for (coord_index, c) in coords.iter().enumerate() {
        if !c.x.is_finite() {
            return Err(GeometryError::NonFiniteCoordinate {
                coord_index,
                axis: 0,
            });
        }
        if !c.y.is_finite() {
            return Err(GeometryError::NonFiniteCoordinate {
                coord_index,
                axis: 1,
            });
        }
    }
    Ok(())
}

/// Whether `p` lies on the closed segment `[a, b]`, endpoints included.
///
/// Exact: a bounding-box test on stored coordinates, then
/// `orient2d(a, b, p) = 0`. A degenerate segment (`a == b`) contains exactly
/// the point `a`.
///
/// # Errors
///
/// [`GeometryError::NonFiniteCoordinate`] if any of `p`, `a`, `b` has a
/// non-finite component (positions 0, 1, 2). No answer is given for a point
/// that is not a point.
pub fn on_segment(p: Coord2d, a: Coord2d, b: Coord2d) -> Result<bool, GeometryError> {
    require_finite(&[p, a, b])?;
    Ok(kernel::on_segment(p, a, b))
}

/// The first pair of edges of `ring` that touch anywhere other than at the
/// vertex two consecutive edges share, by edge index; `None` if the ring is
/// simple in the OGC sense.
///
/// Edge `i` runs from vertex `i` to vertex `i + 1` (the last edge closes the
/// ring). Pairs are found in order of the first edge, and a crossing, a
/// pinch (the ring passing through one point twice), a vertex on a
/// non-adjacent edge and a fold-back of consecutive edges are all contacts:
/// the absence of proper crossings alone is not simplicity. Exact: every test
/// is `orient2d` on the ring's stored coordinates or a comparison of them.
/// `O(n²)` in the vertex count.
///
/// Takes a validated ring, so its coordinates are finite and the answer is
/// always defined.
pub fn first_self_contact(ring: &LinearRing) -> Option<(usize, usize)> {
    kernel::first_self_contact(ring.coords().as_slice())
}

/// The exact location of `p` relative to `polygon`: interior, boundary
/// (on an edge of any ring, vertices included) or exterior.
///
/// The same function the DE-9IM matrix uses for a point against a polygon:
/// an exact on-edge test over every ring, then an exact crossing-number test
/// over the exterior and each hole.
///
/// # Errors
///
/// [`GeometryError::NonFiniteCoordinate`] if `p` has a non-finite component
/// (position 0).
pub fn locate(p: Coord2d, polygon: &Polygon) -> Result<Location, GeometryError> {
    require_finite(&[p])?;
    Ok(line_analysis::locate(p, polygon))
}

/// The exact location of `p` relative to an areal operand: interior,
/// boundary or exterior.
///
/// A point on the boundary of any part is on the boundary, including a point
/// where parts touch: the boundary is topological.
///
/// # Errors
///
/// [`GeometryError::NonFiniteCoordinate`] if `p` has a non-finite component
/// (position 0).
pub fn locate_in_areal(p: Coord2d, area: &Areal<'_>) -> Result<Location, GeometryError> {
    require_finite(&[p])?;
    Ok(area.locate(p))
}

/// Every event on the segment `p → q` against `area`, each with the location
/// of the open stretch that follows it.
///
/// The events are: the start `p` (always first); every vertex of `area`
/// strictly inside the segment, once each; and every proper crossing of the
/// segment with a ring edge, **named by that edge** ([`EventAt::Crossing`])
/// and not constructed, except a crossing that coincides with a listed
/// vertex. The events are distinct points. They are **not** ordered along
/// the segment: ordering needs the crossings' coordinates, which this crate
/// does not construct. A caller that needs them constructs them exactly.
///
/// Every location is decided by `orient2d` on input coordinates.
///
/// # Errors
///
/// [`GeometryError::NonFiniteCoordinate`] if `p` or `q` has a non-finite
/// component (positions 0 and 1).
pub fn segment_events(
    p: Coord2d,
    q: Coord2d,
    area: &Areal<'_>,
) -> Result<Vec<SegmentEvent>, GeometryError> {
    require_finite(&[p, q])?;
    Ok(line_analysis::segment_events(p, q, area))
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::types::CoordSeq;

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    fn square() -> Polygon {
        Polygon::new_simple(
            LinearRing::new(CoordSeq::new(vec![
                c(0.0, 0.0),
                c(4.0, 0.0),
                c(4.0, 4.0),
                c(0.0, 4.0),
            ]))
            .unwrap(),
        )
        .unwrap()
    }

    /// Every checked item refuses a non-finite coordinate in every position
    /// and on both axes, naming the position, and answers as the unchecked
    /// kernel does on finite input.
    #[test]
    fn checked_items_refuse_non_finite_coordinates_by_position() {
        let ok = c(1.0, 0.0);
        for bad in [c(f64::NAN, 0.0), c(0.0, f64::INFINITY)] {
            // The x component is checked first, so a finite x means y is the bad axis.
            let axis = u8::from(bad.x.is_finite());
            for pos in 0..3 {
                let mut args = [ok, c(0.0, 0.0), c(4.0, 0.0)];
                args[pos] = bad;
                assert_eq!(
                    on_segment(args[0], args[1], args[2]),
                    Err(GeometryError::NonFiniteCoordinate {
                        coord_index: pos,
                        axis
                    }),
                    "on_segment, position {pos}"
                );
            }
            let sq = square();
            let area = Areal::new(&sq);
            let at = |coord_index| GeometryError::NonFiniteCoordinate { coord_index, axis };
            assert_eq!(locate(bad, &sq), Err(at(0)));
            assert_eq!(locate_in_areal(bad, &area), Err(at(0)));
            assert_eq!(segment_events(bad, ok, &area).err(), Some(at(0)));
            assert_eq!(segment_events(ok, bad, &area).err(), Some(at(1)));
        }
    }

    #[test]
    fn checked_items_agree_with_the_kernel_on_finite_input() {
        let sq = square();
        let area = Areal::new(&sq);
        assert_eq!(on_segment(c(2.0, 0.0), c(0.0, 0.0), c(4.0, 0.0)), Ok(true));
        assert_eq!(on_segment(c(2.0, 1.0), c(0.0, 0.0), c(4.0, 0.0)), Ok(false));
        for (p, want) in [
            (c(2.0, 2.0), Location::Interior),
            (c(4.0, 1.0), Location::Boundary),
            (c(0.0, 0.0), Location::Boundary),
            (c(5.0, 2.0), Location::Exterior),
        ] {
            assert_eq!(locate(p, &sq), Ok(want), "{p:?}");
            assert_eq!(locate_in_areal(p, &area), Ok(want), "{p:?}");
            assert_eq!(locate(p, &sq).unwrap(), line_analysis::locate(p, &sq));
        }
        // A segment from inside to outside crosses one edge, named by that
        // edge and not constructed.
        let ev = segment_events(c(2.0, 2.0), c(6.0, 2.0), &area).unwrap();
        assert_eq!(ev.len(), 2);
        assert_eq!(ev[0].at, EventAt::Point(c(2.0, 2.0)));
        assert_eq!(ev[0].ahead, Location::Interior);
        assert_eq!(
            ev[1].at,
            EventAt::Crossing {
                a: c(4.0, 0.0),
                b: c(4.0, 4.0)
            }
        );
        assert_eq!(ev[1].ahead, Location::Exterior);
        assert_eq!(first_self_contact(sq.exterior()), None);
        assert_eq!(rings(&sq).count(), 1);
    }
}
