// src/topology/exact_point.rs

//! Location of an **exact** point relative to an areal operand, with a
//! direction.
//!
//! The point is a stored coordinate or a constructed rational point such as
//! a crossing, and the direction gives the location of the points
//! immediately beside a ray leaving it.
//!
//! # Exactness
//!
//! Every test is one of three exact primitives of `arithmetica`, on the
//! stored endpoints of the mask's edges and the point itself:
//!
//! - `orient2d` of three stored points, or `orient2d_rational` of two stored
//!   points and the rational one;
//! - `rational_compare_along` on an axis, which compares one coordinate of
//!   the rational point with a stored value, exactly and without division;
//! - `rational_coords_equal`.
//!
//! For a stored point the rational predicates are not used at all: the
//! stored branch is the `f64` code of `line_analysis` and `kernel`, called
//! as it is, so a stored point gets exactly the answer the unprepared
//! analyses give it. The rational branch is tested against the stored one
//! on rational points that happen to be representable (`tests`): the two
//! must agree, and do.
//!
//! # The wedge, at a rational apex
//!
//! `line_analysis::location_ahead_dir` orders the directions leaving a point
//! by `upper_half` (is the direction's polar angle in `[0, π)`) and then by
//! `orient2d` within a half. With a rational apex `x` and stored `u`, `v`:
//! `upper_half(x, w)` compares `w.y` with `x.y`, then `w.x` with `x.x`, both
//! by `rational_compare_along` on an axis; and `orient(x, u, v)` is
//! `orient(u, v, x)` by cyclic symmetry of the determinant, which
//! `orient2d_rational` computes exactly. So the comparator is the same
//! function on the same real numbers.

use core::cmp::Ordering;

use arithmetica::{RationalCoord2d, Sign, orient2d, orient2d_rational, rational_compare_along};

use crate::exact::kernel::{edge_crosses_ray, on_segment};
use crate::topology::areal::{AreaOperand, Window};
use crate::topology::label::Location;
use crate::types::Coord2d;

/// A point with an exact position: stored, or a constructed rational.
#[non_exhaustive]
#[derive(Clone, Copy, Debug)]
pub enum ExactPoint<'a> {
    /// A stored coordinate.
    Stored(Coord2d),
    /// A constructed point, exact (a crossing of two stored segments).
    Rational(&'a RationalCoord2d),
}

const AXIS_X: ([f64; 2], [f64; 2]) = ([0.0, 0.0], [1.0, 0.0]);
const AXIS_Y: ([f64; 2], [f64; 2]) = ([0.0, 0.0], [0.0, 1.0]);

/// `rational_compare_along` on an axis between a rational point and a
/// stored value lifted to rational. The axis endpoints are literals, so
/// the call cannot fail; a failure would be a broken invariant and is
/// reported as `Equal`-free `None`.
fn cmp_coord(r: &RationalCoord2d, axis: ([f64; 2], [f64; 2]), stored: Coord2d) -> Option<Ordering> {
    let s = RationalCoord2d::from_f64([stored.x, stored.y]).ok()?;
    rational_compare_along(axis.0, axis.1, r, &s).ok()
}

impl ExactPoint<'_> {
    /// The closed box the point lies in: itself for a stored point; the
    /// representable neighbors of each coordinate for a rational one.
    pub(crate) fn window(&self) -> Window {
        match self {
            ExactPoint::Stored(c) => (*c, *c),
            ExactPoint::Rational(r) => {
                let [bx, by] = r.f64_bracket();
                (Coord2d::new(bx.lo, by.lo), Coord2d::new(bx.hi, by.hi))
            }
        }
    }

    /// Exact `on_segment` of the point on the stored segment `[a, b]`.
    fn on_segment(&self, a: Coord2d, b: Coord2d) -> bool {
        match self {
            ExactPoint::Stored(p) => on_segment(*p, a, b),
            ExactPoint::Rational(r) => {
                // Collinear, and within the segment's box, both exact.
                let collinear =
                    orient2d_rational([a.x, a.y], [b.x, b.y], r).ok() == Some(Sign::Zero);
                if !collinear {
                    return false;
                }
                // Within the segment's closed box, coordinate by coordinate.
                let within = |axis, lo: Coord2d, hi: Coord2d| {
                    cmp_coord(r, axis, lo) != Some(Ordering::Less)
                        && cmp_coord(r, axis, hi) != Some(Ordering::Greater)
                };
                within(
                    AXIS_X,
                    Coord2d::new(a.x.min(b.x), 0.0),
                    Coord2d::new(a.x.max(b.x), 0.0),
                ) && within(
                    AXIS_Y,
                    Coord2d::new(0.0, a.y.min(b.y)),
                    Coord2d::new(0.0, a.y.max(b.y)),
                )
            }
        }
    }

    /// Whether the stored edge `a → b` crosses the `+x` ray from the point,
    /// under `kernel::edge_crosses_ray`'s rule.
    fn edge_crosses_ray(&self, a: Coord2d, b: Coord2d) -> bool {
        match self {
            ExactPoint::Stored(p) => edge_crosses_ray(a, b, *p),
            ExactPoint::Rational(r) => {
                // `(a.y > p.y) != (b.y > p.y)`, then the side.
                let above = |c: Coord2d| cmp_coord(r, AXIS_Y, c) == Some(Ordering::Less);
                if above(a) == above(b) {
                    return false;
                }
                let s = orient2d_rational([a.x, a.y], [b.x, b.y], r).ok();
                if b.y > a.y {
                    s == Some(Sign::Positive)
                } else {
                    s == Some(Sign::Negative)
                }
            }
        }
    }

    /// `line_analysis::upper_half` at this apex.
    // The stored branch is `line_analysis::upper_half` verbatim: `w.y == x.y`
    // compares two STORED coordinates (so `-0.0 == 0.0`, as there), not two
    // computed values, and the test is exact.
    #[allow(clippy::float_cmp)]
    fn upper_half(&self, w: Coord2d) -> bool {
        match self {
            ExactPoint::Stored(x) => w.y > x.y || (w.y == x.y && w.x > x.x),
            ExactPoint::Rational(r) => match cmp_coord(r, AXIS_Y, w) {
                Some(Ordering::Less) => true,
                Some(Ordering::Equal) => cmp_coord(r, AXIS_X, w) == Some(Ordering::Less),
                _ => false,
            },
        }
    }

    /// `line_analysis::cmp_angle` at this apex.
    fn cmp_angle(&self, u: Coord2d, v: Coord2d) -> Ordering {
        match (self.upper_half(u), self.upper_half(v)) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            _ => {
                let s = match self {
                    ExactPoint::Stored(x) => orient2d([x.x, x.y], [u.x, u.y], [v.x, v.y]).ok(),
                    // orient(x, u, v) = orient(u, v, x): the determinant is
                    // invariant under a cyclic permutation of its rows.
                    ExactPoint::Rational(r) => orient2d_rational([u.x, u.y], [v.x, v.y], r).ok(),
                };
                match s {
                    Some(Sign::Positive) => Ordering::Less,
                    Some(Sign::Negative) => Ordering::Greater,
                    _ => Ordering::Equal,
                }
            }
        }
    }

    /// The point's `x` is beyond `max_x` (so the `+x` ray meets no edge of
    /// an operand whose envelope ends there), or its `y` is outside
    /// `[min_y, max_y]`.
    fn outside_ray_reach(&self, env: (f64, f64, f64, f64)) -> bool {
        match self {
            ExactPoint::Stored(p) => p.x > env.2 || p.y < env.1 || p.y > env.3,
            ExactPoint::Rational(r) => {
                cmp_coord(r, AXIS_X, Coord2d::new(env.2, 0.0)) == Some(Ordering::Greater)
                    || cmp_coord(r, AXIS_Y, Coord2d::new(0.0, env.1)) == Some(Ordering::Less)
                    || cmp_coord(r, AXIS_Y, Coord2d::new(0.0, env.3)) == Some(Ordering::Greater)
            }
        }
    }
}

/// A ray leaving the apex along a ring edge (as `line_analysis::Ray`).
#[derive(Clone, Copy)]
struct Ray {
    to: Coord2d,
    is_out: bool,
}

/// The rays of every ring edge containing the apex.
fn incident_rays<A: AreaOperand + ?Sized>(x: &ExactPoint<'_>, area: &A) -> Vec<Ray> {
    let mut out = Vec::new();
    let stored_eq = |c: Coord2d| match x {
        ExactPoint::Stored(p) => *p == c,
        ExactPoint::Rational(r) => RationalCoord2d::from_f64([c.x, c.y])
            .ok()
            .is_some_and(|s| arithmetica::rational_coords_equal(r, &s)),
    };
    area.for_edges_near(&x.window(), &mut |a, b| {
        if stored_eq(a) {
            out.push(Ray {
                to: b,
                is_out: true,
            });
        } else if stored_eq(b) {
            out.push(Ray {
                to: a,
                is_out: false,
            });
        } else if x.on_segment(a, b) {
            out.push(Ray {
                to: b,
                is_out: true,
            });
            out.push(Ray {
                to: a,
                is_out: false,
            });
        }
        true
    });
    out
}

/// Exact location of the point in the operand: on a ring is `Boundary`;
/// inside some part is `Interior`; otherwise `Exterior`. The ray parity
/// per ring over the edges near the `+x` ray, combined per part as
/// `inside_part_off_boundary` does.
pub(crate) fn locate<A: AreaOperand + ?Sized>(x: &ExactPoint<'_>, area: &A) -> Location {
    let mut on = false;
    area.for_edges_near(&x.window(), &mut |a, b| {
        on = x.on_segment(a, b);
        !on
    });
    if on {
        return Location::Boundary;
    }
    let Some(env) = area.envelope() else {
        return Location::Exterior;
    };
    if x.outside_ray_reach(env) {
        return Location::Exterior;
    }
    let (lo, hi) = x.window();
    let ray: Window = (lo, Coord2d::new(env.2, hi.y));
    let mut parity: Vec<(usize, usize)> = Vec::new();
    area.for_edges_with_ref_near(&ray, &mut |part, ring, a, b| {
        if x.edge_crosses_ray(a, b) {
            parity.push((part, ring));
        }
        true
    });
    parity.sort_unstable();
    let mut odd: Vec<(usize, usize)> = Vec::new();
    let mut i = 0;
    while i < parity.len() {
        let mut j = i;
        while j < parity.len() && parity[j] == parity[i] {
            j += 1;
        }
        if (j - i) % 2 == 1 {
            odd.push(parity[i]);
        }
        i = j;
    }
    let inside = odd
        .iter()
        .filter(|(_, ring)| *ring == 0)
        .any(|&(part, _)| !odd.iter().any(|&(q, ring)| q == part && ring > 0));
    if inside {
        Location::Interior
    } else {
        Location::Exterior
    }
}

/// The location of the points immediately to the **left** of the ray from
/// `x` toward `toward`, just after `x`: the winding side a windowed
/// arrangement asks for (the face to the left of a half-edge leaving `x`).
///
/// - `x` off every ring: the location of `x` itself (its neighborhood is
///   uniform).
/// - `x` on a ring, the ray not along a ring edge: the wedge the ray enters
///   (`location_ahead_dir`), which is the location of a neighborhood of the
///   open stretch, left side included.
/// - `x` on a ring and the ray along a ring edge: `Interior` exactly when
///   that edge runs in the ray's direction, because the region lies to the
///   left of every directed ring edge; `Exterior` otherwise.
pub(crate) fn location_left_of<A: AreaOperand + ?Sized>(
    x: &ExactPoint<'_>,
    toward: Coord2d,
    area: &A,
) -> Location {
    let rays = incident_rays(x, area);
    if rays.is_empty() {
        return locate(x, area);
    }
    if let Some(r) = rays
        .iter()
        .find(|r| x.cmp_angle(r.to, toward) == Ordering::Equal)
    {
        return if r.is_out {
            Location::Interior
        } else {
            Location::Exterior
        };
    }
    // `location_ahead_dir`'s wedge: the ray with the greatest angle below
    // the probe, wrapping to the greatest overall.
    let mut below: Option<&Ray> = None;
    let mut greatest: Option<&Ray> = None;
    for r in &rays {
        if x.cmp_angle(r.to, toward) == Ordering::Less
            && below.is_none_or(|b| x.cmp_angle(b.to, r.to) == Ordering::Less)
        {
            below = Some(r);
        }
        if greatest.is_none_or(|g| x.cmp_angle(g.to, r.to) == Ordering::Less) {
            greatest = Some(r);
        }
    }
    match below.or(greatest) {
        Some(r) if r.is_out => Location::Interior,
        _ => Location::Exterior,
    }
}
