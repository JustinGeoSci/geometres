// src/topology/line_analysis.rs

//! Exact DE-9IM computation for line × polygon and line × line pairs, where a
//! "line" is a `LineString` or a `MultiLineString` ([`LineSet`]).
//!
//! Both entry points ([`lines_polygon_matrix`], [`lines_lines_matrix`]) build an
//! [`IntersectionMatrix`] using only exact geometric tests: every sign
//! decision is made by `arithmetica::orient2d` on the original `f64`
//! coordinates, and every ordering decision compares original coordinates.
//! No intersection point is ever constructed, so no rounding enters any
//! classification.
//!
//! # Why not sample segment midpoints
//!
//! A segment of a line can cross a polygon boundary, so the segment is in the
//! interior over one stretch and in the exterior over another. The midpoint
//! of the segment lies in only one of those stretches, and testing it alone
//! reports the other as empty: a plausible-looking wrong matrix. The midpoint
//! itself is also not representable in `f64` in general. The approach here
//! classifies every sub-interval of every segment instead.
//!
//! # Sub-interval classification without constructing split points
//!
//! Cut a segment `s = (p, q)` at every point where it meets the polygon
//! boundary. Each resulting open sub-interval lies wholly in the interior,
//! the boundary, or the exterior of the polygon. Every sub-interval starts at
//! an *event point*, and the location of the points immediately after an
//! event, in the direction of `q`, is the sub-interval's location. So it is
//! enough to classify "just after" each event point:
//!
//! - **Start `p`.** Its location, or the wedge it faces if `p` is on the
//!   boundary.
//! - **A polygon vertex strictly inside `s`.** Classified by an angular wedge
//!   test at the vertex (see [`location_ahead`]).
//! - **A proper crossing of a ring edge `e` with `s`.** The crossing point is
//!   not representable, but the side after it is: the polygon interior is
//!   always to the left of a directed ring edge, so the stretch toward `q` is
//!   interior exactly when `q` is strictly left of `e`.
//!
//! The order of the events along `s` is never needed, because each
//! sub-interval is reached from its own left event.
//!
//! # Orientation convention relied on
//!
//! `Polygon::new` requires counterclockwise exterior rings and clockwise
//! holes, so the polygon interior is on the left of every directed ring
//! edge. The wedge test and the crossing test both depend on this.
//!
//! # Complexity
//!
//! `O(n·m)` for `n` line segments and `m` ring edges per event, which is
//! adequate for validation-grade use and matches the naive `is_simple` check.

use std::cmp::Ordering;

use arithmetica::Sign;

use crate::exact::kernel::{on_segment, orient};
use crate::topology::areal::{AreaOperand, Areal, point_window, segment_window};
use crate::topology::de9im::IntersectionMatrix;
use crate::topology::label::{Dimension, Location};
use crate::types::{Coord2d, LineString, LinearRing, Polygon};

// ===========================================================================
// Exact primitives
// ===========================================================================
//
// `orient`, `on_segment` and `ring_contains_off_boundary` live in
// `exact::kernel`, so that `Polygon::new` and `LinearRing::is_simple` can use
// them without `types` depending on `topology`. They are imported above.

/// Iterate over the exterior ring followed by every hole.
pub fn rings(poly: &Polygon) -> impl Iterator<Item = &LinearRing> {
    std::iter::once(poly.exterior()).chain(poly.holes().iter())
}

/// Exact location of `p` relative to `poly`.
///
/// On a ring is `Boundary`; inside the exterior and outside every hole is
/// `Interior`.
pub fn locate(p: Coord2d, poly: &Polygon) -> Location {
    Areal::single(poly).locate(p)
}

// ===========================================================================
// Angular wedge test
// ===========================================================================

/// A ray leaving an *event point* along a polygon ring edge.
#[derive(Clone, Copy)]
struct Ray {
    /// A point on the ray other than the event point.
    to: Coord2d,
    /// `true` if the ring edge is directed away from the event point.
    /// The polygon interior is on the left of an *outgoing* ray and on the
    /// right of an *incoming* one.
    is_out: bool,
}

/// All ring-edge rays that leave `x`, for every ring edge containing `x`.
///
/// An edge with `x` strictly inside it contributes two rays (one incoming
/// toward its start, one outgoing toward its end).
fn incident_rays<A: AreaOperand + ?Sized>(x: Coord2d, area: &A) -> Vec<Ray> {
    let mut out = Vec::new();
    // An edge containing `x` has `x` in its box.
    area.for_edges_near(&point_window(x), &mut |a, b| {
        if a == x {
            out.push(Ray {
                to: b,
                is_out: true,
            });
        } else if b == x {
            out.push(Ray {
                to: a,
                is_out: false,
            });
        } else if on_segment(x, a, b) {
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

/// `true` if the direction `x → w` has polar angle in `[0, π)`.
///
/// Compares coordinates only, so it is exact.
// `w.y == x.y` compares two STORED input coordinates, not two computed
// values: two f64s are either the same value or they are not, so the test is
// exact. `float_cmp` cannot tell this apart from comparing arithmetic results,
// which is the hazard it exists for.
#[allow(clippy::float_cmp)]
fn upper_half(x: Coord2d, w: Coord2d) -> bool {
    w.y > x.y || (w.y == x.y && w.x > x.x)
}

/// Exact counterclockwise ordering of the directions `x → u` and `x → v`
/// over `[0, 2π)`. `Equal` means the two rays point in the same direction.
fn cmp_angle(x: Coord2d, u: Coord2d, v: Coord2d) -> Ordering {
    match (upper_half(x, u), upper_half(x, v)) {
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        // Same half-plane: opposite directions cannot occur, so a zero
        // orientation means the same direction.
        _ => match orient(x, u, v) {
            Sign::Positive => Ordering::Less,
            Sign::Negative => Ordering::Greater,
            Sign::Zero => Ordering::Equal,
        },
    }
}

/// Location of the points immediately after `x` in the direction of `toward`.
///
/// `rays` must be the non-empty set from [`incident_rays`] (so `x` is on the
/// polygon boundary). The rays divide the neighborhood of `x` into wedges.
/// The wedge counterclockwise of an *outgoing* ray is interior (interior is
/// on its left); the wedge counterclockwise of an *incoming* ray is exterior.
/// If the probe direction coincides with a ray, the stretch lies along the
/// boundary.
#[cfg(test)]
fn location_ahead(x: Coord2d, toward: Coord2d, rays: &[Ray]) -> Location {
    location_ahead_dir(x, toward, rays).0
}

/// As [`location_ahead`], and when the stretch lies along the boundary also
/// reports whether the polygon edge runs in the probe direction (`true`) or
/// against it (`false`).
fn location_ahead_dir(x: Coord2d, toward: Coord2d, rays: &[Ray]) -> (Location, Option<bool>) {
    if let Some(r) = rays
        .iter()
        .find(|r| cmp_angle(x, r.to, toward) == Ordering::Equal)
    {
        return (Location::Boundary, Some(r.is_out));
    }
    // The wedge containing the probe starts at the ray with the greatest
    // angle below the probe, wrapping to the greatest overall if none is below.
    let mut below: Option<&Ray> = None;
    let mut greatest: Option<&Ray> = None;
    for r in rays {
        if cmp_angle(x, r.to, toward) == Ordering::Less
            && below.is_none_or(|b| cmp_angle(x, b.to, r.to) == Ordering::Less)
        {
            below = Some(r);
        }
        if greatest.is_none_or(|g| cmp_angle(x, g.to, r.to) == Ordering::Less) {
            greatest = Some(r);
        }
    }
    match below.or(greatest) {
        Some(r) if r.is_out => (Location::Interior, None),
        _ => (Location::Exterior, None),
    }
}

// ===========================================================================
// Coverage of a segment by collinear segments
// ===========================================================================

/// `true` if the segment `[a, b]` is covered end to end by the segments in
/// `segments` that are collinear with it.
///
/// Collinearity is exact. Along a non-vertical line the `x` coordinate is a
/// strictly monotone parameter (the `y` coordinate for a vertical line), so
/// the overlap intervals and the gap test compare original coordinates, with
/// no division and no rounding. Segments may come from several line strings
/// and may overlap one another; the interval sweep handles both.
// `a.x == b.x` is a verticality test on stored endpoint coordinates: the
// segment is vertical exactly when the two x values are the same f64. Exact,
// not an approximate comparison.
#[allow(clippy::float_cmp)]
fn covered_by_segments(
    a: Coord2d,
    b: Coord2d,
    segments: impl Iterator<Item = (Coord2d, Coord2d)>,
) -> bool {
    let along = |c: Coord2d| if a.x == b.x { c.y } else { c.x };
    let (lo, hi) = (along(a).min(along(b)), along(a).max(along(b)));

    let mut intervals: Vec<(f64, f64)> = Vec::new();
    for (c, d) in segments {
        if orient(a, b, c) != Sign::Zero || orient(a, b, d) != Sign::Zero {
            continue;
        }
        let ilo = along(c).min(along(d)).max(lo);
        let ihi = along(c).max(along(d)).min(hi);
        if ilo < ihi {
            intervals.push((ilo, ihi));
        }
    }
    intervals.sort_by(|x, y| x.0.total_cmp(&y.0));

    let mut reach = lo;
    for (ilo, ihi) in intervals {
        if ilo > reach {
            return false;
        }
        reach = reach.max(ihi);
        if reach >= hi {
            return true;
        }
    }
    reach >= hi
}

// ===========================================================================
// Linear geometry as a set of line strings
// ===========================================================================

/// The boundary of a set of line strings under the OGC mod-2 rule.
///
/// Every component contributes its start and its end. A point is in the
/// boundary exactly when it occurs an odd number of times among all those
/// endpoints. For one open line this gives its two endpoints; for a closed
/// line the start and end coincide, count twice, and cancel, so the boundary
/// is empty; for several lines, an endpoint shared by two of them
/// is an interior point of the union, not a boundary point.
///
/// Sorting uses `total_cmp` on `x + 0.0` and `y + 0.0`. The `+ 0.0` maps
/// `-0.0` to `+0.0`, so coordinates that compare equal under `==` are adjacent
/// after the sort and are counted together.
fn mod2_boundary(parts: &[&LineString]) -> Vec<Coord2d> {
    let mut ends: Vec<Coord2d> = parts.iter().flat_map(|l| [l.start(), l.end()]).collect();
    ends.sort_by(|a, b| {
        (a.x + 0.0)
            .total_cmp(&(b.x + 0.0))
            .then((a.y + 0.0).total_cmp(&(b.y + 0.0)))
    });
    let mut boundary = Vec::new();
    let mut i = 0;
    while i < ends.len() {
        let mut j = i + 1;
        while j < ends.len() && ends[j] == ends[i] {
            j += 1;
        }
        if (j - i) % 2 == 1 {
            boundary.push(ends[i]);
        }
        i = j;
    }
    boundary
}

/// A linear geometry (`LineString` or non-empty `MultiLineString`) viewed as
/// a set of line strings plus its OGC mod-2 boundary.
///
/// The point set is the union of the components. Components may touch,
/// cross, or overlap one another; no collection-level invariant is assumed,
/// because `MultiLineString::new` enforces none. Its interior is that union
/// minus the boundary points.
pub(crate) struct LineSet<'a> {
    parts: Vec<&'a LineString>,
    boundary: Vec<Coord2d>,
}

impl<'a> LineSet<'a> {
    /// A set holding one line string.
    pub(crate) fn single(line: &'a LineString) -> Self {
        Self::new(vec![line])
    }

    /// A set of line strings; computes the mod-2 boundary once.
    pub(crate) fn new(parts: Vec<&'a LineString>) -> Self {
        let boundary = mod2_boundary(&parts);
        Self { parts, boundary }
    }

    /// The boundary points (mod-2 rule).
    pub(crate) fn boundary(&self) -> &[Coord2d] {
        &self.boundary
    }

    /// Every segment of every component.
    pub(crate) fn segments(&self) -> impl Iterator<Item = (Coord2d, Coord2d)> + '_ {
        self.parts
            .iter()
            .flat_map(|l| l.coords().as_slice().windows(2).map(|w| (w[0], w[1])))
    }

    /// The box over every part, from the parts' cached boxes; `None` for
    /// no parts.
    pub(crate) fn bbox(&self) -> Option<crate::types::Bbox> {
        crate::types::bbox_fold(self.parts.iter().map(|l| l.bbox()))
    }

    /// Every vertex of every component.
    pub(crate) fn vertices(&self) -> impl Iterator<Item = Coord2d> + '_ {
        self.parts
            .iter()
            .flat_map(|l| l.coords().as_slice().iter().copied())
    }

    /// `true` if `x` lies on any segment (interior or boundary of the set).
    pub(crate) fn contains(&self, x: Coord2d) -> bool {
        self.segments().any(|(a, b)| on_segment(x, a, b))
    }

    /// Exact location of `x` relative to this linear geometry.
    pub(crate) fn locate(&self, x: Coord2d) -> Location {
        if self.boundary.contains(&x) {
            Location::Boundary
        } else if self.contains(x) {
            Location::Interior
        } else {
            Location::Exterior
        }
    }
}

// ===========================================================================
// Line × polygon
// ===========================================================================

/// Which locations the positive-length pieces of a line occupy.
///
/// Positive length only. Isolated contacts (points) are found separately
/// by the matrix functions, which need each contact point's own location.
// A summary of which locations a line's pieces occupy is a flag bundle by
// nature: the flags are independent observations from one traversal, and
// packing them into an enum would lose the fact that several hold at once.
#[allow(clippy::struct_excessive_bools)]
#[derive(Default)]
pub(crate) struct PieceSummary {
    /// Some positive-length piece lies in the polygon interior.
    pub(crate) interior: bool,
    /// Some positive-length piece lies along the polygon boundary.
    pub(crate) boundary: bool,
    /// Some positive-length piece lies in the polygon exterior.
    pub(crate) exterior: bool,
    /// Some piece along the boundary runs in the same direction as the
    /// polygon's ring edge there (so both interiors are on the same side).
    pub(crate) shared_same: bool,
    /// Some piece along the boundary runs against the polygon's ring edge.
    pub(crate) shared_opp: bool,
}

impl PieceSummary {
    fn add_dir(&mut self, (loc, dir): (Location, Option<bool>)) {
        self.add(loc);
        match dir {
            Some(true) => self.shared_same = true,
            Some(false) => self.shared_opp = true,
            None => {}
        }
    }

    fn add(&mut self, loc: Location) {
        match loc {
            Location::Interior => self.interior = true,
            Location::Boundary => self.boundary = true,
            Location::Exterior => self.exterior = true,
        }
    }
}

/// Classify every sub-interval of the segment `(p, q)` against `area`.
///
/// The events are those of [`segment_events`]: the start, every ring vertex
/// strictly inside the segment, and every proper crossing except one at a
/// ring vertex, whose wedge event already stands. Such a crossing is not
/// classified by its one edge's side test, which is wrong where a hole
/// touches the exterior at that point.
pub(crate) fn summarize_segment<A: AreaOperand + ?Sized>(
    p: Coord2d,
    q: Coord2d,
    area: &A,
    out: &mut PieceSummary,
) {
    for ev in segment_events(p, q, area) {
        match ev.at {
            EventAt::Point(x) => {
                if area.locate(x) == Location::Boundary {
                    out.add_dir(location_ahead_dir(x, q, &incident_rays(x, area)));
                } else {
                    out.add(ev.ahead);
                }
            }
            EventAt::Crossing { .. } => out.add(ev.ahead),
        }
    }
}

// ===========================================================================
// Ordered events, for clipping
// ===========================================================================

/// Where an event on a segment lies.
///
/// Either a stored input coordinate (the segment's start, or a polygon
/// vertex lying on the segment) or the proper crossing of the segment with
/// one ring edge, named by that edge. A crossing is named here rather than
/// constructed: classification never needs its coordinates, and the
/// caller that does (clipping, to emit geometry) constructs it exactly.
///
/// The two variants are exhaustive, and the enum is deliberately not
/// `#[non_exhaustive]`: a point where the location along a segment against
/// a polygon's boundary can change is a stored vertex (of the segment or of
/// the boundary, which is where every collinear overlap begins and ends) or a
/// proper crossing with one boundary edge. A caller matches both and needs
/// no arm for a kind of event that cannot exist.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EventAt {
    /// An input coordinate.
    Point(Coord2d),
    /// The proper crossing of the segment with the ring edge `a → b`.
    Crossing {
        /// Start of the ring edge, in ring order.
        a: Coord2d,
        /// End of the ring edge, in ring order.
        b: Coord2d,
    },
}

/// One event on a segment and the location of the stretch that follows it.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SegmentEvent {
    /// Where the event lies.
    pub at: EventAt,
    /// The location of the open stretch from this event to the next one (or
    /// to the segment's end), in the polygon: interior, boundary or exterior.
    pub ahead: Location,
}

/// Every event on the segment `p → q`, each with the location of the
/// stretch after it. These are the three event kinds [`summarize_segment`]
/// classifies, kept as a list instead of folded into flags.
///
/// # What the list is, exactly
///
/// - the start `p`, always first in the list;
/// - every polygon vertex lying strictly inside the segment, once each,
///   however many rings or edges meet there;
/// - every **proper** crossing of the segment with a ring edge, except one
///   that coincides with a polygon vertex already listed.
///
/// The events are distinct points: vertices are deduplicated by exact
/// coordinate equality; a proper crossing excludes `p` and `q` (both are
/// strictly off the edge) and excludes the edge's own endpoints; two proper
/// crossings at one point would be two ring edges meeting at a point interior
/// to both, which a valid polygon does not have. The one coincidence a valid
/// polygon allows, a hole vertex lying on an exterior edge where the
/// segment crosses that edge exactly at the vertex, is detected exactly
/// (`on_segment` of the vertex on both the segment and the edge) and the
/// crossing dropped, because the vertex's wedge test sees every ring at that
/// point and the crossing's side test sees only one edge.
///
/// The list is **not** ordered along the segment; ordering needs the
/// crossing coordinates, which this module does not construct.
///
/// Every location is decided by `orient2d` on input coordinates, as for the
/// DE-9IM matrix. Nothing here is approximate.
pub(crate) fn segment_events<A: AreaOperand + ?Sized>(
    p: Coord2d,
    q: Coord2d,
    area: &A,
) -> Vec<SegmentEvent> {
    let mut out = Vec::new();
    // A vertex on the segment, and an edge crossing it, lie in its box.
    let w = segment_window(p, q);

    // The start: on a ring, the wedge test decides which side the stretch
    // leaving `p` is on.
    let ahead = match area.locate(p) {
        Location::Boundary => location_ahead_dir(p, q, &incident_rays(p, area)).0,
        loc => loc,
    };
    out.push(SegmentEvent {
        at: EventAt::Point(p),
        ahead,
    });

    // Polygon vertices strictly inside the segment.
    let mut seen: Vec<Coord2d> = Vec::new();
    area.for_vertices_near(&w, &mut |v| {
        if v != p && v != q && on_segment(v, p, q) && !seen.contains(&v) {
            seen.push(v);
            let rays = incident_rays(v, area);
            out.push(SegmentEvent {
                at: EventAt::Point(v),
                ahead: location_ahead_dir(v, q, &rays).0,
            });
        }
        true
    });

    // Proper crossings with ring edges.
    area.for_edges_near(&w, &mut |a, b| {
        let (sa, sb) = (orient(p, q, a), orient(p, q, b));
        let (sp, sq) = (orient(a, b, p), orient(a, b, q));
        let opposite = |x: Sign, y: Sign| {
            (x == Sign::Positive && y == Sign::Negative)
                || (x == Sign::Negative && y == Sign::Positive)
        };
        if !(opposite(sa, sb) && opposite(sp, sq)) {
            return true;
        }
        if seen.iter().any(|&v| on_segment(v, a, b)) {
            // A listed vertex lies on this edge and on the segment; the
            // segment and the edge's line meet in one point, so the
            // crossing is that vertex. Its wedge event already stands.
            return true;
        }
        out.push(SegmentEvent {
            at: EventAt::Crossing { a, b },
            // Interior is left of every directed ring edge, so the
            // stretch after the crossing is interior exactly when `q` is
            // strictly left of the edge.
            ahead: if sq == Sign::Positive {
                Location::Interior
            } else {
                Location::Exterior
            },
        });
        true
    });
    out
}

/// Single-line convenience wrapper over [`lines_polygon_matrix`]; used by the unit tests.
#[cfg(test)]
pub(crate) fn line_polygon_matrix(l: &LineString, poly: &Polygon) -> IntersectionMatrix {
    lines_polygon_matrix(&LineSet::single(l), &Areal::single(poly))
}

/// The point cells a set of contact points contributes: `out[k]` is set when
/// some point of the line's interior (`row = 0`) or boundary (`row = 1`) has
/// areal location `k` (interior, boundary, exterior).
fn mark(out: &mut [[bool; 3]; 2], row: usize, loc: Location) {
    let k = match loc {
        Location::Interior => 0,
        Location::Boundary => 1,
        Location::Exterior => 2,
    };
    out[row][k] = true;
}

/// Computes the DE-9IM matrix of a linear geometry (A) against an areal one
/// (B, a polygon or a valid multipolygon), exactly.
///
/// # Cell derivation
///
/// With the OGC definitions (interior of a line is the line minus its mod-2
/// boundary; the boundary of B is every part's rings, topologically):
///
/// - **Stretches.** `II`, `IB`, `IE` are `Line` if a positive-length piece of
///   the line lies in B's interior, boundary or exterior.
/// - **Points.** Otherwise a cell is `Point` if an isolated point of the
///   line's interior has that location in B. The candidates are exactly the
///   points where something happens: a line vertex that is not a boundary
///   point, a ring vertex of B strictly inside a segment, and a proper
///   crossing of a segment with a ring edge (which lies on B's boundary,
///   unless a ring vertex is there, in which case that vertex is the
///   candidate). A line passing through the point where two parts of B
///   touch meets B's boundary there.
/// - `BI`, `BB`, `BE`: `Point` per boundary point of the line.
/// - `EI` and `EE`: always `Area`; a curve cannot fill a region.
/// - `EB`: `Line` unless the line covers every ring of B.
pub(crate) fn lines_polygon_matrix<A: AreaOperand + ?Sized>(
    lines: &LineSet<'_>,
    area: &A,
) -> IntersectionMatrix {
    let mut pieces = PieceSummary::default();
    for (p, q) in lines.segments() {
        summarize_segment(p, q, area, &mut pieces);
    }

    let bnd = lines.boundary();
    let mut pts = [[false; 3]; 2];
    // Line vertices: boundary points by row 1, the rest by row 0.
    for v in lines.vertices() {
        let row = usize::from(bnd.contains(&v));
        mark(&mut pts, row, area.locate(v));
    }
    // Ring vertices strictly inside a segment (they lie in the segment's
    // box). A line boundary point lying inside another component's segment
    // is a boundary point, already marked. Marking is idempotent, so a
    // vertex inside two segments is marked twice to the same effect.
    for (p, q) in lines.segments() {
        area.for_vertices_near(&segment_window(p, q), &mut |v| {
            if !bnd.contains(&v) && v != p && v != q && on_segment(v, p, q) {
                mark(&mut pts, 0, area.locate(v));
            }
            true
        });
    }
    // Proper crossings, unless a ring vertex or a line boundary point is
    // there (those are marked above, with their own locations). A crossing
    // edge, and a vertex on both the segment and the edge, lie in the
    // segment's box.
    for (p, q) in lines.segments() {
        let w = segment_window(p, q);
        area.for_edges_near(&w, &mut |a, b| {
            if !proper_crossing(p, q, a, b) {
                return true;
            }
            let mut at_vertex = bnd
                .iter()
                .any(|&v| on_segment(v, p, q) && on_segment(v, a, b));
            if !at_vertex {
                area.for_vertices_near(&w, &mut |v| {
                    at_vertex = on_segment(v, p, q) && on_segment(v, a, b);
                    !at_vertex
                });
            }
            if !at_vertex {
                mark(&mut pts, 0, Location::Boundary);
            }
            !pts[0][1]
        });
    }

    // The boundary is fully covered only if every ring edge is covered by
    // collinear line segments; an edge outside the lines' box cannot be.
    let lines_window = {
        let mut w: Option<(Coord2d, Coord2d)> = None;
        for v in lines.vertices() {
            w = Some(match w {
                None => (v, v),
                Some((lo, hi)) => (
                    Coord2d::new(lo.x.min(v.x), lo.y.min(v.y)),
                    Coord2d::new(hi.x.max(v.x), hi.y.max(v.y)),
                ),
            });
        }
        w.unwrap_or((Coord2d::new(0.0, 0.0), Coord2d::new(0.0, 0.0)))
    };
    let mut boundary_fully_covered = !area.has_edge_outside(&lines_window);
    if boundary_fully_covered {
        area.for_edges_near(&lines_window, &mut |a, b| {
            boundary_fully_covered = covered_by_segments(a, b, lines.segments());
            boundary_fully_covered
        });
    }

    let cell = |line: bool, point: bool, d: Dimension| {
        if line {
            d
        } else if point {
            Dimension::Point
        } else {
            Dimension::Empty
        }
    };
    let dim = |flag: bool, d: Dimension| if flag { d } else { Dimension::Empty };
    IntersectionMatrix {
        interior_interior: cell(pieces.interior, pts[0][0], Dimension::Line),
        interior_boundary: cell(pieces.boundary, pts[0][1], Dimension::Line),
        interior_exterior: cell(pieces.exterior, pts[0][2], Dimension::Line),
        boundary_interior: dim(pts[1][0], Dimension::Point),
        boundary_boundary: dim(pts[1][1], Dimension::Point),
        boundary_exterior: dim(pts[1][2], Dimension::Point),
        exterior_interior: Dimension::Area,
        exterior_boundary: dim(!boundary_fully_covered, Dimension::Line),
        exterior_exterior: Dimension::Area,
    }
}

// ===========================================================================
// Line × line
// ===========================================================================

/// `true` if the segments `(a, b)` and `(c, d)` are collinear and overlap in
/// a segment of positive length.
// `a.x == b.x` is the same exact verticality test as in `covered_by_segments`.
#[allow(clippy::float_cmp)]
pub(crate) fn collinear_overlap(a: Coord2d, b: Coord2d, c: Coord2d, d: Coord2d) -> bool {
    if orient(a, b, c) != Sign::Zero || orient(a, b, d) != Sign::Zero {
        return false;
    }
    let along = |p: Coord2d| if a.x == b.x { p.y } else { p.x };
    let lo = along(a).min(along(b)).max(along(c).min(along(d)));
    let hi = along(a).max(along(b)).min(along(c).max(along(d)));
    lo < hi
}

/// `true` if segments `(a, b)` and `(c, d)` cross at a single point that is
/// interior to both.
pub(crate) fn proper_crossing(a: Coord2d, b: Coord2d, c: Coord2d, d: Coord2d) -> bool {
    let opposite = |x: Sign, y: Sign| {
        (x == Sign::Positive && y == Sign::Negative) || (x == Sign::Negative && y == Sign::Positive)
    };
    opposite(orient(a, b, c), orient(a, b, d)) && opposite(orient(c, d, a), orient(c, d, b))
}

/// Single-line convenience wrapper over [`lines_lines_matrix`]; used by the unit tests.
#[cfg(test)]
pub(crate) fn line_line_matrix(a: &LineString, b: &LineString) -> IntersectionMatrix {
    lines_lines_matrix(&LineSet::single(a), &LineSet::single(b))
}

/// Computes the DE-9IM matrix of line `a` (geometry A) against line `b`
/// (geometry B), exactly.
///
/// # Cell derivation
///
/// - `II`: `Line` if any pair of segments overlaps over positive length.
///   Otherwise `Point` if the interiors share an isolated point: a proper
///   crossing, or a vertex of one line lying on the other, in each case
///   excluding points that are endpoints of either line.
/// - `IB`, `BI`: `Point` if an endpoint of one line lies on the other line
///   without being an endpoint of that line.
/// - `BB`: `Point` if the lines share an endpoint.
/// - `BE`, `EB`: `Point` if an endpoint of one line is not on the other line.
/// - `IE`, `EI`: `Line` if some segment of one line is not covered by the
///   other. The uncovered part is relatively open, so it has positive length.
/// - `EE`: always `Area`.
///
/// A closed line has an empty boundary under the OGC mod-2 rule, so its
/// `B*` / `*B` cells are empty.
///
/// # Several components
///
/// The derivation applies to the union of each side's components with the
/// mod-2 boundary of the whole side. A vertex that is a component endpoint
/// but not a mod-2 boundary point is an interior point of its side, so it can
/// produce an `II` point contact.
pub(crate) fn lines_lines_matrix(a: &LineSet<'_>, b: &LineSet<'_>) -> IntersectionMatrix {
    let (a_bnd, b_bnd) = (a.boundary(), b.boundary());
    let any_bnd = |x: Coord2d| a_bnd.contains(&x) || b_bnd.contains(&x);

    // ---- II ----
    let mut overlap = false;
    let mut ii_point = false;
    for (a0, a1) in a.segments() {
        for (b0, b1) in b.segments() {
            if collinear_overlap(a0, a1, b0, b1) {
                overlap = true;
            } else if proper_crossing(a0, a1, b0, b1) {
                // The crossing is interior to both segments; it is a line
                // boundary point only if an endpoint lies on both segments.
                let at_bnd = a_bnd
                    .iter()
                    .chain(b_bnd.iter())
                    .any(|&e| on_segment(e, a0, a1) && on_segment(e, b0, b1));
                ii_point |= !at_bnd;
            }
        }
    }
    ii_point |= a.vertices().any(|v| !any_bnd(v) && b.contains(v));
    ii_point |= b.vertices().any(|v| !any_bnd(v) && a.contains(v));

    // ---- Boundary cells ----
    let ib = b_bnd.iter().any(|&e| !a_bnd.contains(&e) && a.contains(e));
    let bi = a_bnd.iter().any(|&e| !b_bnd.contains(&e) && b.contains(e));
    let bb = a_bnd.iter().any(|e| b_bnd.contains(e));
    let be = a_bnd.iter().any(|&e| !b.contains(e));
    let eb = b_bnd.iter().any(|&e| !a.contains(e));

    // ---- Exterior cells ----
    let ie = a
        .segments()
        .any(|(p, q)| !covered_by_segments(p, q, b.segments()));
    let ei = b
        .segments()
        .any(|(p, q)| !covered_by_segments(p, q, a.segments()));

    let dim = |flag: bool, d: Dimension| if flag { d } else { Dimension::Empty };
    IntersectionMatrix {
        interior_interior: if overlap {
            Dimension::Line
        } else {
            dim(ii_point, Dimension::Point)
        },
        interior_boundary: dim(ib, Dimension::Point),
        interior_exterior: dim(ie, Dimension::Line),
        boundary_interior: dim(bi, Dimension::Point),
        boundary_boundary: dim(bb, Dimension::Point),
        boundary_exterior: dim(be, Dimension::Point),
        exterior_interior: dim(ei, Dimension::Line),
        exterior_boundary: dim(eb, Dimension::Point),
        exterior_exterior: Dimension::Area,
    }
}

// ===========================================================================
// Tests
// ===========================================================================

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
    use crate::types::CoordSeq;

    fn c(x: f64, y: f64) -> Coord2d {
        Coord2d::new(x, y)
    }

    fn ring(pts: &[(f64, f64)]) -> LinearRing {
        LinearRing::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
    }

    fn line(pts: &[(f64, f64)]) -> LineString {
        LineString::new(CoordSeq::new(pts.iter().map(|&(x, y)| c(x, y)).collect())).unwrap()
    }

    /// 4×4 square with a 2×2 hole in the middle.
    fn square_with_hole() -> Polygon {
        Polygon::new(
            ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            vec![ring(&[(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)])],
        )
        .unwrap()
    }

    fn square() -> Polygon {
        Polygon::new_simple(ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])).unwrap()
    }

    #[test]
    fn locate_classifies_all_regions() {
        let p = square_with_hole();
        assert_eq!(locate(c(0.5, 0.5), &p), Location::Interior);
        assert_eq!(
            locate(c(2.0, 2.0), &p),
            Location::Exterior,
            "inside the hole"
        );
        assert_eq!(locate(c(5.0, 5.0), &p), Location::Exterior);
        assert_eq!(
            locate(c(0.0, 0.0), &p),
            Location::Boundary,
            "exterior vertex"
        );
        assert_eq!(locate(c(2.0, 0.0), &p), Location::Boundary, "exterior edge");
        assert_eq!(locate(c(1.0, 2.0), &p), Location::Boundary, "hole edge");
        assert_eq!(locate(c(1.0, 1.0), &p), Location::Boundary, "hole vertex");
    }

    #[test]
    fn locate_ray_through_vertex_is_counted_once() {
        // A +x ray from (-1, 2) passes exactly through the vertex (0, 2) of
        // this chevron; the half-open rule must count it once, not twice or
        // zero times.
        let chevron = Polygon::new_simple(ring(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (2.0, 2.0),
        ]))
        .unwrap();
        assert_eq!(
            locate(c(1.0, 2.0), &chevron),
            Location::Exterior,
            "in the notch"
        );
        assert_eq!(locate(c(3.0, 2.0), &chevron), Location::Interior);
        assert_eq!(locate(c(-1.0, 2.0), &chevron), Location::Exterior);
    }

    #[test]
    fn cmp_angle_orders_counterclockwise_from_positive_x() {
        let o = c(0.0, 0.0);
        let dirs = [
            c(1.0, 0.0),
            c(1.0, 1.0),
            c(0.0, 1.0),
            c(-1.0, 0.0),
            c(0.0, -1.0),
            c(1.0, -1.0),
        ];
        for w in dirs.windows(2) {
            assert_eq!(cmp_angle(o, w[0], w[1]), Ordering::Less);
            assert_eq!(cmp_angle(o, w[1], w[0]), Ordering::Greater);
        }
        assert_eq!(cmp_angle(o, c(1.0, 1.0), c(5.0, 5.0)), Ordering::Equal);
    }

    #[test]
    fn location_ahead_at_convex_corner() {
        // CCW square, corner (0,0): out-ray to (4,0), in-ray to (0,4).
        let p = square();
        let rays = incident_rays(c(0.0, 0.0), &Areal::single(&p));
        assert_eq!(rays.len(), 2);
        let at = |t: (f64, f64)| location_ahead(c(0.0, 0.0), c(t.0, t.1), &rays);
        assert_eq!(at((1.0, 1.0)), Location::Interior);
        assert_eq!(at((-1.0, 1.0)), Location::Exterior);
        assert_eq!(at((-1.0, -1.0)), Location::Exterior);
        assert_eq!(at((1.0, -1.0)), Location::Exterior);
        assert_eq!(at((2.0, 0.0)), Location::Boundary, "along the out-edge");
        assert_eq!(at((0.0, 2.0)), Location::Boundary, "along the in-edge");
    }

    #[test]
    fn location_ahead_at_edge_interior_point() {
        let p = square();
        let rays = incident_rays(c(2.0, 0.0), &Areal::single(&p));
        assert_eq!(
            rays.len(),
            2,
            "an edge interior contributes both directions"
        );
        assert_eq!(
            location_ahead(c(2.0, 0.0), c(2.0, 1.0), &rays),
            Location::Interior
        );
        assert_eq!(
            location_ahead(c(2.0, 0.0), c(2.0, -1.0), &rays),
            Location::Exterior
        );
        assert_eq!(
            location_ahead(c(2.0, 0.0), c(3.0, 0.0), &rays),
            Location::Boundary
        );
        assert_eq!(
            location_ahead(c(2.0, 0.0), c(1.0, 0.0), &rays),
            Location::Boundary
        );
    }

    #[test]
    fn covered_by_segments_handles_gaps_and_joins() {
        let cover = |a: Coord2d, b: Coord2d, pl: &[Coord2d]| {
            covered_by_segments(a, b, pl.windows(2).map(|w| (w[0], w[1])))
        };
        let a = c(0.0, 0.0);
        let b = c(4.0, 0.0);
        let pl = |pts: &[(f64, f64)]| pts.iter().map(|&(x, y)| c(x, y)).collect::<Vec<_>>();
        // Two touching collinear pieces cover it.
        assert!(cover(a, b, &pl(&[(-1.0, 0.0), (2.0, 0.0), (5.0, 0.0)])));
        assert!(cover(a, b, &pl(&[(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)])));
        // A gap of 1 is not covered.
        assert!(!cover(
            a,
            b,
            &pl(&[(0.0, 0.0), (1.0, 0.0), (1.0, 5.0), (3.0, 0.0), (4.0, 0.0)])
        ));
        // A parallel line does not cover.
        assert!(!cover(a, b, &pl(&[(0.0, 1.0), (4.0, 1.0)])));
        // Vertical edge uses y as the parameter.
        assert!(cover(
            c(0.0, 0.0),
            c(0.0, 4.0),
            &pl(&[(0.0, 5.0), (0.0, -1.0)])
        ));
    }

    #[test]
    fn segment_crossing_the_boundary_reports_both_sides() {
        // The segment (-1,1)→(5,1) is exterior, interior, exterior. Sampling
        // only its midpoint (2,1) would see just the interior stretch and
        // report IE and IB as empty.
        let m = line_polygon_matrix(&line(&[(-1.0, 1.0), (5.0, 1.0)]), &square());
        assert_eq!(m.to_string(), "101FF0212");
    }

    #[test]
    fn line_through_hole_crosses_boundary_twice() {
        let m = line_polygon_matrix(&line(&[(-1.0, 2.0), (5.0, 2.0)]), &square_with_hole());
        // Enters the interior, crosses the hole boundary, runs through the
        // hole (exterior), re-enters the interior, and leaves: II, IB (points),
        // IE all present.
        assert_eq!(m.to_string(), "101FF0212");
    }

    #[test]
    fn line_along_boundary_edge_is_boundary_overlap() {
        let m = line_polygon_matrix(&line(&[(0.0, 0.0), (4.0, 0.0)]), &square());
        assert_eq!(m.interior_boundary, Dimension::Line);
        assert_eq!(m.boundary_boundary, Dimension::Point);
        assert_eq!(m.interior_interior, Dimension::Empty);
        assert_eq!(m.interior_exterior, Dimension::Empty);
        assert_eq!(
            m.exterior_boundary,
            Dimension::Line,
            "rest of the boundary is uncovered"
        );
    }

    #[test]
    fn closed_line_has_no_boundary_points() {
        let tri = line(&[(1.0, 1.0), (3.0, 1.0), (2.0, 3.0), (1.0, 1.0)]);
        let m = line_polygon_matrix(&tri, &square());
        assert_eq!(m.boundary_interior, Dimension::Empty);
        assert_eq!(m.boundary_boundary, Dimension::Empty);
        assert_eq!(m.boundary_exterior, Dimension::Empty);
        assert_eq!(m.interior_interior, Dimension::Line);
        // Closed line × closed line also has empty boundary rows/columns.
        let other = line(&[(5.0, 5.0), (7.0, 5.0), (6.0, 7.0), (5.0, 5.0)]);
        let mm = line_line_matrix(&tri, &other);
        assert_eq!(mm.to_string(), "FF1FFF1F2");
    }

    #[test]
    fn endpoint_touching_boundary_from_outside() {
        let m = line_polygon_matrix(&line(&[(-2.0, 2.0), (0.0, 2.0)]), &square());
        assert_eq!(m.boundary_boundary, Dimension::Point);
        assert_eq!(m.boundary_exterior, Dimension::Point);
        assert_eq!(
            m.interior_boundary,
            Dimension::Empty,
            "touch is at an endpoint, not interior"
        );
        assert_eq!(m.interior_exterior, Dimension::Line);
    }

    #[test]
    fn interior_vertex_touching_boundary_is_an_interior_point_contact() {
        // The middle vertex (4,2) sits on the square's right edge; the line
        // comes from outside and goes back outside.
        let m = line_polygon_matrix(&line(&[(6.0, 0.0), (4.0, 2.0), (6.0, 4.0)]), &square());
        assert_eq!(m.interior_boundary, Dimension::Point);
        assert_eq!(m.interior_interior, Dimension::Empty);
    }

    #[test]
    fn line_line_t_junction_and_overlap() {
        // Endpoint of b on the interior of a: the contact is interior(a) ∩
        // boundary(b). It is not in interior(b), so II stays empty.
        let a = line(&[(0.0, 0.0), (4.0, 0.0)]);
        let b = line(&[(2.0, 0.0), (2.0, 3.0)]);
        let m = line_line_matrix(&a, &b);
        assert_eq!(m.to_string(), "F01FF0102");
        // Collinear overlap: II = Line.
        let c2 = line(&[(1.0, 0.0), (6.0, 0.0)]);
        let m2 = line_line_matrix(&a, &c2);
        assert_eq!(m2.interior_interior, Dimension::Line);
        assert_eq!(m2.interior_exterior, Dimension::Line);
        assert_eq!(m2.exterior_interior, Dimension::Line);
    }

    #[test]
    fn identical_lines_have_empty_exteriors() {
        let a = line(&[(0.0, 0.0), (2.0, 2.0), (4.0, 0.0)]);
        let m = line_line_matrix(&a, &a.clone());
        assert_eq!(m.to_string(), "1FFF0FFF2");
    }

    #[test]
    fn crossing_at_a_line_endpoint_is_not_an_interior_contact() {
        // b's start (2,2) lies in the interior of a's second segment, and
        // a's own start is (2,2) too (a self-touching line): the point is in
        // the boundary of both lines, never in an interior.
        let a = line(&[(2.0, 2.0), (0.0, 0.0), (4.0, 0.0), (2.0, 2.0), (3.0, 3.0)]);
        let b = line(&[(2.0, 2.0), (2.0, 5.0)]);
        let m = line_line_matrix(&a, &b);
        assert_eq!(m.boundary_boundary, Dimension::Point);
    }
}
