// src/prepared/line.rs

//! [`PreparedPolygonal::locate_on_line`]: many stored points on one
//! horizontal line, located exactly in one pass over the edges that meet
//! the line.

use arithmetica::Sign;

use super::{Kind, PreparedPolygonal};
use crate::error::GeometryError;
use crate::exact::kernel::orient;
use crate::index::EdgeRef;
use crate::topology::label::Location;
use crate::types::{Coord2d, Envelope};

impl PreparedPolygonal {
    /// The exact location of every point `(xs[i], y)`, in input order:
    /// `result[i] == self.locate(Coord2d::new(xs[i], y))` for every `i`.
    ///
    /// Built for many points on one row, such as the centers of a row of
    /// raster cells. The points cost one pass over the edges that meet the
    /// line, not one location each.
    ///
    /// `xs` may be in any order, and may repeat values. Input already in
    /// ascending order is detected and not sorted again.
    ///
    /// # Why the answers are `locate`'s
    ///
    /// No coordinate is computed. Every decision is one of the exact tests
    /// [`locate`](Self::locate) makes, taken for a run of points at a time:
    ///
    /// - **Only the edges whose closed `y` range contains `y` can matter.**
    ///   Any other edge's box misses every point on the line, so it can
    ///   neither hold one nor cross the `+x` ray from one.
    /// - **A non-horizontal edge `a → b`.** Along the line,
    ///   `orient(a, b, (x, y))` is linear in `x` with a nonzero slope. So
    ///   over the points in ascending `x` its sign runs: first the sign of
    ///   "left of the edge" (positive when `b.y > a.y`, else negative), then
    ///   zero, then the opposite sign. The edge meets the line at one
    ///   point, inside the edge's closed `x` range, so the points before
    ///   that range are all left of it and the points after it all right:
    ///   comparisons of stored values place them. Two binary searches with
    ///   the exact sign, over the points inside the range only, find the
    ///   zero run.
    ///   - **A point in the zero run is on the edge,** because its `y` lies
    ///     in the edge's `y` range, so the edge's box contains it. That is
    ///     `on_segment`, and the answer is `Boundary`.
    ///   - **If the edge straddles the line** under the half-open rule,
    ///     `(a.y > y) != (b.y > y)`, every point left of it has the edge
    ///     crossing its ray. That is `edge_crosses_ray`, and it toggles the
    ///     parity of the edge's ring for that prefix of the points.
    /// - **A horizontal edge on the line** holds exactly the points whose
    ///   `x` lies in its closed range, found by comparing stored values. It
    ///   straddles nothing.
    ///
    /// `Boundary` comes first, as in `locate`. Otherwise a point is
    /// `Interior` when some part has odd parity on its exterior and even
    /// parity on every hole. That is `locate`'s combination, computed here
    /// by one sweep over the prefix ends.
    ///
    /// # Cost
    ///
    /// One edge-index query for the line, or a scan for a mask built
    /// without an index. Then, for each of the `k` edges meeting the line,
    /// two binary searches over the points by stored value and about
    /// `log m` exact orientations, for `m` points inside the edge's `x`
    /// range (usually one or two). Then a sort of at most `k` prefix ends
    /// and a pass over the points. Unsorted input adds a sort of the
    /// points.
    ///
    /// # Errors
    ///
    /// [`GeometryError::NonFiniteCoordinate`], checking `y` first:
    ///
    /// - a non-finite `y` gives `coord_index: 0, axis: 1`;
    /// - a non-finite `xs[i]` gives `coord_index: i, axis: 0`.
    pub fn locate_on_line(&self, y: f64, xs: &[f64]) -> Result<Vec<Location>, GeometryError> {
        if !y.is_finite() {
            return Err(GeometryError::NonFiniteCoordinate {
                coord_index: 0,
                axis: 1,
            });
        }
        if let Some(i) = xs.iter().position(|x| !x.is_finite()) {
            return Err(GeometryError::NonFiniteCoordinate {
                coord_index: i,
                axis: 0,
            });
        }
        if xs.windows(2).all(|w| w[0] <= w[1]) {
            return Ok(self.locate_sorted(y, xs));
        }
        // Sort a permutation, locate in ascending order, and put the
        // answers back in input order. The values are finite, so
        // `total_cmp` orders them numerically (it also separates `-0.0`
        // from `0.0`, which is harmless: both give the same answers).
        let mut order: Vec<usize> = (0..xs.len()).collect();
        order.sort_unstable_by(|&i, &j| xs[i].total_cmp(&xs[j]));
        let sorted: Vec<f64> = order.iter().map(|&i| xs[i]).collect();
        let located = self.locate_sorted(y, &sorted);
        let mut out = vec![Location::Exterior; xs.len()];
        for (k, &i) in order.iter().enumerate() {
            out[i] = located[k];
        }
        Ok(out)
    }

    /// [`Self::locate_on_line`] for finite `xs` in ascending order.
    fn locate_sorted(&self, y: f64, xs: &[f64]) -> Vec<Location> {
        let n = xs.len();
        let parts = self.mask.polygons();
        if n == 0 || parts.is_empty() {
            return vec![Location::Exterior; n];
        }
        // A global number for each ring: `base[part] + ring`.
        let mut base = Vec::with_capacity(parts.len() + 1);
        base.push(0usize);
        for p in parts {
            base.push(base[base.len() - 1] + 1 + p.holes().len());
        }

        let mut boundary = vec![false; n];
        // `(end, ring)`: the ring's parity toggles for the points before
        // `end`.
        let mut toggles: Vec<(usize, usize)> = Vec::new();
        self.for_edges_on_row(y, |r, a, b| {
            place_edge(a, b, y, xs, &mut boundary, |end| {
                toggles.push((end, base[r.part] + r.ring));
            });
        });

        // Sweep from the right: before reading point `i`, apply every
        // toggle whose prefix contains `i` (`end > i`).
        toggles.sort_unstable_by_key(|t| core::cmp::Reverse(t.0));
        let mut odd = vec![false; base[parts.len()]];
        let mut ext_odd = vec![false; parts.len()];
        let mut odd_holes = vec![0usize; parts.len()];
        let mut inside_parts = 0usize;
        let part_of = |ring: usize| base.partition_point(|&b| b <= ring) - 1;
        let mut next = 0;
        let mut out = vec![Location::Exterior; n];
        for i in (0..n).rev() {
            while next < toggles.len() && toggles[next].0 > i {
                let ring = toggles[next].1;
                next += 1;
                let part = part_of(ring);
                let was = ext_odd[part] && odd_holes[part] == 0;
                odd[ring] = !odd[ring];
                if ring == base[part] {
                    ext_odd[part] = odd[ring];
                } else if odd[ring] {
                    odd_holes[part] += 1;
                } else {
                    odd_holes[part] -= 1;
                }
                let now = ext_odd[part] && odd_holes[part] == 0;
                if now && !was {
                    inside_parts += 1;
                } else if was && !now {
                    inside_parts -= 1;
                }
            }
            out[i] = if boundary[i] {
                Location::Boundary
            } else if inside_parts > 0 {
                Location::Interior
            } else {
                Location::Exterior
            };
        }
        out
    }

    /// Calls `visit` with every edge whose closed `y` range contains `y`,
    /// each once, in any order: through the edge index for an indexed mask,
    /// by a scan otherwise. The index returns every edge whose box meets
    /// the line's stretch across the mask's box, which is exactly these
    /// edges.
    fn for_edges_on_row(&self, y: f64, mut visit: impl FnMut(EdgeRef, Coord2d, Coord2d)) {
        let Some(env) = self.envelope else {
            return;
        };
        if y < env.min().y || y > env.max().y {
            return;
        }
        match &self.kind {
            Kind::Simple { edge_index, .. } => {
                let row = Envelope::new(Coord2d::new(env.min().x, y), Coord2d::new(env.max().x, y));
                let Ok(row) = row else {
                    return;
                };
                edge_index.candidates(&row, |r| {
                    let (a, b) = self.edge(r);
                    visit(r, a, b);
                });
            }
            Kind::Small(_) | Kind::Delegating(_) => {
                for (part, p) in self.mask.polygons().iter().enumerate() {
                    let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
                    for (ring, rg) in rings.enumerate() {
                        for (edge, (a, b)) in rg.edges().enumerate() {
                            if a.y.min(b.y) <= y && a.y.max(b.y) >= y {
                                visit(EdgeRef { part, ring, edge }, a, b);
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Places the edge `a → b`, whose closed `y` range contains `y`, against the
/// ascending points `xs` on the line: marks the points on it in `boundary`,
/// and, if it straddles the line, reports the end of the prefix of points
/// left of it through `toggle`.
///
/// Every comparison of coordinates here is exact by intent: stored values
/// compared as they are, never within a tolerance.
#[allow(clippy::float_cmp)]
fn place_edge(
    a: Coord2d,
    b: Coord2d,
    y: f64,
    xs: &[f64],
    boundary: &mut [bool],
    mut toggle: impl FnMut(usize),
) {
    if a.y == b.y {
        // Horizontal, so on the line: the points within its closed range.
        let (lo, hi) = if a.x <= b.x { (a.x, b.x) } else { (b.x, a.x) };
        let from = xs.partition_point(|&x| x < lo);
        let to = xs.partition_point(|&x| x <= hi);
        for flag in &mut boundary[from..to] {
            *flag = true;
        }
        return;
    }
    let left = if b.y > a.y {
        Sign::Positive
    } else {
        Sign::Negative
    };
    let sign = |x: f64| orient(a, b, Coord2d::new(x, y));
    // The edge meets the line at one point, whose `x` lies in the edge's
    // closed `x` range: every point before that range is left of it, and
    // every point after it is right, decided by comparing stored values.
    // Only the points inside the range need the sign, which runs `left`,
    // then zero, then the opposite, in ascending `x`.
    let (lo, hi) = if a.x <= b.x { (a.x, b.x) } else { (b.x, a.x) };
    let first = xs.partition_point(|&x| x < lo);
    let last = first + xs[first..].partition_point(|&x| x <= hi);
    let zero_from = first + xs[first..last].partition_point(|&x| sign(x) == left);
    let zero_to = zero_from + xs[zero_from..last].partition_point(|&x| sign(x) == Sign::Zero);
    for flag in &mut boundary[zero_from..zero_to] {
        *flag = true;
    }
    if (a.y > y) != (b.y > y) {
        toggle(zero_from);
    }
}
