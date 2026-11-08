// src/topology/collection.rs

//! Exact DE-9IM computation where either operand is a `GeometryCollection`.
//! The collection operand is restricted, as set out below.
//!
//! # The restriction, and the decomposition
//!
//! A collection is accepted when its polygonal leaves together satisfy the
//! multipolygon part rule: every pair meets at most at finitely many
//! points, decided exactly by the same `place_ring` test `MultiPolygon::new`
//! applies ([`first_overlapping_parts_ref`]). An accepted collection
//! decomposes exactly into three strata:
//!
//! - an **areal part `A`**, the polygonal leaves. The restriction makes them
//!   interior-disjoint, so the invariant of [`AreaOperand`] (the region on
//!   the left of every directed ring edge) holds and every existing areal
//!   proof applies unchanged ([`ArealParts`]);
//! - a **lineal part `L`**, every line leaf in one [`LineSet`], whose mod-2
//!   boundary is computed over the line leaves jointly and only them;
//! - a **punctal part `P`**, the stored coordinates of the point leaves.
//!
//! A collection outside the restriction is refused with
//! [`GeometryError::CollectionElementsOverlap`], naming both leaf paths.
//! Full union semantics over overlapping elements is future state; the
//! explicit route is to union the elements and relate the vouched result.
//!
//! # The location rule (the collection's interior and boundary)
//!
//! The collection's point set is `A ∪ L ∪ P`. The location of a point `x`,
//! in order ([`gc_locate`]):
//!
//! 1. interior of `A` → Interior;  2. boundary of `A` → Boundary;
//! 3. on `L`, not a mod-2 boundary point, outside `A` → Interior;
//! 4. a mod-2 boundary point of `L` outside `A` → Boundary (a `P` point
//!    there stays Boundary, which is a convention rather than a
//!    consequence of the topology);
//! 5. a `P` point outside `A ∪ L` → Interior;  6. otherwise Exterior.
//!
//! Rules 1 to 3, 5 and 6 are forced by point-set topology: a set of dimension
//! at most 1 can neither create nor destroy 2-dimensional interior, so
//! lower-dimensional material inside `A` is absorbed, and `L` off `A` keeps
//! the line's own structure. So
//!
//! ```text
//! I(GC) = Int(A) ⊔ (L ∖ Ā ∖ ∂₂L) ⊔ (P ∖ (Ā ∪ L))
//! B(GC) = Bnd(A) ⊔ (∂₂L ∖ Ā)            (∂₂L: the mod-2 boundary of L)
//! E(GC) = ℝ² ∖ (A ∪ L ∖ P)ᶜ-complement
//! ```
//!
//! # How the matrix is assembled ([`merged_matrix`])
//!
//! Both operands are viewed as strata ([`Strata`]; a simple operand has one
//! stratum). Each cell is the maximum dimension over its stratum-pair
//! contributions, which come from four exact sources:
//!
//! - **Area × area:** [`polygons_matrix`], unchanged. Its `II`, `IB`, `BI`,
//!   `BB` cells are collection cells directly. Its `IE` and `EI` cells are
//!   too, by a proved reduction: `Int(A) ∖ (other side)` and
//!   `Int(other) ∖ (this side)` keep their dimension when the closed,
//!   nowhere-dense `L ∪ P` is removed from an open set, because a 0- or
//!   1-dimensional set can neither empty an open set nor lower its
//!   dimension. Its `BE`/`EB` cells are **not** usable. A boundary stretch
//!   outside the other area can still be covered by the other side's lines,
//!   which is 1-dimensional against 1-dimensional and needs the walk.
//! - **The jointly labeled walk** ([`walk_segment`]): every segment of a
//!   lineal stratum, and every ring edge (for the boundary-coverage cells),
//!   classified against several operands *at once*. `segment_events` gives
//!   each areal operand's events with the location ahead of each;
//!   collinear-overlap sweeps give each lineal operand's. Merging them
//!   needs their order along the segment, decided by predicates on input
//!   coordinates only: stored-coordinate comparison for two input points,
//!   and `arithmetica`'s exact rational comparison
//!   (`rational_compare_along`, `from_segment_intersection`) where a proper
//!   crossing is involved. That is a sign decision, and **no coordinate is
//!   constructed or returned**. Each open sub-interval between consecutive distinct
//!   events has one location per operand (constant on the interval), so a
//!   stretch's joint label is exact; an event point's own label comes from
//!   its event, and a representable candidate point is simply re-located
//!   directly.
//! - **Masking:** a stretch or point of `L` inside its own side's `Ā` is
//!   absorbed. It is skipped because its points belong to the `A` strata,
//!   whose pairs the other sources measure (a subset cannot exceed its
//!   superset's dimension).
//! - **Point strata and mod-2 boundary points:** located directly with
//!   [`gc_locate`] against the other side.
//!
//! The exterior row and column are never computed per part. A stretch
//! counts toward `E` only when its labels against **every** stratum of the
//! other side are exterior, and an isolated candidate additionally checks
//! the other side's `P`. `BE`/`EB` need no point case at all: `Bnd(A)` minus
//! the other side's closed point set is relatively open in the boundary
//! curves, so it is a positive-length stretch or empty. Isolated boundary
//! survivors cannot exist.
//!
//! # Exactness
//!
//! Every decision is `orient2d` on stored coordinates, a stored-coordinate
//! comparison, or an exact rational sign (`arithmetica::rational`). The
//! rational comparisons cannot fail on validated geometry, since every
//! coordinate is finite by construction. The fallbacks that
//! `deny(unwrap_used)` requires for these unreachable error arms return a
//! neutral value rather than panic.

use std::cmp::Ordering;

use arithmetica::rational::{RationalCoord2d, rational_compare_along};

use crate::GeometryError;
use crate::exact::kernel::on_segment;
use crate::topology::areal::{AreaOperand, Window, inside_part_off_boundary, on_part_boundary};
use crate::topology::de9im::{IntersectionMatrix, Operand};
use crate::topology::label::{Dimension, Location};
use crate::topology::line_analysis::{
    EventAt, LineSet, collinear_overlap, proper_crossing, segment_events,
};
use crate::topology::polygon_analysis::polygons_matrix;
use crate::types::{
    Coord2d, Geometry, GeometryCollection, LineString, LinearRing, Polygon,
    first_overlapping_parts_ref,
};

// ===========================================================================
// Decomposition
// ===========================================================================

/// A **restricted** collection operand: the three strata of the module rule.
pub(crate) struct GcParts<'a> {
    /// The polygonal leaves, interior-disjoint by the restriction.
    pub(crate) area: Option<ArealParts<'a>>,
    /// The line leaves, with their joint mod-2 boundary.
    pub(crate) lines: Option<LineSet<'a>>,
    /// The point leaves' coordinates, duplicates kept (they change nothing).
    pub(crate) pts: Vec<Coord2d>,
}

/// The polygonal leaves of a collection as an areal operand.
///
/// The same operand `Areal` presents for a `Polygon` or `MultiPolygon`,
/// over borrowed parts gathered from scattered leaves. The restriction
/// check has established the one invariant every areal analysis relies on:
/// the parts meet at most at finitely many points, so the region lies on
/// the left of every directed ring edge. The accessors visit parts in leaf
/// order, the order the equivalent `MultiPolygon` would hold them in.
pub(crate) struct ArealParts<'a> {
    parts: Vec<&'a Polygon>,
}

impl ArealParts<'_> {
    fn rings(&self) -> impl Iterator<Item = &crate::types::LinearRing> + '_ {
        self.parts
            .iter()
            .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes().iter()))
    }
}

impl AreaOperand for ArealParts<'_> {
    fn for_edges_near(&self, _w: &Window, f: &mut dyn FnMut(Coord2d, Coord2d) -> bool) {
        for ring in self.rings() {
            for (a, b) in ring.edges() {
                if !f(a, b) {
                    return;
                }
            }
        }
    }

    fn for_vertices_near(&self, _w: &Window, f: &mut dyn FnMut(Coord2d) -> bool) {
        for v in self.rings().flat_map(|r| r.coords().as_slice()) {
            if !f(*v) {
                return;
            }
        }
    }

    fn for_edges_with_ref_near(
        &self,
        _w: &Window,
        f: &mut dyn FnMut(usize, usize, Coord2d, Coord2d) -> bool,
    ) {
        for (part, p) in self.parts.iter().enumerate() {
            let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
            for (ring, r) in rings.enumerate() {
                for (a, b) in r.edges() {
                    if !f(part, ring, a, b) {
                        return;
                    }
                }
            }
        }
    }

    fn has_edge_outside(&self, _w: &Window) -> bool {
        false
    }

    fn locate(&self, p: Coord2d) -> Location {
        if self.parts.iter().any(|q| on_part_boundary(p, q)) {
            Location::Boundary
        } else if self.parts.iter().any(|q| inside_part_off_boundary(p, q)) {
            Location::Interior
        } else {
            Location::Exterior
        }
    }

    fn envelope(&self) -> Option<(f64, f64, f64, f64)> {
        crate::types::bbox_fold(self.parts.iter().map(|p| p.exterior().bbox()))
    }
}

/// The operand a `GeometryCollection` reduces to.
///
/// The leaves are gathered depth first with their paths (an explicit stack,
/// never recursion). The polygonal leaves are checked against the
/// multipolygon part rule; a violation is the refusal, with both leaf
/// paths. A collection of one stratum reduces to the existing operand kind:
/// line leaves become the `MultiLineString` operand (the same `LineSet`,
/// the same joint mod-2 boundary), point leaves and the empty collection
/// become the point-set operand. Anything with polygonal leaves, or with
/// several strata, is the collection operand of [`merged_matrix`].
///
/// # Errors
///
/// [`GeometryError::CollectionElementsOverlap`] when two polygonal leaves'
/// interiors meet or their boundaries share a stretch.
pub(crate) fn collection_operand(gc: &GeometryCollection) -> Result<Operand<'_>, GeometryError> {
    let mut polys: Vec<&Polygon> = Vec::new();
    let mut poly_paths: Vec<Vec<usize>> = Vec::new();
    let mut lines: Vec<&LineString> = Vec::new();
    let mut pts: Vec<Coord2d> = Vec::new();

    let mut stack = vec![gc.geometries().iter().enumerate()];
    let mut path: Vec<usize> = Vec::new();
    while let Some(top) = stack.last_mut() {
        let Some((k, g)) = top.next() else {
            stack.pop();
            path.pop();
            continue;
        };
        path.push(k);
        match g {
            Geometry::Polygon(p) => {
                polys.push(p);
                poly_paths.push(path.clone());
            }
            Geometry::MultiPolygon(m) => {
                for (part, p) in m.polygons().iter().enumerate() {
                    polys.push(p);
                    let mut leaf = path.clone();
                    leaf.push(part);
                    poly_paths.push(leaf);
                }
            }
            Geometry::LineString(l) => lines.push(l),
            Geometry::MultiLineString(m) => lines.extend(m.line_strings()),
            Geometry::Point(p) => pts.push(p.coord()),
            Geometry::MultiPoint(m) => pts.extend(m.points().iter().map(|p| p.coord())),
            Geometry::GeometryCollection(inner) => {
                stack.push(inner.geometries().iter().enumerate());
                continue; // the element index stays on the path while inside
            }
        }
        path.pop();
    }

    if let Some((first, second)) = first_overlapping_parts_ref(&polys) {
        return Err(GeometryError::CollectionElementsOverlap {
            first: poly_paths[first].clone(),
            second: poly_paths[second].clone(),
        });
    }
    Ok(if polys.is_empty() && lines.is_empty() {
        Operand::Points(pts)
    } else if polys.is_empty() && pts.is_empty() {
        Operand::Lines(LineSet::new(lines))
    } else {
        Operand::Collection(GcParts {
            area: (!polys.is_empty()).then_some(ArealParts { parts: polys }),
            lines: (!lines.is_empty()).then(|| LineSet::new(lines)),
            pts,
        })
    })
}

/// The boundary of a restricted collection, as `relate` locates it: every
/// ring of the polygonal leaves (leaf order, exterior then holes), and the
/// mod-2 endpoints of the line leaves that lie in the areal part's exterior,
/// in exact lexicographic order.
///
/// By the location rule, a point of the areal part's boundary is a
/// collection boundary point; a point of its interior is an interior point
/// (absorption); outside it, a mod-2 endpoint of the lineal part is a
/// boundary point and nothing else is, since points have no boundary. An
/// endpoint on the areal boundary is already a point of a ring, so it is not
/// listed again.
///
/// # Errors
///
/// [`GeometryError::CollectionElementsOverlap`], as `relate` refuses.
pub(crate) fn collection_boundary(
    gc: &GeometryCollection,
) -> Result<(Vec<&LinearRing>, Vec<Coord2d>), GeometryError> {
    Ok(match collection_operand(gc)? {
        Operand::Lines(l) => (Vec::new(), l.boundary().to_vec()),
        Operand::Collection(parts) => {
            let rings: Vec<&LinearRing> = parts.area.as_ref().map_or_else(Vec::new, |a| {
                a.parts
                    .iter()
                    .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes().iter()))
                    .collect()
            });
            let pts = parts.lines.as_ref().map_or_else(Vec::new, |l| {
                l.boundary()
                    .iter()
                    .copied()
                    .filter(|&x| {
                        parts
                            .area
                            .as_ref()
                            .is_none_or(|a| a.locate(x) == Location::Exterior)
                    })
                    .collect()
            });
            (rings, pts)
        }
        // A collection reduces to points, lines or the mixed operand; points
        // have no boundary.
        _ => (Vec::new(), Vec::new()),
    })
}

// ===========================================================================
// Strata
// ===========================================================================

/// An operand decomposed into its strata; a simple operand has one.
pub(crate) struct Strata<'s> {
    area: Option<&'s dyn AreaOperand>,
    lines: Option<&'s LineSet<'s>>,
    pts: &'s [Coord2d],
}

impl<'s> Strata<'s> {
    /// The strata view of any operand.
    pub(crate) fn of(op: &'s Operand<'s>) -> Self {
        match op {
            Operand::Points(p) => Self {
                area: None,
                lines: None,
                pts: p,
            },
            Operand::Lines(l) => Self {
                area: None,
                lines: Some(l),
                pts: &[],
            },
            Operand::Area(a) => Self {
                area: Some(a),
                lines: None,
                pts: &[],
            },
            Operand::Collection(gc) => Self::of_parts(gc),
        }
    }

    /// The strata of a decomposed collection.
    pub(crate) fn of_parts(gc: &'s GcParts<'s>) -> Self {
        Self {
            area: gc.area.as_ref().map(|a| a as &dyn AreaOperand),
            lines: gc.lines.as_ref(),
            pts: &gc.pts,
        }
    }

    /// An areal operand alone, as on the side of a prepared mask.
    pub(crate) fn of_area(a: &'s dyn AreaOperand) -> Self {
        Self {
            area: Some(a),
            lines: None,
            pts: &[],
        }
    }

    /// The mod-2 boundary of the lineal stratum, empty without one.
    fn mod2(&self) -> &[Coord2d] {
        self.lines.map_or(&[], |l| l.boundary())
    }
}

/// The location of `x` against a strata-decomposed operand: the module's
/// six-step rule, every test exact.
fn gc_locate(x: Coord2d, s: &Strata<'_>) -> Location {
    if let Some(a) = s.area {
        match a.locate(x) {
            Location::Interior => return Location::Interior,
            Location::Boundary => return Location::Boundary,
            Location::Exterior => {}
        }
    }
    if let Some(l) = s.lines
        && l.contains(x)
    {
        return if l.boundary().contains(&x) {
            Location::Boundary
        } else {
            Location::Interior
        };
    }
    if s.pts.contains(&x) {
        return Location::Interior;
    }
    Location::Exterior
}

// ===========================================================================
// The jointly labeled walk
// ===========================================================================

/// One operand a walked segment is classified against.
enum Classifier<'c> {
    /// An areal operand: labels are `Location`s.
    Areal(&'c dyn AreaOperand),
    /// A lineal operand: labels are on (`true`) or off the line set.
    Lineal(&'c LineSet<'c>),
}

/// A classifier's label for a stretch or at a point.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Lab {
    /// Against an areal classifier.
    Loc(Location),
    /// Against a lineal classifier: on the line set or off it.
    On(bool),
}

impl Lab {
    fn initial(c: &Classifier<'_>) -> Self {
        match c {
            // Corrected by the start event before any stretch is emitted.
            Classifier::Areal(_) => Lab::Loc(Location::Exterior),
            Classifier::Lineal(_) => Lab::On(false),
        }
    }
}

/// Where an event lies on the walked segment.
#[derive(Clone, Copy)]
enum Pos {
    /// An input coordinate on the segment.
    At(Coord2d),
    /// The proper crossing of the segment with the edge or line segment
    /// `a → b`, named and never constructed.
    Cross(Coord2d, Coord2d),
}

/// One classifier's event.
struct Ev {
    which: usize,
    pos: Pos,
    /// The classifier's exact label at the event point itself.
    point: Lab,
    /// The label of the open stretch after the event, where it changes.
    ahead: Option<Lab>,
}

/// The exact parameter of an event position along `p → q`, for ordering.
///
/// An input point lifts exactly (`from_f64`); a proper crossing is the
/// rational line–line intersection, which exists because a proper crossing
/// is never parallel. Every coordinate reaching this module is finite by
/// the type invariants, so `None` is unreachable. The comparator orders a
/// `None` as `Equal`, which can only affect unreachable inputs.
fn param(pos: &Pos, p: Coord2d, q: Coord2d) -> Option<RationalCoord2d> {
    match *pos {
        Pos::At(x) => RationalCoord2d::from_f64([x.x, x.y]).ok(),
        Pos::Cross(a, b) => RationalCoord2d::from_segment_intersection(
            [p.x, p.y],
            [q.x, q.y],
            [a.x, a.y],
            [b.x, b.y],
        )
        .ok()
        .flatten(),
    }
}

/// The events of one classifier on the segment `p → q`.
fn classifier_events(which: usize, c: &Classifier<'_>, p: Coord2d, q: Coord2d) -> Vec<Ev> {
    match c {
        Classifier::Areal(a) => segment_events(p, q, *a)
            .into_iter()
            .map(|ev| {
                let (pos, point) = match ev.at {
                    EventAt::Point(x) => (Pos::At(x), Lab::Loc(a.locate(x))),
                    EventAt::Crossing { a: ea, b: eb } => {
                        (Pos::Cross(ea, eb), Lab::Loc(Location::Boundary))
                    }
                };
                Ev {
                    which,
                    pos,
                    point,
                    ahead: Some(Lab::Loc(ev.ahead)),
                }
            })
            .collect(),
        Classifier::Lineal(l) => lineal_events(which, l, p, q),
    }
}

/// The events of a lineal classifier: collinear-overlap toggles, proper
/// crossings, and vertices touching the segment.
// `a.x == b.x` and the parameter comparisons are on stored coordinates of
// collinear points: exact, the hazard `float_cmp` exists for is absent.
#[allow(clippy::float_cmp)]
fn lineal_events(which: usize, l: &LineSet<'_>, p: Coord2d, q: Coord2d) -> Vec<Ev> {
    let along = |c: Coord2d| if p.x == q.x { c.y } else { c.x };
    let forward = along(q) > along(p);
    let pm = |c: Coord2d| if forward { along(c) } else { -along(c) };

    // Union of the collinear overlaps, as parameter toggles.
    let mut toggles: Vec<(f64, Coord2d, i32)> = Vec::new();
    let mut out = Vec::new();
    for (c, d) in l.segments() {
        if collinear_overlap(p, q, c, d) {
            // The overlap endpoints are input coordinates on the segment:
            // the inner two of {p, q, c, d} by parameter.
            let (lo, hi) = if pm(c) <= pm(d) { (c, d) } else { (d, c) };
            let lo = if pm(lo) < pm(p) { p } else { lo };
            let hi = if pm(hi) > pm(q) { q } else { hi };
            if pm(lo) < pm(hi) {
                toggles.push((pm(lo), lo, 1));
                toggles.push((pm(hi), hi, -1));
            }
        } else if proper_crossing(p, q, c, d) {
            out.push(Ev {
                which,
                pos: Pos::Cross(c, d),
                point: Lab::On(true),
                ahead: None,
            });
        }
        // Every vertex on the walked segment is an event, in the overlap
        // branch too: the union sweep suppresses a toggle where overlaps
        // stack (depth 2 → 1 changes nothing), yet the vertex may coincide
        // with another classifier's proper crossing, and the merged group
        // must then carry the input coordinate so the candidate path can
        // test it against the mod-2 boundaries and point strata.
        for v in [c, d] {
            if v != p && v != q && on_segment(v, p, q) {
                out.push(Ev {
                    which,
                    pos: Pos::At(v),
                    point: Lab::On(true),
                    ahead: None,
                });
            }
        }
    }
    toggles.sort_by(|a, b| a.0.total_cmp(&b.0).then(b.2.cmp(&a.2)));
    let mut depth = 0;
    for (_, at, delta) in toggles {
        let was = depth > 0;
        depth += delta;
        let now = depth > 0;
        if was != now {
            out.push(Ev {
                which,
                pos: Pos::At(at),
                point: Lab::On(true),
                ahead: Some(Lab::On(now)),
            });
        }
    }
    out
}

/// What one walk of `p → q` yields, in order along the segment.
enum WalkItem {
    /// A point where some classifier has an event: its input coordinate
    /// when it has one (a pure crossing has none), and the joint labels at
    /// the point itself.
    Point(Option<Coord2d>, Vec<Lab>),
    /// An open sub-interval of positive length between consecutive events
    /// (or an event and `q`), with its constant joint labels.
    Stretch(Vec<Lab>),
}

/// Walks `p → q` against `cls`: each event point's joint labels and each
/// open stretch's joint labels, in order along the segment.
///
/// A classifier with no event at a position contributes its current
/// stretch label there, which is exact: a point with none of that
/// classifier's events is on none of its features, so its location is the
/// stretch's constant one. The segment's endpoints are the caller's to
/// handle. With no classifiers the whole segment is one unlabeled stretch.
fn walk_segment(p: Coord2d, q: Coord2d, cls: &[Classifier<'_>]) -> Vec<WalkItem> {
    let mut cur: Vec<Lab> = cls.iter().map(Lab::initial).collect();
    let mut out = Vec::new();
    if cls.is_empty() {
        out.push(WalkItem::Stretch(cur));
        return out;
    }
    let mut evs: Vec<Ev> = Vec::new();
    for (i, c) in cls.iter().enumerate() {
        evs.extend(classifier_events(i, c, p, q));
    }
    // Order along the segment, exactly. `rational_compare_along` is a
    // strict weak ordering (its doc carries the proof), so the sort is
    // sound; `None` parameters are unreachable on validated geometry.
    let params: Vec<Option<RationalCoord2d>> = evs.iter().map(|e| param(&e.pos, p, q)).collect();
    let q_param = param(&Pos::At(q), p, q);
    let cmp = |x: &Option<RationalCoord2d>, y: &Option<RationalCoord2d>| -> Ordering {
        match (x, y) {
            (Some(a), Some(b)) => {
                rational_compare_along([p.x, p.y], [q.x, q.y], a, b).unwrap_or(Ordering::Equal)
            }
            _ => Ordering::Equal,
        }
    };
    let mut order: Vec<usize> = (0..evs.len()).collect();
    order.sort_by(|&i, &j| cmp(&params[i], &params[j]));

    // The stretch from `p` to the first event. An areal classifier always
    // has a start event at `p` (`segment_events` emits it), and a lineal
    // overlap reaching `p` toggles at `p`, so whenever this fires the
    // initial labels are the true ones: off every lineal classifier, and no
    // areal classifier present.
    let p_param = param(&Pos::At(p), p, q);
    if order
        .first()
        .is_none_or(|&f| cmp(&params[f], &p_param) != Ordering::Equal)
    {
        out.push(WalkItem::Stretch(cur.clone()));
    }

    let mut k = 0;
    while k < order.len() {
        let mut j = k + 1;
        while j < order.len() && cmp(&params[order[k]], &params[order[j]]) == Ordering::Equal {
            j += 1;
        }
        let mut at_point = cur.clone();
        let mut coord: Option<Coord2d> = None;
        for &e in &order[k..j] {
            at_point[evs[e].which] = evs[e].point;
            if let Pos::At(x) = evs[e].pos {
                coord = Some(x);
            }
        }
        out.push(WalkItem::Point(coord, at_point));
        for &e in &order[k..j] {
            if let Some(ahead) = evs[e].ahead {
                cur[evs[e].which] = ahead;
            }
        }
        // The stretch to the next event, or to `q` when it has length.
        if j < order.len() || cmp(&params[order[k]], &q_param) != Ordering::Equal {
            out.push(WalkItem::Stretch(cur.clone()));
        }
        k = j;
    }
    out
}

// ===========================================================================
// Matrix assembly
// ===========================================================================

/// Row and column indices.
const I: usize = 0;
const B: usize = 1;
const E: usize = 2;

/// The nine cells as maxima over contributions.
struct Acc {
    cells: [[Dimension; 3]; 3],
    /// `true` raises `(row, col)` transposed: the walked side is operand B.
    swap: bool,
}

impl Acc {
    fn new() -> Self {
        let mut cells = [[Dimension::Empty; 3]; 3];
        cells[E][E] = Dimension::Area;
        Self { cells, swap: false }
    }

    fn raise(&mut self, row: usize, col: usize, d: Dimension) {
        let (r, c) = if self.swap { (col, row) } else { (row, col) };
        if d > self.cells[r][c] {
            self.cells[r][c] = d;
        }
    }

    fn row_of(loc: Location) -> usize {
        match loc {
            Location::Interior => I,
            Location::Boundary => B,
            Location::Exterior => E,
        }
    }

    fn matrix(&self) -> IntersectionMatrix {
        let c = &self.cells;
        IntersectionMatrix {
            interior_interior: c[I][I],
            interior_boundary: c[I][B],
            interior_exterior: c[I][E],
            boundary_interior: c[B][I],
            boundary_boundary: c[B][B],
            boundary_exterior: c[B][E],
            exterior_interior: c[E][I],
            exterior_boundary: c[E][B],
            exterior_exterior: c[E][E],
        }
    }
}

/// The column of a stretch or unrepresentable point against the other
/// side, from its walk labels: the module's location rule, with the `P`
/// stratum absent (a finite point set contains no stretch, and an
/// unrepresentable crossing is never an input point, since a crossing at
/// an input point merges with that point's own event and takes its labels).
fn col_of(area_lab: Option<Lab>, line_lab: Option<Lab>) -> usize {
    match area_lab {
        Some(Lab::Loc(Location::Interior)) => return I,
        Some(Lab::Loc(Location::Boundary)) => return B,
        _ => {}
    }
    if line_lab == Some(Lab::On(true)) {
        return I;
    }
    E
}

/// The contributions of one side's lineal stratum, walked jointly against
/// its own areal mask and the other side.
fn lineal_contributions(a: &Strata<'_>, b: &Strata<'_>, acc: &mut Acc) {
    let Some(l) = a.lines else { return };
    let mut cls: Vec<Classifier<'_>> = Vec::new();
    // The walked line set itself is a classifier, for one reason: a proper
    // crossing can coincide with one of the walked side's own vertices (a
    // self-overlapping line folds back across its own endpoint), and that
    // vertex may be a mod-2 boundary point. As a classifier, the vertex is
    // a coordinate-bearing event, so the merged group carries its input
    // coordinate and the candidate path below re-locates it exactly. Its
    // labels are never read.
    cls.push(Classifier::Lineal(l));
    let mask = a.area.map(|ar| {
        cls.push(Classifier::Areal(ar));
        cls.len() - 1
    });
    let tgt_a = b.area.map(|ar| {
        cls.push(Classifier::Areal(ar));
        cls.len() - 1
    });
    let tgt_l = b.lines.map(|ln| {
        cls.push(Classifier::Lineal(ln));
        cls.len() - 1
    });
    let masked = |labels: &[Lab]| {
        mask.is_some_and(|m| matches!(labels[m], Lab::Loc(Location::Interior | Location::Boundary)))
    };

    // A representable candidate point, relocated directly: exact, and it
    // sees the point strata and mod-2 boundaries the stretch labels cannot.
    let candidate = |x: Coord2d, acc: &mut Acc| {
        if a.area.is_some_and(|ar| ar.locate(x) != Location::Exterior) {
            return; // absorbed into the areal strata
        }
        let row = if a.mod2().contains(&x) { B } else { I };
        acc.raise(row, Acc::row_of(gc_locate(x, b)), Dimension::Point);
    };

    for (p, q) in l.segments() {
        candidate(p, acc);
        candidate(q, acc);
        for item in walk_segment(p, q, &cls) {
            match item {
                WalkItem::Point(Some(x), _) => candidate(x, acc),
                WalkItem::Point(None, labels) => {
                    // A pure crossing: not an input point, so never a mod-2
                    // boundary point of either side and never in a point
                    // stratum.
                    if !masked(&labels) {
                        let col = col_of(tgt_a.map(|i| labels[i]), tgt_l.map(|i| labels[i]));
                        acc.raise(I, col, Dimension::Point);
                    }
                }
                WalkItem::Stretch(labels) => {
                    if !masked(&labels) {
                        let col = col_of(tgt_a.map(|i| labels[i]), tgt_l.map(|i| labels[i]));
                        acc.raise(I, col, Dimension::Line);
                    }
                }
            }
        }
    }
}

/// The boundary-coverage contributions of one side's areal stratum: each
/// ring edge walked against the other side, for the cells a per-part
/// matrix cannot decide (`B × E` above all: a boundary stretch escapes the
/// other side only when it escapes its area *and* its lines *at once*).
fn boundary_contributions(a: &Strata<'_>, b: &Strata<'_>, acc: &mut Acc) {
    let Some(ar) = a.area else { return };
    let mut cls: Vec<Classifier<'_>> = Vec::new();
    let tgt_a = b.area.map(|x| {
        cls.push(Classifier::Areal(x));
        cls.len() - 1
    });
    let tgt_l = b.lines.map(|x| {
        cls.push(Classifier::Lineal(x));
        cls.len() - 1
    });
    let mut edges: Vec<(Coord2d, Coord2d)> = Vec::new();
    ar.for_edges_near(&ar.all_window(), &mut |p, q| {
        edges.push((p, q));
        true
    });
    for (p, q) in edges {
        for item in walk_segment(p, q, &cls) {
            // Isolated boundary points never decide a cell here: `Bnd(A)`
            // minus the other side's closed point set is relatively open in
            // the boundary curves, so `B × E` has no point case, and every
            // point contact the other cells need is raised by the other
            // side's own candidates and loops.
            if let WalkItem::Stretch(labels) = item {
                let col = col_of(tgt_a.map(|i| labels[i]), tgt_l.map(|i| labels[i]));
                acc.raise(B, col, Dimension::Line);
            }
        }
    }
}

/// The contributions of one side's punctal stratum *and* of its lineal
/// stratum's mod-2 boundary points: direct exact locates.
fn point_contributions(a: &Strata<'_>, b: &Strata<'_>, acc: &mut Acc) {
    let on_own_lines = |x: Coord2d| a.lines.is_some_and(|l| l.contains(x));
    let in_own_area = |x: Coord2d| a.area.is_some_and(|ar| ar.locate(x) != Location::Exterior);
    for &x in a.pts {
        if in_own_area(x) || on_own_lines(x) {
            continue; // absorbed
        }
        acc.raise(I, Acc::row_of(gc_locate(x, b)), Dimension::Point);
    }
    for &x in a.mod2() {
        if in_own_area(x) {
            continue; // absorbed
        }
        acc.raise(B, Acc::row_of(gc_locate(x, b)), Dimension::Point);
    }
}

/// The exact DE-9IM matrix of two strata-decomposed operands. It is the
/// entry point for every pair with a collection in it.
///
/// A pure-stratum collection goes through the same generic analyses the
/// equivalent multi-geometry does (`polygons_matrix` over an areal operand,
/// the same `LineSet`, the same point set), so it agrees with the
/// multi-geometry by construction rather than by coincidence.
pub(crate) fn merged_matrix(a: &Strata<'_>, b: &Strata<'_>) -> IntersectionMatrix {
    let mut acc = Acc::new();

    // Area × area, and the proved open-set reductions for I/E against the
    // other side's lower-dimensional strata.
    match (a.area, b.area) {
        (Some(aa), Some(ba)) => {
            let m = polygons_matrix(aa, ba);
            for (row, col, d) in [
                (I, I, m.interior_interior),
                (I, B, m.interior_boundary),
                (B, I, m.boundary_interior),
                (B, B, m.boundary_boundary),
                (I, E, m.interior_exterior),
                (E, I, m.exterior_interior),
            ] {
                acc.raise(row, col, d);
            }
        }
        // `Int(A)` minus a closed set of dimension at most 1 keeps
        // dimension 2; with no area on the other side nothing else can
        // reach `Area` in these cells.
        (Some(_), None) => acc.raise(I, E, Dimension::Area),
        (None, Some(_)) => acc.raise(E, I, Dimension::Area),
        (None, None) => {}
    }

    // This side's strata against the other, then the other against this.
    lineal_contributions(a, b, &mut acc);
    boundary_contributions(a, b, &mut acc);
    point_contributions(a, b, &mut acc);
    acc.swap = true;
    lineal_contributions(b, a, &mut acc);
    boundary_contributions(b, a, &mut acc);
    point_contributions(b, a, &mut acc);

    acc.matrix()
}

/// The dimension of a collection operand: its highest stratum.
pub(crate) fn parts_dimension(gc: &GcParts<'_>) -> Dimension {
    if gc.area.is_some() {
        Dimension::Area
    } else if gc.lines.is_some() {
        Dimension::Line
    } else if gc.pts.is_empty() {
        Dimension::Empty
    } else {
        Dimension::Point
    }
}
