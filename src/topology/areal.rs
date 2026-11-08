// src/topology/areal.rs

//! An areal operand for the exact analyses: one polygon, or the parts of a
//! valid multipolygon. The edge-filtered form, [`IndexedAreal`], serves
//! prepared geometry.
//!
//! # One ring set, because the parts are interior-disjoint
//!
//! A valid `MultiPolygon`'s parts meet at most at finitely many points
//! (`MultiPolygon::new` enforces it). So the union of their rings is a set of
//! simple closed curves that touch only at isolated points, and **every**
//! property the single-polygon analyses rely on still holds of that set:
//!
//! - the region lies to the **left** of every directed ring edge: each part's
//!   exterior is counterclockwise and its holes clockwise, and no other
//!   part's interior is on the right of any edge, because interiors are
//!   disjoint;
//! - a point off every ring is inside the region exactly when it is inside
//!   some part;
//! - at a point where several rings meet, the wedge test over all incident
//!   rays decides the side exactly as it does for one polygon's rings.
//!
//! The exact classifiers in `line_analysis` therefore run unchanged over all
//! parts' rings at once: one classification across all parts of A against
//! all parts of B, rather than pairwise matrices combined.
//!
//! # The boundary is topological
//!
//! A point on the boundary of **any** part is on the multipolygon's boundary,
//! including a point where two or more parts touch. That is the point-set
//! definition (every neighborhood of a touching point meets the complement)
//! and the OGC definition of a multipolygon's boundary. The linear mod-2 rule,
//! extended to areas, would make a point where an even number of parts touch
//! interior, although that point has no neighborhood inside the set. Mod-2
//! remains the rule for **lines**, where OGC specifies it.
//!
//! So one location function serves points and stretches alike, and a
//! multipolygon's location function is the single-polygon one taken over
//! all parts.
//!
//! # Two operands, one set of analyses
//!
//! The analyses ask an areal operand four things, through [`AreaOperand`]:
//! the ring edges whose closed box meets a window, the ring vertices whose
//! edges do, the exact location of a point, and the envelope. [`Areal`]
//! answers from every edge, ignoring the window. [`IndexedAreal`]
//! answers from an `EdgeIndex`, visiting only the candidate edges. Because
//! every analysis is written once, over the trait, the two cannot diverge
//! in logic. What the indexed operand must demonstrate is that **skipping an
//! edge whose box misses the window never changes an answer**, and the
//! argument is given at each accessor below.

use crate::exact::kernel::on_segment;
use crate::index::EdgeIndex;
use crate::topology::label::Location;
use crate::types::{Coord2d, Envelope, LinearRing, MultiPolygon, Polygon, Polygonal};

/// A closed box, as `(min, max)` corners, that an accessor filters by.
pub(crate) type Window = (Coord2d, Coord2d);

/// The window of a single point.
pub(crate) fn point_window(p: Coord2d) -> Window {
    (p, p)
}

/// The window of the segment `p → q`.
pub(crate) fn segment_window(p: Coord2d, q: Coord2d) -> Window {
    (
        Coord2d::new(p.x.min(q.x), p.y.min(q.y)),
        Coord2d::new(p.x.max(q.x), p.y.max(q.y)),
    )
}

/// The window of an envelope given as `(min_x, min_y, max_x, max_y)`.
pub(crate) fn envelope_window(e: (f64, f64, f64, f64)) -> Window {
    (Coord2d::new(e.0, e.1), Coord2d::new(e.2, e.3))
}

/// What the exact analyses ask of an areal operand.
///
/// Every edge and vertex accessor takes a window and may visit **more** than
/// the edges whose boxes meet it, but never fewer; both implementations
/// visit exactly those edges, the unindexed one by a scan of box
/// comparisons and the prepared one through its index. The callbacks
/// return `true` to continue and `false` to stop, so a search can end at
/// its first hit without visiting the rest.
pub(crate) trait AreaOperand {
    /// Every ring edge `(a, b)` whose closed box meets `w`, in ring
    /// direction, each once.
    fn for_edges_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d, Coord2d) -> bool);

    /// Every ring vertex that starts an edge whose closed box meets `w`;
    /// a vertex may be visited once per ring it belongs to.
    fn for_vertices_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d) -> bool);

    /// As [`for_edges_near`](Self::for_edges_near), with each edge's part
    /// and ring index (`0` the exterior, `k` the hole `k − 1`), for a ray
    /// parity kept per ring.
    fn for_edges_with_ref_near(
        &self,
        w: &Window,
        f: &mut dyn FnMut(usize, usize, Coord2d, Coord2d) -> bool,
    );

    /// Whether some ring edge's closed box does **not** meet `w`.
    fn has_edge_outside(&self, w: &Window) -> bool;

    /// Exact location: on any part's ring is `Boundary`; otherwise inside
    /// some part is `Interior`; otherwise `Exterior`.
    fn locate(&self, p: Coord2d) -> Location;

    /// `(min_x, min_y, max_x, max_y)` over every part's exterior; `None` for
    /// no parts.
    fn envelope(&self) -> Option<(f64, f64, f64, f64)>;

    /// The edges near the whole operand: every edge.
    fn all_window(&self) -> Window {
        self.envelope().map_or(
            (Coord2d::new(0.0, 0.0), Coord2d::new(0.0, 0.0)),
            envelope_window,
        )
    }
}

// ===========================================================================
// WindowedAreal: an unindexed operand's edges near one window, kept
// ===========================================================================

/// An unindexed operand with the edges an analysis against a small operand
/// can ask for collected once: the edges whose box meets the small
/// operand's window `W`, and the edges whose box meets the strip to the
/// right of `W` (the `+x` rays of `locate` from points in `W`).
///
/// Every accessor visits exactly the edges the unindexed operand's accessor
/// would visit, in the same ring order: a query window inside `W` is
/// answered from the first list (an edge meeting the query window meets
/// `W`), a ray window inside the strip from the second, and any other
/// window by the unindexed operand itself. `locate` is the edge-filtered
/// location the prepared operand uses (`exact_point::locate`), over these
/// accessors. So every analysis returns what it returns on the unindexed
/// operand; only the box comparisons against far edges are not repeated.
pub(crate) struct WindowedAreal<'a> {
    base: Base<'a>,
    window: Window,
    near: Vec<(usize, usize, Coord2d, Coord2d)>,
    strip: Vec<(usize, usize, Coord2d, Coord2d)>,
    envelope: Option<(f64, f64, f64, f64)>,
}

/// The operand a [`WindowedAreal`] answers for when a window is outside its
/// lists: unindexed, or a prepared mask's indexed operand.
enum Base<'a> {
    Plain(Areal<'a>),
    Indexed(IndexedAreal<'a>),
}

impl Base<'_> {
    fn op(&self) -> &dyn AreaOperand {
        match self {
            Base::Plain(a) => a,
            Base::Indexed(a) => a,
        }
    }
}

/// Ring edges above which an unindexed operand analyzed against a smaller
/// one is windowed. A performance choice only: both paths give the same
/// answers.
pub(crate) const WINDOW_MIN_EDGES: usize = 256;

impl<'a> WindowedAreal<'a> {
    /// Whether `base` has more ring edges than [`WINDOW_MIN_EDGES`].
    pub(crate) fn worth_it(base: &Areal<'_>) -> bool {
        base.rings().map(LinearRing::num_vertices).sum::<usize>() > WINDOW_MIN_EDGES
    }

    /// `base` windowed to `window`, if it has more than [`WINDOW_MIN_EDGES`]
    /// edges; `None` otherwise. One pass over the edges collects both lists
    /// and the envelope.
    pub(crate) fn of(base: &Areal<'a>, window: Window) -> Option<Self> {
        if !Self::worth_it(base) {
            return None;
        }
        let strip_of = |a: Coord2d, b: Coord2d| {
            a.x.max(b.x) >= window.0.x && a.y.min(b.y) <= window.1.y && a.y.max(b.y) >= window.0.y
        };
        let (mut near, mut strip) = (Vec::new(), Vec::new());
        // The envelope is over exteriors, as `Areal`'s.
        let env = crate::types::bbox_fold(base.parts.iter().map(|p| p.exterior().bbox()));
        for (part, p) in base.parts.iter().enumerate() {
            let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
            for (ring, r) in rings.enumerate() {
                for (a, b) in r.edges() {
                    if edge_meets(&window, a, b) {
                        near.push((part, ring, a, b));
                    }
                    if strip_of(a, b) {
                        strip.push((part, ring, a, b));
                    }
                }
            }
        }
        Some(Self {
            base: Base::Plain(Areal::from_parts(base.parts)),
            window,
            near,
            strip,
            envelope: env,
        })
    }

    /// A prepared mask's indexed operand with its lists for `window`, each
    /// taken from the index by one query: the active edges meeting
    /// `window`, and those meeting the strip from `window`'s left edge to
    /// the operand's right edge at `window`'s height (every `+x` ray from a
    /// point in `window`). The index returns exactly the active edges whose
    /// box meets a query window, in ascending `(part, ring, edge)` order, so
    /// each list filtered to a smaller window is what the indexed operand
    /// itself returns for it, in the same order.
    pub(crate) fn of_indexed(base: IndexedAreal<'a>, window: Window) -> Self {
        let near = base
            .candidates(&window)
            .into_iter()
            .map(|(p, r, e)| {
                let (a, b) = base.edge(p, r, e);
                (p, r, a, b)
            })
            .collect();
        let envelope = base.envelope;
        let strip = match envelope {
            Some(env) if env.2 >= window.0.x => base
                .candidates(&(window.0, Coord2d::new(env.2, window.1.y)))
                .into_iter()
                .map(|(p, r, e)| {
                    let (a, b) = base.edge(p, r, e);
                    (p, r, a, b)
                })
                .collect(),
            _ => Vec::new(),
        };
        Self {
            base: Base::Indexed(base),
            window,
            near,
            strip,
            envelope,
        }
    }

    /// The list that holds every edge meeting `w`, or `None` if neither
    /// does and the unindexed operand must answer.
    fn list_for(&self, w: &Window, rays: bool) -> Option<&[(usize, usize, Coord2d, Coord2d)]> {
        let inside =
            |o: &Window| w.0.x >= o.0.x && w.0.y >= o.0.y && w.1.x <= o.1.x && w.1.y <= o.1.y;
        if inside(&self.window) {
            Some(&self.near)
        } else if rays
            && w.0.x >= self.window.0.x
            && w.0.y >= self.window.0.y
            && w.1.y <= self.window.1.y
        {
            Some(&self.strip)
        } else {
            None
        }
    }
}

impl AreaOperand for WindowedAreal<'_> {
    fn for_edges_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d, Coord2d) -> bool) {
        match self.list_for(w, false) {
            Some(list) => {
                for &(_, _, a, b) in list {
                    if edge_meets(w, a, b) && !f(a, b) {
                        return;
                    }
                }
            }
            None => self.base.op().for_edges_near(w, f),
        }
    }

    fn for_vertices_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d) -> bool) {
        match self.list_for(w, false) {
            Some(list) => {
                for &(_, _, a, b) in list {
                    if edge_meets(w, a, b) && !f(a) {
                        return;
                    }
                }
            }
            None => self.base.op().for_vertices_near(w, f),
        }
    }

    fn for_edges_with_ref_near(
        &self,
        w: &Window,
        f: &mut dyn FnMut(usize, usize, Coord2d, Coord2d) -> bool,
    ) {
        match self.list_for(w, true) {
            Some(list) => {
                for &(part, ring, a, b) in list {
                    if edge_meets(w, a, b) && !f(part, ring, a, b) {
                        return;
                    }
                }
            }
            None => self.base.op().for_edges_with_ref_near(w, f),
        }
    }

    fn has_edge_outside(&self, w: &Window) -> bool {
        self.base.op().has_edge_outside(w)
    }

    fn locate(&self, p: Coord2d) -> Location {
        crate::topology::exact_point::locate(
            &crate::topology::exact_point::ExactPoint::Stored(p),
            self,
        )
    }

    fn envelope(&self) -> Option<(f64, f64, f64, f64)> {
        self.envelope
    }
}

// ===========================================================================
// Areal: every edge, filtered by its box
// ===========================================================================

/// Whether the closed box of the edge `a → b` meets the window `w`.
///
/// Comparisons of stored coordinates only, so the answer is exact; an edge
/// for which it is `false` has no point in `w`, so no exact test an
/// analysis would run on it inside `w` can succeed.
fn edge_meets(w: &Window, a: Coord2d, b: Coord2d) -> bool {
    a.x.min(b.x) <= w.1.x && a.x.max(b.x) >= w.0.x && a.y.min(b.y) <= w.1.y && a.y.max(b.y) >= w.0.y
}

/// One polygon, or the parts of a valid multipolygon.
pub struct Areal<'a> {
    parts: &'a [Polygon],
}

pub(crate) fn on_part_boundary(p: Coord2d, part: &Polygon) -> bool {
    core::iter::once(part.exterior())
        .chain(part.holes())
        .any(|r| r.edges().any(|(a, b)| on_segment(p, a, b)))
}

pub(crate) fn inside_part_off_boundary(p: Coord2d, part: &Polygon) -> bool {
    crate::exact::kernel::ring_contains_off_boundary(p, part.exterior())
        && !part
            .holes()
            .iter()
            .any(|h| crate::exact::kernel::ring_contains_off_boundary(p, h))
}

impl<'a> Areal<'a> {
    /// One polygon.
    pub(crate) fn single(p: &'a Polygon) -> Self {
        Self {
            parts: core::slice::from_ref(p),
        }
    }

    /// The parts of a multipolygon.
    pub(crate) fn multi(mp: &'a MultiPolygon) -> Self {
        Self::from_parts(mp.polygons())
    }

    /// The areal operand a polygon or a multipolygon describes.
    ///
    /// Takes the sealed [`Polygonal`] trait rather than a slice of parts:
    /// every exact classification over an `Areal` relies on its parts
    /// meeting at most at finitely many points (the module documentation
    /// gives the argument), and only a `Polygon` or a `MultiPolygon` that
    /// passed its constructor guarantees that. An arbitrary slice of
    /// polygons could overlap, and the classifications would then answer
    /// for a region that does not exist.
    pub fn new<P: Polygonal>(operand: &'a P) -> Self {
        Self::from_parts(operand.parts())
    }

    /// Parts known to meet at most at finitely many points: a polygon's one
    /// part, or a validated multipolygon's.
    pub(crate) fn from_parts(parts: &'a [Polygon]) -> Self {
        Self { parts }
    }

    /// Every ring of every part.
    pub(crate) fn rings(&self) -> impl Iterator<Item = &'a LinearRing> + '_ {
        self.parts
            .iter()
            .flat_map(|p| core::iter::once(p.exterior()).chain(p.holes().iter()))
    }
}

// The edge accessors visit exactly the edges whose closed box meets the
// window, in ring order: the set a prepared operand's index returns
// (`IndexedAreal`), found here by a scan of box comparisons. The trait
// allows any superset of that set, and every analysis is written for it, so
// visiting fewer edges than "all of them" changes no answer; it only skips
// the exact predicates an edge outside the window could never satisfy.
// `has_edge_outside` then reports the skipped edges, as the indexed operand
// does.
impl AreaOperand for Areal<'_> {
    fn for_edges_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d, Coord2d) -> bool) {
        for ring in self.rings() {
            for (a, b) in ring.edges() {
                if edge_meets(w, a, b) && !f(a, b) {
                    return;
                }
            }
        }
    }

    fn for_vertices_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d) -> bool) {
        for ring in self.rings() {
            for (a, b) in ring.edges() {
                if edge_meets(w, a, b) && !f(a) {
                    return;
                }
            }
        }
    }

    fn for_edges_with_ref_near(
        &self,
        w: &Window,
        f: &mut dyn FnMut(usize, usize, Coord2d, Coord2d) -> bool,
    ) {
        for (part, p) in self.parts.iter().enumerate() {
            let rings = core::iter::once(p.exterior()).chain(p.holes().iter());
            for (ring, r) in rings.enumerate() {
                for (a, b) in r.edges() {
                    if edge_meets(w, a, b) && !f(part, ring, a, b) {
                        return;
                    }
                }
            }
        }
    }

    fn has_edge_outside(&self, w: &Window) -> bool {
        self.rings()
            .any(|r| r.edges().any(|(a, b)| !edge_meets(w, a, b)))
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
        // The exteriors' cached boxes (holes lie inside their exteriors).
        crate::types::bbox_fold(self.parts.iter().map(|p| p.exterior().bbox()))
    }
}

// ===========================================================================
// IndexedAreal: the edges near a window, through an EdgeIndex
// ===========================================================================

/// The parts of a valid multipolygon with an [`EdgeIndex`] over all their
/// edges, and an optional set of **active** parts (a shallow restriction to
/// candidate parts). Edges and vertices of inactive parts
/// are never visited and never located against, exactly as if the operand
/// were the active parts alone.
///
/// # Why skipping far edges changes no answer
///
/// Each accessor visits the edges whose closed box meets the window `w`,
/// and the analyses ask only questions whose answer an edge outside `w`
/// cannot affect:
///
/// - **edges containing a point, or crossing a segment** (`incident_rays`,
///   `segment_events`, the contact and crossing searches): an edge that
///   contains a point of `w`, or meets a segment inside `w`, has a point in
///   `w`, so its box meets `w`;
/// - **vertices on a segment**: such a vertex lies in the segment's box;
/// - **the ray parity of a point** (`locate`): an edge crosses the +x ray
///   from `p` only if its box meets the ray's box, `[p.x, max_x] × [p.y]`
///   (`edge_crosses_ray` documents why), so the parity over the candidate
///   edges is the parity over all;
/// - **edges wholly outside another operand's box** (`has_edge_outside`,
///   used by the ring classification): such an edge lies in that operand's
///   exterior, which is its whole classification.
///
/// The candidate set is exact, precisely the edges whose boxes meet `w`
/// (the `index` module doc), so nothing is visited twice.
pub struct IndexedAreal<'a> {
    parts: &'a [Polygon],
    index: &'a EdgeIndex,
    /// `active[i]` for part `i`; `None` means every part.
    active: Option<&'a [bool]>,
    /// The envelope of the active parts' exteriors.
    envelope: Option<(f64, f64, f64, f64)>,
}

impl<'a> IndexedAreal<'a> {
    /// The operand over `parts`, whose edges `index` indexes (built by
    /// `EdgeIndex::new(parts)`), restricted to the parts `active` marks.
    /// `envelope` is the active parts' exterior envelope, which the caller
    /// has from its part boxes.
    pub(crate) fn new(
        parts: &'a [Polygon],
        index: &'a EdgeIndex,
        active: Option<&'a [bool]>,
        envelope: Option<(f64, f64, f64, f64)>,
    ) -> Self {
        Self {
            parts,
            index,
            active,
            envelope,
        }
    }

    fn is_active(&self, part: usize) -> bool {
        self.active.is_none_or(|a| a[part])
    }

    fn ring(&self, part: usize, ring: usize) -> &'a LinearRing {
        let p = &self.parts[part];
        if ring == 0 {
            p.exterior()
        } else {
            &p.holes()[ring - 1]
        }
    }

    fn edge(&self, part: usize, ring: usize, edge: usize) -> (Coord2d, Coord2d) {
        let c = self.ring(part, ring).coords().as_slice();
        (c[edge], c[(edge + 1) % c.len()])
    }

    fn window_envelope(w: &Window) -> Option<Envelope> {
        // Both corners are stored coordinates and `min <= max` by
        // construction; the constructor refuses only non-finite values,
        // which no validated geometry has.
        Envelope::new(w.0, w.1).ok()
    }

    /// The active edges whose box meets `w`, as `(part, ring, edge)`
    /// triples, in a deterministic order (by triple), each once.
    fn candidates(&self, w: &Window) -> Vec<(usize, usize, usize)> {
        let Some(env) = Self::window_envelope(w) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        self.index.candidates(&env, |r| {
            if self.is_active(r.part) {
                out.push((r.part, r.ring, r.edge));
            }
        });
        out.sort_unstable();
        out
    }

    /// The number of active edges.
    fn active_edge_count(&self) -> usize {
        self.parts
            .iter()
            .enumerate()
            .filter(|(i, _)| self.is_active(*i))
            .map(|(_, p)| {
                core::iter::once(p.exterior())
                    .chain(p.holes())
                    .map(LinearRing::num_vertices)
                    .sum::<usize>()
            })
            .sum()
    }
}

impl AreaOperand for IndexedAreal<'_> {
    fn for_edges_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d, Coord2d) -> bool) {
        for (p, r, e) in self.candidates(w) {
            let (a, b) = self.edge(p, r, e);
            if !f(a, b) {
                return;
            }
        }
    }

    fn for_vertices_near(&self, w: &Window, f: &mut dyn FnMut(Coord2d) -> bool) {
        for (p, r, e) in self.candidates(w) {
            let (a, _) = self.edge(p, r, e);
            if !f(a) {
                return;
            }
        }
    }

    fn for_edges_with_ref_near(
        &self,
        w: &Window,
        f: &mut dyn FnMut(usize, usize, Coord2d, Coord2d) -> bool,
    ) {
        for (p, r, e) in self.candidates(w) {
            let (a, b) = self.edge(p, r, e);
            if !f(p, r, a, b) {
                return;
            }
        }
    }

    fn has_edge_outside(&self, w: &Window) -> bool {
        self.candidates(w).len() < self.active_edge_count()
    }

    fn locate(&self, p: Coord2d) -> Location {
        // The edge-filtered location, shared with the exact-point oracle:
        // on a ring among the edges whose box contains `p`; otherwise the
        // ray parity per ring over the edges whose box meets the ray's box,
        // combined per part as `inside_part_off_boundary` does.
        crate::topology::exact_point::locate(
            &crate::topology::exact_point::ExactPoint::Stored(p),
            self,
        )
    }

    fn envelope(&self) -> Option<(f64, f64, f64, f64)> {
        self.envelope
    }
}
