// src/ops/delaunay.rs

//! [`delaunay`]: the Delaunay triangulation of the distinct vertices of a
//! geometry, exact and canonical.
//!
//! # Why it is exact
//!
//! Every decision is one of two exact predicates on stored input
//! coordinates: `orient2d` (which side of a directed line a point lies on)
//! and `incircle` (whether a point lies inside the circle through three
//! others). Nothing is constructed: no circumcenter, no midpoint, no
//! distance. The result is a set of triangles whose vertices are input
//! points, bit for bit, with the defining property decided exactly: **no
//! input point lies strictly inside any output triangle's circumcircle.**
//!
//! # The tie rule, and why the result is canonical
//!
//! A point set with four or more cocircular points has more than one
//! Delaunay triangulation. The result here is made canonical by a rule a
//! reader can apply by hand:
//!
//! > Group the output triangles across every shared edge whose fourth point
//! > is exactly cocircular (`incircle` is zero). Each group is a face of the
//! > Delaunay subdivision: a convex polygon whose vertices lie on one empty
//! > circle. Cut each face into triangles by the fan from its
//! > lexicographically smallest vertex (smallest `x`, then smallest `y`).
//!
//! Why this makes the output depend only on the point set, and not on the
//! input order or the construction:
//!
//! 1. **The Delaunay subdivision is unique.** Its faces are the convex hulls
//!    of the maximal cocircular sets whose circle is empty (contains no
//!    input point strictly inside). Every Delaunay triangulation refines it:
//!    each triangle's circumcircle is empty, so each triangle lies in
//!    exactly one face, the one inscribed in the same circle.
//! 2. **The faces are recovered exactly.** For an internal edge between
//!    triangles of one face, the opposite vertex lies on the shared
//!    circumcircle, so `incircle` is zero; between triangles of different
//!    faces it lies strictly outside, so `incircle` is strictly negative.
//!    The triangles of one face are connected through shared edges (the
//!    dual of a triangulated convex polygon is a tree), so grouping across
//!    zero edges yields the faces, whichever Delaunay triangulation the
//!    construction produced.
//! 3. **The fan depends only on the face.** A face is a convex polygon; the
//!    fan from one of its vertices is a triangulation of it, and the vertex
//!    is chosen by the lexicographic order of stored coordinates, which does
//!    not depend on input order.
//!
//! # The construction
//!
//! Incremental insertion (Bowyer 1981; Watson 1981), with the unbounded face
//! handled by ghost triangles (one virtual triangle per hull edge, sharing
//! a symbolic vertex "beyond" the hull), so that every edge has two
//! triangles and insertion outside the hull is the same code path as
//! insertion inside. A new point's *cavity* is the set of triangles whose
//! circumdisk contains it (for a ghost triangle: the open half-plane beyond
//! its hull edge, plus the open edge itself); the *cavity* is deleted and
//! re-triangulated as the fan from the new point. Points are inserted in a
//! space-filling-curve order of the stored coordinates, and each search for
//! the cavity walks from the last triangle created (a *visibility walk*,
//! which terminates on a Delaunay triangulation: Devillers, Pion & Teillaud,
//! *Walking in a triangulation*, Int. J. Found. Comput. Sci. 13(2), 2002);
//! a walk that exceeds its budget falls back to a scan of every triangle,
//! so the insertion order and the walk affect speed and nothing else.
//!
//! # Certification
//!
//! The canonicalization pass computes `incircle` for every internal edge to
//! find the tie groups, so the Delaunay property of the
//! whole result is re-checked on every call at no extra cost, and each
//! output triangle's orientation is certified by one `orient2d`. A failure
//! of either check, which the argument above says cannot happen, is
//! returned as [`GeometryError::Unvouched`] with
//! [`UnvouchedReason::NotDelaunay`], never as a plausible wrong
//! triangulation.
//!
//! # Degenerate input
//!
//! Fewer than three distinct points, or all points collinear, admit no
//! triangle. [`Triangulation::triangles`] is then empty, and
//! [`Triangulation::edges`] returns the Delaunay graph those points still
//! have: the edges between lexicographically consecutive points along the
//! line. (For collinear points the lexicographic order is the order along
//! the line, each consecutive pair has an empty circle (the circle on that
//! segment as diameter), and a non-consecutive pair has none, because every
//! circle through its two points contains the points between them.) An
//! input with no coordinates gives the empty triangulation.

use core::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use arithmetica::Sign;

use crate::GeometryError;
use crate::ops::basic::{coords_of, key, lex};
use crate::predicates::{Relatable, borrowed};
use crate::types::{
    Coord2d, CoordSeq, Geometry, GeometryCollection, LineString, LinearRing, MultiLineString,
    Polygon, UnvouchedReason,
};

/// The symbolic vertex shared by the ghost triangles. Never an index into
/// the point array: construction refuses more points than `u32::MAX - 1`.
const GHOST: u32 = u32::MAX;

/// A triangle slot that has no neighbor yet, during wiring only.
const UNWIRED: u32 = u32::MAX;

// ======================================================================
// Exact predicates, counted
// ======================================================================

/// Exact orientation of `c` about the directed line `a → b`.
///
/// Every coordinate comes from a validated geometry, so it is finite and
/// `orient2d` cannot return `Err`; the fallback is unreachable and kept only
/// because the crate denies `unwrap` in production code.
fn orient(a: Coord2d, b: Coord2d, c: Coord2d) -> Sign {
    crate::count!(delaunay_orient_calls);
    arithmetica::orient2d([a.x, a.y], [b.x, b.y], [c.x, c.y]).unwrap_or(Sign::Zero)
}

/// Exact position of `d` relative to the circle through `a`, `b`, `c`
/// (counterclockwise): `Positive` strictly inside, `Zero` on the circle,
/// `Negative` strictly outside. The fallback is unreachable as in
/// [`orient`].
fn in_circle(a: Coord2d, b: Coord2d, c: Coord2d, d: Coord2d) -> Sign {
    crate::count!(delaunay_incircle_calls);
    arithmetica::incircle([a.x, a.y], [b.x, b.y], [c.x, c.y], [d.x, d.y]).unwrap_or(Sign::Zero)
}

// ======================================================================
// Insertion order (speed only; the result is canonical regardless)
// ======================================================================

/// A monotone map of a finite `f64` axis value to 32 bits: the sign-flipped
/// bit pattern's top half, so that larger values get larger keys.
fn axis_key(x: f64) -> u32 {
    let b = x.to_bits();
    let m = if b >> 63 == 1 { !b } else { b | (1 << 63) };
    u32::try_from(m >> 32).unwrap_or(u32::MAX)
}

/// Spreads the 32 bits of `v` into the even bit positions of a `u64`.
fn spread(v: u32) -> u64 {
    let mut x = u64::from(v);
    x = (x | (x << 16)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x << 8)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x << 2)) & 0x3333_3333_3333_3333;
    x = (x | (x << 1)) & 0x5555_5555_5555_5555;
    x
}

/// The Morton (Z-order) key of a point: a space-filling-curve order that
/// keeps nearby points nearby in the insertion sequence, so the walk from
/// the previous insertion is short. A function of the stored coordinates
/// alone, so the insertion sequence of a point set does not depend on the
/// input order.
fn morton(c: Coord2d) -> u64 {
    spread(axis_key(c.x)) | (spread(axis_key(c.y)) << 1)
}

// ======================================================================
// The mesh
// ======================================================================

/// One triangle of the working mesh: vertex indices in coherent cyclic
/// order (counterclockwise for a real triangle), and across the edge
/// opposite `v[i]`, the neighbor `n[i]`.
///
/// A ghost triangle carries [`GHOST`] as one vertex. Coherence means the
/// two triangles of an edge traverse it in opposite directions, so a ghost
/// for the hull edge that the interior triangle traverses `a → b` is
/// cyclically `(b, a, GHOST)`.
#[derive(Clone, Copy)]
struct Tri {
    v: [u32; 3],
    n: [u32; 3],
    alive: bool,
}

/// Construction failed an internal invariant. Mapped to the public refusal
/// by [`delaunay`]; expected never to occur (module documentation,
/// "Certification").
struct Broken;

struct Mesh<'p> {
    pts: &'p [Coord2d],
    tris: Vec<Tri>,
    free: Vec<u32>,
    /// Cavity stamps, parallel to `tris`: `mark[t] == stamp` marks `t` as
    /// in the cavity of the insertion being processed.
    mark: Vec<u32>,
    stamp: u32,
    /// A live triangle the next walk starts from.
    last: u32,
}

impl Mesh<'_> {
    fn coord(&self, v: u32) -> Coord2d {
        self.pts[v as usize]
    }

    fn tri(&self, t: u32) -> &Tri {
        &self.tris[t as usize]
    }

    /// The position of [`GHOST`] in `t`'s vertices, if `t` is a ghost.
    fn ghost_pos(&self, t: u32) -> Option<usize> {
        self.tri(t).v.iter().position(|&v| v == GHOST)
    }

    /// Does `p` lie strictly inside `t`'s circumdisk?
    ///
    /// For a real triangle, exactly the open circumdisk. For a ghost with
    /// directed edge `u → w` (its own traversal of its hull edge), the
    /// region the ghost stands for: the open half-plane strictly left of
    /// `u → w`, plus the open segment `(u, w)`. So a point beyond the hull,
    /// or on a hull edge's interior, lands in a cavity.
    fn violates(&self, t: u32, p: Coord2d) -> bool {
        let tri = self.tri(t);
        if let Some(g) = self.ghost_pos(t) {
            let u = self.coord(tri.v[(g + 1) % 3]);
            let w = self.coord(tri.v[(g + 2) % 3]);
            match orient(u, w, p) {
                Sign::Positive => true,
                Sign::Zero => strictly_between(u, w, p),
                Sign::Negative => false,
            }
        } else {
            let [a, b, c] = tri.v.map(|v| self.coord(v));
            in_circle(a, b, c, p) == Sign::Positive
        }
    }

    /// A triangle whose circumdisk contains `p`: the seed of `p`'s cavity.
    ///
    /// A *visibility walk* from [`Mesh::last`], which terminates on a Delaunay
    /// triangulation (module documentation); a walk that exceeds its budget
    /// falls back to scanning every triangle, which cannot fail, so this
    /// returns an error only for a point that violates no triangle: a
    /// duplicate, which deduplication has excluded.
    fn find_violated(&self, p: Coord2d) -> Result<u32, Broken> {
        let mut t = self.last;
        if let Some(g) = self.ghost_pos(t) {
            if self.violates(t, p) {
                return Ok(t);
            }
            t = self.tri(t).n[g];
        }
        let mut budget = 2 * self.tris.len() + 16;
        'walk: while budget > 0 {
            budget -= 1;
            let tri = self.tri(t);
            for i in 0..3 {
                let a = self.coord(tri.v[(i + 1) % 3]);
                let b = self.coord(tri.v[(i + 2) % 3]);
                if orient(a, b, p) == Sign::Negative {
                    let nb = tri.n[i];
                    if self.ghost_pos(nb).is_some() {
                        // p is strictly right of the hull edge a → b, which
                        // is the ghost's region.
                        return Ok(nb);
                    }
                    t = nb;
                    continue 'walk;
                }
            }
            // No edge has p strictly on the outer side: the closed triangle
            // contains p, and p is not a vertex, so p is strictly inside the
            // circumdisk (a point of a closed triangle other than a corner
            // is strictly inside the circle through the corners).
            return Ok(t);
        }
        // Budget exhausted: scan. Deterministic, and correct regardless of
        // the walk: every non-duplicate point violates some triangle.
        for t in 0..self.tris.len() {
            let ti = u32::try_from(t).map_err(|_| Broken)?;
            if self.tris[t].alive && self.violates(ti, p) {
                return Ok(ti);
            }
        }
        Err(Broken)
    }

    fn alloc(&mut self, v: [u32; 3]) -> Result<u32, Broken> {
        if let Some(t) = self.free.pop() {
            self.tris[t as usize] = Tri {
                v,
                n: [UNWIRED; 3],
                alive: true,
            };
            return Ok(t);
        }
        let t = u32::try_from(self.tris.len()).map_err(|_| Broken)?;
        if t == UNWIRED {
            return Err(Broken);
        }
        self.tris.push(Tri {
            v,
            n: [UNWIRED; 3],
            alive: true,
        });
        self.mark.push(0);
        Ok(t)
    }

    /// Inserts point `p_ix`: finds the cavity, deletes it, and fans the new
    /// point to the cavity boundary (module documentation).
    fn insert(&mut self, p_ix: u32) -> Result<(), Broken> {
        let p = self.coord(p_ix);
        let seed = self.find_violated(p)?;
        self.stamp = self.stamp.checked_add(1).ok_or(Broken)?;

        // The cavity: every triangle whose circumdisk contains p, found by
        // search from the seed (the set is edge-connected).
        let mut cavity = vec![seed];
        self.mark[seed as usize] = self.stamp;
        let mut head = 0;
        while head < cavity.len() {
            let t = cavity[head];
            head += 1;
            for i in 0..3 {
                let nb = self.tri(t).n[i];
                if self.mark[nb as usize] != self.stamp && self.violates(nb, p) {
                    self.mark[nb as usize] = self.stamp;
                    cavity.push(nb);
                }
            }
        }

        // The boundary: each directed edge a cavity triangle traverses whose
        // neighbor is outside the cavity, keyed by its first vertex. The
        // cavity is a disk, so the boundary is one directed cycle and each
        // key occurs once; a collision reports the invariant broken.
        let mut boundary: HashMap<u32, (u32, u32)> = HashMap::new();
        let mut start = GHOST;
        for &t in &cavity {
            for i in 0..3 {
                let nb = self.tri(t).n[i];
                if self.mark[nb as usize] != self.stamp {
                    let from = self.tri(t).v[(i + 1) % 3];
                    let to = self.tri(t).v[(i + 2) % 3];
                    if boundary.insert(from, (to, nb)).is_some() {
                        return Err(Broken);
                    }
                    if start == GHOST || from == GHOST {
                        // Any deterministic start works; GHOST never keys
                        // two edges, so prefer it for a stable choice.
                        start = from;
                    }
                }
            }
        }

        // The fan: one new triangle (from, to, p) per boundary edge, which
        // keeps the boundary edge's direction and so the mesh's coherent
        // orientation; p is on the cavity side, strictly left of every
        // boundary edge, so a real new triangle is counterclockwise.
        let m = boundary.len();
        if m < 3 {
            return Err(Broken);
        }
        let mut created = Vec::with_capacity(m);
        let mut from = start;
        for _ in 0..m {
            let &(to, outside) = boundary.get(&from).ok_or(Broken)?;
            let nt = self.alloc([from, to, p_ix])?;
            created.push((nt, outside));
            from = to;
        }
        if from != start {
            return Err(Broken);
        }
        for (k, &(nt, outside)) in created.iter().enumerate() {
            let prev = created[(k + m - 1) % m].0;
            let next = created[(k + 1) % m].0;
            // Edges of (from, to, p): opposite v[2] = p is (from, to), the
            // boundary edge; opposite v[0] = from is (to, p), shared with
            // the next fan triangle; opposite v[1] = to is (p, from),
            // shared with the previous.
            self.tris[nt as usize].n = [next, prev, outside];
            let [ef, et, _] = self.tris[nt as usize].v;
            // The outside triangle traverses the shared edge the other way.
            let j = (0..3)
                .find(|&j| {
                    let o = self.tri(outside);
                    o.v[(j + 1) % 3] == et && o.v[(j + 2) % 3] == ef
                })
                .ok_or(Broken)?;
            self.tris[outside as usize].n[j] = nt;
        }
        for &t in &cavity {
            self.tris[t as usize].alive = false;
            self.free.push(t);
        }
        // Start the next walk at a real new triangle; at least one exists
        // (a boundary edge between two ghosts of one straight hull run is
        // impossible, since the cavity of a point that sees one sees both).
        self.last = created
            .iter()
            .map(|&(t, _)| t)
            .find(|&t| self.ghost_pos(t).is_none())
            .ok_or(Broken)?;
        Ok(())
    }
}

/// Is `p` strictly inside the open segment `(u, w)`? All three are
/// collinear; the comparison is on stored coordinates of an axis the
/// segment spans.
#[expect(
    clippy::float_cmp,
    reason = "exact comparison of stored input coordinates"
)]
fn strictly_between(u: Coord2d, w: Coord2d, p: Coord2d) -> bool {
    if u.x == w.x {
        (u.y < p.y && p.y < w.y) || (w.y < p.y && p.y < u.y)
    } else {
        (u.x < p.x && p.x < w.x) || (w.x < p.x && p.x < u.x)
    }
}

/// The first triangulation: a chain of collinear points `q[0..]` (in
/// lexicographic order, which is the order along the line) fanned to the
/// first non-collinear point `p`, plus the ghost ring around the hull.
///
/// `side` is `orient(q[0], q[1], p)`, nonzero. Every `orient(q[i], q[i+1],
/// p)` has that same sign, because the chain is monotone along one line, so
/// with the chain read in the direction that puts `p` on the left, every
/// fan triangle is counterclockwise. This is the one triangulation of
/// these points: an edge between non-consecutive chain points would pass
/// through the points between them.
fn first_fan<'p>(
    pts: &'p [Coord2d],
    chain: &[u32],
    p_ix: u32,
    side: Sign,
) -> Result<Mesh<'p>, Broken> {
    let q: Vec<u32> = if side == Sign::Positive {
        chain.to_vec()
    } else {
        chain.iter().rev().copied().collect()
    };
    let m = q.len() - 1;
    let mut tris = Vec::with_capacity(4 * m + 4);
    for i in 0..m {
        tris.push(Tri {
            v: [q[i], q[i + 1], p_ix],
            n: [UNWIRED; 3],
            alive: true,
        });
    }
    // The hull, counterclockwise: q[0] … q[m], then p. One ghost per hull
    // edge, traversing it in the opposite direction (coherence).
    for i in 0..m {
        tris.push(Tri {
            v: [q[i + 1], q[i], GHOST],
            n: [UNWIRED; 3],
            alive: true,
        });
    }
    tris.push(Tri {
        v: [p_ix, q[m], GHOST],
        n: [UNWIRED; 3],
        alive: true,
    });
    tris.push(Tri {
        v: [q[0], p_ix, GHOST],
        n: [UNWIRED; 3],
        alive: true,
    });
    wire_all(&mut tris)?;
    let mark = vec![0; tris.len()];
    Ok(Mesh {
        pts,
        tris,
        free: Vec::new(),
        mark,
        stamp: 0,
        last: 0,
    })
}

/// Wires every triangle's neighbors by pairing the two traversals of each
/// undirected edge. Fails if any edge is traversed other than exactly
/// twice, which the construction above cannot produce.
fn wire_all(tris: &mut [Tri]) -> Result<(), Broken> {
    let mut open: HashMap<(u32, u32), (u32, usize)> = HashMap::new();
    for t in 0..tris.len() {
        let ti = u32::try_from(t).map_err(|_| Broken)?;
        for i in 0..3 {
            let a = tris[t].v[(i + 1) % 3];
            let b = tris[t].v[(i + 2) % 3];
            let k = (a.min(b), a.max(b));
            match open.remove(&k) {
                Some((t2, j)) => {
                    tris[t].n[i] = t2;
                    tris[t2 as usize].n[j] = ti;
                }
                None => {
                    if open.insert(k, (ti, i)).is_some() {
                        return Err(Broken);
                    }
                }
            }
        }
    }
    if open.is_empty() { Ok(()) } else { Err(Broken) }
}

// ======================================================================
// Canonicalization: the tie rule
// ======================================================================

/// Root with path halving.
fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        let g = parent[parent[x as usize] as usize];
        parent[x as usize] = g;
        x = g;
    }
    x
}

/// The canonical triangles of the finished mesh: each cocircular face cut
/// by the fan from its lexicographically smallest vertex (module
/// documentation), singleton faces kept as built. Also the certification:
/// a strictly positive `incircle` on any internal edge, which a Delaunay
/// triangulation cannot have, is reported as broken.
fn canonical_triangles(mesh: &Mesh<'_>) -> Result<Vec<[u32; 3]>, Broken> {
    let n = mesh.tris.len();
    let mut parent: Vec<u32> = (0..u32::try_from(n).map_err(|_| Broken)?).collect();
    let real = |t: u32| mesh.tris[t as usize].alive && mesh.ghost_pos(t).is_none();
    for t in 0..parent.len() {
        let ti = u32::try_from(t).map_err(|_| Broken)?;
        if !real(ti) {
            continue;
        }
        let tri = mesh.tri(ti);
        let [a, b, c] = tri.v.map(|v| mesh.coord(v));
        for i in 0..3 {
            let nb = tri.n[i];
            if nb <= ti || !real(nb) {
                continue;
            }
            // The vertex of nb opposite the shared edge: nb's own index of
            // the edge back to ti.
            let j = (0..3).find(|&j| mesh.tri(nb).n[j] == ti).ok_or(Broken)?;
            let d = mesh.coord(mesh.tri(nb).v[j]);
            match in_circle(a, b, c, d) {
                Sign::Positive => return Err(Broken),
                Sign::Zero => {
                    let (ra, rb) = (find(&mut parent, ti), find(&mut parent, nb));
                    parent[ra.max(rb) as usize] = ra.min(rb);
                }
                Sign::Negative => {}
            }
        }
    }
    // Singleton faces pass through; a face of several triangles is re-cut
    // by the fan rule. Member lists are gathered only for the latter.
    let mut size = vec![0u32; n];
    let mut reals: Vec<u32> = Vec::new();
    for t in 0..n {
        let ti = u32::try_from(t).map_err(|_| Broken)?;
        if real(ti) {
            reals.push(ti);
            size[find(&mut parent, ti) as usize] += 1;
        }
    }
    let mut out = Vec::with_capacity(reals.len());
    let mut multi: HashMap<u32, Vec<u32>> = HashMap::new();
    for &ti in &reals {
        let r = find(&mut parent, ti);
        if size[r as usize] == 1 {
            out.push(mesh.tri(ti).v);
        } else {
            multi.entry(r).or_default().push(ti);
        }
    }
    let mut roots: Vec<u32> = multi.keys().copied().collect();
    roots.sort_unstable();
    for r in roots {
        fan_face(mesh, &multi[&r], &mut out)?;
    }
    Ok(out)
}

/// Replaces a cocircular face's triangles by the fan from the face's
/// lexicographically smallest vertex. The face is a convex polygon
/// inscribed in an empty circle, so its boundary is one directed cycle and
/// the fan triangles are counterclockwise.
fn fan_face(mesh: &Mesh<'_>, members: &[u32], out: &mut Vec<[u32; 3]>) -> Result<(), Broken> {
    let set: HashSet<u32> = members.iter().copied().collect();
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &t in members {
        let tri = mesh.tri(t);
        for i in 0..3 {
            if !set.contains(&tri.n[i])
                && next
                    .insert(tri.v[(i + 1) % 3], tri.v[(i + 2) % 3])
                    .is_some()
            {
                return Err(Broken);
            }
        }
    }
    let v0 = next
        .keys()
        .copied()
        .min_by(|&a, &b| lex(mesh.coord(a), mesh.coord(b)))
        .ok_or(Broken)?;
    let mut cycle = Vec::with_capacity(next.len());
    let mut cur = v0;
    for _ in 0..next.len() {
        cycle.push(cur);
        cur = *next.get(&cur).ok_or(Broken)?;
    }
    if cur != v0 || cycle.len() != next.len() {
        return Err(Broken);
    }
    for k in 1..cycle.len() - 1 {
        out.push([v0, cycle[k], cycle[k + 1]]);
    }
    Ok(())
}

// ======================================================================
// The public type and entry point
// ======================================================================

/// A Delaunay triangulation: distinct input points, and triangles over them
/// whose circumcircles contain no input point strictly inside.
///
/// Canonical: the triangulation, and the order everything is reported in,
/// depend only on the set of points, not on their order in the input (the
/// documentation of src/ops/delaunay.rs gives the cocircular tie rule). Each
/// triangle is counterclockwise and starts at its lexicographically smallest
/// vertex; triangles are sorted by their vertex coordinates; edges are reported
/// with their lexicographically smaller endpoint first, sorted.
#[derive(Clone, Debug)]
pub struct Triangulation {
    points: Vec<Coord2d>,
    triangles: Vec<[u32; 3]>,
    edges: Vec<[u32; 2]>,
}

impl Triangulation {
    /// The distinct input points, in order of first occurrence (as
    /// [`unique_points`](crate::ops::unique_points) reports them).
    #[must_use]
    pub fn points(&self) -> &[Coord2d] {
        &self.points
    }

    /// The triangles, as indices into [`Self::points`], each
    /// counterclockwise, in the canonical order above.
    #[must_use]
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// The triangles as polygons, in [`Self::triangles`]' order.
    ///
    /// A [`GeometryCollection`], because neighboring triangles share edges,
    /// which a valid [`MultiPolygon`](crate::types::MultiPolygon) does not
    /// allow. Each triangle is built through the crate-internal validated
    /// constructors: its counterclockwise orientation, with nonzero area and
    /// therefore three distinct vertices, is certified by an exact
    /// `orient2d` when the triangulation is built.
    #[must_use]
    pub fn to_polygons(&self) -> GeometryCollection {
        let polys = self
            .triangles
            .iter()
            .map(|t| {
                let ring = LinearRing::from_valid(CoordSeq::new(
                    t.iter().map(|&v| self.points[v as usize]).collect(),
                ));
                Geometry::Polygon(Polygon::from_valid(ring, Vec::new()))
            })
            .collect();
        // Depth 1: no element is a collection, so the nesting limit cannot
        // be exceeded.
        GeometryCollection::from_valid(polys, 1)
    }

    /// Every edge of the triangulation exactly once, as two-point
    /// [`LineString`]s: the triangles' edges or, for degenerate input (see
    /// the module documentation), the chain of collinear points. Each edge
    /// runs from its lexicographically smaller endpoint; two distinct
    /// points are never equal, so each is a valid line string.
    #[must_use]
    pub fn edges(&self) -> MultiLineString {
        MultiLineString::new(
            self.edges
                .iter()
                .map(|e| {
                    LineString::from_valid(CoordSeq::new(vec![
                        self.points[e[0] as usize],
                        self.points[e[1] as usize],
                    ]))
                })
                .collect(),
        )
    }
}

/// The Delaunay triangulation of the distinct vertices of `g`.
///
/// Exact and canonical; the documentation of src/ops/delaunay.rs gives the
/// predicates, the cocircular tie rule and the degenerate cases. Duplicate
/// points are used once (`-0.0` and `0.0` are the same point, as everywhere
/// in this crate); an input with no coordinates gives the empty triangulation.
///
/// # Errors
///
/// - [`GeometryError::TooManyEntries`] for more than `u32::MAX - 1`
///   distinct points, never a truncation.
/// - [`GeometryError::Unvouched`] with [`UnvouchedReason::NotDelaunay`] if
///   the certification re-check fails, which is expected never to occur
///   (see src/ops/delaunay.rs).
pub fn delaunay(g: &impl Relatable) -> Result<Triangulation, GeometryError> {
    let mut seen: HashSet<(u64, u64)> = HashSet::new();
    let mut points: Vec<Coord2d> = Vec::new();
    for c in coords_of(borrowed(g)) {
        if seen.insert(key(c)) {
            points.push(c);
        }
    }
    if u32::try_from(points.len()).unwrap_or(GHOST) == GHOST {
        return Err(GeometryError::TooManyEntries {
            count: points.len(),
        });
    }
    build(points).map_err(|Broken| GeometryError::Unvouched {
        reason: UnvouchedReason::NotDelaunay,
        inputs: Vec::new(),
    })
}

/// Builds the triangulation of deduplicated `points` (the module
/// documentation's construction, canonicalization and certification).
fn build(points: Vec<Coord2d>) -> Result<Triangulation, Broken> {
    let mut order: Vec<u32> = (0..u32::try_from(points.len()).map_err(|_| Broken)?).collect();
    order.sort_unstable_by(|&i, &j| {
        let (a, b) = (points[i as usize], points[j as usize]);
        morton(a).cmp(&morton(b)).then_with(|| lex(a, b))
    });

    // Collinear prefix, in lexicographic order; then the mesh.
    let mut chain: Vec<u32> = Vec::new();
    let mut mesh: Option<Mesh<'_>> = None;
    for &ix in &order {
        let p = points[ix as usize];
        if let Some(m) = mesh.as_mut() {
            m.insert(ix)?;
            continue;
        }
        if chain.len() >= 2 {
            let side = orient(points[chain[0] as usize], points[chain[1] as usize], p);
            if side != Sign::Zero {
                mesh = Some(first_fan(&points, &chain, ix, side)?);
                continue;
            }
        }
        let at = chain.partition_point(|&q| lex(points[q as usize], p) == Ordering::Less);
        chain.insert(at, ix);
    }

    let triangles = match &mesh {
        Some(m) => canonical_triangles(m)?,
        None => Vec::new(),
    };
    finish(points, triangles, &chain)
}

/// Certifies, canonically orders, and assembles the output.
fn finish(
    points: Vec<Coord2d>,
    mut triangles: Vec<[u32; 3]>,
    chain: &[u32],
) -> Result<Triangulation, Broken> {
    for t in &mut triangles {
        let [a, b, c] = t.map(|v| points[v as usize]);
        if orient(a, b, c) != Sign::Positive {
            return Err(Broken);
        }
        // Rotate the lexicographically smallest vertex first; a rotation
        // keeps the orientation.
        let k = (0..3)
            .min_by(|&i, &j| lex(points[t[i] as usize], points[t[j] as usize]))
            .unwrap_or(0);
        *t = [t[k], t[(k + 1) % 3], t[(k + 2) % 3]];
    }
    triangles.sort_unstable_by(|s, t| {
        for i in 0..3 {
            let o = lex(points[s[i] as usize], points[t[i] as usize]);
            if o != Ordering::Equal {
                return o;
            }
        }
        Ordering::Equal
    });

    let mut edges: Vec<[u32; 2]> = if triangles.is_empty() {
        chain.windows(2).map(|w| [w[0], w[1]]).collect()
    } else {
        let mut es = Vec::with_capacity(3 * triangles.len());
        for t in &triangles {
            for i in 0..3 {
                let (u, w) = (t[i], t[(i + 1) % 3]);
                let e = if lex(points[u as usize], points[w as usize]) == Ordering::Less {
                    [u, w]
                } else {
                    [w, u]
                };
                es.push(e);
            }
        }
        es
    };
    edges.sort_unstable_by(|a, b| {
        lex(points[a[0] as usize], points[b[0] as usize])
            .then_with(|| lex(points[a[1] as usize], points[b[1] as usize]))
    });
    edges.dedup();
    Ok(Triangulation {
        points,
        triangles,
        edges,
    })
}
