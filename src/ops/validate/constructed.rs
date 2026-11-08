// src/ops/validate/constructed.rs

//! Validity of constructed geometries: what the constructors do not check.

use crate::exact::kernel::self_contacts;
use crate::predicates::{Borrowed, Relatable, borrowed};
use crate::types::{Geometry, GeometryCollection, LinearRing, MultiPolygon, Polygon};

use super::{ValidityIssue, disconnected_interior};

/// A validity issue of a constructed geometry, with where it is.
///
/// `path` leads from the geometry passed to the polygon the issue is about:
/// for a collection, the element index at each level of nesting; then, for a
/// multipolygon, the part index. A `Polygon` passed directly has an empty
/// path; part 2 of a multipolygon has `[2]`; part 0 of a multipolygon that
/// is element 1 of a collection has `[1, 0]`. The ring and any vertex or edge
/// indices are in the issue itself, counted in the polygon's stored rings.
#[non_exhaustive]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocatedIssue {
    /// Element indices through the collections, then the part index.
    pub path: Vec<usize>,
    /// The issue.
    pub issue: ValidityIssue,
}

/// Every validity issue of `g` that its construction does not check: ring
/// self-intersection ([`ValidityIssue::SelfIntersection`], which
/// construction leaves to preparation) and a disconnected interior
/// ([`ValidityIssue::DisconnectedInterior`], which `Polygon::new` admits).
/// Everything else OGC requires, construction has already refused.
///
/// Points and lines have no such issue (OGC lets a line cross itself). A
/// polygon reports each ring's self-contacts by edge index and, only if no
/// ring has one, a disconnected interior. A multipolygon reports each part's
/// issues; its parts already meet at points only. A collection reports each
/// element's, at every depth, walked without recursion. Issues come in the
/// order the polygons are met, each with its [`LocatedIssue::path`].
///
/// Exact: the self-contact test and the touch graph (`ops::validate`) are
/// exact on the stored coordinates.
#[must_use]
pub fn validity_issues(g: &impl Relatable) -> Vec<LocatedIssue> {
    let mut out = Vec::new();
    match borrowed(g) {
        Borrowed::Polygon(p) => polygon_into(p, &[], &mut out),
        Borrowed::MultiPolygon(m) => multi_polygon_into(m, &[], &mut out),
        Borrowed::GeometryCollection(c) => collection_into(c, &mut out),
        Borrowed::Point(_)
        | Borrowed::MultiPoint(_)
        | Borrowed::LineString(_)
        | Borrowed::MultiLineString(_) => {}
    }
    out
}

/// `true` if `g` is valid in the OGC sense: [`validity_issues`] is empty.
#[must_use]
pub fn is_valid(g: &impl Relatable) -> bool {
    validity_issues(g).is_empty()
}

/// The issues of one polygon, located at `path`.
fn polygon_into(p: &Polygon, path: &[usize], out: &mut Vec<LocatedIssue>) {
    let rings: Vec<&LinearRing> = core::iter::once(p.exterior()).chain(p.holes()).collect();
    let before = out.len();
    for (ring, r) in rings.iter().enumerate() {
        for edges in self_contacts(r.coords().as_slice(), false) {
            out.push(LocatedIssue {
                path: path.to_vec(),
                issue: ValidityIssue::SelfIntersection { ring, edges },
            });
        }
    }
    // The connectivity theorem needs simple rings; placement is already
    // guaranteed by construction.
    if out.len() == before
        && let Some((ring, vertex)) = disconnected_interior(&rings)
    {
        out.push(LocatedIssue {
            path: path.to_vec(),
            issue: ValidityIssue::DisconnectedInterior { ring, vertex },
        });
    }
}

fn multi_polygon_into(m: &MultiPolygon, path: &[usize], out: &mut Vec<LocatedIssue>) {
    let mut sub = path.to_vec();
    for (part, p) in m.polygons().iter().enumerate() {
        sub.push(part);
        polygon_into(p, &sub, out);
        sub.pop();
    }
}

/// Every polygon of a collection, at every depth, with an explicit stack of
/// element iterators and the path to each.
fn collection_into(c: &GeometryCollection, out: &mut Vec<LocatedIssue>) {
    let mut stack = vec![c.geometries().iter().enumerate()];
    let mut path: Vec<usize> = Vec::new();
    while let Some(top) = stack.last_mut() {
        let Some((k, g)) = top.next() else {
            stack.pop();
            path.pop();
            continue;
        };
        path.push(k);
        match g {
            Geometry::Polygon(p) => polygon_into(p, &path, out),
            Geometry::MultiPolygon(m) => multi_polygon_into(m, &path, out),
            Geometry::GeometryCollection(inner) => {
                stack.push(inner.geometries().iter().enumerate());
                // The element index stays on the path while inside it.
                continue;
            }
            Geometry::Point(_)
            | Geometry::MultiPoint(_)
            | Geometry::LineString(_)
            | Geometry::MultiLineString(_) => {}
        }
        path.pop();
    }
}
