# geometres

Exact, non-constructive computational geometry over `f64` coordinates.

`geometres` provides validated geometry types, the complete DE-9IM
topological predicates, validation with every invalidity reported, and
measurements whose sign is certified and whose magnitude carries a documented
bound. Every sign decision is an exact predicate from
[`arithmetica`](https://crates.io/crates/arithmetica), so a predicate never
returns a plausible wrong answer on near-degenerate input.

## The contract

Every public function does one of three things:

- it validates;
- it decides a predicate exactly;
- it returns a quantity whose sign is certified and whose magnitude has a
  documented bound (`area`, `length`, `distance`).

**No public function returns a coordinate that is not an input coordinate.**
Selecting among input values is allowed (an envelope takes the least and
greatest of its coordinates), and computing new ones is not. 

## Geometry types

`Point`, `LineString`, `LinearRing`, `Polygon`, `MultiPoint`,
`MultiLineString`, `MultiPolygon`, `GeometryCollection`, the `Geometry` enum
over all of them, and the axis-aligned `Envelope`. Every constructor
validates, so a value that exists has passed its checks: coordinates are
finite, a ring has at least three distinct vertices, a polygon's exterior runs
counterclockwise and its holes clockwise, each hole lies inside the exterior,
and the parts of a `MultiPolygon` meet at most at finitely many points. The
two checks construction leaves to later are described under Validity.

The `3d` feature adds `Coord3d` and the 3D counterparts of each type. Their
topology operates on the XY projection, and their winding is certified by the
same routine as in 2D.

## Predicates

`relate` computes one exact DE-9IM matrix for any two geometries, in any
combination of types, collections included. The free functions `equals`,
`disjoint`, `intersects`, `touches`, `crosses`, `within`, `contains`,
`overlaps`, `covers` and `covered_by` each read that matrix one way, and each
is also a method on `IntersectionMatrix`. A collection is accepted when its
polygonal elements together form a valid `MultiPolygon`, and refused with
`GeometryError::CollectionElementsOverlap` otherwise.

`PreparedPolygonal` prepares a polygonal mask once for repeated queries
(`relate`, `intersects`, `contains`, `covers`, `within`, `locate`,
`locate_on_line`, `dwithin`, `distance`). Every prepared answer is **identical**
to the answer of the unprepared operation. Preparation changes the cost of a
query rather than its result.

## Validity

`validate`, `validate_line` and `validate_multipolygon` take raw
coordinates and return every reason they cannot form a valid geometry, by
ring and vertex index. A constructed geometry has passed every check except
two, which construction leaves to later: ring self-intersection and a
disconnected interior. `validity_issues` reports those, at any depth of a
collection, and `is_valid` answers yes or no. Orientation is decided by a
certified area sign, and self-contact, ring placement and the connectivity
of the interior by exact tests. No issue is reported by constructing a
point.

**An operation does not refuse its input because the input is invalid.**
`validate` is the tool for detecting invalid geometry, and a caller who
operates on invalid geometry receives a result that is not vouched for.

## Operations

Each operation states what it guarantees.

| Operation | Guarantee |
|---|---|
| `area`, `signed_area` | Certified sign; computed as if in twice the working precision, with the bound documented |
| `length`, `perimeter` | Compensated sum of edge lengths |
| `distance` | Exactly `0.0` if and only if the operands intersect; otherwise positive and within a documented relative bound of the true distance |
| `dwithin` | Exact |
| `hausdorff_within`, `frechet_within` | Exact |
| `hausdorff_distance`, `frechet_distance` | Exactly `0.0` when the exact zero test holds; otherwise within a documented relative bound |
| `is_simple` | Exact, by the OGC definition for each type |
| `boundary`, `normalize`, `reverse`, `unique_points` | Exact; every output coordinate is an input coordinate |
| `line_merge`, `line_merge_directed`, `line_sequence`, `shared_paths`, `split_by_point` | Exact; every output coordinate is an input coordinate |
| `convex_hull` | Exact; the hull's vertices are input points |
| `concave_hull` | The threshold is a selection rule in `f64`; the result, a valid polygon of input vertices that contains every input point, is certified exactly or refused |
| `delaunay` | Exact and canonical: every decision is `orient2d` or `incircle`, and cocircular ties are broken by a fixed rule, so the triangulation depends only on the point set |
| `simplify`, `simplify_vw` | Vertex choice in `f64`; the result is vouched exactly or refused |
| `simplify_preserve_topology` | Every topological check is exact, and the result is vouched exactly |
| `split_pinches` | Exact; the loops are runs of the ring's own edges |
| `coverage_validate`, `coverage_validate_with_gap` | Exact, including the gap-width test |
| `coverage_simplify` | Every removal is decided by exact topology checks, shared edges are simplified once, and the result is checked exactly |
| `coverage_union`, `polygonize_noded` | Exact; every output vertex is an input vertex |

The `codec` module is a deterministic binary encoding of every type,
round-trip exact and validating on decode. The `index` module answers
candidate queries over polygon parts and edges, with no false negatives.

## Dependencies

`geometres` depends on [`arithmetica`](https://crates.io/crates/arithmetica)
for exact predicates and on [`situs`](https://crates.io/crates/situs) for
spatial indexing. No `situs` type appears in its public API. It contains no
`unsafe` code.

## Minimum supported Rust version

Rust 1.88.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
