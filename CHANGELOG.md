# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crate follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-11-08

First release.

### Added

#### Geometry types

- `Coord2d`, `CoordSeq`, `Point`, `LineString`, `LinearRing`, `Polygon`,
  `MultiPoint`, `MultiLineString`, `MultiPolygon` and `GeometryCollection`,
  the `Geometry` enum over them (non-exhaustive), and the closed-interval
  `Envelope`. Every constructor validates. Coordinates are finite, a ring
  has at least three distinct vertices, an exterior runs counterclockwise
  and its holes clockwise, each hole lies inside its exterior, and the parts
  of a `MultiPolygon` meet at most at finitely many points.
- Under the `3d` feature, `Coord3d`, `CoordSeq3d` and the 3D counterpart of
  each type. Their topology operates on the XY projection, and their winding
  is certified by the same routine as in 2D.
- `Unvouched` and `MaybeVouched`, which carry a result that is not vouched
  for together with the reason (`UnvouchedReason`).

#### Predicates

- `relate`, the exact DE-9IM matrix (`IntersectionMatrix`) of any two
  geometries, collections included. `IntersectionMatrix::from_pattern` reads
  a nine-character pattern (`PatternError`).
- `equals`, `disjoint`, `intersects`, `touches`, `crosses`, `within`,
  `contains`, `overlaps`, `covers` and `covered_by`, as free functions and as
  methods on `IntersectionMatrix`.
- `PreparedPolygonal`, a polygonal mask prepared once for repeated `relate`,
  `intersects`, `contains`, `covers`, `within`, `locate`, `locate_on_line`,
  `dwithin` and `distance` queries. Every prepared answer is identical to
  the unprepared one.

#### Validity

- `validate`, `validate_line` and `validate_multipolygon`, which return
  every reason raw coordinates cannot form a valid geometry
  (`ValidityIssue`, `LocatedIssue`).
- `validity_issues` and `is_valid` for the two checks construction leaves
  to later, ring self-intersection and a disconnected interior.

#### Measurement

- `area` and `signed_area`, with a certified sign and a documented bound;
  `length` and `perimeter`, by compensated summation (`Accum`).
- `distance`, exactly `0.0` if and only if the operands intersect and
  otherwise within a documented relative bound; `dwithin`, exact.
- `hausdorff_distance` and `frechet_distance`, discrete, with an exact zero
  and a documented bound, and the exact predicates `hausdorff_within` and
  `frechet_within`.

#### Operations

- `is_simple`, `boundary`, `normalize`, `reverse` and `unique_points`.
- `line_merge`, `line_merge_directed`, `line_sequence`, `shared_paths` and
  `split_by_point`.
- `convex_hull`; `concave_hull` (`ConcaveHullParameter`); `delaunay`
  (`Triangulation`), canonical under cocircular ties.
- `simplify` (Douglas–Peucker), `simplify_vw` (Visvalingam–Whyatt) and
  `simplify_preserve_topology`, each vouched exactly or refused.
- `split_pinches`, the split of a ring at its pinch points.
- `coverage_validate` and `coverage_validate_with_gap` (`CoverageIssue`),
  `coverage_simplify`, and, on noded input, `coverage_union` and
  `polygonize_noded`.

#### Encoding and indexing

- `codec`, a deterministic binary encoding of every type, round-trip exact
  and validating on decode (`CodecError`). Tags `0x40` to `0x7F` are
  reserved for crates that encode their own types through its extension API.
- `index`, candidate queries over polygon parts and edges with no false
  negatives (`PartIndex`, `EdgeIndex`).

#### Properties

- No public function returns a coordinate that is not an input coordinate.
- Every sign decision is an exact predicate from `arithmetica`.
- Depends on `arithmetica` 0.1.1 and `situs` 0.1.0. No `situs` type appears
  in the public API. No `unsafe` code.
- Minimum supported Rust version 1.88.
- The `unstable-counters` feature exposes internal counters that are outside
  semver.

[0.1.0]: https://github.com/JustinGeoSci/geometres/releases/tag/v0.1.0
