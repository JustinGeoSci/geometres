// src/ops/area.rs

//! Signed and unsigned area of rings and polygons.
//!
//! The shoelace sum is evaluated with the compensated accumulator in
//! [`super::accum`], so it is as accurate as a computation in twice the
//! working precision (bound and conditions in that module).
//!
//! # Certified sign
//!
//! `signed_area_ring` (and so `LinearRing::signed_area`, which `Polygon::new`
//! uses for its winding and degeneracy checks) and `signed_area` of a whole
//! polygon also guarantee the **sign**. The compensated result `r` satisfies
//! `|r - exact| <= eps * |exact| + gamma_m^2 * S` with `S = sum(|terms|)`
//! (Ogita, Rump, Oishi 2005). Solving for the error gives
//! `|r - exact| <= (eps * |r| + gamma_m^2 * S) / (1 - eps)`, so when `|r|`
//! exceeds that bound the sign of `r` is the sign of the exact area. The code
//! uses four times the bound (covering the rounding in computing `S` and the
//! paper's constants) because a too-large bound only costs speed, never
//! correctness. When `|r|` does not exceed it, the sum is redone in exact
//! expansion arithmetic, whose largest component carries the exact sign.

use crate::accum::Accum;
use crate::types::{Coord2d, LinearRing, Polygon};
use arithmetica::expansion::{grow_expansion, two_product};

/// Twice the signed shoelace sum over `edges`, with a **certified sign**.
///
/// One certified routine governs every winding and area decision in the
/// crate, in both dimensions. An edge is a pair
/// of 2D coordinates: a `LinearRing`'s own edges for the 2D types, and the XY
/// projection of a `LinearRing3d`'s edges for `Polygon3d`, whose Newell `z`
/// is algebraically this same sum.
///
/// The compensated sum is accepted when [`sign_is_certified`] proves its sign;
/// otherwise the sum is redone in exact expansion arithmetic over the same
/// edges, which is why the iterator must be `Clone`.
///
/// # Term order
///
/// For each edge `(a, b)` the two products `a.x·b.y` and `−b.x·a.y` enter the
/// accumulator in that order, edge after edge. That order is part of the
/// contract: the compensated value, and so whether the exact fallback fires,
/// depends on it, and the 2D results are pinned byte for byte by
/// `the_edge_form_reproduces_the_ring_form_bit_for_bit` below.
pub(crate) fn certified_twice_area_edges<I>(edges: I) -> f64
where
    I: Iterator<Item = (Coord2d, Coord2d)> + Clone,
{
    let mut acc = Accum::new();
    let mut abs_sum = 0.0_f64;
    let mut terms = 0.0_f64;
    for (a, b) in edges.clone() {
        acc.add_product(a.x, b.y);
        acc.add_product(-b.x, a.y);
        abs_sum += (a.x * b.y).abs() + (b.x * a.y).abs();
        terms += 2.0;
    }
    let r = acc.value();
    if !r.is_finite() || !abs_sum.is_finite() {
        return r;
    }
    if sign_is_certified(r, abs_sum, terms) {
        return r;
    }
    exact_twice_area_edges(edges)
}

/// Twice the signed shoelace sum over the edges of `rings`, with a
/// **certified sign**: [`certified_twice_area_edges`] over the rings' edges in
/// order.
///
/// Both [`signed_area_ring`] and [`signed_area`] are this function, over one
/// ring and over all of a polygon's rings respectively, so the polygon-level
/// net area has an exact sign too. That matters for the thin-annulus case,
/// where a large exterior and a nearly equal hole cancel to a net area below
/// the rounding of the compensated sum; `area` asserts on the sign of that
/// result.
fn certified_twice_area<'a, I>(rings: I) -> f64
where
    I: Iterator<Item = &'a LinearRing> + Clone,
{
    certified_twice_area_edges(rings.flat_map(LinearRing::edges))
}

/// Signed area of a `LinearRing`: positive for counterclockwise winding,
/// negative for clockwise.
///
/// The magnitude is computed with the shoelace formula in compensated
/// arithmetic (see the module docs). The **sign is exact**: whenever the
/// compensated result cannot be proven to have the right sign, the sum is
/// redone in exact expansion arithmetic. Multiplying by `0.5` at the end is
/// exact. If an intermediate product overflows, the result is non-finite.
pub fn signed_area_ring(ring: &LinearRing) -> f64 {
    certified_twice_area(std::iter::once(ring)) * 0.5
}

/// `true` if the compensated sum `r` of `terms` products with absolute sum
/// `abs_sum` provably has the sign of the exact value (see the module docs).
/// Both inputs must be finite.
fn sign_is_certified(r: f64, abs_sum: f64, terms: f64) -> bool {
    let eps = f64::EPSILON * 0.5;
    let gamma = terms * eps / (1.0 - terms * eps);
    let bound = 4.0 * (eps * r.abs() + gamma * gamma * abs_sum);
    r.abs() > bound
}

/// Twice the signed area over `edges`, summed exactly in expansion arithmetic.
///
/// Every `two_product` is exact and `grow_expansion` is exact, so the
/// expansion equals the real-number shoelace sum over every edge. A
/// non-overlapping expansion has its sign in its largest (last) component;
/// the returned `f64` is the components summed smallest first, which is
/// within an ulp of the exact value. The term order is the one
/// [`certified_twice_area_edges`] uses.
pub(crate) fn exact_twice_area_edges<I>(edges: I) -> f64
where
    I: Iterator<Item = (Coord2d, Coord2d)>,
{
    let mut e: Vec<f64> = Vec::new();
    let mut h: Vec<f64> = Vec::new();
    for (a, b) in edges {
        for (x, y) in [(a.x, b.y), (-b.x, a.y)] {
            let (hi, lo) = two_product(x, y);
            for t in [lo, hi] {
                h.resize(e.len() + 1, 0.0);
                let len = grow_expansion(&e, t, &mut h);
                h.truncate(len);
                std::mem::swap(&mut e, &mut h);
            }
        }
    }
    e.iter().sum()
}

/// Signed area of a `Polygon`, accounting for holes.
///
/// The exterior ring is counterclockwise (positive) and each hole is
/// clockwise (negative), so the sum is the net area of the region. All rings
/// are accumulated in one compensated sum instead of rounding each ring's
/// area first, so cancellation between a large exterior and a large hole does
/// not lose digits, and the **sign is certified** as for a single
/// ring (see `certified_twice_area` in this module).
///
/// `Polygon::new` verifies that every hole lies inside the exterior and that
/// holes are interior-disjoint, so the exact net area of any valid `Polygon`
/// is strictly positive, and with the certified sign this function returns a
/// non-negative value for every valid polygon. A polygon with a
/// self-crossing ring (constructible, invalid) can have a negative net sum
/// ([`area`]). It can be non-finite on overflow.
pub fn signed_area(polygon: &Polygon) -> f64 {
    certified_twice_area(std::iter::once(polygon.exterior()).chain(polygon.holes().iter())) * 0.5
}

/// The area of a `Polygon`: [`signed_area`] with no change of sign.
///
/// Hole placement is validated at construction and the sign of the sum is
/// certified, so a **valid** polygon never yields a negative value. Taking the
/// absolute value would only conceal an invalid polygon, whose holes can
/// exceed its exterior and produce a negative net area. Such a polygon
/// reports that negative value here rather than a plausible positive one.
///
/// An invalid polygon can still be constructed: `Polygon::new` does not
/// check ring simplicity, so a hole that crosses itself (a spiral, say) can
/// enclose more signed area than the exterior. Its net signed sum is then
/// negative, and that sum is returned rather than made positive or panicked
/// on. A negative area is a visible sign of an invalid input, which
/// [`validate()`](super::validate()) names. An exact
/// zero is reachable too, through the one documented limit of the placement
/// check (holes that tile the exterior along collinear shared edges).
pub fn area(polygon: &Polygon) -> f64 {
    signed_area(polygon)
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
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

    fn ring(pts: &[(f64, f64)]) -> LinearRing {
        let v: Vec<Coord2d> = pts.iter().map(|&(x, y)| Coord2d::new(x, y)).collect();
        LinearRing::new(CoordSeq::new(v)).expect("valid ring")
    }

    /// A polygon whose hole crosses itself (a spiral) is constructible:
    /// `Polygon::new` does not check ring simplicity, and `validate` reports
    /// it. Its net signed sum can be negative, and `area` returns that sum
    /// rather than panicking in a debug build.
    #[test]
    fn area_of_a_polygon_with_a_self_crossing_hole_is_its_signed_sum() {
        let ext = ring(&[(0.0, 0.0), (100.0, 0.0), (100.0, 100.0), (0.0, 100.0)]);
        let spiral = ring(&[
            (10.0, 10.0),
            (10.0, 90.0),
            (90.0, 90.0),
            (90.0, 20.0),
            (20.0, 20.0),
            (20.0, 80.0),
            (80.0, 80.0),
            (80.0, 30.0),
            (30.0, 30.0),
            (30.0, 70.0),
            (70.0, 70.0),
            (70.0, 40.0),
            (40.0, 40.0),
            (40.0, 60.0),
            (60.0, 60.0),
            (60.0, 10.0),
        ]);
        let p = Polygon::new(ext, vec![spiral]).expect("constructible: simplicity is not checked");
        assert!(!crate::ops::is_valid(&p));
        assert_eq!(area(&p), -1400.0);
    }

    #[test]
    fn exact_path_gives_the_exact_sign_of_a_sliver_the_compensated_sum_cannot_certify() {
        // Triangle (0,0), (1e8+1, 1e8) , (1e8, 1e8-1) shifted far away: true
        // area is tiny relative to the terms, so certification must fall back.
        let o = 1.0e9_f64;
        let ccw = ring(&[(o, o), (o + 4.0, o + 1.0), (o + 3.0, o + 1.0)]);
        let cw = ring(&[(o, o), (o + 3.0, o + 1.0), (o + 4.0, o + 1.0)]);
        assert_eq!(signed_area_ring(&ccw), 0.5);
        assert_eq!(signed_area_ring(&cw), -0.5);
    }

    #[test]
    fn exactly_collinear_ring_far_from_the_origin_has_exactly_zero_area() {
        let o = 1.0e9_f64;
        let r = ring(&[(o, o), (o + 1.0, o + 1.0), (o + 2.0, o + 2.0)]);
        assert_eq!(signed_area_ring(&r), 0.0);
    }

    #[test]
    fn certification_accepts_a_clear_sign_and_refuses_one_inside_the_error_bound() {
        // Unit square: r = 2, terms sum to 4 -> far outside any bound.
        assert!(sign_is_certified(2.0, 4.0, 8.0));
        assert!(sign_is_certified(-2.0, 4.0, 8.0));
        // r no larger than gamma^2 * S: the sign cannot be trusted.
        let eps = f64::EPSILON * 0.5;
        let gamma = 8.0 * eps / (1.0 - 8.0 * eps);
        assert!(!sign_is_certified(gamma * gamma * 1.0e18, 1.0e18, 8.0));
        assert!(!sign_is_certified(0.0, 1.0, 8.0));
        assert!(!sign_is_certified(1.0e-30, 1.0e18, 8.0));
        // The bound grows with the term count.
        assert!(sign_is_certified(1.0e-3, 1.0e18, 8.0));
        assert!(!sign_is_certified(1.0e-3, 1.0e18, 8.0e6));
    }

    #[test]
    fn a_far_square_is_still_certified_without_the_exact_path() {
        let near = ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let far = ring(&[
            (1e9, 1e9),
            (1e9 + 1.0, 1e9),
            (1e9 + 1.0, 1e9 + 1.0),
            (1e9, 1e9 + 1.0),
        ]);
        let cert = |r: &LinearRing| {
            let (mut acc, mut s, mut m) = (Accum::new(), 0.0_f64, 0.0_f64);
            for (a, b) in r.edges() {
                acc.add_product(a.x, b.y);
                acc.add_product(-b.x, a.y);
                s += (a.x * b.y).abs() + (b.x * a.y).abs();
                m += 2.0;
            }
            sign_is_certified(acc.value(), s, m)
        };
        assert!(cert(&near));
        assert!(
            cert(&far),
            "area 1 vs S ~ 1e18: 4 gamma^2 S ~ 1e-12, still certified"
        );
    }

    #[test]
    fn exact_twice_area_matches_a_hand_computed_value() {
        let r = ring(&[(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (0.0, 2.0)]);
        assert_eq!(exact_twice_area_edges(r.edges()), 12.0);
    }

    /// The certified sum operates on edges so that `Polygon3d` can use it on
    /// projected edges, and the ring form is a wrapper over it. These bit
    /// patterns are the results of the earlier ring-only form, recorded on
    /// fixtures that cover both paths: the sliver at `1e9`, the far square
    /// and the thin annulus are certified by the compensated sum, and the
    /// sliver at `1e15` is not, so it takes the exact fallback. The test
    /// checks which path each takes, so the coverage claim is verified rather
    /// than assumed. A reordered summation would change the compensated
    /// value and could change which path fires; equality of bits proves the
    /// generalization reordered nothing.
    #[test]
    fn the_edge_form_reproduces_the_ring_form_bit_for_bit() {
        let o = 1.0e9_f64;
        let sliver = ring(&[(o, o), (o + 4.0, o + 1.0), (o + 3.0, o + 1.0)]);
        let far = ring(&[(o, o), (o + 1.0, o), (o + 1.0, o + 1.0), (o, o + 1.0)]);
        let ext = ring(&[
            (o, o),
            (o + 1000.0, o),
            (o + 1000.0, o + 1000.0),
            (o, o + 1000.0),
        ]);
        let d = 2f64.powi(-20);
        let hole = ring(&[
            (o + d, o + d),
            (o + d, o + 1000.0 - d),
            (o + 1000.0 - d, o + 1000.0 - d),
            (o + 1000.0 - d, o + d),
        ]);
        let p = 1.0e15_f64;
        let far_sliver = ring(&[(p, p), (p + 4.0, p + 1.0), (p + 3.0, p + 1.0)]);

        // Whether the compensated sum alone certifies the sign.
        let certified = |edges: Vec<(Coord2d, Coord2d)>| {
            let (mut acc, mut s, mut m) = (Accum::new(), 0.0_f64, 0.0_f64);
            for (a, b) in edges {
                acc.add_product(a.x, b.y);
                acc.add_product(-b.x, a.y);
                s += (a.x * b.y).abs() + (b.x * a.y).abs();
                m += 2.0;
            }
            sign_is_certified(acc.value(), s, m)
        };
        assert!(certified(sliver.edges().collect()));
        assert!(certified(far.edges().collect()));
        assert!(certified(ext.edges().chain(hole.edges()).collect()));
        assert!(
            !certified(far_sliver.edges().collect()),
            "takes the exact fallback"
        );

        let before: [(&str, f64, u64); 4] = [
            (
                "sliver",
                certified_twice_area(std::iter::once(&sliver)),
                0x3ff0_0000_0000_0000,
            ),
            (
                "far square",
                certified_twice_area(std::iter::once(&far)),
                0x4000_0000_0000_0000,
            ),
            (
                "thin annulus",
                certified_twice_area([&ext, &hole].into_iter()),
                0x3f7f_3fff_ff80_0000,
            ),
            (
                "sliver at 1e15",
                certified_twice_area(std::iter::once(&far_sliver)),
                0x3ff0_0000_0000_0000,
            ),
        ];
        for (name, now, bits) in before {
            assert_eq!(
                now.to_bits(),
                bits,
                "{name}: {now:e} vs {:e}",
                f64::from_bits(bits)
            );
        }
        // The edge form, called directly on the same edges, gives the same bits.
        assert_eq!(
            certified_twice_area_edges(ext.edges().chain(hole.edges())).to_bits(),
            0x3f7f_3fff_ff80_0000
        );
        assert_eq!(
            certified_twice_area_edges(far_sliver.edges()).to_bits(),
            0x3ff0_0000_0000_0000
        );
    }

    #[test]
    fn overflow_is_reported_as_non_finite_not_a_made_up_sign() {
        let r = ring(&[(0.0, 0.0), (1e200, 0.0), (1e200, 1e200), (0.0, 1e200)]);
        assert!(!signed_area_ring(&r).is_finite());
    }
}
