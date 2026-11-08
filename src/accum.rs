// src/accum.rs

//! Compensated accumulation of sums and dot products.
//!
//! # Numerical argument
//!
//! [`Accum`] implements `Dot2` (and its special case `Sum2`) from
//! T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product",
//! SIAM J. Sci. Comput. 26(6), 2005, Algorithm 5.3. Each product is split
//! into the rounded product and its exact rounding error with `two_product`,
//! and each running-sum step is split the same way with `two_sum`; the
//! rounding errors are accumulated separately and added once at the end.
//!
//! The paper proves that, barring overflow and underflow, the result `r` of
//! summing `m` terms satisfies
//!
//! ```text
//! |r - exact| <= eps * |exact| + gamma_m^2 * sum(|terms|)
//! ```
//!
//! with `eps = 2^-53` and `gamma_m = m * eps / (1 - m * eps)`. That is the
//! accuracy of a computation carried out in twice the working precision and
//! then rounded once. The second term is of order `m^2 * eps^2`, which is
//! negligible next to `eps * sum(|terms|)` (the bound for a plain sum) for any
//! realistic `m`. The shoelace area of a polygon far from the origin is the
//! motivating case: the individual products are huge, the true area is small,
//! and a plain `f64` sum loses most or all of its significant digits to
//! cancellation, while `Dot2` does not.
//!
//! The bound does not hold if a `two_product` overflows (its Veltkamp split
//! multiplies by `2^27 + 1`, so operands above about `1.3e300` overflow) or
//! if a product's rounding error underflows (products below about `1e-292`).
//! In the overflow case the result is non-finite, not a finite wrong value:
//! an infinity can only become finite again through `inf - inf`, which is NaN.
//! In the underflow case the error terms are lost, so the result degrades to
//! a plain sum's accuracy at an absolute scale below `1e-292`.

use arithmetica::expansion::{two_product, two_sum};

/// Running compensated sum of products and values: Ogita–Rump–Oishi `Dot2`.
///
/// # What it guarantees
///
/// **An approximate magnitude with a proven bound. It certifies no sign by
/// itself.** T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product",
/// SIAM J. Sci. Comput. 26(6), 2005, Algorithm 5.3, proves that after `m`
/// terms (each [`add`](Accum::add) or [`add_product`](Accum::add_product)
/// is one term) the result `r` of [`value`](Accum::value) satisfies
///
/// ```text
/// |r − exact| ≤ ε·|exact| + γ_m²·Σ|terms|,   ε = 2⁻⁵³,  γ_m = m·ε / (1 − m·ε)
/// ```
///
/// which is the accuracy of a computation in twice the working precision,
/// rounded once. When `|exact|` is far smaller than `Σ|terms|`, which is
/// the cancellation this type exists for, the second term can exceed
/// `|exact|`, and then the sign of `r` is not the sign of the exact sum. A caller that
/// needs the sign must prove it separately. `ops::area` does: its
/// certification test compares `|r|` against four times the bound above and
/// falls back to exact expansion arithmetic when the test fails. That test is
/// private to `ops::area`; this type does not perform it.
///
/// # Where the bound fails
///
/// - **Overflow.** `two_product`'s Veltkamp split multiplies by `2²⁷ + 1`, so
///   operands above about `1.3e300` overflow. The result is then non-finite,
///   **never a finite wrong value**. An infinity becomes finite again only
///   through `inf − inf`, which is NaN.
/// - **Underflow.** When a product's rounding error underflows (products
///   below about `1e-292`) the error term is lost, and the result degrades
///   to a plain sum's accuracy at that absolute scale.
///
/// The module documentation carries the full argument.
#[derive(Clone, Copy, Debug, Default)]
pub struct Accum {
    /// Running rounded sum.
    p: f64,
    /// Running sum of the exact rounding errors.
    s: f64,
}

impl Accum {
    /// An empty accumulator: [`value`](Accum::value) is `0.0`.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds the plain value `a`.
    pub fn add(&mut self, a: f64) {
        let (p, q) = two_sum(self.p, a);
        self.p = p;
        self.s += q;
    }

    /// Adds the product `a * b`, whose rounding error is captured exactly.
    pub fn add_product(&mut self, a: f64, b: f64) {
        let (h, r) = two_product(a, b);
        let (p, q) = two_sum(self.p, h);
        self.p = p;
        self.s += q + r;
    }

    /// The accumulated value, rounded once. Within the bound stated on
    /// [`Accum`]; non-finite if any term overflowed.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.p + self.s
    }

    /// The accumulated value as an unevaluated pair `(p, s)` whose exact sum
    /// is the compensated value, with `s` the accumulated rounding errors.
    /// Using the pair instead of [`Accum::value`] keeps the extra precision
    /// for a later product.
    #[must_use]
    pub fn parts(&self) -> (f64, f64) {
        (self.p, self.s)
    }
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

    #[test]
    fn recovers_a_sum_that_a_plain_sum_loses() {
        // 1e16 + 1 - 1e16 is 0 in plain f64 (the 1 is absorbed) but 1 exactly.
        let mut acc = Accum::new();
        acc.add(1e16);
        acc.add(1.0);
        acc.add(-1e16);
        assert_eq!(acc.value(), 1.0);
        assert_eq!((1e16_f64 + 1.0) - 1e16, 0.0);
    }

    #[test]
    fn recovers_a_dot_product_that_a_plain_one_loses() {
        // (2^27+1)^2 = 2^54 + 2^28 + 1; the +1 is below the rounded product.
        let a = 134_217_729.0_f64;
        let mut acc = Accum::new();
        acc.add_product(a, a);
        acc.add_product(-134_217_728.0, 134_217_730.0); // 2^27 * (2^27 + 2)
        // exact: (2^54 + 2^28 + 1) - (2^54 + 2^28) = 1
        assert_eq!(acc.value(), 1.0);
    }

    #[test]
    fn overflow_gives_a_non_finite_value() {
        let mut acc = Accum::new();
        acc.add_product(1e300, 1e300);
        assert!(!acc.value().is_finite());
    }
}
