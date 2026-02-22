//! PMNS mixing angles from G₂ Coxeter geometry.
//!
//! G₂ exponents: m₁=1, m₂=5; |W|=12; h=6; rank=2
//!
//! sin²θ₁₃ = tan(m₁π/|W|)/|W| = tan(π/12)/12 = (2-√3)/12
//! sin²θ₁₂ = tan(m₂π/|W|)/|W| = tan(5π/12)/12 = (2+√3)/12
//! sin²θ₂₃ = rank × tan(m₁π/|W|) = 2tan(π/12) = 4-2√3

use crate::algebra::groups::G2;
use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;

/// PMNS sin²θ₁₃ with overrides.
pub fn sin2_theta13_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let w = ctx.get("pmns_weyl_order", G2.weyl_order as f64);
    let m1 = ctx.get("pmns_exponent_m1", G2.exponents[0] as f64);
    let angle = S::from_f64(m1) * S::pi() / S::from_f64(w);
    angle.tan() / S::from_f64(w)
}

/// PMNS sin²θ₁₃ = tan(m₁π/|W(G₂)|)/|W(G₂)| = (2-√3)/12 ≈ 0.02233.
pub fn sin2_theta13<S: Scalar>() -> S {
    let w = G2.weyl_order; // |W(G₂)| = 12
    let m1 = G2.exponents[0]; // m₁ = 1
    let angle = S::from_u64(m1 as u64) * S::pi() / S::from_u64(w);
    angle.tan() / S::from_u64(w)
}

/// PMNS sin²θ₁₂ with overrides.
pub fn sin2_theta12_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let w = ctx.get("pmns_weyl_order", G2.weyl_order as f64);
    let m2 = ctx.get("pmns_exponent_m2", G2.exponents[1] as f64);
    let angle = S::from_f64(m2) * S::pi() / S::from_f64(w);
    angle.tan() / S::from_f64(w)
}

/// PMNS sin²θ₁₂ = tan(m₂π/|W(G₂)|)/|W(G₂)| = (2+√3)/12 ≈ 0.3110.
pub fn sin2_theta12<S: Scalar>() -> S {
    let w = G2.weyl_order; // |W(G₂)| = 12
    let m2 = G2.exponents[1]; // m₂ = 5
    let angle = S::from_u64(m2 as u64) * S::pi() / S::from_u64(w);
    angle.tan() / S::from_u64(w)
}

/// PMNS sin²θ₂₃ with overrides.
pub fn sin2_theta23_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let w = ctx.get("pmns_weyl_order", G2.weyl_order as f64);
    let m1 = ctx.get("pmns_exponent_m1", G2.exponents[0] as f64);
    let rank = ctx.get("pmns_rank", G2.rank as f64);
    let angle = S::from_f64(m1) * S::pi() / S::from_f64(w);
    S::from_f64(rank) * angle.tan()
}

/// PMNS sin²θ₂₃ = rank(G₂) × tan(m₁π/|W(G₂)|) = 2(2-√3) = 4-2√3 ≈ 0.5359.
pub fn sin2_theta23<S: Scalar>() -> S {
    let w = G2.weyl_order; // |W(G₂)| = 12
    let m1 = G2.exponents[0]; // m₁ = 1
    let angle = S::from_u64(m1 as u64) * S::pi() / S::from_u64(w);
    S::from_u64(G2.rank as u64) * angle.tan()
}

/// Verify the G₂ Coxeter rules:
/// 1. sum(s₁₂² + s₁₃²) = rank/h = 2/6 = 1/3
/// 2. prod(s₁₂² × s₁₃²) = 1/|W|² = 1/144
/// 3. Discriminant 48² - 576 = 1728 = |W|³ = 12³
pub fn verify_coxeter_rules<S: Scalar>() -> (S, S) {
    let s12: S = sin2_theta12();
    let s13: S = sin2_theta13();

    let sum = s12.clone() + s13.clone();
    let prod = s12 * s13;

    (sum, prod)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_pmns_angles() {
        set_precision(50);

        let s13: DefaultScalar = sin2_theta13();
        let s12: DefaultScalar = sin2_theta12();
        let s23: DefaultScalar = sin2_theta23();

        let s13 = s13.to_f64();
        let s12 = s12.to_f64();
        let s23 = s23.to_f64();

        assert!((s13 - 0.02233).abs() < 0.001, "sin²θ₁₃ = {}", s13);
        assert!((s12 - 0.3110).abs() < 0.01, "sin²θ₁₂ = {}", s12);
        assert!((s23 - 0.536).abs() < 0.01, "sin²θ₂₃ = {}", s23);
    }

    #[test]
    fn test_coxeter_rules() {
        set_precision(50);
        let (sum, prod): (DefaultScalar, DefaultScalar) = verify_coxeter_rules();

        let sum_err = (sum.to_f64() - 1.0 / 3.0).abs();
        // prod = 1/|W|² = 1/144
        let prod_err = (prod.to_f64() - 1.0 / 144.0).abs();

        #[cfg(feature = "arbitrary-precision")]
        {
            assert!(sum_err < 1e-40, "sum = {}", sum);
            assert!(prod_err < 1e-40, "prod = {}", prod);
        }
        #[cfg(not(feature = "arbitrary-precision"))]
        {
            assert!(sum_err < 1e-14, "sum = {}", sum);
            assert!(prod_err < 1e-14, "prod = {}", prod);
        }
    }
}
