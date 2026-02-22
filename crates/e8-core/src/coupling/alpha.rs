//! Fine structure constant from the E8 continued fraction tower.
//!
//! 1/α = [244; 14, 13, 193] × e^{-γ} = (44665/183) × e^{-γ}
//!
//! The CF coefficients are Lie algebra invariants:
//! - a₀ = 244 = |Φ(E8)| + rank(E8)/2 (Killing form)
//! - a₁ = 14 = dim(G₂)
//! - a₂ = 13 = |W(G₂)| + 1 (prime)
//! - a₃ = 193 = |W(D₄)| + 1 (prime)

use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;
use crate::special::cf::{cf_to_float, ALPHA_CF_COEFFS};

/// Compute 1/α with overridable CF coefficients.
pub fn alpha_inverse_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let eng = S::exp_neg_gamma();

    if ctx.is_empty() {
        let cf_value: S = cf_to_float(&ALPHA_CF_COEFFS);
        return cf_value * eng;
    }

    // Allow overriding individual CF coefficients
    let a0 = ctx.get("alpha_cf_a0", ALPHA_CF_COEFFS[0] as f64) as u64;
    let a1 = ctx.get("alpha_cf_a1", ALPHA_CF_COEFFS[1] as f64) as u64;
    let a2 = ctx.get("alpha_cf_a2", ALPHA_CF_COEFFS[2] as f64) as u64;
    let a3 = ctx.get("alpha_cf_a3", ALPHA_CF_COEFFS[3] as f64) as u64;

    let cf_value: S = cf_to_float(&[a0, a1, a2, a3]);
    cf_value * eng
}

/// Compute 1/α at full precision.
pub fn alpha_inverse<S: Scalar>() -> S {
    alpha_inverse_with_ctx(&OverrideContext::new())
}

/// Compute α with overrides.
pub fn alpha_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    S::one() / alpha_inverse_with_ctx::<S>(ctx)
}

/// Compute α at full precision.
pub fn alpha<S: Scalar>() -> S {
    S::one() / alpha_inverse::<S>()
}

/// Verify the continued fraction decomposition.
pub fn verify_cf_decomposition() -> bool {
    let (p2, q2) = crate::special::cf::alpha_cf_rational_level2();
    let (p4, q4) = crate::special::cf::alpha_cf_rational();
    p2 == 44665 && q2 == 183 && p4 == 8623762 && q4 == 35333
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_alpha_inverse() {
        set_precision(50);
        let ai: DefaultScalar = alpha_inverse();
        let val = ai.to_f64();
        assert!(
            (val - 137.035999177).abs() < 0.001,
            "1/α = {}",
            val
        );
    }

    #[test]
    fn test_cf_decomposition() {
        assert!(verify_cf_decomposition());
    }
}
