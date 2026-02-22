//! Weinberg angle at GUT scale and M_Z.
//!
//! GUT: sin²θ_W = 3/8 = Tr(T₃²)/Tr(Q²) (THEOREM)
//! M_Z: sin²θ_W = (3/(|W(G₂)|+1)) × (1 + (h-1)α/(hπ)) (DERIVED)

use crate::algebra::groups::identities::W_G2_PLUS_1;
use crate::algebra::groups::{E8, G2};
use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;

use super::alpha::alpha_with_ctx;

/// sin²θ_W at the GUT scale = Tr(T₃²)/Tr(Q²) = 30/80 = 3/8 (THEOREM).
pub fn sin2_theta_w_gut<S: Scalar>() -> S {
    S::from_u64(3) / S::from_u64(E8.rank as u64)
}

/// sin²θ_W at M_Z with overrides.
pub fn sin2_theta_w_mz_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let a: S = alpha_with_ctx(ctx);

    let tree_num = ctx.get("weinberg_tree_num", 3.0);
    let tree_den = ctx.get("weinberg_tree_den", W_G2_PLUS_1 as f64);
    let corr_num = ctx.get("weinberg_correction_num", (G2.coxeter_number - 1) as f64);
    let corr_den = ctx.get("weinberg_correction_den", G2.coxeter_number as f64);

    let correction = S::from_f64(corr_num) * a
        / (S::from_f64(corr_den) * S::pi());

    let tree = S::from_f64(tree_num) / S::from_f64(tree_den);

    tree * (S::one() + correction)
}

/// sin²θ_W at M_Z scale = (3/(|W(G₂)|+1)) × (1 + (h(G₂)-1)α/(h(G₂)π)).
///
/// Tree level: 30/130 = 3/13 from trace doubling at M_Z.
///   130 = 10 × (|W(G₂)|+1), 80 = 10 × rank(E₈)
/// Radiative correction: (h(G₂)-1)/h(G₂) = 5/6.
pub fn sin2_theta_w_mz<S: Scalar>() -> S {
    sin2_theta_w_mz_with_ctx(&OverrideContext::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_sin2_gut() {
        set_precision(50);
        let s: DefaultScalar = sin2_theta_w_gut();
        assert!((s.to_f64() - 0.375).abs() < 1e-15);
    }

    #[test]
    fn test_sin2_mz() {
        set_precision(50);
        let s: DefaultScalar = sin2_theta_w_mz();
        let val = s.to_f64();
        assert!(
            (val - 0.23122).abs() < 0.001,
            "sin²θ_W(M_Z) = {}",
            val
        );
    }
}
