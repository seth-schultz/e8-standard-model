//! CP phases from octonionic non-associativity.
//!
//! δ_CKM = m₂π/dim(G₂) = 5π/14 (from arg(C_Fritz) = π/dim(Im(O)) in Fano plane)
//! δ_PMNS = (N_gen × m₂)π/dim(G₂) = 15π/14 (CP complementarity)

use crate::algebra::groups::identities::{DIM_G2, DIM_IM_OCTONIONS, N_GEN};
use crate::algebra::groups::G2;
use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;

/// CKM CP phase with overrides: δ_CKM = m₂ × fritzsch / 2 = m₂π/(2×dim(Im(O))).
/// Default: 5π/14 ≈ 64.3° (since dim(G₂) = 2 × dim(Im(O))).
pub fn delta_ckm_rad_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let m2 = G2.exponents[1]; // m₂ = 5
    let fritzsch: S = fritzsch_phase_with_ctx(ctx);
    S::from_u64(m2 as u64) * fritzsch / S::from_u64(2)
}

/// CKM CP phase: δ_CKM = m₂π/dim(G₂) = 5π/14 ≈ 64.3°.
/// Origin: arg(C_Fritz) = π/dim(Im(O)) = π/7 from octonionic associator `[e₆,e₃,e₁]`.
/// Mapped to PDG convention: π/7 → 5π/14.
pub fn delta_ckm_rad<S: Scalar>() -> S {
    let m2 = G2.exponents[1]; // m₂ = 5
    S::from_u64(m2 as u64) * S::pi() / S::from_u64(DIM_G2 as u64)
}

/// CKM CP phase in degrees.
pub fn delta_ckm_deg<S: Scalar>() -> S {
    delta_ckm_rad::<S>() * S::from_u64(180) / S::pi()
}

/// CKM CP phase in degrees with overrides.
pub fn delta_ckm_deg_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    delta_ckm_rad_with_ctx::<S>(ctx) * S::from_u64(180) / S::pi()
}

/// PMNS CP phase: δ_PMNS = (N_gen × m₂)π/dim(G₂) = 15π/14 ≈ 192.9°.
/// From CP complementarity: sin(δ_CKM) = cos(π/7).
pub fn delta_pmns_rad<S: Scalar>() -> S {
    let m2 = G2.exponents[1]; // m₂ = 5
    let coeff = N_GEN * m2;   // 3 × 5 = 15
    S::from_u64(coeff as u64) * S::pi() / S::from_u64(DIM_G2 as u64)
}

/// PMNS CP phase in degrees.
pub fn delta_pmns_deg<S: Scalar>() -> S {
    delta_pmns_rad::<S>() * S::from_u64(180) / S::pi()
}

/// The Fritzsch phase: arg(C) = π/dim(Im(O)) = π/7 from the octonionic cross product.
/// 2π/dim(G₂) = 2π/14 = π/7.
pub fn fritzsch_phase<S: Scalar>() -> S {
    S::pi() / S::from_u64(DIM_IM_OCTONIONS as u64)
}

/// Fritzsch phase with overrides.
pub fn fritzsch_phase_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let cp_num = ctx.get("cp_phase_numerator", 1.0);
    let cp_den = ctx.get("cp_phase_denominator", DIM_IM_OCTONIONS as f64);
    S::from_f64(cp_num) * S::pi() / S::from_f64(cp_den)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_delta_ckm() {
        set_precision(50);
        let d: DefaultScalar = delta_ckm_deg();
        let d = d.to_f64();
        // PDG: 65.5 ± 2.8°
        assert!(
            (d - 64.286).abs() < 0.01,
            "δ_CKM = {:.3}°",
            d
        );
    }

    #[test]
    fn test_delta_pmns() {
        set_precision(50);
        let d: DefaultScalar = delta_pmns_deg();
        let d = d.to_f64();
        // Measured: 197 ± 30°
        assert!(
            (d - 192.857).abs() < 0.01,
            "δ_PMNS = {:.3}°",
            d
        );
    }

    #[test]
    fn test_cp_complementarity() {
        set_precision(50);
        // sin(δ_CKM) = cos(π/7)
        let sin_ckm = delta_ckm_rad::<DefaultScalar>().sin();
        let cos_pi7 = fritzsch_phase::<DefaultScalar>().cos();
        let diff = (sin_ckm - cos_pi7).abs().to_f64();

        #[cfg(feature = "arbitrary-precision")]
        assert!(diff < 1e-40, "CP complementarity fails: diff = {}", diff);
        #[cfg(not(feature = "arbitrary-precision"))]
        assert!(diff < 1e-14, "CP complementarity fails: diff = {}", diff);
    }
}
