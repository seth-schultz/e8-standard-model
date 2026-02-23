//! The E8 mass formula: Σ = f × m_P × exp(-((A + ΔA)×R + δ)/dim(so(8))).
//!
//! QCD corrections: ΔA shifts the effective A-value for colored sectors.
//! - ΔA_leptons = 0 (no color charge)
//! - ΔA_neutrinos = 0 (no color charge)
//! - ΔA_up = -(α_s × C_F) / (π × 3 × |Φ(E₆)|) = -(α_s × 4/3) / (π × 216)
//! - ΔA_down = ΔA_up × (9 - 2/7) = ΔA_up × 61/7

use crate::algebra::groups::identities::N_C;
use crate::algebra::groups::{G2, SU3};
use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;

use super::constants::{delta, m_planck_mev, mertens_r, norm_factor, CF_FUNDAMENTAL, E6_SU3_FACTOR};

/// All four sector sums from the E8 mass formula.
pub struct SectorSums<S: Scalar> {
    pub leptons: S,  // A = dim(su(3))+1 = 9, f = 1
    pub up: S,       // A = dim(su(3)) = 8, f = 1/C₂(SU3,fund)
    pub down: S,     // A = dim(su(3))+1 = 9, f = N_c/C₂(SU3,fund)
    pub neutrino: S, // A = dim(G₂) = 14, f = √((|W(G₂)|-rank(G₂))/(|W(G₂)|+1))
}

/// Compute sector mass sum with f64 parameters for override support.
///
/// # Arguments
/// * `a` - The representation quantum number (A-value) as f64 for override support
/// * `f` - The lattice gauge dressing factor
/// * `norm` - Normalization factor (default 28 = dim(so(8)))
pub fn sector_sum_f64<S: Scalar>(a: f64, f: &S, norm: f64) -> S {
    let r: S = mertens_r();
    let d: S = delta();
    let n = S::from_f64(norm);
    let mp: S = m_planck_mev();

    let exponent = (S::from_f64(a) * r + d).neg() / n;
    let boltzmann = exponent.exp();

    f.clone() * mp * boltzmann
}

/// Compute sector mass sum: Σ = f × m_P × exp(-(A × R + δ)/28).
///
/// # Arguments
/// * `a` - The representation quantum number (A-value)
/// * `f` - The lattice gauge dressing factor
pub fn sector_sum<S: Scalar>(a: u32, f: &S) -> S {
    let r: S = mertens_r();
    let d: S = delta();
    let n28: S = norm_factor();
    let mp: S = m_planck_mev();

    let exponent = (S::from_u64(a as u64) * r + d).neg() / n28;
    let boltzmann = exponent.exp();

    f.clone() * mp * boltzmann
}

/// Compute sector sum with precise f = p/q rational factor.
pub fn sector_sum_rational<S: Scalar>(a: u32, f_num: u64, f_den: u64) -> S {
    let f = S::from_u64(f_num) / S::from_u64(f_den);
    sector_sum(a, &f)
}

/// α_s(M_Z) = 0.11794 — the E8-derived strong coupling (not a free parameter).
const ALPHA_S_MZ: f64 = 0.11794;

/// ΔA for the up quark sector: -(α_s × C_F) / (π × E6_SU3_FACTOR).
///
/// C_F = 4/3, E6_SU3_FACTOR = 216 = 3 × |Φ(E₆)|.
/// Physical origin: gluon dressing shifts the effective lattice distance.
pub fn delta_a_up() -> f64 {
    -(ALPHA_S_MZ * CF_FUNDAMENTAL) / (std::f64::consts::PI * E6_SU3_FACTOR)
}

/// ΔA for the down quark sector: ΔA_up × 61/7.
///
/// 61/7 = 9 - 2/7 where:
/// - 9 = dim(u(3)) = A_down
/// - 2/7 = rank(G₂)/dim(Im(O))
/// Enhanced correction because 5̄ representation couples more strongly.
pub fn delta_a_down() -> f64 {
    delta_a_up() * 61.0 / 7.0
}

/// Compute all four sector mass sums with overrides.
pub fn compute_all_sector_sums_with_ctx<S: Scalar>(ctx: &OverrideContext) -> SectorSums<S> {
    let (c2_num, c2_den) = SU3.c2_fundamental;

    let norm = ctx.get("norm_factor", 28.0);

    // QCD ΔA corrections (overridable)
    let da_up = ctx.get("delta_a_up", delta_a_up());
    let da_down = ctx.get("delta_a_down", delta_a_down());

    // A-values with QCD corrections (overridable)
    let a_lep = ctx.get("a_lepton", (SU3.dimension + 1) as f64);
    let a_up = ctx.get("a_up", SU3.dimension as f64 + da_up);
    let a_down = ctx.get("a_down", (SU3.dimension + 1) as f64 + da_down);
    let a_nu = ctx.get("a_neutrino", G2.dimension as f64);

    // f-factors (overridable)
    let f_lep_default = 1.0;
    let f_up_default = c2_den as f64 / c2_num as f64; // 3/4
    let f_down_default = (N_C * c2_den) as f64 / c2_num as f64; // 9/4
    let w_minus_r = G2.weyl_order - G2.rank as u64;
    let w_plus_1 = G2.weyl_order + 1;
    let f_nu_default = (w_minus_r as f64 / w_plus_1 as f64).sqrt(); // sqrt(10/13)

    let f_lep: S = S::from_f64(ctx.get("f_lepton", f_lep_default));
    let f_up: S = S::from_f64(ctx.get("f_up", f_up_default));
    let f_down: S = S::from_f64(ctx.get("f_down", f_down_default));
    let f_nu_val = ctx.get("f_neutrino", f_nu_default);
    let f_nu: S = if ctx.overrides().contains_key("f_neutrino") {
        S::from_f64(f_nu_val)
    } else {
        // Use high-precision computation for default
        (S::from_u64(w_minus_r) / S::from_u64(w_plus_1)).sqrt()
    };

    SectorSums {
        leptons: sector_sum_f64(a_lep, &f_lep, norm),
        up: sector_sum_f64(a_up, &f_up, norm),
        down: sector_sum_f64(a_down, &f_down, norm),
        neutrino: sector_sum_f64(a_nu, &f_nu, norm),
    }
}

/// Compute all four sector mass sums.
///
/// A-values from representation theory (with QCD ΔA corrections for quarks):
///   leptons (1 of SU5): A = dim(u(3)) = dim(su(3)) + 1 = 9
///   up (10 of SU5):     A = dim(su(3)) + ΔA_up ≈ 8 - 0.000232
///   down (5̄ of SU5):   A = dim(u(3)) + ΔA_down ≈ 9 - 0.00202
///   neutrinos:          A = dim(G₂) = 14
///
/// ΔA_up = -(α_s × C_F) / (π × 216),  ΔA_down = ΔA_up × 61/7
///
/// f-factors from color Casimir scaling:
///   f_lep  = 1 (color singlet)
///   f_up   = 1/C₂(SU3,fund) = 3/4
///   f_down = N_c/C₂(SU3,fund) = 9/4
///   f_ν    = √((|W(G₂)| - rank(G₂)) / (|W(G₂)| + 1)) = √(10/13)
pub fn compute_all_sector_sums<S: Scalar>() -> SectorSums<S> {
    compute_all_sector_sums_with_ctx(&OverrideContext::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::override_context::OverrideContext;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_lepton_sum() {
        crate::precision::set_precision(50);
        let sigma: DefaultScalar = sector_sum_rational(9, 1, 1);
        let val = sigma.to_f64();
        assert!(val > 1880.0 && val < 1890.0, "Σ_lep = {}", val);
    }

    #[test]
    fn test_up_sum() {
        crate::precision::set_precision(50);
        let sigma: DefaultScalar = sector_sum_rational(8, 3, 4);
        let val = sigma.to_f64();
        assert!(val > 170_000.0 && val < 180_000.0, "Σ_up = {}", val);
    }

    #[test]
    fn test_down_sum() {
        crate::precision::set_precision(50);
        let sigma: DefaultScalar = sector_sum_rational(9, 9, 4);
        let val = sigma.to_f64();
        assert!(val > 4200.0 && val < 4300.0, "Σ_down = {}", val);
    }

    #[test]
    fn test_a_values() {
        assert_eq!(SU3.dimension + 1, 9);
        assert_eq!(SU3.dimension, 8);
        assert_eq!(G2.dimension, 14);
    }

    #[test]
    fn test_f_factors() {
        let (num, den) = SU3.c2_fundamental;
        assert_eq!(den, 3);
        assert_eq!(num, 4);
        assert_eq!(N_C * den, 9);
    }

    #[test]
    fn test_override_changes_result() {
        crate::precision::set_precision(50);
        let default_sums = compute_all_sector_sums::<DefaultScalar>();
        let default_lep = default_sums.leptons.to_f64();

        let mut params = std::collections::HashMap::new();
        params.insert("a_lepton".to_string(), 10.0); // different A
        let ctx = OverrideContext::from_params(params).unwrap();
        let override_sums = compute_all_sector_sums_with_ctx::<DefaultScalar>(&ctx);
        let override_lep = override_sums.leptons.to_f64();

        assert!(
            (default_lep - override_lep).abs() > 1.0,
            "Override should change result: default={}, override={}",
            default_lep,
            override_lep
        );
    }

    #[test]
    fn test_delta_a_values() {
        let da_up = delta_a_up();
        // ΔA_up = -(0.11794 × 4/3) / (π × 216) ≈ -0.000232
        assert!(da_up < 0.0, "ΔA_up should be negative");
        assert!((da_up - (-0.000232)).abs() < 0.0001, "ΔA_up = {}", da_up);

        let da_down = delta_a_down();
        // ΔA_down = ΔA_up × 61/7 ≈ -0.00202
        assert!(da_down < 0.0, "ΔA_down should be negative");
        assert!((da_down / da_up - 61.0 / 7.0).abs() < 1e-10, "ratio = {}", da_down / da_up);
    }

    #[test]
    fn test_qcd_corrections_increase_quark_sums() {
        // QCD corrections (ΔA < 0) should slightly increase sector sums
        // because exp(-(A+ΔA)R/28) > exp(-AR/28) when ΔA < 0
        crate::precision::set_precision(50);

        // Without QCD corrections
        let mut params = std::collections::HashMap::new();
        params.insert("delta_a_up".to_string(), 0.0);
        params.insert("delta_a_down".to_string(), 0.0);
        let ctx_no_qcd = OverrideContext::from_params(params).unwrap();
        let sums_no_qcd = compute_all_sector_sums_with_ctx::<DefaultScalar>(&ctx_no_qcd);

        // With default QCD corrections
        let sums_qcd = compute_all_sector_sums::<DefaultScalar>();

        assert!(
            sums_qcd.up.to_f64() > sums_no_qcd.up.to_f64(),
            "QCD correction should increase Σ_up"
        );
        assert!(
            sums_qcd.down.to_f64() > sums_no_qcd.down.to_f64(),
            "QCD correction should increase Σ_down"
        );
        // Leptons should be identical
        assert!(
            (sums_qcd.leptons.to_f64() - sums_no_qcd.leptons.to_f64()).abs() < 1e-10,
            "Lepton sum should be unchanged"
        );
    }
}
