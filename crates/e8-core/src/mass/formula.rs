//! The E8 mass formula: Σ = f × m_P × exp(-(A×R + δ)/dim(so(8))).

use crate::algebra::groups::identities::N_C;
use crate::algebra::groups::{G2, SU3};
use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;

use super::constants::{delta, m_planck_mev, mertens_r, norm_factor};

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

/// Compute all four sector mass sums with overrides.
pub fn compute_all_sector_sums_with_ctx<S: Scalar>(ctx: &OverrideContext) -> SectorSums<S> {
    let (c2_num, c2_den) = SU3.c2_fundamental;

    let norm = ctx.get("norm_factor", 28.0);

    // A-values (overridable)
    let a_lep = ctx.get("a_lepton", (SU3.dimension + 1) as f64);
    let a_up = ctx.get("a_up", SU3.dimension as f64);
    let a_down = ctx.get("a_down", (SU3.dimension + 1) as f64);
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
/// A-values from representation theory:
///   leptons (1 of SU5): A = dim(u(3)) = dim(su(3)) + 1 = 9
///   up (10 of SU5):     A = dim(su(3)) = 8
///   down (5̄ of SU5):   A = dim(u(3)) = dim(su(3)) + 1 = 9
///   neutrinos:          A = dim(G₂) = 14
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
}
