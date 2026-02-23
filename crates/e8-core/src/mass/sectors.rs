//! Sector configurations for all fermion families.
//!
//! Each sector's Koide parameters are derived from E8 representation theory:
//! - r⁴ from SU(5) Yukawa structure (+ QCD corrections for quarks)
//! - φ from G₂ Coxeter geometry (+ QCD corrections for quarks)
//! - Σ from the mass formula with group-theoretic A and f
//!
//! QCD corrections to Koide parameters:
//! - Up quarks:   Δr⁴ = -α_s/(4π) × 15/14 × (1 - α_s/34)
//!                Δφ  = +α_s/(400π) × 57/56
//! - Down quarks: Δr⁴ = Δr⁴_up / 14
//!                Δφ  = -Δφ_up / √2

use crate::algebra::groups::{G2, SU5};
use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;

use super::formula::compute_all_sector_sums_with_ctx;
use super::koide::{koide_masses, KoideParams};

/// All 9 charged fermion masses + 3 neutrino masses.
#[derive(Debug, Clone)]
pub struct AllMasses<S: Scalar> {
    // Charged leptons (MeV)
    pub electron: S,
    pub muon: S,
    pub tau: S,

    // Up-type quarks (MeV)
    pub up: S,
    pub charm: S,
    pub top: S,

    // Down-type quarks (MeV)
    pub down: S,
    pub strange: S,
    pub bottom: S,

    // Neutrinos (meV)
    pub nu1: S,
    pub nu2: S,
    pub nu3: S,

    // Sector sums (MeV, except neutrinos in MeV)
    pub sigma_lep: S,
    pub sigma_up: S,
    pub sigma_down: S,
    pub sigma_nu: S,
}

/// α_s(M_Z) = 0.11794 — the E8-derived strong coupling.
const ALPHA_S_MZ: f64 = 0.11794;

/// QCD correction to r⁴ for up quarks.
///
/// Δr⁴_up = -α_s/(4π) × 15/14 × (1 - α_s/34)
///
/// E8 origins of coefficients:
/// - 15/14 = (dim(G₂)+1)/dim(G₂) = 15/14
/// - 34 = h(E₈) + rank(D₄) = 30 + 4
pub fn delta_r4_up() -> f64 {
    let alpha_over_4pi = ALPHA_S_MZ / (4.0 * std::f64::consts::PI);
    let coeff = 15.0 / 14.0; // (dim(G₂)+1)/dim(G₂)
    let higher_order = 1.0 - ALPHA_S_MZ / 34.0; // 34 = h(E₈) + rank(D₄)
    -alpha_over_4pi * coeff * higher_order
}

/// QCD correction to r⁴ for down quarks.
///
/// Δr⁴_down = Δr⁴_up / 14
///
/// 14 = dim(G₂) = dim(Aut(O))
pub fn delta_r4_down() -> f64 {
    delta_r4_up() / 14.0
}

/// QCD correction to φ for up quarks.
///
/// Δφ_up = +α_s/(400π) × 57/56
///
/// E8 origins of coefficients:
/// - 57/56 = (rank(D₄)×dim(G₂)+1)/(rank(D₄)×dim(G₂)) = (4×14+1)/(4×14)
/// - 400 = 20² where 20 = |Φ(SU5)|
pub fn delta_phi_up() -> f64 {
    let coeff = 57.0 / 56.0; // (rank(D₄)×dim(G₂)+1)/(rank(D₄)×dim(G₂))
    ALPHA_S_MZ / (400.0 * std::f64::consts::PI) * coeff
}

/// QCD correction to φ for down quarks.
///
/// Δφ_down = -Δφ_up / √2
///
/// √2 = minimal vector norm in E8 lattice.
pub fn delta_phi_down() -> f64 {
    -delta_phi_up() / std::f64::consts::SQRT_2
}

/// Compute all masses with overrides.
pub fn compute_all_masses_with_ctx<S: Scalar>(ctx: &OverrideContext) -> AllMasses<S> {
    let sums = compute_all_sector_sums_with_ctx::<S>(ctx);

    let h = G2.coxeter_number; // h(G₂) = 6
    let w = G2.weyl_order;     // |W(G₂)| = 12

    // QCD corrections to Koide parameters
    let dr4_up = delta_r4_up();
    let dr4_down = delta_r4_down();
    let dphi_up = delta_phi_up();
    let dphi_down = delta_phi_down();

    // ═══════════════════════════════════════════════════════════
    // Charged leptons: r⁴ = 4 = (√2)⁴, φ = 2/9 (Z₃ variational)
    // No QCD corrections (color singlet)
    // ═══════════════════════════════════════════════════════════
    let r4_lep = ctx.get("r4_lepton", 4.0);
    let phi_lep = ctx.get("phi_lepton", 2.0 / 9.0);
    let lep_params = KoideParams {
        r_fourth: S::from_f64(r4_lep),
        phi: S::from_f64(phi_lep),
        sigma: sums.leptons.clone(),
    };
    let lep = koide_masses(&lep_params);

    // ═══════════════════════════════════════════════════════════
    // Up quarks: r⁴ = dim(∧²(5)) + Δr⁴_up = 10 + Δr⁴_up
    //            φ = (h-1)⁴/h⁵ + Δφ_up = 625/7776 + Δφ_up
    // QCD corrections from gluon dressing of Yukawa couplings
    // ═══════════════════════════════════════════════════════════
    let n_su5 = SU5.rank + 1; // 5
    let dim_antisym = n_su5 * (n_su5 - 1) / 2; // C(5,2) = 10
    let h_minus_1 = h - 1; // 5
    let phi_up_default = (h_minus_1 as f64).powi(4) / (h as f64).powi(5) + dphi_up;

    let r4_up = ctx.get("r4_up", dim_antisym as f64 + dr4_up);
    let phi_up_val = ctx.get("phi_up", phi_up_default);
    let up_params = KoideParams {
        r_fourth: S::from_f64(r4_up),
        phi: S::from_f64(phi_up_val),
        sigma: sums.up.clone(),
    };
    let up = koide_masses(&up_params);

    // ═══════════════════════════════════════════════════════════
    // Down quarks: r⁴ = (10 - √2) + Δr⁴_down
    //              φ = 1/h(G₂) + Δφ_down = 1/6 + Δφ_down
    // QCD corrections suppressed by 1/dim(G₂) relative to up sector
    // ═══════════════════════════════════════════════════════════
    let r4_down_default = dim_antisym as f64 - std::f64::consts::SQRT_2 + dr4_down;
    let r4_down = ctx.get("r4_down", r4_down_default);
    let phi_down_val = ctx.get("phi_down", 1.0 / h as f64 + dphi_down);

    let down_params = KoideParams {
        r_fourth: S::from_f64(r4_down),
        phi: S::from_f64(phi_down_val),
        sigma: sums.down.clone(),
    };
    let down = koide_masses(&down_params);

    // ═══════════════════════════════════════════════════════════
    // Neutrinos: r⁴ = 4 (same as leptons)
    //            φ = 2/9 + π/|W(G₂)| (Majorana G₂ Weyl shift)
    // ═══════════════════════════════════════════════════════════
    let r4_nu = ctx.get("r4_neutrino", 4.0);
    let phi_nu_base = ctx.get("phi_neutrino_base", 2.0 / 9.0);
    let phi_nu_shift_den = ctx.get("phi_neutrino_shift_den", w as f64);
    let phi_nu = S::from_f64(phi_nu_base)
        + S::pi() / S::from_f64(phi_nu_shift_den);

    let sigma_nu_mev = sums.neutrino.clone() * S::from_f64(1e9);
    let nu_params = KoideParams {
        r_fourth: S::from_f64(r4_nu),
        phi: phi_nu,
        sigma: sigma_nu_mev,
    };
    let nu = koide_masses(&nu_params);

    AllMasses {
        tau: lep[0].clone(),
        electron: lep[1].clone(),
        muon: lep[2].clone(),

        top: up[0].clone(),
        up: up[1].clone(),
        charm: up[2].clone(),

        bottom: down[0].clone(),
        down: down[1].clone(),
        strange: down[2].clone(),

        nu3: nu[0].clone(),
        nu1: nu[1].clone(),
        nu2: nu[2].clone(),

        sigma_lep: sums.leptons,
        sigma_up: sums.up,
        sigma_down: sums.down,
        sigma_nu: sums.neutrino,
    }
}

/// Compute all masses from zero free parameters.
pub fn compute_all_masses<S: Scalar>() -> AllMasses<S> {
    compute_all_masses_with_ctx(&OverrideContext::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_all_masses() {
        set_precision(50);
        let m: AllMasses<DefaultScalar> = compute_all_masses();

        let e = m.electron.to_f64();
        let mu = m.muon.to_f64();
        let tau = m.tau.to_f64();
        assert!((e - 0.511).abs() / 0.511 < 0.01, "e = {}", e);
        assert!((mu - 105.658).abs() / 105.658 < 0.01, "mu = {}", mu);
        assert!((tau - 1776.86).abs() / 1776.86 < 0.01, "tau = {}", tau);

        let u = m.up.to_f64();
        let c = m.charm.to_f64();
        let t = m.top.to_f64();
        assert!(u > 1.5 && u < 3.5, "u = {}", u);
        assert!(c > 1200.0 && c < 1350.0, "c = {}", c);
        assert!(t > 170_000.0 && t < 176_000.0, "t = {}", t);

        let d = m.down.to_f64();
        let s = m.strange.to_f64();
        let b = m.bottom.to_f64();
        assert!(d > 3.0 && d < 6.0, "d = {}", d);
        assert!(s > 85.0 && s < 105.0, "s = {}", s);
        assert!(b > 4100.0 && b < 4300.0, "b = {}", b);
    }

    #[test]
    fn test_qcd_koide_corrections() {
        // Verify QCD corrections have the expected signs and magnitudes
        let dr4_up = delta_r4_up();
        let dr4_down = delta_r4_down();
        let dphi_up = delta_phi_up();
        let dphi_down = delta_phi_down();

        // Δr⁴_up < 0 (gluon dressing reduces effective coupling)
        assert!(dr4_up < 0.0, "Δr⁴_up = {}", dr4_up);
        // Δr⁴_down = Δr⁴_up / 14, same sign
        assert!(dr4_down < 0.0, "Δr⁴_down = {}", dr4_down);
        assert!((dr4_down * 14.0 - dr4_up).abs() < 1e-15, "ratio check");

        // Δφ_up > 0
        assert!(dphi_up > 0.0, "Δφ_up = {}", dphi_up);
        // Δφ_down = -Δφ_up/√2 < 0
        assert!(dphi_down < 0.0, "Δφ_down = {}", dphi_down);
        assert!(
            (dphi_down + dphi_up / std::f64::consts::SQRT_2).abs() < 1e-15,
            "Δφ_down relation"
        );

        // Corrections are small (perturbative: proportional to α_s/(4π))
        assert!(dr4_up.abs() < 0.02, "|Δr⁴_up| should be < 0.02");
        assert!(dphi_up.abs() < 0.001, "|Δφ_up| should be < 0.001");
    }
}
