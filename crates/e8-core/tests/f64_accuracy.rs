//! Verify f64 mode produces results within 1e-8 relative of MPFR reference.
//!
//! These reference values were computed with 50-digit MPFR precision.
//! The f64 path should agree to ~15 significant digits for most quantities,
//! except CKM elements (Newton solver is precision-sensitive).

use e8_core::precision::scalar::Scalar;
use e8_core::precision::DefaultScalar;

fn check_rel(name: &str, got: f64, reference: f64, tol: f64) {
    if reference == 0.0 {
        assert!(got.abs() < tol, "{name}: expected ~0, got {got:.6e}");
        return;
    }
    let rel = ((got - reference) / reference).abs();
    assert!(
        rel < tol,
        "{name}: relative error {rel:.3e} exceeds {tol:.0e}\n  got: {got:.15e}\n  ref: {reference:.15e}"
    );
}

#[test]
fn f64_gauge_couplings() {
    e8_core::precision::set_precision(50);

    let alpha_inv: DefaultScalar = e8_core::coupling::alpha::alpha_inverse();
    check_rel("1/α(0)", alpha_inv.to_f64(), 137.0359991770789, 1e-8);

    let s2w_gut: DefaultScalar = e8_core::coupling::weinberg::sin2_theta_w_gut();
    check_rel("sin²θ_W(GUT)", s2w_gut.to_f64(), 0.375, 1e-12);

    let s2w_mz: DefaultScalar = e8_core::coupling::weinberg::sin2_theta_w_mz();
    check_rel("sin²θ_W(M_Z)", s2w_mz.to_f64(), 0.2312159268200373, 1e-8);

    let alpha_s: DefaultScalar = e8_core::coupling::alpha_s::alpha_s_mz();
    check_rel("α_s(M_Z)", alpha_s.to_f64(), 0.1179413093345547, 1e-6);
}

#[test]
fn f64_masses() {
    e8_core::precision::set_precision(50);
    let masses = e8_core::mass::sectors::compute_all_masses::<DefaultScalar>();

    check_rel("m_e", masses.electron.to_f64(), 0.5109740402604931, 1e-8);
    check_rel("m_μ", masses.muon.to_f64(), 105.6542637585243, 1e-8);
    check_rel("m_τ", masses.tau.to_f64(), 1776.898345471456, 1e-8);
    check_rel("m_t", masses.top.to_f64(), 172502.4484437401, 1e-8);
    check_rel("Σ_lep", masses.sigma_lep.to_f64(), 1883.063583270241, 1e-8);
}

#[test]
fn f64_mixing_angles() {
    e8_core::precision::set_precision(50);

    let s12: DefaultScalar = e8_core::mixing::pmns::sin2_theta12();
    check_rel("sin²θ₁₂", s12.to_f64(), 0.3110042339640731, 1e-10);

    let s23: DefaultScalar = e8_core::mixing::pmns::sin2_theta23();
    check_rel("sin²θ₂₃", s23.to_f64(), 0.5358983848622454, 1e-10);

    let d_ckm: DefaultScalar = e8_core::mixing::cp_phase::delta_ckm_deg();
    check_rel("δ_CKM", d_ckm.to_f64(), 64.28571428571429, 1e-12);
}

#[test]
fn f64_higgs() {
    e8_core::precision::set_precision(50);

    let lambda: DefaultScalar = e8_core::higgs::quartic::higgs_quartic();
    check_rel("λ_H", lambda.to_f64(), 0.1315323374301730, 1e-8);

    let m_h: DefaultScalar = e8_core::higgs::mass::higgs_mass_default();
    check_rel("m_H", m_h.to_f64(), 125.1242618125252, 1e-6);
}
