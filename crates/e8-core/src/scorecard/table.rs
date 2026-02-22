//! Full scorecard computation: all 49 quantities from E8.
//!
//! Matches the paper's Table 9 exactly:
//! - #1-4: Gauge couplings
//! - #5-13: Fermion masses (9)
//! - #14-17: Sector mass sums (4, including Σ_ν)
//! - #18-28: CKM matrix (9 elements + δ_CKM + J)
//! - #29-32: PMNS mixing (3 angles + δ_PMNS)
//! - #33-37: Neutrino masses (3) + Δm² (2)
//! - #38-41: Higgs, second scalar, QCD vacuum
//! - #42-49: Structural predictions

use super::prediction::{Category, Prediction, Status};
use crate::coupling::alpha::alpha_inverse_with_ctx;
use crate::coupling::alpha_s::{alpha_s_mz, m_gut_gev};
use crate::coupling::weinberg::{sin2_theta_w_gut, sin2_theta_w_mz_with_ctx};
use crate::higgs::mass::higgs_mass_with_ctx;
use crate::higgs::quartic::{higgs_quartic_with_ctx, second_scalar_mass};
use crate::mass::sectors::compute_all_masses_with_ctx;
use crate::mixing::ckm::build_ckm_with_ctx;
use crate::mixing::cp_phase::{delta_ckm_deg_with_ctx, delta_pmns_deg};
use crate::mixing::pmns::{sin2_theta12_with_ctx, sin2_theta13_with_ctx, sin2_theta23_with_ctx};
use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;
use crate::precision::{set_precision, DefaultScalar};

/// Compute the full scorecard with parameter overrides.
pub fn compute_scorecard_with_ctx(digits: u32, ctx: &OverrideContext) -> Vec<Prediction> {
    set_precision(digits);
    let mut results = Vec::new();

    // ═══════════════════════════════════════════════════════════
    // #1-4: GAUGE COUPLINGS
    // ═══════════════════════════════════════════════════════════

    let alpha_inv = alpha_inverse_with_ctx::<DefaultScalar>(ctx).to_f64();
    results.push(Prediction::new(
        "1/α(0)",
        Category::Gauge,
        "[244;14,13,193] × e^{-γ}",
        alpha_inv,
        Some(137.035999177),
        Some(0.000000021),
        "",
        Status::DerivedStar,
    ));

    let s2w_mz = sin2_theta_w_mz_with_ctx::<DefaultScalar>(ctx).to_f64();
    results.push(Prediction::new(
        "sin²θ_W(M_Z)",
        Category::Gauge,
        "(3/13)×[1 + 5α/(6π)]",
        s2w_mz,
        Some(0.23122),
        Some(0.00003),
        "",
        Status::Derived,
    ));

    let as_val = alpha_s_mz::<DefaultScalar>().to_f64();
    results.push(Prediction::new(
        "α_s(M_Z)",
        Category::Gauge,
        "RGE from α_GUT with m_t threshold",
        as_val,
        Some(0.1180),
        Some(0.0009),
        "",
        Status::Derived,
    ));

    let s2w_gut = sin2_theta_w_gut::<DefaultScalar>().to_f64();
    results.push(Prediction::new(
        "sin²θ_W(GUT)",
        Category::Gauge,
        "Tr(T₃²)/Tr(Q²) = 30/80 = 3/8",
        s2w_gut,
        Some(0.375),
        None,
        "",
        Status::Theorem,
    ));

    // ═══════════════════════════════════════════════════════════
    // #5-13: FERMION MASSES (MeV)
    // ═══════════════════════════════════════════════════════════
    let masses = compute_all_masses_with_ctx::<DefaultScalar>(ctx);

    let m_e = masses.electron.to_f64();
    results.push(Prediction::new("m_e", Category::LeptonMass, "Koide(Σ_ℓ, √2, 2/9)", m_e, Some(0.51099895), Some(0.00000003), "MeV", Status::Derived));

    let m_mu = masses.muon.to_f64();
    results.push(Prediction::new("m_μ", Category::LeptonMass, "Koide(Σ_ℓ, √2, 2/9)", m_mu, Some(105.6583755), Some(0.0000023), "MeV", Status::Derived));

    let m_tau = masses.tau.to_f64();
    results.push(Prediction::new("m_τ", Category::LeptonMass, "Koide(Σ_ℓ, √2, 2/9)", m_tau, Some(1776.86), Some(0.12), "MeV", Status::Derived));

    let m_u = masses.up.to_f64();
    results.push(Prediction::new("m_u", Category::UpQuarkMass, "Koide(Σ_u, 10^{1/4}, 5⁴/6⁵)", m_u, Some(2.16), Some(0.49), "MeV", Status::Derived));

    let m_c = masses.charm.to_f64();
    results.push(Prediction::new("m_c", Category::UpQuarkMass, "Koide(Σ_u, 10^{1/4}, 5⁴/6⁵)", m_c, Some(1270.0), Some(20.0), "MeV", Status::Derived));

    let m_t = masses.top.to_f64();
    results.push(Prediction::new("m_t", Category::UpQuarkMass, "Koide(Σ_u, 10^{1/4}, 5⁴/6⁵)", m_t, Some(172760.0), Some(300.0), "MeV", Status::Derived));

    let m_d = masses.down.to_f64();
    results.push(Prediction::new("m_d", Category::DownQuarkMass, "Koide(Σ_d, (10-√2)^{1/4}, 1/6)", m_d, Some(4.67), Some(0.48), "MeV", Status::Derived));

    let m_s = masses.strange.to_f64();
    results.push(Prediction::new("m_s", Category::DownQuarkMass, "Koide(Σ_d, (10-√2)^{1/4}, 1/6)", m_s, Some(93.4), Some(8.6), "MeV", Status::Derived));

    let m_b = masses.bottom.to_f64();
    results.push(Prediction::new("m_b", Category::DownQuarkMass, "Koide(Σ_d, (10-√2)^{1/4}, 1/6)", m_b, Some(4180.0), Some(30.0), "MeV", Status::Derived));

    // ═══════════════════════════════════════════════════════════
    // #14-17: SECTOR MASS SUMS (MeV)
    // ═══════════════════════════════════════════════════════════

    let sigma_lep = masses.sigma_lep.to_f64();
    let sigma_lep_exp = 0.51099895 + 105.6583755 + 1776.86;
    results.push(Prediction::new("Σ_lep", Category::MassSum, "m_P × exp(-(9R+δ)/28)", sigma_lep, Some(sigma_lep_exp), Some(sigma_lep_exp * 0.0001), "MeV", Status::Derived));

    let sigma_up = masses.sigma_up.to_f64();
    let sigma_up_exp = 2.16 + 1270.0 + 172760.0;
    results.push(Prediction::new("Σ_up", Category::MassSum, "(3/4)m_P × exp(-(8R+δ)/28)", sigma_up, Some(sigma_up_exp), Some(sigma_up_exp * 0.002), "MeV", Status::Derived));

    let sigma_down = masses.sigma_down.to_f64();
    let sigma_down_exp = 4.67 + 93.4 + 4180.0;
    results.push(Prediction::new("Σ_down", Category::MassSum, "(9/4)m_P × exp(-(9R+δ)/28)", sigma_down, Some(sigma_down_exp), Some(sigma_down_exp * 0.01), "MeV", Status::Derived));

    let sigma_nu = masses.nu1.to_f64() + masses.nu2.to_f64() + masses.nu3.to_f64();
    results.push(Prediction::new("Σ_ν", Category::MassSum, "√(10/13)·m_P × exp(-(14R+δ)/28)", sigma_nu, None, None, "meV", Status::Derived));

    // ═══════════════════════════════════════════════════════════
    // #18-28: CKM MATRIX
    // ═══════════════════════════════════════════════════════════
    let ckm = build_ckm_with_ctx(&masses, ctx);

    let ckm_data = [
        ("V_ud", 0.97373, 0.00031),
        ("V_us", 0.2243, 0.0005),
        ("V_ub", 0.00382, 0.00020),
        ("V_cd", 0.221, 0.004),
        ("V_cs", 0.975, 0.006),
        ("V_cb", 0.0408, 0.0014),
        ("V_td", 0.0086, 0.0002),
        ("V_ts", 0.0415, 0.0009),
        ("V_tb", 1.014, 0.029),
    ];

    for (i, &(name, exp_val, exp_err)) in ckm_data.iter().enumerate() {
        results.push(Prediction::new(
            name, Category::CKM, "Fritzsch + octonionic CP",
            ckm.magnitudes[i].to_f64(), Some(exp_val), Some(exp_err), "", Status::Derived,
        ));
    }

    let d_ckm = delta_ckm_deg_with_ctx::<DefaultScalar>(ctx).to_f64();
    results.push(Prediction::new("δ_CKM", Category::CKM, "5π/14 from [e₆,e₃,e₁]", d_ckm, Some(65.5), Some(2.8), "deg", Status::Derived));

    results.push(Prediction::new(
        "J (Jarlskog)", Category::CKM, "Im(V_ud V_cs V*_us V*_cd)",
        ckm.jarlskog.to_f64(), Some(3.08e-5), Some(0.15e-5), "", Status::Derived,
    ));

    // ═══════════════════════════════════════════════════════════
    // #29-32: PMNS MIXING
    // ═══════════════════════════════════════════════════════════
    let s12 = sin2_theta12_with_ctx::<DefaultScalar>(ctx).to_f64();
    let s23 = sin2_theta23_with_ctx::<DefaultScalar>(ctx).to_f64();
    let s13 = sin2_theta13_with_ctx::<DefaultScalar>(ctx).to_f64();

    results.push(Prediction::new("sin²θ₁₂", Category::PMNS, "(2+√3)/12 from G₂ Coxeter", s12, Some(0.307), Some(0.013), "", Status::Derived));
    results.push(Prediction::new("sin²θ₂₃", Category::PMNS, "4-2√3 from G₂ rank×tan", s23, Some(0.546), Some(0.021), "", Status::Derived));
    results.push(Prediction::new("sin²θ₁₃", Category::PMNS, "(2-√3)/12 from G₂ Coxeter", s13, Some(0.02220), Some(0.00068), "", Status::Derived));

    let d_pmns = delta_pmns_deg::<DefaultScalar>().to_f64();
    results.push(Prediction::new("δ_PMNS", Category::PMNS, "15π/14", d_pmns, Some(197.0), Some(30.0), "deg", Status::Derived));

    // ═══════════════════════════════════════════════════════════
    // #33-37: NEUTRINO MASSES
    // ═══════════════════════════════════════════════════════════
    let nu1 = masses.nu1.to_f64();
    let nu2 = masses.nu2.to_f64();
    let nu3 = masses.nu3.to_f64();

    results.push(Prediction::new("m_ν₁", Category::Neutrino, "Koide(Σ_ν, √2, 2/9+π/12)", nu1, None, None, "meV", Status::Derived));
    results.push(Prediction::new("m_ν₂", Category::Neutrino, "Koide(Σ_ν, √2, 2/9+π/12)", nu2, None, None, "meV", Status::Derived));
    results.push(Prediction::new("m_ν₃", Category::Neutrino, "Koide(Σ_ν, √2, 2/9+π/12)", nu3, None, None, "meV", Status::Derived));

    let dm21 = (nu2 * nu2 - nu1 * nu1) * 1e-6;
    results.push(Prediction::new("Δm²₂₁", Category::Neutrino, "m₂² - m₁²", dm21, Some(7.53e-5), Some(0.18e-5), "eV²", Status::Derived));

    let dm31 = (nu3 * nu3 - nu1 * nu1) * 1e-6;
    results.push(Prediction::new("Δm²₃₁", Category::Neutrino, "m₃² - m₁²", dm31, Some(2.453e-3), Some(0.033e-3), "eV²", Status::Derived));

    // ═══════════════════════════════════════════════════════════
    // #38-41: HIGGS, SECOND SCALAR, QCD VACUUM
    // ═══════════════════════════════════════════════════════════

    let m_h = higgs_mass_with_ctx::<DefaultScalar>(ctx).to_f64();
    results.push(Prediction::new("m_H", Category::Higgs, "v×√(2λ), v=√2·m_t", m_h, Some(125.25), Some(0.17), "GeV", Status::Derived));

    let lambda = higgs_quartic_with_ctx::<DefaultScalar>(ctx).to_f64();
    results.push(Prediction::new("λ_H", Category::Higgs, "7π⁴/72²", lambda, Some(0.1315), Some(0.001), "", Status::Derived));

    let m_h_float: DefaultScalar = higgs_mass_with_ctx(ctx);
    let m_s_val = second_scalar_mass(&m_h_float).to_f64();
    results.push(Prediction::new("m_S", Category::Higgs, "m_H × √(7/12)", m_s_val, Some(95.4), Some(2.0), "GeV", Status::Derived));

    results.push(Prediction::new("θ̄_QCD", Category::Higgs, "0 (E8 parity + Fritzsch det)", 0.0, Some(0.0), None, "", Status::Theorem));

    // ═══════════════════════════════════════════════════════════
    // #42-49: STRUCTURAL PREDICTIONS
    // ═══════════════════════════════════════════════════════════

    results.push(Prediction::new("|Φ_E8|", Category::Structural, "E8 root count", 240.0, None, None, "", Status::Theorem));
    results.push(Prediction::new("Plaquettes", Category::Structural, "Triangular plaquettes in E8", 2240.0, None, None, "", Status::Theorem));
    results.push(Prediction::new("Per root", Category::Structural, "28 = dim(so(8))", 28.0, None, None, "", Status::Theorem));
    results.push(Prediction::new("Generations", Category::Structural, "3 from E8 ⊃ SU(3)_gen", 3.0, Some(3.0), None, "", Status::Theorem));
    results.push(Prediction::new("d = 8", Category::Structural, "Hurwitz ∩ Milnor = {8}", 8.0, None, None, "", Status::Theorem));

    // #47: y_t = √2 × m_t / v — top Yukawa coupling
    // v = 246.22 GeV (experimental VEV), m_t from E8 prediction
    let v_mev = DefaultScalar::from_f64(246220.0); // 246.22 GeV in MeV
    let sqrt2 = DefaultScalar::from_u64(2).sqrt();
    let y_t = (sqrt2 * masses.top.clone() / v_mev).to_f64();
    results.push(Prediction::new("y_t", Category::Structural, "√2·m_t/v ≈ 1 (E8 critical)", y_t, Some(0.991), None, "", Status::Derived));

    // #48: M_GUT
    let m_gut = m_gut_gev::<DefaultScalar>().to_f64();
    results.push(Prediction::new("M_GUT", Category::Structural, "(3/8) × m_P", m_gut, None, None, "GeV", Status::Derived));

    results.push(Prediction::new("Confinement", Category::Structural, "yes (d > 4 lattice theorem)", 1.0, Some(1.0), None, "", Status::Theorem));

    results
}

/// Compute the full scorecard of all 49 quantities (default E8 parameters).
pub fn compute_scorecard(digits: u32) -> Vec<Prediction> {
    compute_scorecard_with_ctx(digits, &OverrideContext::new())
}

/// Compute total chi-squared for predictions with experimental data.
pub fn chi_squared(predictions: &[Prediction]) -> f64 {
    predictions
        .iter()
        .filter_map(|p| {
            match (p.experimental, p.exp_uncertainty) {
                (Some(exp), Some(err)) if err > 0.0 => {
                    let pull = (p.predicted - exp) / err;
                    Some(pull * pull)
                }
                _ => None,
            }
        })
        .sum()
}

/// Print scorecard summary statistics.
pub fn scorecard_summary(predictions: &[Prediction]) -> ScorecardSummary {
    let total = predictions.len();
    let theorems = predictions.iter().filter(|p| p.status == Status::Theorem).count();
    let derived = predictions.iter().filter(|p| p.status == Status::Derived).count();
    let derived_star = predictions.iter().filter(|p| p.status == Status::DerivedStar).count();

    let with_pull: Vec<&Prediction> = predictions
        .iter()
        .filter(|p| p.pull_sigma.is_some())
        .collect();

    let within_1sigma = with_pull
        .iter()
        .filter(|p| p.pull_sigma.unwrap_or(0.0).abs() < 1.0)
        .count();
    let within_2sigma = with_pull
        .iter()
        .filter(|p| p.pull_sigma.unwrap_or(0.0).abs() < 2.0)
        .count();

    ScorecardSummary {
        total,
        theorems,
        derived,
        derived_star,
        with_experimental: with_pull.len(),
        within_1sigma,
        within_2sigma,
        chi_squared: chi_squared(predictions),
    }
}

/// Summary statistics for the scorecard.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ScorecardSummary {
    pub total: usize,
    pub theorems: usize,
    pub derived: usize,
    pub derived_star: usize,
    pub with_experimental: usize,
    pub within_1sigma: usize,
    pub within_2sigma: usize,
    pub chi_squared: f64,
}
