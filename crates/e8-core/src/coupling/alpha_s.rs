//! Strong coupling constant α_s(M_Z) from RGE running.
//!
//! M_GUT = sin²θ_W × m_P = (3/8) × m_P
//! 1/α_GUT = sin²θ_W × (44665/183) × e^{-γ} = 51.389
//!
//! Two-step one-loop RGE with flavor threshold at m_t:
//!   1/α_s(M_Z) = 1/α_GUT + b₃⁽⁶⁾/(2π) ln(m_t/M_GUT)
//!                          + b₃⁽⁵⁾/(2π) ln(M_Z/m_t)
//! where b₃⁽ⁿᶠ⁾ = (11N_c - 2n_f)/3:  b₃⁽⁶⁾ = 7, b₃⁽⁵⁾ = 23/3.

use rug::Float;

use crate::algebra::groups::identities::N_C;
use crate::algebra::groups::E8;
use crate::precision::{exp_neg_gamma, ln, pi, planck_mass_gev, precision_bits};
use crate::special::cf::cf_to_float;

/// QCD beta function coefficient: b₃⁽ⁿᶠ⁾ = (11N_c - 2n_f)/3.
fn beta_qcd(n_f: u32) -> (u32, u32) {
    let numerator = 11 * N_C - 2 * n_f;
    (numerator, 3)
}

/// M_GUT = sin²θ_W(GUT) × m_P = (3/rank(E8)) × m_P in GeV.
pub fn m_gut_gev() -> Float {
    let prec = precision_bits();
    let sin2 = Float::with_val(prec, 3) / Float::with_val(prec, E8.rank);
    sin2 * planck_mass_gev()
}

/// 1/α_GUT = sin²θ_W(GUT) × (44665/183) × e^{-γ} = 51.389.
pub fn alpha_gut_inverse() -> Float {
    let prec = precision_bits();
    let sin2 = Float::with_val(prec, 3) / Float::with_val(prec, E8.rank);
    let cf = cf_to_float::<Float>(&crate::special::cf::ALPHA_CF_COEFFS);
    let eng = exp_neg_gamma();
    Float::with_val(prec, Float::with_val(prec, sin2 * cf) * eng)
}

/// α_s(M_Z) via two-step one-loop RGE from the GUT scale.
///
/// Uses the E8-predicted top mass from Koide for the flavor threshold.
/// Two thresholds: n_f=6 above m_t, n_f=5 below m_t.
///
/// 1/α_s(M_Z) = 1/α_GUT + b₃⁽⁶⁾/(2π) ln(m_t/M_GUT)
///                        + b₃⁽⁵⁾/(2π) ln(M_Z/m_t)
pub fn alpha_s_mz() -> Float {
    // Use the E8-predicted top mass
    let masses = crate::mass::sectors::compute_all_masses::<Float>();
    let m_t_mev = &masses.top;
    let m_t_gev = Float::with_val(precision_bits(), m_t_mev / Float::with_val(precision_bits(), 1000));
    alpha_s_mz_with_mt(&m_t_gev)
}

/// α_s(M_Z) with explicit top mass (in GeV) for the flavor threshold.
pub fn alpha_s_mz_with_mt(m_t_gev: &Float) -> Float {
    let prec = precision_bits();
    let m_z = Float::with_val(prec, 91.1876f64); // M_Z in GeV (experimental input)
    let two_pi = Float::with_val(prec, Float::with_val(prec, 2) * pi());

    // Beta coefficients: b₃ = (11N_c - 2n_f)/3
    let (b6_num, b6_den) = beta_qcd(6); // n_f=6 above m_t: (33-12)/3 = 7
    let (b5_num, b5_den) = beta_qcd(5); // n_f=5 below m_t: (33-10)/3 = 23/3
    let b3_6 = Float::with_val(prec, b6_num) / Float::with_val(prec, b6_den);
    let b3_5 = Float::with_val(prec, b5_num) / Float::with_val(prec, b5_den);

    let alpha_gut_inv = alpha_gut_inverse();
    let m_gut = m_gut_gev();

    // Step 1: GUT → m_t (6 active flavors)
    let log_mt_mgut = ln(&(Float::with_val(prec, m_t_gev / &m_gut)));
    let step1 = Float::with_val(prec, &b3_6 / &two_pi) * log_mt_mgut;

    // Step 2: m_t → M_Z (5 active flavors)
    let log_mz_mt = ln(&(Float::with_val(prec, &m_z / m_t_gev)));
    let step2 = Float::with_val(prec, &b3_5 / &two_pi) * log_mz_mt;

    let alpha_s_inv = Float::with_val(prec, &alpha_gut_inv + &step1) + step2;
    Float::with_val(prec, 1) / alpha_s_inv
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_m_gut() {
        set_precision(50);
        let mg = m_gut_gev();
        let val = mg.to_f64();
        // M_GUT ≈ 4.58 × 10¹⁸ GeV
        assert!(val > 4e18 && val < 5e18, "M_GUT = {:.3e}", val);
    }

    #[test]
    fn test_alpha_s() {
        set_precision(50);
        let as_mz = alpha_s_mz();
        let val = as_mz.to_f64();
        // Paper: α_s(M_Z) = 0.11794 (PDG: 0.1180 ± 0.0009)
        assert!(
            (val - 0.1179).abs() < 0.002,
            "α_s(M_Z) = {}",
            val
        );
    }

    #[test]
    fn test_beta_coefficients() {
        // b₃⁽⁶⁾ = (33-12)/3 = 21/3 = 7
        let (num6, den6) = beta_qcd(6);
        assert_eq!(num6, 21); // 11×3 - 2×6 = 33 - 12 = 21
        assert_eq!(den6, 3);
        assert_eq!(num6 / den6, 7); // integer result

        // b₃⁽⁵⁾ = (33-10)/3 = 23/3
        let (num5, den5) = beta_qcd(5);
        assert_eq!(num5, 23);
        assert_eq!(den5, 3);
    }
}
