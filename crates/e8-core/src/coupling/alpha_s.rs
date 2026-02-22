//! Strong coupling constant α_s(M_Z) from RGE running.
//!
//! M_GUT = sin²θ_W × m_P = (3/8) × m_P
//! 1/α_GUT = sin²θ_W × (44665/183) × e^{-γ} = 51.389
//!
//! Two-step one-loop RGE with flavor threshold at m_t:
//!   1/α_s(M_Z) = 1/α_GUT + b₃⁽⁶⁾/(2π) ln(m_t/M_GUT)
//!                          + b₃⁽⁵⁾/(2π) ln(M_Z/m_t)
//! where b₃⁽ⁿᶠ⁾ = (11N_c - 2n_f)/3:  b₃⁽⁶⁾ = 7, b₃⁽⁵⁾ = 23/3.

use crate::algebra::groups::identities::N_C;
use crate::algebra::groups::E8;
use crate::precision::scalar::Scalar;
use crate::special::cf::cf_to_float;

/// QCD beta function coefficient: b₃⁽ⁿᶠ⁾ = (11N_c - 2n_f)/3.
fn beta_qcd(n_f: u32) -> (u32, u32) {
    let numerator = 11 * N_C - 2 * n_f;
    (numerator, 3)
}

/// M_GUT = sin²θ_W(GUT) × m_P = (3/rank(E8)) × m_P in GeV.
pub fn m_gut_gev<S: Scalar>() -> S {
    let sin2 = S::from_u64(3) / S::from_u64(E8.rank as u64);
    sin2 * S::planck_mass_gev()
}

/// 1/α_GUT = sin²θ_W(GUT) × (44665/183) × e^{-γ} = 51.389.
pub fn alpha_gut_inverse<S: Scalar>() -> S {
    let sin2 = S::from_u64(3) / S::from_u64(E8.rank as u64);
    let cf: S = cf_to_float(&crate::special::cf::ALPHA_CF_COEFFS);
    let eng = S::exp_neg_gamma();
    sin2 * cf * eng
}

/// α_s(M_Z) via two-step one-loop RGE from the GUT scale.
///
/// Uses the E8-predicted top mass from Koide for the flavor threshold.
/// Two thresholds: n_f=6 above m_t, n_f=5 below m_t.
///
/// 1/α_s(M_Z) = 1/α_GUT + b₃⁽⁶⁾/(2π) ln(m_t/M_GUT)
///                        + b₃⁽⁵⁾/(2π) ln(M_Z/m_t)
pub fn alpha_s_mz<S: Scalar>() -> S {
    let masses = crate::mass::sectors::compute_all_masses::<S>();
    let m_t_gev = masses.top / S::from_u64(1000);
    alpha_s_mz_with_mt(&m_t_gev)
}

/// α_s(M_Z) with explicit top mass (in GeV) for the flavor threshold.
pub fn alpha_s_mz_with_mt<S: Scalar>(m_t_gev: &S) -> S {
    let m_z = S::from_f64(91.1876); // M_Z in GeV (experimental input)
    let two_pi = S::from_u64(2) * S::pi();

    // Beta coefficients: b₃ = (11N_c - 2n_f)/3
    let (b6_num, b6_den) = beta_qcd(6); // n_f=6 above m_t: (33-12)/3 = 7
    let (b5_num, b5_den) = beta_qcd(5); // n_f=5 below m_t: (33-10)/3 = 23/3
    let b3_6 = S::from_u64(b6_num as u64) / S::from_u64(b6_den as u64);
    let b3_5 = S::from_u64(b5_num as u64) / S::from_u64(b5_den as u64);

    let alpha_gut_inv: S = alpha_gut_inverse();
    let m_gut: S = m_gut_gev();

    // Step 1: GUT → m_t (6 active flavors)
    let log_mt_mgut = (m_t_gev.clone() / m_gut).ln();
    let step1 = b3_6 / two_pi.clone() * log_mt_mgut;

    // Step 2: m_t → M_Z (5 active flavors)
    let log_mz_mt = (m_z / m_t_gev.clone()).ln();
    let step2 = b3_5 / two_pi * log_mz_mt;

    let alpha_s_inv = alpha_gut_inv + step1 + step2;
    S::one() / alpha_s_inv
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_m_gut() {
        set_precision(50);
        let mg: rug::Float = m_gut_gev();
        let val = mg.to_f64();
        // M_GUT ≈ 4.58 × 10¹⁸ GeV
        assert!(val > 4e18 && val < 5e18, "M_GUT = {:.3e}", val);
    }

    #[test]
    fn test_alpha_s() {
        set_precision(50);
        let as_mz: rug::Float = alpha_s_mz();
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
