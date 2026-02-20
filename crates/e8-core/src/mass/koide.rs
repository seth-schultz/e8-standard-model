//! Koide parametrization: √m_k = M × (1 + r × cos(2πk/3 + φ)).
//!
//! CRITICAL: m_k = M² × val² even when val < 0 (sign-flip for quarks).

use rug::Float;

use crate::precision::{cos, pi, precision_bits, sqrt};

/// Koide sector parameters.
#[derive(Debug, Clone)]
pub struct KoideParams {
    /// r⁴ value (determines the Koide ratio)
    pub r_fourth: Float,
    /// Koide phase φ
    pub phi: Float,
    /// Sector mass sum Σ in MeV
    pub sigma: Float,
}

/// Compute three masses from Koide parametrization.
/// Returns [m_heavy, m_light, m_middle] for k=0, k=1, k=2.
///
/// Formula:
///   M² = 2Σ/(6 + 3r²)
///   √m_k = M × (1 + r × cos(2πk/3 + φ))
///   m_k = M² × val_k²  (val_k can be negative!)
pub fn koide_masses(params: &KoideParams) -> [Float; 3] {
    let prec = precision_bits();
    let two = Float::with_val(prec, 2);
    let three = Float::with_val(prec, 3);
    let six = Float::with_val(prec, 6);
    let two_pi_3 = &two * pi() / &three;

    // r = (r⁴)^{1/4}
    let r = sqrt(&sqrt(&params.r_fourth));

    // M² = 2Σ/(6 + 3r²)
    let r_sq = Float::with_val(prec, &r * &r);
    let denominator = Float::with_val(prec, &six + Float::with_val(prec, &three * &r_sq));
    let m_sq = Float::with_val(prec, &two * &params.sigma) / denominator;

    let mut masses = [
        Float::with_val(prec, 0),
        Float::with_val(prec, 0),
        Float::with_val(prec, 0),
    ];

    for k in 0..3 {
        let angle = Float::with_val(prec, &two_pi_3 * k as u32) + &params.phi;
        let val = Float::with_val(prec, 1) + &r * cos(&angle);
        // m_k = M² × val² — CRITICAL: this is val², so negative val gives positive mass
        masses[k] = Float::with_val(prec, Float::with_val(prec, &m_sq * &val) * &val);
    }

    masses
}

/// Compute the Koide quality parameter Q = (Σ√m)²/(3Σm).
/// For leptons: Q = 2/3 exactly (theorem).
pub fn koide_q(m1: &Float, m2: &Float, m3: &Float) -> Float {
    let prec = precision_bits();
    let s1 = sqrt(m1);
    let s2 = sqrt(m2);
    let s3 = sqrt(m3);

    let sum_sqrt = Float::with_val(prec, Float::with_val(prec, &s1 + &s2) + &s3);
    let sum_sqrt_sq = Float::with_val(prec, &sum_sqrt * &sum_sqrt);

    let sum_m = Float::with_val(prec, Float::with_val(prec, m1 + m2) + m3);
    let three = Float::with_val(prec, 3);

    sum_sqrt_sq / (three * sum_m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::{mpf_u64, set_precision};

    #[test]
    fn test_lepton_koide() {
        set_precision(50);
        let prec = precision_bits();

        // Experimental lepton masses
        let sigma = Float::with_val(prec, 0.51099895f64)
            + Float::with_val(prec, 105.6583755f64)
            + Float::with_val(prec, 1776.86f64);

        let params = KoideParams {
            r_fourth: mpf_u64(4), // r⁴ = 4 → r = √2
            phi: Float::with_val(prec, 2) / Float::with_val(prec, 9), // φ = 2/9
            sigma,
        };

        let masses = koide_masses(&params);
        // k=0 → tau, k=1 → electron, k=2 → muon
        let tau = masses[0].to_f64();
        let electron = masses[1].to_f64();
        let muon = masses[2].to_f64();

        // Check approximate values
        assert!(
            (tau - 1776.86).abs() / 1776.86 < 0.001,
            "tau = {}",
            tau
        );
        assert!(
            (muon - 105.658).abs() / 105.658 < 0.001,
            "muon = {}",
            muon
        );
        assert!(
            (electron - 0.511).abs() / 0.511 < 0.01,
            "electron = {}",
            electron
        );
    }
}
