//! Koide parametrization: √m_k = M × (1 + r × cos(2πk/3 + φ)).
//!
//! CRITICAL: m_k = M² × val² even when val < 0 (sign-flip for quarks).

use crate::precision::scalar::Scalar;

/// Koide sector parameters.
#[derive(Debug, Clone)]
pub struct KoideParams<S: Scalar> {
    /// r⁴ value (determines the Koide ratio)
    pub r_fourth: S,
    /// Koide phase φ
    pub phi: S,
    /// Sector mass sum Σ in MeV
    pub sigma: S,
}

/// Compute three masses from Koide parametrization.
/// Returns [m_heavy, m_light, m_middle] for k=0, k=1, k=2.
///
/// Formula:
///   M² = 2Σ/(6 + 3r²)
///   √m_k = M × (1 + r × cos(2πk/3 + φ))
///   m_k = M² × val_k²  (val_k can be negative!)
pub fn koide_masses<S: Scalar>(params: &KoideParams<S>) -> [S; 3] {
    let two = S::from_u64(2);
    let three = S::from_u64(3);
    let six = S::from_u64(6);
    let two_pi_3 = two.clone() * S::pi() / three.clone();

    // r = (r⁴)^{1/4}
    let r = params.r_fourth.sqrt().sqrt();

    // M² = 2Σ/(6 + 3r²)
    let r_sq = r.clone() * r.clone();
    let denominator = six + three * r_sq;
    let m_sq = two * params.sigma.clone() / denominator;

    // Compute masses for k=0,1,2
    let m0 = {
        let angle = params.phi.clone();
        let val = S::one() + r.clone() * angle.cos();
        m_sq.clone() * val.clone() * val
    };
    let m1 = {
        let angle = two_pi_3.clone() + params.phi.clone();
        let val = S::one() + r.clone() * angle.cos();
        m_sq.clone() * val.clone() * val
    };
    let m2 = {
        let angle = two_pi_3.clone() + two_pi_3 + params.phi.clone();
        let val = S::one() + r * angle.cos();
        m_sq * val.clone() * val
    };

    [m0, m1, m2]
}

/// Compute the Koide quality parameter Q = (Σ√m)²/(3Σm).
/// For leptons: Q = 2/3 exactly (theorem).
pub fn koide_q<S: Scalar>(m1: &S, m2: &S, m3: &S) -> S {
    let s1 = m1.sqrt();
    let s2 = m2.sqrt();
    let s3 = m3.sqrt();

    let sum_sqrt = s1 + s2 + s3;
    let sum_sqrt_sq = sum_sqrt.clone() * sum_sqrt;

    let sum_m = m1.clone() + m2.clone() + m3.clone();

    sum_sqrt_sq / (S::from_u64(3) * sum_m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::scalar::Scalar;

    #[test]
    fn test_lepton_koide() {
        set_precision(50);

        // Experimental lepton masses
        let sigma = rug::Float::from_f64(0.51099895)
            + rug::Float::from_f64(105.6583755)
            + rug::Float::from_f64(1776.86);

        let params = KoideParams {
            r_fourth: rug::Float::from_u64(4), // r⁴ = 4 → r = √2
            phi: rug::Float::from_u64(2) / rug::Float::from_u64(9), // φ = 2/9
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
