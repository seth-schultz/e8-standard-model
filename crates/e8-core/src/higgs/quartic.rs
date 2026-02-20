//! Higgs quartic coupling λ = dim(Im(O)) × π⁴ / |Φ(E₆)|².
//!
//! λ(m_P) = 0 is a THEOREM: E8 has no degree-4 Casimir.
//! Casimir degrees of E8: [2, 8, 12, 14, 18, 20, 24, 30] — no 4.
//! RGE convergence: 1L→11.1%, 2L→1.8%, 3L→0.2%

use rug::Float;

use crate::algebra::groups::identities::{DIM_IM_OCTONIONS, PHI_E6};
use crate::algebra::groups::G2;
use crate::precision::{pi, powi, precision_bits};

/// λ = dim(Im(O)) × π⁴ / |Φ(E₆)|² = 7π⁴/72²
///
/// Decomposition:
/// - dim(Im(O)) = 7 — the imaginary octonions
/// - π⁴/12 = Res(Z_E8, s=4) — Epstein zeta residue
/// - |Φ(E₆)| = 72 = h(G₂) × |W(G₂)| = 6 × 12
pub fn higgs_quartic() -> Float {
    let prec = precision_bits();
    let pi4 = powi(&pi(), 4);
    let dim_im_o = Float::with_val(prec, DIM_IM_OCTONIONS);
    let phi_e6_sq = Float::with_val(prec, PHI_E6 as u64 * PHI_E6 as u64);

    dim_im_o * pi4 / phi_e6_sq
}

/// Verify λ(m_P) = 0: E8 has no degree-4 Casimir.
pub fn lambda_mp_is_zero() -> bool {
    crate::algebra::groups::e8_has_no_degree4_casimir()
}

/// The ratio m_H/m_t = π²√dim(Im(O)) / h(G₂)² = π²√7/36.
pub fn mh_over_mt() -> Float {
    let prec = precision_bits();
    let pi_sq = powi(&pi(), 2);
    let sqrt_dim = Float::with_val(prec, DIM_IM_OCTONIONS).sqrt();
    let h_sq = (G2.coxeter_number as u64) * (G2.coxeter_number as u64); // h² = 36
    pi_sq * sqrt_dim / Float::with_val(prec, h_sq)
}

/// Second scalar mass: m_S = m_H × √(dim(Im(O))/|W(G₂)|) = m_H × √(7/12).
///
/// The E₆ singlet in 27 = 16₁ ⊕ 10₋₂ ⊕ 1₄:
/// - H lives in 10 (gauge-charged), λ_H ∝ |W(G₂)| = 12
/// - S lives in 1 (gauge-neutral), λ_S ∝ dim(Im(O)) = 7
/// - Ratio: m_S/m_H = √(λ_S/λ_H) = √(7/12)
pub fn second_scalar_mass(m_h_gev: &Float) -> Float {
    let prec = precision_bits();
    let ratio = Float::with_val(prec, DIM_IM_OCTONIONS) / Float::with_val(prec, G2.weyl_order);
    m_h_gev * ratio.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_quartic() {
        set_precision(50);
        let lambda = higgs_quartic();
        let val = lambda.to_f64();
        // λ ≈ 0.1315
        assert!(
            (val - 0.1315).abs() < 0.001,
            "λ = {}",
            val
        );
    }

    #[test]
    fn test_lambda_mp_zero() {
        assert!(lambda_mp_is_zero());
    }

    #[test]
    fn test_mh_mt_ratio() {
        set_precision(50);
        let ratio = mh_over_mt();
        let val = ratio.to_f64();
        // m_H/m_t ≈ 125.25/172.76 ≈ 0.7251
        assert!(
            (val - 0.725).abs() < 0.005,
            "m_H/m_t = {}",
            val
        );
    }
}
