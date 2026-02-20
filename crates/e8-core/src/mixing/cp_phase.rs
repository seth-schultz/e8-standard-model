//! CP phases from octonionic non-associativity.
//!
//! δ_CKM = m₂π/dim(G₂) = 5π/14 (from arg(C_Fritz) = π/dim(Im(O)) in Fano plane)
//! δ_PMNS = (N_gen × m₂)π/dim(G₂) = 15π/14 (CP complementarity)

use rug::Float;

use crate::algebra::groups::identities::{DIM_G2, DIM_IM_OCTONIONS, N_GEN};
use crate::algebra::groups::G2;
use crate::precision::{pi, precision_bits};

/// CKM CP phase: δ_CKM = m₂π/dim(G₂) = 5π/14 ≈ 64.3°.
/// Origin: arg(C_Fritz) = π/dim(Im(O)) = π/7 from octonionic associator [e₆,e₃,e₁].
/// Mapped to PDG convention: π/7 → 5π/14.
pub fn delta_ckm_rad() -> Float {
    let prec = precision_bits();
    let m2 = G2.exponents[1]; // m₂ = 5
    Float::with_val(prec, m2) * pi() / Float::with_val(prec, DIM_G2)
}

/// CKM CP phase in degrees.
pub fn delta_ckm_deg() -> Float {
    let prec = precision_bits();
    delta_ckm_rad() * Float::with_val(prec, 180) / pi()
}

/// PMNS CP phase: δ_PMNS = (N_gen × m₂)π/dim(G₂) = 15π/14 ≈ 192.9°.
/// From CP complementarity: sin(δ_CKM) = cos(π/7).
pub fn delta_pmns_rad() -> Float {
    let prec = precision_bits();
    let m2 = G2.exponents[1]; // m₂ = 5
    let coeff = N_GEN * m2;   // 3 × 5 = 15
    Float::with_val(prec, coeff) * pi() / Float::with_val(prec, DIM_G2)
}

/// PMNS CP phase in degrees.
pub fn delta_pmns_deg() -> Float {
    let prec = precision_bits();
    delta_pmns_rad() * Float::with_val(prec, 180) / pi()
}

/// The Fritzsch phase: arg(C) = π/dim(Im(O)) = π/7 from the octonionic cross product.
/// 2π/dim(G₂) = 2π/14 = π/7.
pub fn fritzsch_phase() -> Float {
    let prec = precision_bits();
    pi() / Float::with_val(prec, DIM_IM_OCTONIONS)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_delta_ckm() {
        set_precision(50);
        let d = delta_ckm_deg().to_f64();
        // PDG: 65.5 ± 2.8°
        assert!(
            (d - 64.286).abs() < 0.01,
            "δ_CKM = {:.3}°",
            d
        );
    }

    #[test]
    fn test_delta_pmns() {
        set_precision(50);
        let d = delta_pmns_deg().to_f64();
        // Measured: 197 ± 30°
        assert!(
            (d - 192.857).abs() < 0.01,
            "δ_PMNS = {:.3}°",
            d
        );
    }

    #[test]
    fn test_cp_complementarity() {
        set_precision(50);
        // sin(δ_CKM) = cos(π/7)
        let sin_ckm = crate::precision::sin(&delta_ckm_rad());
        let cos_pi7 = crate::precision::cos(&fritzsch_phase());
        let diff = (sin_ckm - cos_pi7).abs().to_f64();
        assert!(diff < 1e-40, "CP complementarity fails: diff = {}", diff);
    }
}
