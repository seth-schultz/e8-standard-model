//! Fine structure constant from the E8 continued fraction tower.
//!
//! 1/α = [244; 14, 13, 193] × e^{-γ} = (44665/183) × e^{-γ}
//!
//! The CF coefficients are Lie algebra invariants:
//! - a₀ = 244 = |Φ(E8)| + rank(E8)/2 (Killing form)
//! - a₁ = 14 = dim(G₂)
//! - a₂ = 13 = |W(G₂)| + 1 (prime)
//! - a₃ = 193 = |W(D₄)| + 1 (prime)

use rug::Float;

use crate::precision::{exp_neg_gamma, precision_bits};
use crate::special::cf::{cf_to_float, ALPHA_CF_COEFFS};

/// Compute 1/α at full precision.
pub fn alpha_inverse() -> Float {
    let cf_value = cf_to_float(&ALPHA_CF_COEFFS); // 44665/183
    let eng = exp_neg_gamma();
    Float::with_val(precision_bits(), cf_value * eng)
}

/// Compute α at full precision.
pub fn alpha() -> Float {
    let prec = precision_bits();
    Float::with_val(prec, 1) / alpha_inverse()
}

/// Verify the continued fraction decomposition:
/// Level-2: `[244;14,13]` = 44665/183
/// (44665 = 244 × 183 + 13, 183 = 14 × 13 + 1)
/// Full: `[244;14,13,193]` = 8623762/35333
pub fn verify_cf_decomposition() -> bool {
    let (p2, q2) = crate::special::cf::alpha_cf_rational_level2();
    let (p4, q4) = crate::special::cf::alpha_cf_rational();
    p2 == 44665 && q2 == 183 && p4 == 8623762 && q4 == 35333
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_alpha_inverse() {
        set_precision(50);
        let ai = alpha_inverse();
        let val = ai.to_f64();
        // Should be very close to 137.035999177
        assert!(
            (val - 137.035999177).abs() < 0.001,
            "1/α = {}",
            val
        );
    }

    #[test]
    fn test_cf_decomposition() {
        assert!(verify_cf_decomposition());
    }
}
