//! Weinberg angle at GUT scale and M_Z.
//!
//! GUT: sin²θ_W = 3/8 = Tr(T₃²)/Tr(Q²) (THEOREM)
//! M_Z: sin²θ_W = (3/(|W(G₂)|+1)) × (1 + (h-1)α/(hπ)) (DERIVED)

use rug::Float;

use crate::algebra::groups::identities::W_G2_PLUS_1;
use crate::algebra::groups::{E8, G2};
use crate::precision::{pi, precision_bits};

use super::alpha::alpha;

/// sin²θ_W at the GUT scale = Tr(T₃²)/Tr(Q²) = 30/80 = 3/8 (THEOREM).
pub fn sin2_theta_w_gut() -> Float {
    let prec = precision_bits();
    // Tr(T₃²) = 30, Tr(Q²) = 80 = rank(E8) × 10
    Float::with_val(prec, 3) / Float::with_val(prec, E8.rank)
}

/// sin²θ_W at M_Z scale = (3/(|W(G₂)|+1)) × (1 + (h(G₂)-1)α/(h(G₂)π)).
///
/// Tree level: 30/130 = 3/13 from trace doubling at M_Z.
///   130 = 10 × (|W(G₂)|+1), 80 = 10 × rank(E₈)
/// Radiative correction: (h(G₂)-1)/h(G₂) = 5/6.
pub fn sin2_theta_w_mz() -> Float {
    let prec = precision_bits();
    let a = alpha();
    let h = G2.coxeter_number as u64; // h(G₂) = 6
    let h_minus_1 = h - 1;            // 5

    // Correction = (h-1)/h × α/π = 5α/(6π)
    let correction = Float::with_val(prec, h_minus_1) * &a
        / (Float::with_val(prec, h) * pi());

    // Tree = 3/(|W(G₂)|+1) = 3/13
    let tree = Float::with_val(prec, 3) / Float::with_val(prec, W_G2_PLUS_1);

    tree * (Float::with_val(prec, 1) + correction)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_sin2_gut() {
        set_precision(50);
        let s = sin2_theta_w_gut();
        assert!((s.to_f64() - 0.375).abs() < 1e-15);
    }

    #[test]
    fn test_sin2_mz() {
        set_precision(50);
        let s = sin2_theta_w_mz();
        let val = s.to_f64();
        // Experimental: 0.23122 ± 0.00003
        assert!(
            (val - 0.23122).abs() < 0.001,
            "sin²θ_W(M_Z) = {}",
            val
        );
    }
}
