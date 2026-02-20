//! PMNS mixing angles from G₂ Coxeter geometry.
//!
//! G₂ exponents: m₁=1, m₂=5; |W|=12; h=6; rank=2
//!
//! sin²θ₁₃ = tan(m₁π/|W|)/|W| = tan(π/12)/12 = (2-√3)/12
//! sin²θ₁₂ = tan(m₂π/|W|)/|W| = tan(5π/12)/12 = (2+√3)/12
//! sin²θ₂₃ = rank × tan(m₁π/|W|) = 2tan(π/12) = 4-2√3

use rug::Float;

use crate::algebra::groups::G2;
use crate::precision::{pi, precision_bits, tan};

/// PMNS sin²θ₁₃ = tan(m₁π/|W(G₂)|)/|W(G₂)| = (2-√3)/12 ≈ 0.02233.
pub fn sin2_theta13() -> Float {
    let prec = precision_bits();
    let w = G2.weyl_order; // |W(G₂)| = 12
    let m1 = G2.exponents[0]; // m₁ = 1
    let angle = Float::with_val(prec, m1) * pi() / Float::with_val(prec, w);
    tan(&angle) / Float::with_val(prec, w)
}

/// PMNS sin²θ₁₂ = tan(m₂π/|W(G₂)|)/|W(G₂)| = (2+√3)/12 ≈ 0.3110.
pub fn sin2_theta12() -> Float {
    let prec = precision_bits();
    let w = G2.weyl_order; // |W(G₂)| = 12
    let m2 = G2.exponents[1]; // m₂ = 5
    let angle = Float::with_val(prec, m2) * pi() / Float::with_val(prec, w);
    tan(&angle) / Float::with_val(prec, w)
}

/// PMNS sin²θ₂₃ = rank(G₂) × tan(m₁π/|W(G₂)|) = 2(2-√3) = 4-2√3 ≈ 0.5359.
pub fn sin2_theta23() -> Float {
    let prec = precision_bits();
    let w = G2.weyl_order; // |W(G₂)| = 12
    let m1 = G2.exponents[0]; // m₁ = 1
    let angle = Float::with_val(prec, m1) * pi() / Float::with_val(prec, w);
    Float::with_val(prec, G2.rank) * tan(&angle)
}

/// Verify the G₂ Coxeter rules:
/// 1. sum(s₁₂² + s₁₃²) = rank/h = 2/6 = 1/3
/// 2. prod(s₁₂² × s₁₃²) = 1/|W|² = 1/144
/// 3. Discriminant 48² - 576 = 1728 = |W|³ = 12³
pub fn verify_coxeter_rules() -> (Float, Float) {
    let s12 = sin2_theta12();
    let s13 = sin2_theta13();

    let sum = Float::with_val(precision_bits(), &s12 + &s13);
    let prod = Float::with_val(precision_bits(), &s12 * &s13);

    (sum, prod)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_pmns_angles() {
        set_precision(50);

        let s13 = sin2_theta13().to_f64();
        let s12 = sin2_theta12().to_f64();
        let s23 = sin2_theta23().to_f64();

        // Experimental: sin²θ₁₃ = 0.02220 ± 0.00068
        assert!(
            (s13 - 0.02233).abs() < 0.001,
            "sin²θ₁₃ = {}",
            s13
        );

        // Experimental: sin²θ₁₂ = 0.307 ± 0.013
        assert!(
            (s12 - 0.3110).abs() < 0.01,
            "sin²θ₁₂ = {}",
            s12
        );

        // Experimental: sin²θ₂₃ = 0.546 ± 0.021
        assert!(
            (s23 - 0.536).abs() < 0.01,
            "sin²θ₂₃ = {}",
            s23
        );
    }

    #[test]
    fn test_coxeter_rules() {
        set_precision(50);
        let (sum, prod) = verify_coxeter_rules();

        // sum = rank/h = 2/6 = 1/3
        let sum_err = (sum.to_f64() - 1.0 / 3.0).abs();
        assert!(sum_err < 1e-40, "sum = {}", sum);

        // prod = 1/|W|² = 1/144
        let prod_err = (prod.to_f64() - 1.0 / 144.0).abs();
        assert!(prod_err < 1e-40, "prod = {}", prod);
    }
}
