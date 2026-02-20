//! Continued fraction expansion and evaluation.

use rug::Float;

use crate::precision::precision_bits;

/// Convert continued fraction coefficients [a0; a1, a2, ...] to rational p/q.
/// Returns (numerator, denominator).
pub fn cf_to_rational(coeffs: &[u64]) -> (u64, u64) {
    assert!(!coeffs.is_empty());
    let mut p_prev: u64 = 1;
    let mut p_curr: u64 = coeffs[0];
    let mut q_prev: u64 = 0;
    let mut q_curr: u64 = 1;
    for &a in &coeffs[1..] {
        let p_next = a * p_curr + p_prev;
        let q_next = a * q_curr + q_prev;
        p_prev = p_curr;
        p_curr = p_next;
        q_prev = q_curr;
        q_curr = q_next;
    }
    (p_curr, q_curr)
}

/// Evaluate continued fraction [a0; a1, a2, ...] as a high-precision float ratio.
pub fn cf_to_float(coeffs: &[u64]) -> Float {
    let (p, q) = cf_to_rational(coeffs);
    let prec = precision_bits();
    let pf = Float::with_val(prec, p);
    let qf = Float::with_val(prec, q);
    pf / qf
}

/// The alpha continued fraction tower: [244; 14, 13, 193].
/// These are Lie algebra invariants from the E8 ⊃ D4 ⊃ G2 chain.
pub const ALPHA_CF_COEFFS: [u64; 4] = [244, 14, 13, 193];

/// Evaluate 1/alpha CF rational at full depth (all 4 coefficients).
/// [244;14,13,193] = 8623762/35333
pub fn alpha_cf_rational() -> (u64, u64) {
    cf_to_rational(&ALPHA_CF_COEFFS)
}

/// The level-2 truncation [244;14,13] = 44665/183, used for α_GUT.
pub fn alpha_cf_rational_level2() -> (u64, u64) {
    cf_to_rational(&ALPHA_CF_COEFFS[..3])
}

/// Verify the Euclidean algorithm structure:
/// 44665 = 244*183 + 13, 183 = 14*13 + 1
pub fn verify_euclidean() -> bool {
    let (p, q) = alpha_cf_rational_level2();
    p == 44665
        && q == 183
        && 44665 == 244 * 183 + 13
        && 183 == 14 * 13 + 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cf_to_rational_full() {
        // Full 4-coefficient: [244;14,13,193] = 8623762/35333
        let (p, q) = cf_to_rational(&ALPHA_CF_COEFFS);
        assert_eq!(p, 8623762);
        assert_eq!(q, 35333);
    }

    #[test]
    fn test_cf_to_rational_level2() {
        // 3-coefficient truncation: [244;14,13] = 44665/183
        let (p, q) = alpha_cf_rational_level2();
        assert_eq!(p, 44665);
        assert_eq!(q, 183);
    }

    #[test]
    fn test_euclidean_structure() {
        assert!(verify_euclidean());
    }

    #[test]
    fn test_cf_to_float() {
        crate::precision::set_precision(50);
        let val = cf_to_float(&ALPHA_CF_COEFFS);
        let ratio = 8623762.0 / 35333.0;
        assert!((val.to_f64() - ratio).abs() < 1e-10);
    }
}
