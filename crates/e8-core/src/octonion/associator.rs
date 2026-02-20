//! Octonionic associator: [eₐ, e_b, e_c] = (eₐ e_b) e_c - eₐ (e_b e_c).

use super::multiply::{is_fano_triple, octonion_multiply};

/// Compute the associator [eₐ, e_b, e_c] for imaginary octonion units.
/// Returns a vector of (coefficient, index) pairs.
/// The associator is zero if and only if {a, b, c} is a Fano triple.
pub fn associator(a: u8, b: u8, c: u8) -> Vec<(i8, u8)> {
    // (eₐ eᵦ) eᵧ
    let (s1, r1) = octonion_multiply(a, b);
    let (s2, r2) = octonion_multiply(r1, c);
    let lhs_sign = s1 * s2;
    let lhs_idx = r2;

    // eₐ (eᵦ eᵧ)
    let (s3, r3) = octonion_multiply(b, c);
    let (s4, r4) = octonion_multiply(a, r3);
    let rhs_sign = s3 * s4;
    let rhs_idx = r4;

    // [eₐ, eᵦ, eᵧ] = lhs - rhs
    if lhs_idx == rhs_idx {
        let net = lhs_sign - rhs_sign;
        if net == 0 {
            vec![]
        } else {
            vec![(net, lhs_idx)]
        }
    } else {
        vec![(lhs_sign, lhs_idx), (-rhs_sign, rhs_idx)]
    }
}

/// Compute |[eₐ, e_b, e_c]|² for imaginary units.
/// Returns 0 for Fano triples, 4 for non-Fano triples.
pub fn associator_norm_squared(a: u8, b: u8, c: u8) -> u32 {
    let assoc = associator(a, b, c);
    let mut norm_sq: i32 = 0;
    for &(coeff, _idx) in &assoc {
        norm_sq += (coeff as i32) * (coeff as i32);
    }
    norm_sq as u32
}

/// Count all non-Fano triples from {1,...,7}.
/// There are C(7,3) = 35 unordered triples.
/// 7 are Fano triples, 28 are non-Fano.
/// With orderings: 168 non-Fano = 7 × 24 = dim(Im(O)) × |Φ(D₄)|
pub fn count_non_fano_triples() -> (u32, u32) {
    let mut fano_count = 0u32;
    let mut non_fano_count = 0u32;

    for a in 1..=7u8 {
        for b in (a + 1)..=7 {
            for c in (b + 1)..=7 {
                if is_fano_triple(a, b, c) {
                    fano_count += 1;
                } else {
                    non_fano_count += 1;
                }
            }
        }
    }
    (fano_count, non_fano_count)
}

/// Compute ||T||² = sum over all ordered triples of |[eₐ,eᵦ,eᵧ]|².
/// Should equal 672 = 7 × 96.
pub fn total_associator_norm_squared() -> u64 {
    let mut total = 0u64;
    for a in 1..=7u8 {
        for b in 1..=7 {
            if b == a {
                continue;
            }
            for c in 1..=7 {
                if c == a || c == b {
                    continue;
                }
                total += associator_norm_squared(a, b, c) as u64;
            }
        }
    }
    total
}

/// Verify that each output direction α receives S_α = 96 (perfectly uniform).
pub fn verify_uniform_distribution() -> bool {
    let mut s_alpha = [0u64; 8]; // index 1-7

    for a in 1..=7u8 {
        for b in 1..=7 {
            if b == a {
                continue;
            }
            for c in 1..=7 {
                if c == a || c == b {
                    continue;
                }
                let assoc = associator(a, b, c);
                for &(coeff, idx) in &assoc {
                    s_alpha[idx as usize] += (coeff as i32 * coeff as i32) as u64;
                }
            }
        }
    }

    // Each direction should have S_α = 96
    (1..=7).all(|i| s_alpha[i] == 96)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fano_triple_vanishes() {
        // [e₁, e₂, e₄] should be 0 (Fano triple)
        assert_eq!(associator_norm_squared(1, 2, 4), 0);
    }

    #[test]
    fn test_non_fano_nonzero() {
        // [e₆, e₃, e₁] should be nonzero (the generation assignment triple)
        let norm = associator_norm_squared(6, 3, 1);
        assert_eq!(norm, 4); // |assoc|² = 4, |assoc| = 2
    }

    #[test]
    fn test_triple_counts() {
        let (fano, non_fano) = count_non_fano_triples();
        assert_eq!(fano, 7);
        assert_eq!(non_fano, 28);
    }

    #[test]
    fn test_total_norm_squared() {
        let total = total_associator_norm_squared();
        assert_eq!(total, 672); // 7 × 96
    }

    #[test]
    fn test_uniform_distribution() {
        assert!(verify_uniform_distribution());
    }
}
