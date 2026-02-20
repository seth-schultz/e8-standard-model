//! Octonionic multiplication via the Fano plane.

/// The 7 Fano plane triples (i, j, k) where eᵢ × eⱼ = eₖ.
/// Convention: indices 1-7 for imaginary units e₁..e₇.
/// Using the standard Fano plane: (1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3)
pub const FANO_TRIPLES: [(u8, u8, u8); 7] = [
    (1, 2, 4),
    (2, 3, 5),
    (3, 4, 6),
    (4, 5, 7),
    (5, 6, 1),
    (6, 7, 2),
    (7, 1, 3),
];

/// Multiplication table: octonion_multiply(a, b) returns (sign, index).
/// For eₐ × e_b = sign * e_index.
/// Returns (0, 0) if a=0 or b=0 (real unit).
/// Returns (+1 or -1, index) for imaginary units.
pub fn octonion_multiply(a: u8, b: u8) -> (i8, u8) {
    if a == 0 {
        return (1, b);
    }
    if b == 0 {
        return (1, a);
    }
    if a == b {
        return (-1, 0); // eᵢ² = -1
    }

    // Search Fano triples
    for &(i, j, k) in &FANO_TRIPLES {
        if a == i && b == j {
            return (1, k);
        }
        if a == j && b == i {
            return (-1, k);
        }
        if a == j && b == k {
            return (1, i);
        }
        if a == k && b == j {
            return (-1, i);
        }
        if a == k && b == i {
            return (1, j);
        }
        if a == i && b == k {
            return (-1, j);
        }
    }

    unreachable!("Invalid octonion indices: {} {}", a, b);
}

/// Check if (a, b, c) is a Fano triple (cyclic order).
pub fn is_fano_triple(a: u8, b: u8, c: u8) -> bool {
    for &(i, j, k) in &FANO_TRIPLES {
        if (a == i && b == j && c == k)
            || (a == j && b == k && c == i)
            || (a == k && b == i && c == j)
        {
            return true;
        }
    }
    false
}

/// Get the sign of the Fano triple: eₐ × e_b = ±e_c.
/// Returns Some(+1) or Some(-1) if it's a Fano triple, None otherwise.
pub fn fano_sign(a: u8, b: u8, c: u8) -> Option<i8> {
    let (sign, result) = octonion_multiply(a, b);
    if result == c {
        Some(sign)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fano_triples() {
        // e₁ × e₂ = e₄
        assert_eq!(octonion_multiply(1, 2), (1, 4));
        // e₂ × e₁ = -e₄
        assert_eq!(octonion_multiply(2, 1), (-1, 4));
        // eᵢ² = -1
        for i in 1..=7 {
            assert_eq!(octonion_multiply(i, i), (-1, 0));
        }
    }

    #[test]
    fn test_is_fano_triple() {
        assert!(is_fano_triple(1, 2, 4));
        assert!(is_fano_triple(2, 4, 1)); // cyclic
        assert!(is_fano_triple(4, 1, 2)); // cyclic
        assert!(!is_fano_triple(1, 3, 4)); // not a Fano triple
    }

    #[test]
    fn test_all_products_defined() {
        for a in 1..=7u8 {
            for b in 1..=7u8 {
                if a != b {
                    let (sign, idx) = octonion_multiply(a, b);
                    assert!(sign == 1 || sign == -1);
                    assert!(idx >= 1 && idx <= 7);
                }
            }
        }
    }
}
