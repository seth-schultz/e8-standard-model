//! Divisor sum functions.

/// Compute σ₃(k) = sum of cubes of divisors of k.
pub fn sigma3(k: u64) -> u64 {
    if k == 0 {
        return 0;
    }
    let mut sum = 0u64;
    let mut d = 1;
    while d * d <= k {
        if k % d == 0 {
            sum += d * d * d;
            let other = k / d;
            if other != d {
                sum += other * other * other;
            }
        }
        d += 1;
    }
    sum
}

/// Shell population for E8 lattice: N_k = 240 * σ₃(k).
pub fn shell_population(k: u64) -> u64 {
    240 * sigma3(k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sigma3() {
        assert_eq!(sigma3(1), 1);
        assert_eq!(sigma3(2), 1 + 8); // 1³ + 2³ = 9
        assert_eq!(sigma3(3), 1 + 27); // 1³ + 3³ = 28
        assert_eq!(sigma3(4), 1 + 8 + 64); // 1³ + 2³ + 4³ = 73
    }

    #[test]
    fn test_shell_population() {
        // Shell 1: N_1 = 240 * σ₃(1) = 240
        assert_eq!(shell_population(1), 240);
        // Shell 2: N_2 = 240 * σ₃(2) = 240 * 9 = 2160
        assert_eq!(shell_population(2), 2160);
    }
}
