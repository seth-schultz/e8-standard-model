//! E8 root system: 240 roots in 8 dimensions.

/// An E8 root as 8 rational coordinates (stored as 2× to avoid fractions).
/// Each coordinate is either 0, ±2 (integer type) or ±1 (half-integer type).
/// Actual coordinates: root[i] / 2.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Root {
    /// Coordinates multiplied by 2 to keep everything integer.
    pub coords: [i8; 8],
}

impl Root {
    /// The actual coordinate value at index i.
    pub fn coord(&self, i: usize) -> f64 {
        self.coords[i] as f64 / 2.0
    }

    /// Norm squared: |α|² = sum(coords²)/4. Should be 2 for all E8 roots.
    pub fn norm_squared_times_4(&self) -> i32 {
        self.coords.iter().map(|&c| (c as i32) * (c as i32)).sum()
    }

    /// Inner product with another root, times 4.
    pub fn inner_product_times_4(&self, other: &Root) -> i32 {
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(&a, &b)| (a as i32) * (b as i32))
            .sum()
    }
}

/// Generate all 240 roots of the E8 lattice.
pub fn generate_e8_roots() -> Vec<Root> {
    let mut roots = Vec::with_capacity(240);

    // Type I: ±eᵢ ± eⱼ for i < j → 4 * C(8,2) = 4 * 28 = 112 roots
    // In doubled coordinates: two entries are ±2, rest are 0
    for i in 0..8 {
        for j in (i + 1)..8 {
            for &si in &[2i8, -2] {
                for &sj in &[2i8, -2] {
                    let mut coords = [0i8; 8];
                    coords[i] = si;
                    coords[j] = sj;
                    roots.push(Root { coords });
                }
            }
        }
    }

    // Type II: (±½, ±½, ..., ±½) with even number of minus signs → 128 roots
    // In doubled coordinates: all entries are ±1, with even number of -1s
    for mask in 0u16..256 {
        let minus_count = mask.count_ones();
        if minus_count % 2 != 0 {
            continue;
        }
        let mut coords = [1i8; 8];
        for bit in 0..8 {
            if mask & (1 << bit) != 0 {
                coords[bit] = -1;
            }
        }
        roots.push(Root { coords });
    }

    assert_eq!(roots.len(), 240);
    roots
}

/// Verify all roots have |α|² = 2.
pub fn verify_root_norms(roots: &[Root]) -> bool {
    roots
        .iter()
        .all(|r| r.norm_squared_times_4() == 8) // 8/4 = 2
}

/// Compute the number of root pairs with inner product = v/4.
/// For E8 shell 1:
///   <α,β> = 2 → 1 pair per root (itself)
///   <α,β> = 1 → 56 pairs per root
///   <α,β> = 0 → 126 pairs per root
///   <α,β> = -1 → 56 pairs per root
///   <α,β> = -2 → 1 pair per root (negative)
pub fn inner_product_distribution(roots: &[Root]) -> [(i32, usize); 5] {
    let n = roots.len();
    let mut counts = [0usize; 5]; // indices for ip = -2, -1, 0, 1, 2

    for i in 0..n {
        for j in (i + 1)..n {
            let ip4 = roots[i].inner_product_times_4(&roots[j]);
            let ip = ip4 / 4; // actual inner product (integer for E8)
            match ip {
                -2 => counts[0] += 1,
                -1 => counts[1] += 1,
                0 => counts[2] += 1,
                1 => counts[3] += 1,
                2 => counts[4] += 1,
                _ => {} // shouldn't happen for E8
            }
        }
    }

    [
        (-2, counts[0]),
        (-1, counts[1]),
        (0, counts[2]),
        (1, counts[3]),
        (2, counts[4]),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_root_count() {
        let roots = generate_e8_roots();
        assert_eq!(roots.len(), 240);
    }

    #[test]
    fn test_root_norms() {
        let roots = generate_e8_roots();
        assert!(verify_root_norms(&roots));
    }

    #[test]
    fn test_type_counts() {
        let roots = generate_e8_roots();
        let integer_count = roots
            .iter()
            .filter(|r| r.coords.iter().any(|&c| c.abs() == 2))
            .count();
        let half_int_count = roots.len() - integer_count;
        assert_eq!(integer_count, 112);
        assert_eq!(half_int_count, 128);
    }

    #[test]
    fn test_inner_product_distribution() {
        let roots = generate_e8_roots();
        let dist = inner_product_distribution(&roots);
        // For 240 roots: each root has 1 negative, 56 with ip=-1, 126 with ip=0,
        // 56 with ip=+1, 1 positive (itself is excluded since we count pairs).
        // Total pairs = C(240,2) = 28680
        // ip=-2: 120 (each root has 1 opposite = 240/2 pairs)
        assert_eq!(dist[0], (-2, 120));
        // ip=+1: 240*56/2 = 6720
        assert_eq!(dist[3], (1, 6720));
    }
}
