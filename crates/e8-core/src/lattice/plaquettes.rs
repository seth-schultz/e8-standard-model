//! Plaquette geometry: triangular plaquettes in the E8 root system.
//!
//! A plaquette is a triple (α, β, γ) of roots with α + β + γ = 0.
//! Equivalently, γ = -(α + β), and all three are roots.

use super::roots::Root;

/// Find all triangular plaquettes {α, β, γ} with α + β + γ = 0.
/// Returns indices into the root array (unordered triples).
/// Expected: 2,240 plaquettes = 240 × 28 / 3.
pub fn find_plaquettes(roots: &[Root]) -> Vec<[usize; 3]> {
    use std::collections::HashMap;

    // Build lookup: coords → index
    let mut lookup: HashMap<[i8; 8], usize> = HashMap::with_capacity(roots.len());
    for (idx, root) in roots.iter().enumerate() {
        lookup.insert(root.coords, idx);
    }

    let mut plaquettes = Vec::new();
    let n = roots.len();

    for i in 0..n {
        for j in (i + 1)..n {
            // Check if -(α + β) is also a root
            let mut neg_sum = [0i8; 8];
            for (k, val) in neg_sum.iter_mut().enumerate() {
                *val = -(roots[i].coords[k].wrapping_add(roots[j].coords[k]));
            }

            if let Some(&k) = lookup.get(&neg_sum) {
                if k > j {
                    // Only count each unordered triple once (i < j < k)
                    plaquettes.push([i, j, k]);
                }
            }
        }
    }

    plaquettes
}

/// Verify plaquette properties:
/// 1. Each plaquette has all inner products = -1 between its roots
/// 2. Count = 2,240
/// 3. Each root participates in 28 plaquettes
pub fn verify_plaquettes(roots: &[Root], plaquettes: &[[usize; 3]]) -> PlaquetteStats {
    let mut per_root = vec![0u32; roots.len()];
    let mut all_ip_minus_one = true;

    for &[i, j, k] in plaquettes {
        per_root[i] += 1;
        per_root[j] += 1;
        per_root[k] += 1;

        // Check inner products (times 4): should be -4 (i.e., ip = -1)
        let ip_ij = roots[i].inner_product_times_4(&roots[j]);
        let ip_jk = roots[j].inner_product_times_4(&roots[k]);
        let ip_ik = roots[i].inner_product_times_4(&roots[k]);

        if ip_ij != -4 || ip_jk != -4 || ip_ik != -4 {
            all_ip_minus_one = false;
        }
    }

    PlaquetteStats {
        count: plaquettes.len(),
        per_root_count: per_root[0], // should be 28 for each root
        all_per_root_equal: per_root.iter().all(|&c| c == per_root[0]),
        all_ip_minus_one,
    }
}

#[derive(Debug)]
pub struct PlaquetteStats {
    pub count: usize,
    pub per_root_count: u32,
    pub all_per_root_equal: bool,
    pub all_ip_minus_one: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::roots::generate_e8_roots;

    #[test]
    fn test_plaquette_count() {
        let roots = generate_e8_roots();
        let plaquettes = find_plaquettes(&roots);
        assert_eq!(plaquettes.len(), 2240);
    }

    #[test]
    fn test_plaquette_properties() {
        let roots = generate_e8_roots();
        let plaquettes = find_plaquettes(&roots);
        let stats = verify_plaquettes(&roots, &plaquettes);

        assert_eq!(stats.count, 2240);
        assert_eq!(stats.per_root_count, 28);
        assert!(stats.all_per_root_equal);
        assert!(stats.all_ip_minus_one);
    }
}
