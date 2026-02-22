//! E8 lattice root system, shells, traces, plaquette geometry, and D₄ coset decomposition.

pub mod roots;
pub mod plaquettes;
pub mod d4_coset;

use roots::Root;

// ═══════════════════════════════════════════════════════════════
// RootSystem trait — extensible root system abstraction
// ═══════════════════════════════════════════════════════════════

/// A root system of a simple Lie algebra.
///
/// Implementors provide the roots, dimension, and geometric properties.
/// The E8 root system is the default implementation used throughout this crate.
pub trait RootSystem {
    /// Spatial dimension of the root system.
    fn dimension(&self) -> usize;

    /// Generate all roots.
    fn roots(&self) -> Vec<Root>;

    /// Number of roots |Φ|.
    fn root_count(&self) -> usize {
        self.roots().len()
    }

    /// Verify all roots have the expected norm squared (times 4, for integer arithmetic).
    /// For simply-laced algebras (ADE), all roots have |α|² = 2, so norm_sq_times_4 = 8.
    fn verify_norms(&self) -> bool {
        self.roots().iter().all(|r| r.norm_squared_times_4() == 8)
    }

    /// Count triangular plaquettes {α, β, γ} with α + β + γ = 0.
    fn plaquette_count(&self) -> usize {
        plaquettes::find_plaquettes(&self.roots()).len()
    }

    /// Plaquettes per root.
    fn plaquettes_per_root(&self) -> usize {
        let n = self.root_count();
        if n == 0 { return 0; }
        self.plaquette_count() * 3 / n
    }
}

/// The E8 root system: 240 roots in 8 dimensions.
pub struct E8RootSystem;

impl RootSystem for E8RootSystem {
    fn dimension(&self) -> usize {
        8
    }

    fn roots(&self) -> Vec<Root> {
        roots::generate_e8_roots()
    }

    fn root_count(&self) -> usize {
        240
    }

    fn plaquette_count(&self) -> usize {
        2240
    }

    fn plaquettes_per_root(&self) -> usize {
        28 // dim(so(8))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_root_system_trait() {
        let rs = E8RootSystem;
        assert_eq!(rs.dimension(), 8);
        assert_eq!(rs.root_count(), 240);
        assert!(rs.verify_norms());
        assert_eq!(rs.plaquette_count(), 2240);
        assert_eq!(rs.plaquettes_per_root(), 28);
    }
}
