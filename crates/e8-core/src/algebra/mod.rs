//! Lie algebra invariants: groups, Casimirs, Weyl orders, embeddings.

pub mod groups;
pub mod embedding;

use groups::LieGroup;

// ═══════════════════════════════════════════════════════════════
// LieAlgebra trait — extensible Lie algebra abstraction
// ═══════════════════════════════════════════════════════════════

/// Invariants of a simple Lie algebra.
///
/// Provides access to the standard group-theoretic quantities
/// used throughout the E8 framework: rank, dimension, Coxeter number,
/// Weyl group order, exponents, and Casimir degrees.
pub trait LieAlgebra {
    /// Name of the algebra (e.g., "E8", "G2").
    fn name(&self) -> &str;

    /// Rank (dimension of Cartan subalgebra).
    fn rank(&self) -> u32;

    /// Dimension of the algebra.
    fn dimension(&self) -> u32;

    /// Number of roots |Φ|.
    fn num_roots(&self) -> u32;

    /// Coxeter number h.
    fn coxeter_number(&self) -> u32;

    /// Dual Coxeter number h∨.
    fn dual_coxeter_number(&self) -> u32;

    /// Order of the Weyl group |W|.
    fn weyl_order(&self) -> u64;

    /// Exponents m₁, m₂, ..., mᵣ.
    fn exponents(&self) -> &[u32];

    /// Casimir degrees d₁, d₂, ..., dᵣ (= exponents + 1).
    fn casimir_degrees(&self) -> &[u32];

    /// Quadratic Casimir C₂(fundamental) as (numerator, denominator).
    fn c2_fundamental(&self) -> (u32, u32);

    /// Whether the algebra has a Casimir of the given degree.
    fn has_casimir_degree(&self, d: u32) -> bool {
        self.casimir_degrees().contains(&d)
    }

    /// Dimension formula check: dim = rank + |Φ|.
    fn verify_dimension(&self) -> bool {
        self.dimension() == self.rank() + self.num_roots()
    }
}

impl LieAlgebra for LieGroup {
    fn name(&self) -> &str { self.name }
    fn rank(&self) -> u32 { self.rank }
    fn dimension(&self) -> u32 { self.dimension }
    fn num_roots(&self) -> u32 { self.num_roots }
    fn coxeter_number(&self) -> u32 { self.coxeter_number }
    fn dual_coxeter_number(&self) -> u32 { self.dual_coxeter }
    fn weyl_order(&self) -> u64 { self.weyl_order }
    fn exponents(&self) -> &[u32] { self.exponents }
    fn casimir_degrees(&self) -> &[u32] { self.casimir_degrees }
    fn c2_fundamental(&self) -> (u32, u32) { self.c2_fundamental }
}

#[cfg(test)]
mod tests {
    use super::*;
    use groups::{E8, G2, SU3};

    #[test]
    fn test_lie_algebra_trait() {
        let e8: &dyn LieAlgebra = &E8;
        assert_eq!(e8.name(), "E8");
        assert_eq!(e8.rank(), 8);
        assert_eq!(e8.dimension(), 248);
        assert!(e8.verify_dimension());
        assert!(!e8.has_casimir_degree(4));
        assert!(e8.has_casimir_degree(30));
    }

    #[test]
    fn test_g2_via_trait() {
        let g2: &dyn LieAlgebra = &G2;
        assert_eq!(g2.coxeter_number(), 6);
        assert_eq!(g2.weyl_order(), 12);
        assert_eq!(g2.exponents(), &[1, 5]);
    }

    #[test]
    fn test_su3_via_trait() {
        let su3: &dyn LieAlgebra = &SU3;
        assert_eq!(su3.c2_fundamental(), (4, 3));
        assert!(su3.verify_dimension());
    }
}
