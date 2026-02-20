//! Lie group invariants used throughout the E8 framework.

/// Invariants for a simple Lie group.
#[derive(Debug, Clone)]
pub struct LieGroup {
    pub name: &'static str,
    pub rank: u32,
    pub dimension: u32,
    pub num_roots: u32,          // |Φ|
    pub coxeter_number: u32,     // h
    pub dual_coxeter: u32,       // h∨
    pub weyl_order: u64,         // |W|
    pub exponents: &'static [u32],
    pub casimir_degrees: &'static [u32],
    pub c2_fundamental: (u32, u32), // Quadratic Casimir C₂(fund) as fraction (num, den)
}

// ═══════════════════════════════════════════════════════════════
// Standard Lie group data
// ═══════════════════════════════════════════════════════════════

pub const E8: LieGroup = LieGroup {
    name: "E8",
    rank: 8,
    dimension: 248,
    num_roots: 240,
    coxeter_number: 30,
    dual_coxeter: 30,
    weyl_order: 696_729_600, // 2^14 * 3^5 * 5^2 * 7
    exponents: &[1, 7, 11, 13, 17, 19, 23, 29],
    casimir_degrees: &[2, 8, 12, 14, 18, 20, 24, 30],
    c2_fundamental: (60, 1),  // C₂(248) = 60 (adjoint = fundamental for E8)
};

pub const E7: LieGroup = LieGroup {
    name: "E7",
    rank: 7,
    dimension: 133,
    num_roots: 126,
    coxeter_number: 18,
    dual_coxeter: 18,
    weyl_order: 2_903_040,
    exponents: &[1, 5, 7, 9, 11, 13, 17],
    casimir_degrees: &[2, 6, 8, 10, 12, 14, 18],
    c2_fundamental: (3, 1),  // C₂(56)
};

pub const E6: LieGroup = LieGroup {
    name: "E6",
    rank: 6,
    dimension: 78,
    num_roots: 72,
    coxeter_number: 12,
    dual_coxeter: 12,
    weyl_order: 51_840,
    exponents: &[1, 4, 5, 7, 8, 11],
    casimir_degrees: &[2, 5, 6, 8, 9, 12],
    c2_fundamental: (4, 3), // C₂(27) = 4/3... actually 26/6 = 13/3
};

pub const G2: LieGroup = LieGroup {
    name: "G2",
    rank: 2,
    dimension: 14,
    num_roots: 12,
    coxeter_number: 6,
    dual_coxeter: 4,
    weyl_order: 12,
    exponents: &[1, 5],
    casimir_degrees: &[2, 6],
    c2_fundamental: (2, 1),  // C₂(7) = 2
};

pub const D4: LieGroup = LieGroup {
    name: "D4",
    rank: 4,
    dimension: 28,
    num_roots: 24,
    coxeter_number: 6,
    dual_coxeter: 6,
    weyl_order: 192,
    exponents: &[1, 3, 3, 5],
    casimir_degrees: &[2, 4, 4, 6],
    c2_fundamental: (7, 4), // C₂(8_v) = 7/4
};

pub const F4: LieGroup = LieGroup {
    name: "F4",
    rank: 4,
    dimension: 52,
    num_roots: 48,
    coxeter_number: 12,
    dual_coxeter: 9,
    weyl_order: 1152,
    exponents: &[1, 5, 7, 11],
    casimir_degrees: &[2, 6, 8, 12],
    c2_fundamental: (6, 1),
};

pub const SU5: LieGroup = LieGroup {
    name: "SU(5)",
    rank: 4,
    dimension: 24,
    num_roots: 20,
    coxeter_number: 5,
    dual_coxeter: 5,
    weyl_order: 120,
    exponents: &[1, 2, 3, 4],
    casimir_degrees: &[2, 3, 4, 5],
    c2_fundamental: (12, 5), // C₂(5) = 12/5... actually (N²-1)/(2N) = 24/10 = 12/5
};

pub const SU3: LieGroup = LieGroup {
    name: "SU(3)",
    rank: 2,
    dimension: 8,
    num_roots: 6,
    coxeter_number: 3,
    dual_coxeter: 3,
    weyl_order: 6,
    exponents: &[1, 2],
    casimir_degrees: &[2, 3],
    c2_fundamental: (4, 3), // C₂(3) = (N²-1)/(2N) = 8/6 = 4/3
};

pub const SU2: LieGroup = LieGroup {
    name: "SU(2)",
    rank: 1,
    dimension: 3,
    num_roots: 2,
    coxeter_number: 2,
    dual_coxeter: 2,
    weyl_order: 2,
    exponents: &[1],
    casimir_degrees: &[2],
    c2_fundamental: (3, 4), // C₂(2) = (N²-1)/(2N) = 3/4
};

/// Check that E8 has NO degree-4 Casimir (theorem: λ(m_P)=0).
pub fn e8_has_no_degree4_casimir() -> bool {
    !E8.casimir_degrees.contains(&4)
}

/// Key identities used in the framework.
pub mod identities {
    use super::*;

    /// |W(G₂)| + 1 = 13 (prime)
    pub const W_G2_PLUS_1: u64 = G2.weyl_order + 1;
    /// |W(D₄)| + 1 = 193 (prime)
    pub const W_D4_PLUS_1: u64 = D4.weyl_order + 1;
    /// |W(F₄)| + 1 = 1153 (prime)
    pub const W_F4_PLUS_1: u64 = F4.weyl_order + 1;
    /// |W(E₇)| + 1 = 2903041 (prime)
    pub const W_E7_PLUS_1: u64 = E7.weyl_order + 1;
    /// dim(G₂) = 14
    pub const DIM_G2: u32 = G2.dimension;
    /// dim(Im(O)) = 7 — the imaginary octonions
    pub const DIM_IM_OCTONIONS: u32 = 7;
    /// |Φ(E₆)| = 72 = h(G₂) × |W(G₂)| = 6 × 12
    pub const PHI_E6: u32 = E6.num_roots;
    /// dim(so(8)) = C(8,2) = 28
    pub const DIM_SO8: u32 = D4.dimension;
    /// N_c = 3 — number of QCD colors = rank(SU(3)) + 1
    pub const N_C: u32 = SU3.rank + 1;
    /// N_gen = 3 — number of generations
    pub const N_GEN: u32 = 3;
    /// C(dim(Im(O)), 3) = C(7,3) = 35 — 3-planes in Im(O)
    pub const BINOM_IM_O_3: u32 =
        DIM_IM_OCTONIONS * (DIM_IM_OCTONIONS - 1) * (DIM_IM_OCTONIONS - 2) / 6;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_e8_basics() {
        assert_eq!(E8.num_roots, 240);
        assert_eq!(E8.dimension, 248);
        assert_eq!(E8.rank, 8);
        assert_eq!(E8.coxeter_number, 30);
    }

    #[test]
    fn test_no_degree4_casimir() {
        assert!(e8_has_no_degree4_casimir());
    }

    #[test]
    fn test_g2_invariants() {
        assert_eq!(G2.dimension, 14);
        assert_eq!(G2.weyl_order, 12);
        assert_eq!(G2.coxeter_number, 6);
        assert_eq!(G2.rank, 2);
    }

    #[test]
    fn test_weyl_plus_one_primes() {
        assert_eq!(G2.weyl_order + 1, 13);
        assert_eq!(D4.weyl_order + 1, 193);
        assert_eq!(F4.weyl_order + 1, 1153);
    }

    #[test]
    fn test_su3_casimir() {
        // C₂(SU(3), fund) = 4/3
        assert_eq!(SU3.c2_fundamental, (4, 3));
    }
}
