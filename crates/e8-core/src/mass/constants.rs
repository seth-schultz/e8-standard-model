//! Physical constants derived from E8 lattice geometry.

use crate::algebra::groups::identities::{BINOM_IM_O_3, DIM_SO8};
use crate::algebra::groups::E8;
use crate::precision::scalar::Scalar;

// ═══════════════════════════════════════════════════════════════
// QCD correction constants from E8 subgroup geometry
// ═══════════════════════════════════════════════════════════════

/// C_F = 4/3 — color Casimir for fundamental representation of SU(3).
pub const CF_FUNDAMENTAL: f64 = 4.0 / 3.0;

/// |Φ(E₆)| = 72 — number of roots of E₆.
pub const PHI_E6_ROOTS: u32 = 72;

/// 3 × |Φ(E₆)| = 216 — from E8 ⊃ E6 × SU(3) decomposition.
/// n_colors × |Φ(E₆)| = 3 × 72.
pub const E6_SU3_FACTOR: f64 = 216.0;

/// R = |Φ(E8)| × e^{-γ} — the Mertens-regularized coordination number.
/// 240 = |Φ(E8)| (number of roots), e^{-γ} from Epstein zeta regularization.
pub fn mertens_r<S: Scalar>() -> S {
    S::from_u64(E8.num_roots as u64) * S::exp_neg_gamma()
}

/// δ = C(dim(Im(O)),3) / (4π⁴) — the 8D nearest-neighbor gravity correction.
/// C(7,3) = 35 = 7(associative) + 28(non-associative) 3-planes in Im(O).
pub fn delta<S: Scalar>() -> S {
    let pi4 = S::pi().powi(4);
    S::from_u64(BINOM_IM_O_3 as u64) / (S::from_u64(4) * pi4)
}

/// The normalization factor: dim(so(8)) = C(8,2) = 28.
pub fn norm_factor<S: Scalar>() -> S {
    S::from_u64(DIM_SO8 as u64)
}

/// Planck mass in MeV (input constant).
pub fn m_planck_mev<S: Scalar>() -> S {
    S::planck_mass_mev()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_mertens_r() {
        crate::precision::set_precision(50);
        let r: DefaultScalar = mertens_r();
        // R ≈ 240 × 0.5615 ≈ 134.76
        assert!(r > 134.0 && r < 136.0);
    }

    #[test]
    fn test_delta() {
        crate::precision::set_precision(50);
        let d: DefaultScalar = delta();
        // δ ≈ 35/(4 × 97.409) ≈ 0.0899
        assert!(d > 0.089 && d < 0.091);
    }

    #[test]
    fn test_constants_match() {
        // Verify the named constants produce expected values
        assert_eq!(E8.num_roots, 240);
        assert_eq!(DIM_SO8, 28);
        assert_eq!(BINOM_IM_O_3, 35);
    }

    #[test]
    fn test_qcd_constants() {
        assert!((CF_FUNDAMENTAL - 4.0 / 3.0).abs() < 1e-15);
        assert_eq!(PHI_E6_ROOTS, 72);
        assert!((E6_SU3_FACTOR - 216.0).abs() < 1e-15);
        // Verify decomposition: 3 × 72 = 216
        assert_eq!(3 * PHI_E6_ROOTS as u64, E6_SU3_FACTOR as u64);
    }
}
