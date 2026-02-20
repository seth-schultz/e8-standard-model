//! Physical constants derived from E8 lattice geometry.

use rug::Float;

use crate::algebra::groups::identities::{BINOM_IM_O_3, DIM_SO8};
use crate::algebra::groups::E8;
use crate::precision::{exp_neg_gamma, mpf_u64, pi, planck_mass_mev, powi, precision_bits};

/// R = |Φ(E8)| × e^{-γ} — the Mertens-regularized coordination number.
/// 240 = |Φ(E8)| (number of roots), e^{-γ} from Epstein zeta regularization.
pub fn mertens_r() -> Float {
    Float::with_val(precision_bits(), E8.num_roots) * exp_neg_gamma()
}

/// δ = C(dim(Im(O)),3) / (4π⁴) — the 8D nearest-neighbor gravity correction.
/// C(7,3) = 35 = 7(associative) + 28(non-associative) 3-planes in Im(O).
pub fn delta() -> Float {
    let pi4 = powi(&pi(), 4);
    Float::with_val(precision_bits(), BINOM_IM_O_3) / (Float::with_val(precision_bits(), 4) * pi4)
}

/// The normalization factor: dim(so(8)) = C(8,2) = 28.
pub fn norm_factor() -> Float {
    mpf_u64(DIM_SO8 as u64)
}

/// Planck mass in MeV (input constant).
pub fn m_planck_mev() -> Float {
    planck_mass_mev()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mertens_r() {
        crate::precision::set_precision(50);
        let r = mertens_r();
        // R ≈ 240 × 0.5615 ≈ 134.76
        assert!(r > 134.0 && r < 136.0);
    }

    #[test]
    fn test_delta() {
        crate::precision::set_precision(50);
        let d = delta();
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
}
