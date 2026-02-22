//! CKM and PMNS mixing matrices from E8 structure.

pub mod ckm;
pub mod pmns;
pub mod cp_phase;

use crate::mass::sectors::AllMasses;
use crate::precision::scalar::Scalar;
use ckm::CKMResult;

// ═══════════════════════════════════════════════════════════════
// MixingMatrix trait — extensible mixing matrix abstraction
// ═══════════════════════════════════════════════════════════════

/// A CKM-like quark mixing matrix computed from masses.
pub trait CKMMixing {
    /// Build the CKM matrix from quark masses.
    fn ckm<S: Scalar>(&self, masses: &AllMasses<S>) -> CKMResult<S>;
}

/// PMNS-like lepton mixing angles.
pub trait PMNSMixing {
    /// sin²θ₁₂ (solar angle).
    fn sin2_theta12<S: Scalar>(&self) -> S;
    /// sin²θ₂₃ (atmospheric angle).
    fn sin2_theta23<S: Scalar>(&self) -> S;
    /// sin²θ₁₃ (reactor angle).
    fn sin2_theta13<S: Scalar>(&self) -> S;
    /// CP phase δ_PMNS in radians.
    fn delta_pmns_rad<S: Scalar>(&self) -> S;
}

/// Fritzsch texture CKM with octonionic CP phase (E8 default).
pub struct FritzschCKM;

impl CKMMixing for FritzschCKM {
    fn ckm<S: Scalar>(&self, masses: &AllMasses<S>) -> CKMResult<S> {
        ckm::build_ckm(masses)
    }
}

/// G₂ Coxeter geometry PMNS mixing (E8 default).
pub struct G2PMNS;

impl PMNSMixing for G2PMNS {
    fn sin2_theta12<S: Scalar>(&self) -> S {
        pmns::sin2_theta12()
    }

    fn sin2_theta23<S: Scalar>(&self) -> S {
        pmns::sin2_theta23()
    }

    fn sin2_theta13<S: Scalar>(&self) -> S {
        pmns::sin2_theta13()
    }

    fn delta_pmns_rad<S: Scalar>(&self) -> S {
        cp_phase::delta_pmns_rad()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_ckm_trait() {
        set_precision(50);
        let masses = crate::mass::sectors::compute_all_masses::<DefaultScalar>();
        let ckm_impl = FritzschCKM;
        let result = ckm_impl.ckm(&masses);
        let v_ud = result.magnitudes[0].to_f64();
        assert!((v_ud - 0.974).abs() < 0.01, "|V_ud| = {}", v_ud);
    }

    #[test]
    fn test_pmns_trait() {
        set_precision(50);
        let pmns = G2PMNS;
        let s12: DefaultScalar = pmns.sin2_theta12();
        assert!((s12.to_f64() - 0.311).abs() < 0.01);
    }
}
