//! Theory composition: combine physics traits into a complete theory.
//!
//! The `Theory` struct composes a root system, mass formula, mass splitting,
//! mixing matrices, gauge couplings, and Higgs sector into a single object
//! that can compute the full scorecard.

use crate::coupling::{E8GaugeCouplings, GaugeCouplings};
use crate::higgs::{E8HiggsSector, HiggsSector};
use crate::lattice::{E8RootSystem, RootSystem};
use crate::mass::{E8MassFormula, KoideSplitting, MassFormula, MassSplitting};
use crate::mass::sectors::AllMasses;
use crate::mixing::{CKMMixing, FritzschCKM, G2PMNS, PMNSMixing};
use crate::mixing::ckm::CKMResult;
use crate::precision::scalar::Scalar;

/// A complete physics theory composed from modular traits.
///
/// Each component can be swapped independently:
/// - `R`: Root system (e.g., E8RootSystem)
/// - `MF`: Mass formula (e.g., E8MassFormula)
/// - `MS`: Mass splitting (e.g., KoideSplitting)
/// - `CKM`: Quark mixing (e.g., FritzschCKM)
/// - `PMNS`: Lepton mixing (e.g., G2PMNS)
/// - `GC`: Gauge couplings (e.g., E8GaugeCouplings)
/// - `HS`: Higgs sector (e.g., E8HiggsSector)
pub struct Theory<R, MF, MS, CKM, PMNS, GC, HS> {
    pub root_system: R,
    pub mass_formula: MF,
    pub mass_splitting: MS,
    pub ckm_mixing: CKM,
    pub pmns_mixing: PMNS,
    pub gauge_couplings: GC,
    pub higgs_sector: HS,
}

impl<R, MF, MS, CKM, PMNS, GC, HS> Theory<R, MF, MS, CKM, PMNS, GC, HS>
where
    R: RootSystem,
    MF: MassFormula,
    MS: MassSplitting,
    CKM: CKMMixing,
    PMNS: PMNSMixing,
    GC: GaugeCouplings,
    HS: HiggsSector,
{
    /// Compute all fermion masses.
    pub fn masses<S: Scalar>(&self) -> AllMasses<S> {
        crate::mass::compute_masses_from(&self.mass_formula, &self.mass_splitting)
    }

    /// Compute the CKM matrix from computed masses.
    pub fn ckm<S: Scalar>(&self, masses: &AllMasses<S>) -> CKMResult<S> {
        self.ckm_mixing.ckm(masses)
    }

    /// Number of roots in the root system.
    pub fn root_count(&self) -> usize {
        self.root_system.root_count()
    }

    /// Verify the root system norms.
    pub fn verify_roots(&self) -> bool {
        self.root_system.verify_norms()
    }
}

/// The E8 Standard Model with all default components.
pub type E8StandardModel = Theory<
    E8RootSystem,
    E8MassFormula,
    KoideSplitting,
    FritzschCKM,
    G2PMNS,
    E8GaugeCouplings,
    E8HiggsSector,
>;

/// Construct the E8 Standard Model (zero free parameters).
pub fn e8_standard_model() -> E8StandardModel {
    Theory {
        root_system: E8RootSystem,
        mass_formula: E8MassFormula,
        mass_splitting: KoideSplitting,
        ckm_mixing: FritzschCKM,
        pmns_mixing: G2PMNS,
        gauge_couplings: E8GaugeCouplings,
        higgs_sector: E8HiggsSector,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_e8_standard_model() {
        set_precision(50);
        let theory = e8_standard_model();

        // Root system
        assert_eq!(theory.root_count(), 240);
        assert!(theory.verify_roots());

        // Masses
        let masses: AllMasses<DefaultScalar> = theory.masses();
        let m_e = masses.electron.to_f64();
        assert!((m_e - 0.511).abs() < 0.01, "m_e = {}", m_e);

        // Gauge couplings
        let ai: DefaultScalar = theory.gauge_couplings.alpha_inverse();
        assert!((ai.to_f64() - 137.036).abs() < 0.001);

        // Higgs
        let m_h: DefaultScalar = theory.higgs_sector.higgs_mass_gev();
        assert!((m_h.to_f64() - 125.3).abs() < 1.0);

        // CKM
        let ckm = theory.ckm(&masses);
        let v_ud = ckm.magnitudes[0].to_f64();
        assert!((v_ud - 0.974).abs() < 0.01, "|V_ud| = {}", v_ud);

        // PMNS
        let s12: DefaultScalar = theory.pmns_mixing.sin2_theta12();
        assert!((s12.to_f64() - 0.311).abs() < 0.01);
    }
}
