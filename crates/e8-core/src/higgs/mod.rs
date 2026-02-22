//! Higgs sector: quartic coupling, mass, and θ_QCD = 0.

pub mod quartic;
pub mod mass;
pub mod theta_qcd;

use crate::precision::scalar::Scalar;

// ═══════════════════════════════════════════════════════════════
// HiggsSector trait — extensible Higgs sector abstraction
// ═══════════════════════════════════════════════════════════════

/// Higgs sector predictions: quartic coupling, masses, θ_QCD.
pub trait HiggsSector {
    /// Higgs quartic coupling λ.
    fn quartic<S: Scalar>(&self) -> S;

    /// Higgs mass m_H in GeV (computed from top mass).
    fn higgs_mass_gev<S: Scalar>(&self) -> S;

    /// Second scalar mass m_S in GeV.
    fn second_scalar_mass_gev<S: Scalar>(&self) -> S;

    /// θ̄_QCD — should be exactly 0 in E8.
    fn theta_qcd(&self) -> f64;
}

/// E8 Higgs sector: λ = 7π⁴/72², m_H from y_t = 1, θ_QCD = 0.
pub struct E8HiggsSector;

impl HiggsSector for E8HiggsSector {
    fn quartic<S: Scalar>(&self) -> S {
        quartic::higgs_quartic()
    }

    fn higgs_mass_gev<S: Scalar>(&self) -> S {
        mass::higgs_mass_default()
    }

    fn second_scalar_mass_gev<S: Scalar>(&self) -> S {
        let m_h: S = mass::higgs_mass_default();
        quartic::second_scalar_mass(&m_h)
    }

    fn theta_qcd(&self) -> f64 {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_higgs_sector_trait() {
        set_precision(50);
        let hs = E8HiggsSector;
        let lambda: rug::Float = hs.quartic();
        assert!((lambda.to_f64() - 0.1315).abs() < 0.001);
        let m_h: rug::Float = hs.higgs_mass_gev();
        assert!((m_h.to_f64() - 125.3).abs() < 1.0);
        assert_eq!(hs.theta_qcd(), 0.0);
    }
}
