//! Gauge coupling constants derived from E8.

pub mod alpha;
pub mod weinberg;
pub mod alpha_s;

use crate::precision::scalar::Scalar;

// ═══════════════════════════════════════════════════════════════
// GaugeCouplings trait — extensible gauge coupling abstraction
// ═══════════════════════════════════════════════════════════════

/// Gauge coupling constants at various scales.
pub trait GaugeCouplings {
    /// 1/α(0) — fine structure constant inverse.
    fn alpha_inverse<S: Scalar>(&self) -> S;

    /// sin²θ_W at GUT scale.
    fn sin2_theta_w_gut<S: Scalar>(&self) -> S;

    /// sin²θ_W at M_Z scale.
    fn sin2_theta_w_mz<S: Scalar>(&self) -> S;

    /// α_s(M_Z) — strong coupling at Z mass.
    fn alpha_s_mz<S: Scalar>(&self) -> S;

    /// M_GUT in GeV.
    fn m_gut_gev<S: Scalar>(&self) -> S;
}

/// E8 gauge couplings: continued fraction α, GUT normalization.
pub struct E8GaugeCouplings;

impl GaugeCouplings for E8GaugeCouplings {
    fn alpha_inverse<S: Scalar>(&self) -> S {
        alpha::alpha_inverse()
    }

    fn sin2_theta_w_gut<S: Scalar>(&self) -> S {
        weinberg::sin2_theta_w_gut()
    }

    fn sin2_theta_w_mz<S: Scalar>(&self) -> S {
        weinberg::sin2_theta_w_mz()
    }

    fn alpha_s_mz<S: Scalar>(&self) -> S {
        alpha_s::alpha_s_mz()
    }

    fn m_gut_gev<S: Scalar>(&self) -> S {
        alpha_s::m_gut_gev()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_gauge_couplings_trait() {
        set_precision(50);
        let gc = E8GaugeCouplings;
        let ai: DefaultScalar = gc.alpha_inverse();
        assert!((ai.to_f64() - 137.036).abs() < 0.001);
        let s2w: DefaultScalar = gc.sin2_theta_w_gut();
        assert!((s2w.to_f64() - 0.375).abs() < 1e-10);
    }
}
