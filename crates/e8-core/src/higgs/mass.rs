//! Higgs boson mass: m_H = v × √(2λ) with v = √2 × m_t (y_t = 1).

use crate::override_context::OverrideContext;
use crate::precision::scalar::Scalar;

use super::quartic::higgs_quartic_with_ctx;

/// Higgs VEV: v = √2 × m_t (consequence of y_t = 1 from E8).
/// m_t in GeV.
pub fn higgs_vev_gev<S: Scalar>(m_t_gev: &S) -> S {
    S::from_u64(2).sqrt() * m_t_gev.clone()
}

/// Higgs mass with overridable quartic coupling.
pub fn higgs_mass_gev_with_ctx<S: Scalar>(m_t_gev: &S, ctx: &OverrideContext) -> S {
    let v = higgs_vev_gev(m_t_gev);
    let lambda: S = higgs_quartic_with_ctx(ctx);
    let two_lambda = S::from_u64(2) * lambda;
    v * two_lambda.sqrt()
}

/// Higgs mass: m_H = v × √(2λ) in GeV.
pub fn higgs_mass_gev<S: Scalar>(m_t_gev: &S) -> S {
    higgs_mass_gev_with_ctx(m_t_gev, &OverrideContext::new())
}

/// Compute m_H using E8-predicted top mass, with overrides.
pub fn higgs_mass_with_ctx<S: Scalar>(ctx: &OverrideContext) -> S {
    let masses = crate::mass::sectors::compute_all_masses_with_ctx::<S>(ctx);
    let m_t = masses.top / S::from_u64(1000);
    higgs_mass_gev_with_ctx(&m_t, ctx)
}

/// Compute m_H using the E8-predicted top mass from Koide.
pub fn higgs_mass_default<S: Scalar>() -> S {
    higgs_mass_with_ctx(&OverrideContext::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::scalar::Scalar;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_higgs_mass() {
        set_precision(50);
        let mh: DefaultScalar = higgs_mass_default();
        let val = mh.to_f64();
        assert!((val - 125.3).abs() < 1.0, "m_H = {} GeV", val);
    }

    #[test]
    fn test_higgs_vev() {
        set_precision(50);
        let m_t = DefaultScalar::from_f64(172.76);
        let v = higgs_vev_gev(&m_t);
        let val = v.to_f64();
        assert!((val - 244.3).abs() < 1.0, "v = {} GeV", val);
    }
}
