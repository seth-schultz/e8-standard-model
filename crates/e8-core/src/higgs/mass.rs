//! Higgs boson mass: m_H = v × √(2λ) with v = √2 × m_t (y_t = 1).

use crate::precision::scalar::Scalar;

use super::quartic::higgs_quartic;

/// Higgs VEV: v = √2 × m_t (consequence of y_t = 1 from E8).
/// m_t in GeV.
pub fn higgs_vev_gev<S: Scalar>(m_t_gev: &S) -> S {
    S::from_u64(2).sqrt() * m_t_gev.clone()
}

/// Higgs mass: m_H = v × √(2λ) in GeV.
pub fn higgs_mass_gev<S: Scalar>(m_t_gev: &S) -> S {
    let v = higgs_vev_gev(m_t_gev);
    let lambda: S = higgs_quartic();
    let two_lambda = S::from_u64(2) * lambda;
    v * two_lambda.sqrt()
}

/// Compute m_H using the E8-predicted top mass from Koide.
pub fn higgs_mass_default<S: Scalar>() -> S {
    let masses = crate::mass::sectors::compute_all_masses::<S>();
    let m_t = masses.top / S::from_u64(1000);
    higgs_mass_gev(&m_t)
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
        // Experimental: 125.25 ± 0.17 GeV
        assert!(
            (val - 125.3).abs() < 1.0,
            "m_H = {} GeV",
            val
        );
    }

    #[test]
    fn test_higgs_vev() {
        set_precision(50);
        let m_t = DefaultScalar::from_f64(172.76);
        let v = higgs_vev_gev(&m_t);
        let val = v.to_f64();
        // v ≈ 244 GeV (≈ √2 × 172.76)
        assert!(
            (val - 244.3).abs() < 1.0,
            "v = {} GeV",
            val
        );
    }
}
