//! Higgs boson mass: m_H = v × √(2λ) with v = √2 × m_t (y_t = 1).

use rug::Float;

use crate::precision::{precision_bits, sqrt};

use super::quartic::higgs_quartic;

/// Higgs VEV: v = √2 × m_t (consequence of y_t = 1 from E8).
/// m_t in GeV.
pub fn higgs_vev_gev(m_t_gev: &Float) -> Float {
    let prec = precision_bits();
    sqrt(&Float::with_val(prec, 2)) * m_t_gev
}

/// Higgs mass: m_H = v × √(2λ) in GeV.
pub fn higgs_mass_gev(m_t_gev: &Float) -> Float {
    let prec = precision_bits();
    let v = higgs_vev_gev(m_t_gev);
    let lambda = higgs_quartic();
    let two_lambda = Float::with_val(prec, 2) * lambda;
    v * sqrt(&two_lambda)
}

/// Compute m_H using the E8-predicted top mass from Koide.
pub fn higgs_mass_default() -> Float {
    let prec = precision_bits();
    let masses = crate::mass::sectors::compute_all_masses::<Float>();
    // Convert predicted m_t from MeV to GeV
    let m_t = Float::with_val(prec, &masses.top / Float::with_val(prec, 1000));
    higgs_mass_gev(&m_t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn test_higgs_mass() {
        set_precision(50);
        let mh = higgs_mass_default();
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
        let prec = precision_bits();
        let m_t = Float::with_val(prec, 172.76f64);
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
