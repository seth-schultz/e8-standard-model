//! Mass formula, Koide parametrization, and sector configurations.

pub mod constants;
pub mod formula;
pub mod koide;
pub mod sectors;

use crate::precision::scalar::Scalar;
use sectors::AllMasses;
use formula::SectorSums;

// ═══════════════════════════════════════════════════════════════
// MassFormula trait — sector mass sums from the lattice
// ═══════════════════════════════════════════════════════════════

/// A formula that computes the four sector mass sums (Σ_lep, Σ_up, Σ_down, Σ_ν)
/// from a lattice theory.
///
/// The E8 implementation uses:
///   Σ = f × m_P × exp(-(A × R + δ)/dim(so(8)))
pub trait MassFormula {
    /// Compute all four sector mass sums.
    fn sector_sums<S: Scalar>(&self) -> SectorSums<S>;
}

/// The E8 mass formula: Σ = f × m_P × exp(-(A × R + δ)/28).
pub struct E8MassFormula;

impl MassFormula for E8MassFormula {
    fn sector_sums<S: Scalar>(&self) -> SectorSums<S> {
        formula::compute_all_sector_sums()
    }
}

// ═══════════════════════════════════════════════════════════════
// MassSplitting trait — sector sum → individual masses
// ═══════════════════════════════════════════════════════════════

/// A splitting rule that produces three masses from a sector sum.
///
/// The E8 implementation uses the Koide parametrization:
///   √m_k = M × (1 + r × cos(2πk/3 + φ))
pub trait MassSplitting {
    /// Split a sector sum into three masses [m₁, m₂, m₃] (heaviest first).
    fn split<S: Scalar>(&self, sigma: &S, r_fourth: &S, phi: &S) -> [S; 3];
}

/// Koide parametrization: √m_k = M × (1 + r × cos(2πk/3 + φ)).
pub struct KoideSplitting;

impl MassSplitting for KoideSplitting {
    fn split<S: Scalar>(&self, sigma: &S, r_fourth: &S, phi: &S) -> [S; 3] {
        koide::koide_masses(&koide::KoideParams {
            r_fourth: r_fourth.clone(),
            phi: phi.clone(),
            sigma: sigma.clone(),
        })
    }
}

/// Compute all 12 fermion masses using a MassFormula + MassSplitting.
///
/// This is the trait-based equivalent of `sectors::compute_all_masses`.
pub fn compute_masses_from<S: Scalar, F: MassFormula, K: MassSplitting>(
    formula: &F,
    splitting: &K,
) -> AllMasses<S> {
    let sums = formula.sector_sums::<S>();

    // Koide parameters from E8 group theory (must match sectors.rs exactly)
    let h = crate::algebra::groups::G2.coxeter_number; // 6
    let w = crate::algebra::groups::G2.weyl_order;     // 12
    let n_su5 = crate::algebra::groups::SU5.rank + 1;  // 5

    // Leptons: r⁴ = 4 = (√2)⁴, φ = 2/9
    let r4_lep = S::from_u64(4);
    let phi_lep = S::from_u64(2) / S::from_u64(9);

    // Up quarks: r⁴ = C(5,2) = 10, φ = 5⁴/6⁵
    let dim_antisym = (n_su5 * (n_su5 - 1) / 2) as u64; // 10
    let r4_up = S::from_u64(dim_antisym);
    let phi_up = S::from_u64(((h - 1) as u64).pow(4))
        / S::from_u64((h as u64).pow(5));

    // Down quarks: r⁴ = 10 - √2, φ = 1/6
    let r4_down = S::from_u64(dim_antisym) - S::from_u64(2).sqrt();
    let phi_down = S::one() / S::from_u64(h as u64);

    // Neutrinos: r⁴ = 4, φ = 2/9 + π/12, sigma in meV (×10⁹)
    let r4_nu = S::from_u64(4);
    let phi_nu = S::from_u64(2) / S::from_u64(9)
        + S::pi() / S::from_u64(w);
    let sigma_nu_mev = sums.neutrino.clone() * S::from_f64(1e9);

    let lep = splitting.split::<S>(&sums.leptons, &r4_lep, &phi_lep);
    let up = splitting.split::<S>(&sums.up, &r4_up, &phi_up);
    let down = splitting.split::<S>(&sums.down, &r4_down, &phi_down);
    let nu = splitting.split::<S>(&sigma_nu_mev, &r4_nu, &phi_nu);

    // k=0 → heaviest, k=1 → lightest, k=2 → middle (Koide convention)
    AllMasses {
        tau: lep[0].clone(),
        electron: lep[1].clone(),
        muon: lep[2].clone(),

        top: up[0].clone(),
        up: up[1].clone(),
        charm: up[2].clone(),

        bottom: down[0].clone(),
        down: down[1].clone(),
        strange: down[2].clone(),

        nu3: nu[0].clone(),
        nu1: nu[1].clone(),
        nu2: nu[2].clone(),

        sigma_lep: sums.leptons,
        sigma_up: sums.up,
        sigma_down: sums.down,
        sigma_nu: sums.neutrino,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_mass_formula_trait() {
        set_precision(50);
        let formula = E8MassFormula;
        let sums: SectorSums<DefaultScalar> = formula.sector_sums();
        let sigma_lep = sums.leptons.to_f64();
        assert!((sigma_lep - 1883.0).abs() < 1.0, "Σ_lep = {}", sigma_lep);
    }

    #[test]
    fn test_trait_based_masses() {
        set_precision(50);
        let masses: AllMasses<DefaultScalar> =
            compute_masses_from(&E8MassFormula, &KoideSplitting);
        let m_e = masses.electron.to_f64();
        assert!((m_e - 0.511).abs() < 0.01, "m_e = {}", m_e);
        let m_t = masses.top.to_f64();
        assert!((m_t - 172760.0).abs() < 500.0, "m_t = {}", m_t);
    }
}
