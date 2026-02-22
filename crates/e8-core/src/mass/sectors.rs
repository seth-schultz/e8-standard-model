//! Sector configurations for all fermion families.
//!
//! Each sector's Koide parameters are derived from E8 representation theory:
//! - r⁴ from SU(5) Yukawa structure
//! - φ from G₂ Coxeter geometry
//! - Σ from the mass formula with group-theoretic A and f

use crate::algebra::groups::{G2, SU5};
use crate::precision::scalar::Scalar;

use super::formula::compute_all_sector_sums;
use super::koide::{koide_masses, KoideParams};

/// All 9 charged fermion masses + 3 neutrino masses.
#[derive(Debug, Clone)]
pub struct AllMasses<S: Scalar> {
    // Charged leptons (MeV)
    pub electron: S,
    pub muon: S,
    pub tau: S,

    // Up-type quarks (MeV)
    pub up: S,
    pub charm: S,
    pub top: S,

    // Down-type quarks (MeV)
    pub down: S,
    pub strange: S,
    pub bottom: S,

    // Neutrinos (meV)
    pub nu1: S,
    pub nu2: S,
    pub nu3: S,

    // Sector sums (MeV, except neutrinos in MeV)
    pub sigma_lep: S,
    pub sigma_up: S,
    pub sigma_down: S,
    pub sigma_nu: S,
}

/// Compute all masses from zero free parameters.
pub fn compute_all_masses<S: Scalar>() -> AllMasses<S> {
    let sums = compute_all_sector_sums::<S>();

    let h = G2.coxeter_number; // h(G₂) = 6
    let w = G2.weyl_order;     // |W(G₂)| = 12

    // ═══════════════════════════════════════════════════════════
    // Charged leptons: r⁴ = 4 = (√2)⁴, φ = 2/9 (Z₃ variational)
    // ═══════════════════════════════════════════════════════════
    let lep_params = KoideParams {
        r_fourth: S::from_u64(4),
        phi: S::from_u64(2) / S::from_u64(9),
        sigma: sums.leptons.clone(),
    };
    let lep = koide_masses(&lep_params);
    // k=0 → tau, k=1 → electron, k=2 → muon

    // ═══════════════════════════════════════════════════════════
    // Up quarks: r⁴ = dim(∧²(5)) = C(5,2) = 10
    //            φ = (h-1)⁴/h⁵ = 5⁴/6⁵ = 625/7776
    //            (from antisymmetric Yukawa eigenvalue ratio)
    // ═══════════════════════════════════════════════════════════
    let h_minus_1 = h - 1; // 5
    let phi_up_num = (h_minus_1 as u64).pow(4);      // 5⁴ = 625
    let phi_up_den = (h as u64).pow(5);               // 6⁵ = 7776
    let phi_up = S::from_u64(phi_up_num) / S::from_u64(phi_up_den);
    // dim(∧²(fund)) = C(N,2) where N = rank(SU5)+1 = 5
    let n_su5 = SU5.rank + 1; // 5
    let dim_antisym = n_su5 * (n_su5 - 1) / 2; // C(5,2) = 10
    let up_params = KoideParams {
        r_fourth: S::from_u64(dim_antisym as u64),
        phi: phi_up,
        sigma: sums.up.clone(),
    };
    let up = koide_masses(&up_params);
    // k=0 → top, k=1 → up, k=2 → charm

    // ═══════════════════════════════════════════════════════════
    // Down quarks: r⁴ = dim(∧²(5)) - |α|_E8 = 10 - √2
    //              (cross-rep Yukawa 10×5̄ loses one root norm)
    //              φ = 1/h(G₂) = 1/6
    // ═══════════════════════════════════════════════════════════
    let r4_down = S::from_u64(dim_antisym as u64) - S::from_u64(2).sqrt();
    let phi_down = S::one() / S::from_u64(h as u64);
    let down_params = KoideParams {
        r_fourth: r4_down,
        phi: phi_down,
        sigma: sums.down.clone(),
    };
    let down = koide_masses(&down_params);
    // k=0 → bottom, k=1 → down, k=2 → strange

    // ═══════════════════════════════════════════════════════════
    // Neutrinos: r⁴ = 4 (same as leptons)
    //            φ = 2/9 + π/|W(G₂)| (Majorana G₂ Weyl shift)
    // Masses in meV (multiply Σ by 10⁹ to convert MeV → meV)
    // ═══════════════════════════════════════════════════════════
    let phi_nu = S::from_u64(2) / S::from_u64(9)
        + S::pi() / S::from_u64(w);
    let sigma_nu_mev = sums.neutrino.clone() * S::from_f64(1e9);
    let nu_params = KoideParams {
        r_fourth: S::from_u64(4),
        phi: phi_nu,
        sigma: sigma_nu_mev,
    };
    let nu = koide_masses(&nu_params);
    // k=0 → ν₃ (heaviest), k=1 → ν₁ (lightest), k=2 → ν₂ (middle)

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
    fn test_all_masses() {
        set_precision(50);
        let m: AllMasses<DefaultScalar> = compute_all_masses();

        // Leptons
        let e = m.electron.to_f64();
        let mu = m.muon.to_f64();
        let tau = m.tau.to_f64();
        assert!((e - 0.511).abs() / 0.511 < 0.01, "e = {}", e);
        assert!((mu - 105.658).abs() / 105.658 < 0.01, "mu = {}", mu);
        assert!((tau - 1776.86).abs() / 1776.86 < 0.01, "tau = {}", tau);

        // Up quarks
        let u = m.up.to_f64();
        let c = m.charm.to_f64();
        let t = m.top.to_f64();
        assert!(u > 1.5 && u < 3.5, "u = {}", u);
        assert!(c > 1200.0 && c < 1350.0, "c = {}", c);
        assert!(t > 170_000.0 && t < 176_000.0, "t = {}", t);

        // Down quarks
        let d = m.down.to_f64();
        let s = m.strange.to_f64();
        let b = m.bottom.to_f64();
        assert!(d > 3.0 && d < 6.0, "d = {}", d);
        assert!(s > 85.0 && s < 105.0, "s = {}", s);
        assert!(b > 4100.0 && b < 4300.0, "b = {}", b);
    }
}
