//! The E8 mass formula: Σ = f × m_P × exp(-(A×R + δ)/dim(so(8))).

use crate::algebra::groups::identities::N_C;
use crate::algebra::groups::{G2, SU3};
use crate::precision::scalar::Scalar;

use super::constants::{delta, m_planck_mev, mertens_r, norm_factor};

/// All four sector sums from the E8 mass formula.
pub struct SectorSums<S: Scalar> {
    pub leptons: S,  // A = dim(su(3))+1 = 9, f = 1
    pub up: S,       // A = dim(su(3)) = 8, f = 1/C₂(SU3,fund)
    pub down: S,     // A = dim(su(3))+1 = 9, f = N_c/C₂(SU3,fund)
    pub neutrino: S, // A = dim(G₂) = 14, f = √((|W(G₂)|-rank(G₂))/(|W(G₂)|+1))
}

/// Compute sector mass sum: Σ = f × m_P × exp(-(A × R + δ)/28).
///
/// # Arguments
/// * `a` - The representation quantum number (A-value)
/// * `f` - The lattice gauge dressing factor
pub fn sector_sum<S: Scalar>(a: u32, f: &S) -> S {
    let r: S = mertens_r();
    let d: S = delta();
    let n28: S = norm_factor();
    let mp: S = m_planck_mev();

    let exponent = (S::from_u64(a as u64) * r + d).neg() / n28;
    let boltzmann = exponent.exp();

    f.clone() * mp * boltzmann
}

/// Compute sector sum with precise f = p/q rational factor.
pub fn sector_sum_rational<S: Scalar>(a: u32, f_num: u64, f_den: u64) -> S {
    let f = S::from_u64(f_num) / S::from_u64(f_den);
    sector_sum(a, &f)
}

/// Compute all four sector mass sums.
///
/// A-values from representation theory:
///   leptons (1 of SU5): A = dim(u(3)) = dim(su(3)) + 1 = 9
///   up (10 of SU5):     A = dim(su(3)) = 8
///   down (5̄ of SU5):   A = dim(u(3)) = dim(su(3)) + 1 = 9
///   neutrinos:          A = dim(G₂) = 14
///
/// f-factors from color Casimir scaling:
///   f_lep  = 1 (color singlet)
///   f_up   = 1/C₂(SU3,fund) = 3/4
///   f_down = N_c/C₂(SU3,fund) = 9/4
///   f_ν    = √((|W(G₂)| - rank(G₂)) / (|W(G₂)| + 1)) = √(10/13)
pub fn compute_all_sector_sums<S: Scalar>() -> SectorSums<S> {
    let (c2_num, c2_den) = SU3.c2_fundamental; // C₂(SU3,fund) = 4/3

    // A-values
    let a_lep = SU3.dimension + 1;  // 9 = dim(u(3))
    let a_up = SU3.dimension;       // 8 = dim(su(3))
    let a_down = SU3.dimension + 1; // 9 = dim(u(3))
    let a_nu = G2.dimension;        // 14 = dim(G₂)

    // f-factors
    let f_lep = S::one();

    // f_up = 1/C₂ = den/num = 3/4
    let f_up = S::from_u64(c2_den as u64) / S::from_u64(c2_num as u64);

    // f_down = N_c/C₂ = N_c × den/num = 3 × 3/4 = 9/4
    let f_down = S::from_u64((N_C * c2_den) as u64) / S::from_u64(c2_num as u64);

    // f_ν = √((|W(G₂)| - rank(G₂)) / (|W(G₂)| + 1)) = √(10/13)
    let w_minus_r = G2.weyl_order - G2.rank as u64; // 12 - 2 = 10
    let w_plus_1 = G2.weyl_order + 1;                // 12 + 1 = 13
    let f_nu = (S::from_u64(w_minus_r) / S::from_u64(w_plus_1)).sqrt();

    SectorSums {
        leptons: sector_sum(a_lep, &f_lep),
        up: sector_sum(a_up, &f_up),
        down: sector_sum(a_down, &f_down),
        neutrino: sector_sum(a_nu, &f_nu),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::DefaultScalar;

    #[test]
    fn test_lepton_sum() {
        crate::precision::set_precision(50);
        let sigma: DefaultScalar = sector_sum_rational(9, 1, 1);
        let val = sigma.to_f64();
        // Σ_lep ≈ 1883 MeV (sum of e + μ + τ)
        assert!(val > 1880.0 && val < 1890.0, "Σ_lep = {}", val);
    }

    #[test]
    fn test_up_sum() {
        crate::precision::set_precision(50);
        let sigma: DefaultScalar = sector_sum_rational(8, 3, 4);
        let val = sigma.to_f64();
        // Σ_up ≈ 174,000 MeV (sum of u + c + t)
        assert!(val > 170_000.0 && val < 180_000.0, "Σ_up = {}", val);
    }

    #[test]
    fn test_down_sum() {
        crate::precision::set_precision(50);
        let sigma: DefaultScalar = sector_sum_rational(9, 9, 4);
        let val = sigma.to_f64();
        // Σ_down ≈ 4237 MeV (sum of d + s + b)
        assert!(val > 4200.0 && val < 4300.0, "Σ_down = {}", val);
    }

    #[test]
    fn test_a_values() {
        // Verify A-values match expected
        assert_eq!(SU3.dimension + 1, 9); // leptons, down
        assert_eq!(SU3.dimension, 8);     // up
        assert_eq!(G2.dimension, 14);     // neutrinos
    }

    #[test]
    fn test_f_factors() {
        // Verify f-factors: 1/C₂ = 3/4, N_c/C₂ = 9/4
        let (num, den) = SU3.c2_fundamental;
        assert_eq!(den, 3); // f_up numerator
        assert_eq!(num, 4); // f_up denominator
        assert_eq!(N_C * den, 9); // f_down numerator
        assert_eq!(num, 4);       // f_down denominator
    }
}
