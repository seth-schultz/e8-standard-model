//! Core numeric abstraction: the Scalar trait.
//!
//! All physics computations are generic over `S: Scalar`, allowing the
//! same code to run with `f64` (fast, ~15 digits) or `rug::Float`
//! (arbitrary precision, default).

use std::fmt::{Debug, Display};
use std::ops::{Add, Sub, Mul, Div, Neg};

/// A scalar number type suitable for physics computations.
///
/// Implementations must provide owned-value arithmetic (clone operands
/// rather than requiring `&Self` operator impls). This keeps generic
/// code straightforward at the cost of occasional extra clones.
pub trait Scalar:
    Clone + Debug + Display + PartialOrd + Sized
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    fn from_u64(v: u64) -> Self;
    fn from_i64(v: i64) -> Self;
    fn from_f64(v: f64) -> Self;
    fn to_f64(&self) -> f64;

    fn zero() -> Self;
    fn one() -> Self;

    fn pi() -> Self;
    fn euler_gamma() -> Self;

    /// e^{-γ} (Mertens constant).
    fn exp_neg_gamma() -> Self {
        Self::euler_gamma().neg().exp()
    }

    /// Planck mass in MeV: 1.220890 × 10²² MeV.
    fn planck_mass_mev() -> Self;

    /// Planck mass in GeV: 1.220890 × 10¹⁹ GeV.
    fn planck_mass_gev() -> Self;

    fn sqrt(&self) -> Self;
    fn exp(&self) -> Self;
    fn ln(&self) -> Self;
    fn sin(&self) -> Self;
    fn cos(&self) -> Self;
    fn tan(&self) -> Self;
    fn acos(&self) -> Self;
    fn abs(&self) -> Self;

    /// self^exp via exp(exp * ln(self)).
    fn pow(&self, exp: &Self) -> Self;

    /// self^n for integer exponent.
    fn powi(&self, n: i64) -> Self;
}

// ─────────────────────────────────────────────────────────────────────
// f64 implementation
// ─────────────────────────────────────────────────────────────────────

impl Scalar for f64 {
    #[inline]
    fn from_u64(v: u64) -> Self { v as f64 }
    #[inline]
    fn from_i64(v: i64) -> Self { v as f64 }
    #[inline]
    fn from_f64(v: f64) -> Self { v }
    #[inline]
    fn to_f64(&self) -> f64 { *self }

    #[inline]
    fn zero() -> Self { 0.0 }
    #[inline]
    fn one() -> Self { 1.0 }

    #[inline]
    fn pi() -> Self { std::f64::consts::PI }
    #[inline]
    fn euler_gamma() -> Self {
        // Euler–Mascheroni constant to full f64 precision.
        0.5772156649015328606065120900824024310421_f64
    }

    fn planck_mass_mev() -> Self { 1.220890e22 }
    fn planck_mass_gev() -> Self { 1.220890e19 }

    #[inline]
    fn sqrt(&self) -> Self { f64::sqrt(*self) }
    #[inline]
    fn exp(&self) -> Self { f64::exp(*self) }
    #[inline]
    fn ln(&self) -> Self { f64::ln(*self) }
    #[inline]
    fn sin(&self) -> Self { f64::sin(*self) }
    #[inline]
    fn cos(&self) -> Self { f64::cos(*self) }
    #[inline]
    fn tan(&self) -> Self { f64::tan(*self) }
    #[inline]
    fn acos(&self) -> Self { f64::acos(*self) }
    #[inline]
    fn abs(&self) -> Self { f64::abs(*self) }

    #[inline]
    fn pow(&self, exp: &Self) -> Self { f64::powf(*self, *exp) }
    #[inline]
    fn powi(&self, n: i64) -> Self { f64::powi(*self, n as i32) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f64_basic_ops() {
        let a = f64::from_u64(3);
        let b = f64::from_u64(4);
        assert_eq!((a + b).to_f64(), 7.0);
        assert_eq!((a.clone() * b.clone()).to_f64(), 12.0);
    }

    #[test]
    fn f64_transcendentals() {
        let pi = f64::pi();
        assert!((pi - std::f64::consts::PI).abs() < 1e-15);
        assert!((pi.sin()).abs() < 1e-15);
        assert!((pi.cos() + 1.0).abs() < 1e-15);
    }

    #[test]
    fn f64_euler_gamma() {
        let g = f64::euler_gamma();
        assert!((g - 0.5772156649015329).abs() < 1e-15);
        let eng = f64::exp_neg_gamma();
        assert!((eng - 0.56145948356).abs() < 1e-8);
    }

    #[test]
    fn f64_pow() {
        let two = f64::from_u64(2);
        let ten = two.powi(10);
        assert_eq!(ten, 1024.0);
    }
}
