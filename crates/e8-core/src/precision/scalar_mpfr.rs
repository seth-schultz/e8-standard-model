//! Scalar implementation for `rug::Float` (arbitrary-precision via MPFR).
//!
//! Gated behind `#[cfg(feature = "arbitrary-precision")]`.

use rug::float::Round;
use rug::ops::{AssignRound, Pow};
use rug::Float;

use super::precision_bits;
use super::scalar::Scalar;

impl Scalar for Float {
    fn from_u64(v: u64) -> Self {
        Float::with_val(precision_bits(), v)
    }

    fn from_i64(v: i64) -> Self {
        Float::with_val(precision_bits(), v)
    }

    fn from_f64(v: f64) -> Self {
        Float::with_val(precision_bits(), v)
    }

    fn to_f64(&self) -> f64 {
        self.to_f64()
    }

    fn zero() -> Self {
        Float::with_val(precision_bits(), 0)
    }

    fn one() -> Self {
        Float::with_val(precision_bits(), 1)
    }

    fn pi() -> Self {
        Float::with_val(precision_bits(), rug::float::Constant::Pi)
    }

    fn euler_gamma() -> Self {
        Float::with_val(precision_bits(), rug::float::Constant::Euler)
    }

    fn exp_neg_gamma() -> Self {
        let g = <Float as Scalar>::euler_gamma();
        (-g).exp()
    }

    fn planck_mass_mev() -> Self {
        let prec = precision_bits();
        let mut f = Float::new(prec);
        f.assign_round(Float::parse("1.220890e22").unwrap(), Round::Nearest);
        f
    }

    fn planck_mass_gev() -> Self {
        let prec = precision_bits();
        let mut f = Float::new(prec);
        f.assign_round(Float::parse("1.220890e19").unwrap(), Round::Nearest);
        f
    }

    fn sqrt(&self) -> Self {
        self.clone().sqrt()
    }

    fn exp(&self) -> Self {
        self.clone().exp()
    }

    fn ln(&self) -> Self {
        self.clone().ln()
    }

    fn sin(&self) -> Self {
        self.clone().sin()
    }

    fn cos(&self) -> Self {
        self.clone().cos()
    }

    fn tan(&self) -> Self {
        self.clone().tan()
    }

    fn acos(&self) -> Self {
        self.clone().acos()
    }

    fn abs(&self) -> Self {
        self.clone().abs()
    }

    fn pow(&self, exp: &Self) -> Self {
        let prec = precision_bits();
        Float::with_val(prec, (exp * self.clone().ln()).exp())
    }

    fn powi(&self, n: i64) -> Self {
        let prec = precision_bits();
        if n >= 0 {
            Float::with_val(prec, Pow::pow(self, n as u32))
        } else {
            Float::with_val(prec, Pow::pow(self, (-n) as u32)).recip()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precision::set_precision;

    #[test]
    fn mpfr_basic_ops() {
        set_precision(50);
        let a = Float::from_u64(3);
        let b = Float::from_u64(4);
        assert_eq!((a.clone() + b.clone()).to_f64(), 7.0);
        assert_eq!((a * b).to_f64(), 12.0);
    }

    #[test]
    fn mpfr_transcendentals() {
        set_precision(50);
        let pi = <Float as Scalar>::pi();
        assert!((pi.to_f64() - std::f64::consts::PI).abs() < 1e-15);
    }

    #[test]
    fn mpfr_matches_f64() {
        set_precision(50);
        let mpfr_gamma = <Float as Scalar>::euler_gamma().to_f64();
        let f64_gamma = <f64 as Scalar>::euler_gamma();
        assert!((mpfr_gamma - f64_gamma).abs() < 1e-15);
    }

    #[test]
    fn mpfr_pow() {
        set_precision(50);
        let two = Float::from_u64(2);
        let ten = Scalar::powi(&two, 10);
        assert_eq!(ten.to_f64(), 1024.0);
    }

    #[test]
    fn mpfr_planck() {
        set_precision(50);
        let mp = <Float as Scalar>::planck_mass_mev();
        assert!(mp.to_f64() > 1.22e22 && mp.to_f64() < 1.23e22);
    }
}
