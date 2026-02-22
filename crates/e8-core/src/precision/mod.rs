//! Adaptive precision arithmetic backend.
//!
//! The [`Scalar`] trait abstracts over number types. With the default
//! `arbitrary-precision` feature, `rug::Float` is the `DefaultScalar`.
//! Without it, `f64` is used.

pub mod scalar;

#[cfg(feature = "arbitrary-precision")]
mod scalar_mpfr;

pub use scalar::Scalar;

/// The default scalar type used throughout E8 computations.
#[cfg(feature = "arbitrary-precision")]
pub type DefaultScalar = rug::Float;

/// The default scalar type used throughout E8 computations (f64 mode).
#[cfg(not(feature = "arbitrary-precision"))]
pub type DefaultScalar = f64;

// ─────────────────────────────────────────────────────────────────────
// Precision control
// ─────────────────────────────────────────────────────────────────────

#[cfg(feature = "arbitrary-precision")]
use std::cell::Cell;

#[cfg(feature = "arbitrary-precision")]
thread_local! {
    static PRECISION_BITS: Cell<u32> = const { Cell::new(832) }; // ~250 decimal digits
}

/// Set working precision in decimal digits.
///
/// With `arbitrary-precision`, this configures the MPFR bit width.
/// Without it, this is a no-op (f64 has fixed ~15-digit precision).
#[cfg(feature = "arbitrary-precision")]
pub fn set_precision(digits: u32) {
    let bits = (digits as f64 * std::f64::consts::LN_10 / std::f64::consts::LN_2).ceil() as u32 + 10;
    PRECISION_BITS.with(|p| p.set(bits));
}

/// No-op precision setter for f64 mode.
#[cfg(not(feature = "arbitrary-precision"))]
pub fn set_precision(_digits: u32) {}

/// Get current working precision in bits.
#[cfg(feature = "arbitrary-precision")]
pub fn precision_bits() -> u32 {
    PRECISION_BITS.with(|p| p.get())
}

/// Get current working precision in decimal digits.
#[cfg(feature = "arbitrary-precision")]
pub fn precision_digits() -> u32 {
    (precision_bits() as f64 * std::f64::consts::LN_2 / std::f64::consts::LN_10).floor() as u32
}

// ─────────────────────────────────────────────────────────────────────
// Legacy MPFR helpers (gated behind arbitrary-precision)
// ─────────────────────────────────────────────────────────────────────

#[cfg(feature = "arbitrary-precision")]
mod legacy {
    use rug::float::Round;
    use rug::ops::{AssignRound, Pow};
    use rug::Float;
    use super::precision_bits;

    /// Create a high-precision float from u64.
    pub fn mpf_u64(n: u64) -> Float {
        Float::with_val(precision_bits(), n)
    }

    /// Create a high-precision float from i64.
    pub fn mpf_i64(n: i64) -> Float {
        Float::with_val(precision_bits(), n)
    }

    /// Create a high-precision float from f64.
    pub fn mpf_f64(x: f64) -> Float {
        Float::with_val(precision_bits(), x)
    }

    /// Create a high-precision float from a decimal string.
    pub fn mpf_str(s: &str) -> Float {
        let prec = precision_bits();
        let mut f = Float::new(prec);
        f.assign_round(Float::parse(s).unwrap(), Round::Nearest);
        f
    }

    /// Pi at current precision.
    pub fn pi() -> Float {
        Float::with_val(precision_bits(), rug::float::Constant::Pi)
    }

    /// Euler-Mascheroni gamma at current precision.
    pub fn euler_gamma() -> Float {
        Float::with_val(precision_bits(), rug::float::Constant::Euler)
    }

    /// e^{-gamma} at current precision (Mertens constant).
    pub fn exp_neg_gamma() -> Float {
        let g = euler_gamma();
        (-g).exp()
    }

    /// Planck mass in MeV.
    pub fn planck_mass_mev() -> Float {
        mpf_str("1.220890e22")
    }

    /// Planck mass in GeV.
    pub fn planck_mass_gev() -> Float {
        mpf_str("1.220890e19")
    }

    /// Create a zero at current precision.
    pub fn zero() -> Float {
        Float::with_val(precision_bits(), 0)
    }

    /// Create a one at current precision.
    pub fn one() -> Float {
        Float::with_val(precision_bits(), 1)
    }

    /// Square root of a high-precision float.
    pub fn sqrt(x: &Float) -> Float {
        x.clone().sqrt()
    }

    /// Absolute value.
    pub fn abs(x: &Float) -> Float {
        x.clone().abs()
    }

    /// Natural logarithm.
    pub fn ln(x: &Float) -> Float {
        x.clone().ln()
    }

    /// Exponential.
    pub fn exp(x: &Float) -> Float {
        x.clone().exp()
    }

    /// Cosine.
    pub fn cos(x: &Float) -> Float {
        x.clone().cos()
    }

    /// Sine.
    pub fn sin(x: &Float) -> Float {
        x.clone().sin()
    }

    /// Tangent.
    pub fn tan(x: &Float) -> Float {
        x.clone().tan()
    }

    /// Arccosine.
    pub fn acos(x: &Float) -> Float {
        x.clone().acos()
    }

    /// Power: base^exp.
    pub fn pow(base: &Float, exp: &Float) -> Float {
        let prec = precision_bits();
        Float::with_val(prec, (exp * base.clone().ln()).exp())
    }

    /// Integer power.
    pub fn powi(base: &Float, n: i32) -> Float {
        let prec = precision_bits();
        if n >= 0 {
            Float::with_val(prec, Pow::pow(base, n as u32))
        } else {
            Float::with_val(prec, Pow::pow(base, (-n) as u32)).recip()
        }
    }
}

#[cfg(feature = "arbitrary-precision")]
pub use legacy::*;

#[cfg(test)]
#[cfg(feature = "arbitrary-precision")]
mod tests {
    use super::*;

    #[test]
    fn test_precision_setup() {
        set_precision(250);
        assert!(precision_bits() >= 830);
        assert!(precision_digits() >= 249);
    }

    #[test]
    fn test_constants() {
        set_precision(50);
        let p = pi();
        let p_f64 = p.to_f64();
        assert!(p_f64 > 3.0 && p_f64 < 4.0);
        assert!((p_f64 - std::f64::consts::PI).abs() < 1e-10);

        let g = euler_gamma();
        assert!(g > 0.5772 && g < 0.5773);

        let eng = exp_neg_gamma();
        assert!(eng > 0.5614 && eng < 0.5615);
    }

    #[test]
    fn test_planck_mass() {
        set_precision(50);
        let mp = planck_mass_mev();
        assert!(mp > 1.22e22 && mp < 1.23e22);
    }
}
