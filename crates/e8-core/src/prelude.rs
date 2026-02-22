//! Convenience re-exports for common usage.
//!
//! ```rust,ignore
//! use e8_core::prelude::*;
//!
//! let theory = e8_standard_model();
//! let masses: AllMasses<f64> = theory.masses();
//! ```

// Scalar abstraction
pub use crate::precision::scalar::Scalar;
pub use crate::precision::DefaultScalar;

// Theory composition
pub use crate::theory::{e8_standard_model, E8StandardModel, Theory};

// Core data types
pub use crate::mass::sectors::AllMasses;
pub use crate::mixing::ckm::CKMResult;
pub use crate::scorecard::prediction::Prediction;
pub use crate::scorecard::table::compute_scorecard;

// Physics traits
pub use crate::coupling::GaugeCouplings;
pub use crate::higgs::HiggsSector;
pub use crate::lattice::RootSystem;
pub use crate::mass::{MassFormula, MassSplitting};
pub use crate::mixing::{CKMMixing, PMNSMixing};
pub use crate::algebra::LieAlgebra;
