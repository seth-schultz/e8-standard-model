//! # E8 Standard Model
//!
//! Compute all 49 Standard Model parameters from the E8 lattice axiom
//! with zero free parameters.
//!
//! ## Quick Start
//!
//! ```rust
//! use e8_core::prelude::*;
//!
//! let predictions = compute_scorecard(50); // 50-digit precision
//! for p in &predictions {
//!     println!("{}: {:.6}", p.name, p.predicted);
//! }
//! ```
//!
//! ## SDK Usage
//!
//! ```rust,ignore
//! use e8_core::prelude::*;
//!
//! let theory = e8_standard_model();
//! let masses: AllMasses<DefaultScalar> = theory.masses();
//! println!("m_e = {} MeV", masses.electron.to_f64());
//! ```

pub mod precision;
pub mod special;
pub mod algebra;
pub mod octonion;
pub mod lattice;
pub mod mass;
pub mod coupling;
pub mod mixing;
pub mod higgs;
pub mod theory;
pub mod scorecard;
pub mod prelude;
