//! Prediction structure for scorecard entries.

use serde::{Deserialize, Serialize};

/// The status of a derived quantity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Status {
    /// Mathematically proven from axioms alone.
    Theorem,
    /// Derived through a complete chain from axioms.
    Derived,
    /// Derived, but one step in the chain is identified (not computed from scratch).
    DerivedStar,
    /// Structurally identified but not fully derived.
    Identified,
}

impl std::fmt::Display for Status {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Status::Theorem => write!(f, "THEOREM"),
            Status::Derived => write!(f, "DERIVED"),
            Status::DerivedStar => write!(f, "DERIVED*"),
            Status::Identified => write!(f, "IDENTIFIED"),
        }
    }
}

/// Category of a prediction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Category {
    Gauge,
    MassSum,
    Koide,
    LeptonMass,
    UpQuarkMass,
    DownQuarkMass,
    CKM,
    Neutrino,
    PMNS,
    Higgs,
    Structural,
}

impl std::fmt::Display for Category {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Category::Gauge => write!(f, "Gauge"),
            Category::MassSum => write!(f, "Mass sums"),
            Category::Koide => write!(f, "Koide"),
            Category::LeptonMass => write!(f, "Lepton masses"),
            Category::UpQuarkMass => write!(f, "Up quark masses"),
            Category::DownQuarkMass => write!(f, "Down quark masses"),
            Category::CKM => write!(f, "CKM"),
            Category::Neutrino => write!(f, "Neutrino"),
            Category::PMNS => write!(f, "PMNS"),
            Category::Higgs => write!(f, "Higgs"),
            Category::Structural => write!(f, "Structural"),
        }
    }
}

/// A single scorecard prediction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Prediction {
    pub name: String,
    pub category: Category,
    pub formula: String,
    pub predicted: f64,
    pub experimental: Option<f64>,
    pub exp_uncertainty: Option<f64>,
    pub unit: String,
    pub status: Status,
    pub pull_sigma: Option<f64>,
    pub pct_error: Option<f64>,
}

impl Prediction {
    /// Create a new prediction and compute pull/error.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: &str,
        category: Category,
        formula: &str,
        predicted: f64,
        experimental: Option<f64>,
        exp_uncertainty: Option<f64>,
        unit: &str,
        status: Status,
    ) -> Self {
        let pull_sigma = match (experimental, exp_uncertainty) {
            (Some(exp), Some(err)) if err > 0.0 => Some((predicted - exp) / err),
            _ => None,
        };

        let pct_error = experimental.map(|exp| {
            if exp.abs() > 0.0 {
                (predicted - exp) / exp * 100.0
            } else {
                0.0
            }
        });

        Self {
            name: name.to_string(),
            category,
            formula: formula.to_string(),
            predicted,
            experimental,
            exp_uncertainty,
            unit: unit.to_string(),
            status,
            pull_sigma,
            pct_error,
        }
    }
}
