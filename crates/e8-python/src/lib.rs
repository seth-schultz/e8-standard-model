//! Python bindings for the E8 Standard Model SDK.

use pyo3::prelude::*;

use e8_core::precision::DefaultScalar;

/// A single scorecard prediction exposed to Python.
#[pyclass]
#[derive(Clone)]
struct Prediction {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    category: String,
    #[pyo3(get)]
    formula: String,
    #[pyo3(get)]
    predicted: f64,
    #[pyo3(get)]
    experimental: Option<f64>,
    #[pyo3(get)]
    exp_uncertainty: Option<f64>,
    #[pyo3(get)]
    unit: String,
    #[pyo3(get)]
    status: String,
    #[pyo3(get)]
    pull_sigma: Option<f64>,
    #[pyo3(get)]
    pct_error: Option<f64>,
}

#[pymethods]
impl Prediction {
    fn __repr__(&self) -> String {
        format!("Prediction({}: {:.6} {})", self.name, self.predicted, self.unit)
    }
}

/// All fermion masses from the E8 theory.
#[pyclass]
#[derive(Clone)]
struct Masses {
    #[pyo3(get)]
    electron: f64,
    #[pyo3(get)]
    muon: f64,
    #[pyo3(get)]
    tau: f64,
    #[pyo3(get)]
    up: f64,
    #[pyo3(get)]
    charm: f64,
    #[pyo3(get)]
    top: f64,
    #[pyo3(get)]
    down: f64,
    #[pyo3(get)]
    strange: f64,
    #[pyo3(get)]
    bottom: f64,
    #[pyo3(get)]
    nu1: f64,
    #[pyo3(get)]
    nu2: f64,
    #[pyo3(get)]
    nu3: f64,
}

#[pymethods]
impl Masses {
    fn __repr__(&self) -> String {
        format!(
            "Masses(e={:.4} MeV, μ={:.2} MeV, τ={:.1} MeV, t={:.0} MeV)",
            self.electron, self.muon, self.tau, self.top
        )
    }
}

/// The E8 Standard Model: 49 predictions from zero free parameters.
#[pyclass]
struct E8StandardModel {
    digits: u32,
}

#[pymethods]
impl E8StandardModel {
    #[new]
    #[pyo3(signature = (digits=50))]
    fn new(digits: u32) -> Self {
        e8_core::precision::set_precision(digits);
        E8StandardModel { digits }
    }

    /// Compute the full scorecard of 49 predictions.
    fn scorecard(&self) -> Vec<Prediction> {
        let preds = e8_core::scorecard::table::compute_scorecard(self.digits);
        preds
            .into_iter()
            .map(|p| Prediction {
                name: p.name,
                category: p.category.to_string(),
                formula: p.formula,
                predicted: p.predicted,
                experimental: p.experimental,
                exp_uncertainty: p.exp_uncertainty,
                unit: p.unit,
                status: p.status.to_string(),
                pull_sigma: p.pull_sigma,
                pct_error: p.pct_error,
            })
            .collect()
    }

    /// Compute all 12 fermion masses.
    fn masses(&self) -> Masses {
        e8_core::precision::set_precision(self.digits);
        let m = e8_core::mass::sectors::compute_all_masses::<DefaultScalar>();
        Masses {
            electron: m.electron.to_f64(),
            muon: m.muon.to_f64(),
            tau: m.tau.to_f64(),
            up: m.up.to_f64(),
            charm: m.charm.to_f64(),
            top: m.top.to_f64(),
            down: m.down.to_f64(),
            strange: m.strange.to_f64(),
            bottom: m.bottom.to_f64(),
            nu1: m.nu1.to_f64(),
            nu2: m.nu2.to_f64(),
            nu3: m.nu3.to_f64(),
        }
    }

    /// 1/α(0) — fine structure constant inverse.
    fn alpha_inverse(&self) -> f64 {
        e8_core::precision::set_precision(self.digits);
        e8_core::coupling::alpha::alpha_inverse::<DefaultScalar>().to_f64()
    }

    /// sin²θ_W at M_Z scale.
    fn sin2_theta_w(&self) -> f64 {
        e8_core::precision::set_precision(self.digits);
        e8_core::coupling::weinberg::sin2_theta_w_mz::<DefaultScalar>().to_f64()
    }

    /// α_s(M_Z) — strong coupling.
    fn alpha_s(&self) -> f64 {
        e8_core::precision::set_precision(self.digits);
        e8_core::coupling::alpha_s::alpha_s_mz::<DefaultScalar>().to_f64()
    }

    /// Higgs mass in GeV.
    fn higgs_mass(&self) -> f64 {
        e8_core::precision::set_precision(self.digits);
        e8_core::higgs::mass::higgs_mass_default::<DefaultScalar>().to_f64()
    }

    /// Higgs quartic coupling λ.
    fn higgs_quartic(&self) -> f64 {
        e8_core::precision::set_precision(self.digits);
        e8_core::higgs::quartic::higgs_quartic::<DefaultScalar>().to_f64()
    }

    fn __repr__(&self) -> String {
        format!("E8StandardModel(digits={})", self.digits)
    }
}

/// E8 Standard Model Python bindings.
#[pymodule]
fn e8(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<E8StandardModel>()?;
    m.add_class::<Prediction>()?;
    m.add_class::<Masses>()?;
    Ok(())
}
