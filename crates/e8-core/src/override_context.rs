//! Parameter override system for theory testing.
//!
//! Scientists can override specific derivable parameters to test alternative
//! theories. All parameters have strict bounds and type enforcement.
//! No user code runs — only numeric overrides from a fixed catalog.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// A single parameter in the catalog.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamSpec {
    pub key: &'static str,
    pub description: &'static str,
    pub category: ParamCategory,
    pub default_value: f64,
    pub min: f64,
    pub max: f64,
}

/// Parameter categories for the theory wizard.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParamCategory {
    Mass,
    Coupling,
    Mixing,
    Higgs,
    Neutrino,
}

impl std::fmt::Display for ParamCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamCategory::Mass => write!(f, "Mass"),
            ParamCategory::Coupling => write!(f, "Coupling"),
            ParamCategory::Mixing => write!(f, "Mixing"),
            ParamCategory::Higgs => write!(f, "Higgs"),
            ParamCategory::Neutrino => write!(f, "Neutrino"),
        }
    }
}

/// The full parameter catalog (~40 overridable parameters).
pub fn parameter_catalog() -> Vec<ParamSpec> {
    vec![
        // ═══════════════════════════════════════════
        // Mass sector: A-values
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "a_lepton",
            description: "A-value for leptons (default: dim(u(3)) = 9)",
            category: ParamCategory::Mass,
            default_value: 9.0,
            min: 1.0,
            max: 30.0,
        },
        ParamSpec {
            key: "a_up",
            description: "A-value for up quarks (default: dim(su(3)) = 8)",
            category: ParamCategory::Mass,
            default_value: 8.0,
            min: 1.0,
            max: 30.0,
        },
        ParamSpec {
            key: "a_down",
            description: "A-value for down quarks (default: dim(u(3)) = 9)",
            category: ParamCategory::Mass,
            default_value: 9.0,
            min: 1.0,
            max: 30.0,
        },
        ParamSpec {
            key: "a_neutrino",
            description: "A-value for neutrinos (default: dim(G2) = 14)",
            category: ParamCategory::Mass,
            default_value: 14.0,
            min: 1.0,
            max: 30.0,
        },

        // ═══════════════════════════════════════════
        // Mass sector: QCD ΔA corrections
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "delta_a_up",
            description: "QCD ΔA for up quarks: -(α_s × C_F)/(π × 216) ≈ -0.000232",
            category: ParamCategory::Mass,
            default_value: -0.000232, // approximate; exact value computed at runtime
            min: -0.01,
            max: 0.01,
        },
        ParamSpec {
            key: "delta_a_down",
            description: "QCD ΔA for down quarks: ΔA_up × 61/7 ≈ -0.00202",
            category: ParamCategory::Mass,
            default_value: -0.00202, // approximate; exact value computed at runtime
            min: -0.1,
            max: 0.1,
        },

        // ═══════════════════════════════════════════
        // Mass sector: f-factors
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "f_lepton",
            description: "Lattice dressing factor for leptons (default: 1)",
            category: ParamCategory::Mass,
            default_value: 1.0,
            min: 0.01,
            max: 10.0,
        },
        ParamSpec {
            key: "f_up",
            description: "Lattice dressing factor for up quarks (default: 3/4)",
            category: ParamCategory::Mass,
            default_value: 0.75,
            min: 0.01,
            max: 10.0,
        },
        ParamSpec {
            key: "f_down",
            description: "Lattice dressing factor for down quarks (default: 9/4)",
            category: ParamCategory::Mass,
            default_value: 2.25,
            min: 0.01,
            max: 10.0,
        },
        ParamSpec {
            key: "f_neutrino",
            description: "Lattice dressing factor for neutrinos (default: sqrt(10/13))",
            category: ParamCategory::Neutrino,
            default_value: 0.8770580193070293, // sqrt(10/13)
            min: 0.01,
            max: 10.0,
        },

        // ═══════════════════════════════════════════
        // Koide parameters: r^4
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "r4_lepton",
            description: "Koide r^4 for leptons (default: 4 = (sqrt2)^4)",
            category: ParamCategory::Mass,
            default_value: 4.0,
            min: 0.1,
            max: 50.0,
        },
        ParamSpec {
            key: "r4_up",
            description: "Koide r^4 for up quarks (default: 10 = dim(antisym 5))",
            category: ParamCategory::Mass,
            default_value: 10.0,
            min: 0.1,
            max: 50.0,
        },
        ParamSpec {
            key: "r4_down",
            description: "Koide r^4 for down quarks (default: 10 - sqrt2)",
            category: ParamCategory::Mass,
            default_value: 8.585786437626905, // 10 - sqrt(2)
            min: 0.1,
            max: 50.0,
        },

        // ═══════════════════════════════════════════
        // Koide parameters: phi
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "phi_lepton",
            description: "Koide phase for leptons (default: 2/9)",
            category: ParamCategory::Mass,
            default_value: 0.2222222222222222, // 2/9
            min: -1.0,
            max: 1.0,
        },
        ParamSpec {
            key: "phi_up",
            description: "Koide phase for up quarks (default: 5^4/6^5 = 625/7776)",
            category: ParamCategory::Mass,
            default_value: 0.08037551440329218, // 625/7776
            min: -1.0,
            max: 1.0,
        },
        ParamSpec {
            key: "phi_down",
            description: "Koide phase for down quarks (default: 1/6)",
            category: ParamCategory::Mass,
            default_value: 0.16666666666666666, // 1/6
            min: -1.0,
            max: 1.0,
        },

        // ═══════════════════════════════════════════
        // Alpha CF coefficients
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "alpha_cf_a0",
            description: "CF coefficient a0 (default: 244 = |Phi(E8)| + rank/2)",
            category: ParamCategory::Coupling,
            default_value: 244.0,
            min: 100.0,
            max: 500.0,
        },
        ParamSpec {
            key: "alpha_cf_a1",
            description: "CF coefficient a1 (default: 14 = dim(G2))",
            category: ParamCategory::Coupling,
            default_value: 14.0,
            min: 1.0,
            max: 100.0,
        },
        ParamSpec {
            key: "alpha_cf_a2",
            description: "CF coefficient a2 (default: 13 = |W(G2)|+1)",
            category: ParamCategory::Coupling,
            default_value: 13.0,
            min: 1.0,
            max: 100.0,
        },
        ParamSpec {
            key: "alpha_cf_a3",
            description: "CF coefficient a3 (default: 193 = |W(D4)|+1)",
            category: ParamCategory::Coupling,
            default_value: 193.0,
            min: 1.0,
            max: 1000.0,
        },

        // ═══════════════════════════════════════════
        // Weinberg angle
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "weinberg_tree_num",
            description: "Weinberg tree-level numerator (default: 3)",
            category: ParamCategory::Coupling,
            default_value: 3.0,
            min: 0.1,
            max: 10.0,
        },
        ParamSpec {
            key: "weinberg_tree_den",
            description: "Weinberg tree-level denominator (default: 13 = |W(G2)|+1)",
            category: ParamCategory::Coupling,
            default_value: 13.0,
            min: 1.0,
            max: 100.0,
        },
        ParamSpec {
            key: "weinberg_correction_num",
            description: "Weinberg correction numerator (default: 5 = h(G2)-1)",
            category: ParamCategory::Coupling,
            default_value: 5.0,
            min: 0.0,
            max: 30.0,
        },
        ParamSpec {
            key: "weinberg_correction_den",
            description: "Weinberg correction denominator (default: 6 = h(G2))",
            category: ParamCategory::Coupling,
            default_value: 6.0,
            min: 1.0,
            max: 30.0,
        },

        // ═══════════════════════════════════════════
        // CKM mixing: self-energy corrections
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "ckm_d1_down_coeff",
            description: "D1 down self-energy coefficient (default: -1, i.e. -m_u)",
            category: ParamCategory::Mixing,
            default_value: -1.0,
            min: -20.0,
            max: 20.0,
        },
        ParamSpec {
            key: "ckm_d2_down_coeff",
            description: "D2 down self-energy coefficient (default: -8 = -dim(su3))",
            category: ParamCategory::Mixing,
            default_value: -8.0,
            min: -50.0,
            max: 50.0,
        },
        ParamSpec {
            key: "ckm_d1_up_coeff",
            description: "D1 up self-energy coefficient (default: 4/3 = C2(SU3,fund))",
            category: ParamCategory::Mixing,
            default_value: 1.3333333333333333,
            min: -20.0,
            max: 20.0,
        },
        ParamSpec {
            key: "ckm_d2_up_coeff",
            description: "D2 up self-energy coefficient (default: 26/9)",
            category: ParamCategory::Mixing,
            default_value: 2.888888888888889,
            min: -20.0,
            max: 20.0,
        },

        // ═══════════════════════════════════════════
        // CKM CP phase
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "cp_phase_numerator",
            description: "CP phase = num*pi/den (default num: 1, giving pi/7)",
            category: ParamCategory::Mixing,
            default_value: 1.0,
            min: 0.0,
            max: 14.0,
        },
        ParamSpec {
            key: "cp_phase_denominator",
            description: "CP phase = num*pi/den (default den: 7 = dim(Im(O)))",
            category: ParamCategory::Mixing,
            default_value: 7.0,
            min: 1.0,
            max: 100.0,
        },

        // ═══════════════════════════════════════════
        // PMNS mixing
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "pmns_weyl_order",
            description: "Weyl order for PMNS (default: 12 = |W(G2)|)",
            category: ParamCategory::Mixing,
            default_value: 12.0,
            min: 2.0,
            max: 200.0,
        },
        ParamSpec {
            key: "pmns_exponent_m1",
            description: "G2 exponent m1 for PMNS (default: 1)",
            category: ParamCategory::Mixing,
            default_value: 1.0,
            min: 1.0,
            max: 30.0,
        },
        ParamSpec {
            key: "pmns_exponent_m2",
            description: "G2 exponent m2 for PMNS (default: 5)",
            category: ParamCategory::Mixing,
            default_value: 5.0,
            min: 1.0,
            max: 30.0,
        },
        ParamSpec {
            key: "pmns_rank",
            description: "Rank for PMNS theta23 (default: 2 = rank(G2))",
            category: ParamCategory::Mixing,
            default_value: 2.0,
            min: 1.0,
            max: 20.0,
        },

        // ═══════════════════════════════════════════
        // Higgs quartic
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "higgs_lambda_num",
            description: "Lambda numerator dim (default: 7 = dim(Im(O)))",
            category: ParamCategory::Higgs,
            default_value: 7.0,
            min: 1.0,
            max: 50.0,
        },
        ParamSpec {
            key: "higgs_lambda_den",
            description: "Lambda denominator |Phi|^2 (default: 5184 = 72^2)",
            category: ParamCategory::Higgs,
            default_value: 5184.0,
            min: 1.0,
            max: 100000.0,
        },

        // ═══════════════════════════════════════════
        // Neutrino-specific
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "r4_neutrino",
            description: "Koide r^4 for neutrinos (default: 4, same as leptons)",
            category: ParamCategory::Neutrino,
            default_value: 4.0,
            min: 0.1,
            max: 50.0,
        },
        ParamSpec {
            key: "phi_neutrino_base",
            description: "Neutrino Koide phase base (default: 2/9)",
            category: ParamCategory::Neutrino,
            default_value: 0.2222222222222222,
            min: -1.0,
            max: 1.0,
        },
        ParamSpec {
            key: "phi_neutrino_shift_den",
            description: "Neutrino phase shift: pi/den (default: 12 = |W(G2)|)",
            category: ParamCategory::Neutrino,
            default_value: 12.0,
            min: 1.0,
            max: 200.0,
        },

        // ═══════════════════════════════════════════
        // Normalization
        // ═══════════════════════════════════════════
        ParamSpec {
            key: "norm_factor",
            description: "Mass formula normalization (default: 28 = dim(so(8)))",
            category: ParamCategory::Mass,
            default_value: 28.0,
            min: 1.0,
            max: 248.0,
        },
    ]
}

/// Holds user-provided parameter overrides.
/// All computation functions call `get()` to check for overrides.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct OverrideContext {
    overrides: HashMap<String, f64>,
}

impl OverrideContext {
    /// Create an empty context (uses all E8 defaults).
    pub fn new() -> Self {
        Self {
            overrides: HashMap::new(),
        }
    }

    /// Create from a map of key-value pairs.
    /// Returns error for unknown keys or out-of-bounds values.
    pub fn from_params(params: HashMap<String, f64>) -> Result<Self, String> {
        let catalog = parameter_catalog();
        let valid_keys: HashMap<&str, &ParamSpec> =
            catalog.iter().map(|p| (p.key, p)).collect();

        for (key, value) in &params {
            let spec = valid_keys
                .get(key.as_str())
                .ok_or_else(|| format!("Unknown parameter: {key}"))?;

            if !value.is_finite() {
                return Err(format!("Parameter {key} must be finite, got {value}"));
            }
            if *value < spec.min || *value > spec.max {
                return Err(format!(
                    "Parameter {key} = {value} out of range [{}, {}]",
                    spec.min, spec.max
                ));
            }
        }

        Ok(Self { overrides: params })
    }

    /// Get a parameter value: override if present, otherwise default.
    pub fn get(&self, key: &str, default: f64) -> f64 {
        self.overrides.get(key).copied().unwrap_or(default)
    }

    /// Check if any overrides are set.
    pub fn is_empty(&self) -> bool {
        self.overrides.is_empty()
    }

    /// Number of overrides.
    pub fn len(&self) -> usize {
        self.overrides.len()
    }

    /// Get all overrides.
    pub fn overrides(&self) -> &HashMap<String, f64> {
        &self.overrides
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_context() {
        let ctx = OverrideContext::new();
        assert!(ctx.is_empty());
        assert_eq!(ctx.get("a_lepton", 9.0), 9.0);
    }

    #[test]
    fn test_override() {
        let mut params = HashMap::new();
        params.insert("a_lepton".to_string(), 10.0);
        let ctx = OverrideContext::from_params(params).unwrap();
        assert_eq!(ctx.get("a_lepton", 9.0), 10.0);
        assert_eq!(ctx.get("a_up", 8.0), 8.0); // not overridden
    }

    #[test]
    fn test_unknown_key() {
        let mut params = HashMap::new();
        params.insert("bogus_key".to_string(), 1.0);
        let result = OverrideContext::from_params(params);
        assert!(result.is_err());
    }

    #[test]
    fn test_out_of_bounds() {
        let mut params = HashMap::new();
        params.insert("a_lepton".to_string(), 999.0); // max is 30
        let result = OverrideContext::from_params(params);
        assert!(result.is_err());
    }

    #[test]
    fn test_nan_rejected() {
        let mut params = HashMap::new();
        params.insert("a_lepton".to_string(), f64::NAN);
        let result = OverrideContext::from_params(params);
        assert!(result.is_err());
    }

    #[test]
    fn test_catalog_size() {
        let catalog = parameter_catalog();
        assert!(catalog.len() >= 35, "Expected ~40 params, got {}", catalog.len());
    }
}
