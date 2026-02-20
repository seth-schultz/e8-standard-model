//! PDG 2024 experimental values for comparison.

/// PDG experimental value with uncertainty.
#[derive(Debug, Clone, Copy)]
pub struct PDGValue {
    pub name: &'static str,
    pub value: f64,
    pub uncertainty: f64,
    pub unit: &'static str,
}

/// All PDG values used in the scorecard.
pub const PDG_VALUES: &[PDGValue] = &[
    // Gauge couplings
    PDGValue { name: "1/alpha", value: 137.035999177, uncertainty: 0.000000021, unit: "" },
    PDGValue { name: "sin2_theta_w_mz", value: 0.23122, uncertainty: 0.00003, unit: "" },
    PDGValue { name: "alpha_s_mz", value: 0.1180, uncertainty: 0.0009, unit: "" },

    // Charged lepton masses (MeV)
    PDGValue { name: "m_e", value: 0.51099895, uncertainty: 0.00000003, unit: "MeV" },
    PDGValue { name: "m_mu", value: 105.6583755, uncertainty: 0.0000023, unit: "MeV" },
    PDGValue { name: "m_tau", value: 1776.86, uncertainty: 0.12, unit: "MeV" },

    // Up-type quark masses (MeV)
    PDGValue { name: "m_u", value: 2.16, uncertainty: 0.49, unit: "MeV" },
    PDGValue { name: "m_c", value: 1270.0, uncertainty: 20.0, unit: "MeV" },
    PDGValue { name: "m_t", value: 172760.0, uncertainty: 300.0, unit: "MeV" },

    // Down-type quark masses (MeV)
    PDGValue { name: "m_d", value: 4.67, uncertainty: 0.48, unit: "MeV" },
    PDGValue { name: "m_s", value: 93.4, uncertainty: 8.6, unit: "MeV" },
    PDGValue { name: "m_b", value: 4180.0, uncertainty: 30.0, unit: "MeV" },

    // CKM matrix elements
    PDGValue { name: "V_ud", value: 0.97373, uncertainty: 0.00031, unit: "" },
    PDGValue { name: "V_us", value: 0.2243, uncertainty: 0.0008, unit: "" },
    PDGValue { name: "V_ub", value: 0.00382, uncertainty: 0.00020, unit: "" },
    PDGValue { name: "V_cd", value: 0.221, uncertainty: 0.004, unit: "" },
    PDGValue { name: "V_cs", value: 0.975, uncertainty: 0.006, unit: "" },
    PDGValue { name: "V_cb", value: 0.0408, uncertainty: 0.0014, unit: "" },
    PDGValue { name: "V_td", value: 0.0086, uncertainty: 0.0002, unit: "" },
    PDGValue { name: "V_ts", value: 0.0415, uncertainty: 0.0009, unit: "" },
    PDGValue { name: "V_tb", value: 1.014, uncertainty: 0.029, unit: "" },

    // Jarlskog invariant
    PDGValue { name: "J", value: 3.08e-5, uncertainty: 0.15e-5, unit: "" },

    // CP phases (degrees)
    PDGValue { name: "delta_CKM", value: 65.5, uncertainty: 2.8, unit: "deg" },
    PDGValue { name: "delta_PMNS", value: 197.0, uncertainty: 30.0, unit: "deg" },

    // Neutrino mass-squared differences
    PDGValue { name: "Dm21_sq", value: 7.53e-5, uncertainty: 0.18e-5, unit: "eV^2" },
    PDGValue { name: "Dm31_sq", value: 2.453e-3, uncertainty: 0.033e-3, unit: "eV^2" },

    // PMNS mixing angles
    PDGValue { name: "sin2_theta13_pmns", value: 0.02220, uncertainty: 0.00068, unit: "" },
    PDGValue { name: "sin2_theta12_pmns", value: 0.307, uncertainty: 0.013, unit: "" },
    PDGValue { name: "sin2_theta23_pmns", value: 0.546, uncertainty: 0.021, unit: "" },

    // Higgs
    PDGValue { name: "m_H", value: 125.25, uncertainty: 0.17, unit: "GeV" },
    PDGValue { name: "lambda_H", value: 0.1315, uncertainty: 0.001, unit: "" },
];
