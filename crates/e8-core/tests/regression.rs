//! Immovable regression gate: all 49 E8 Standard Model predictions.
//!
//! These values are computed at 50-digit precision and must remain
//! bit-identical through any refactoring. If this test fails, the
//! refactoring introduced a numerical regression.

use e8_core::scorecard::table::compute_scorecard;

/// Relative tolerance: 1e-12 (well below the precision floor).
const REL_TOL: f64 = 1e-12;

fn check(name: &str, got: f64, expected: f64) {
    if expected == 0.0 {
        assert!(
            got.abs() < 1e-15,
            "{name}: expected 0.0, got {got:.15e}"
        );
        return;
    }
    let rel = ((got - expected) / expected).abs();
    assert!(
        rel < REL_TOL,
        "{name}: relative error {rel:.3e} exceeds {REL_TOL:.0e}\n  got:      {got:.15e}\n  expected: {expected:.15e}"
    );
}

#[test]
#[cfg(feature = "arbitrary-precision")]
fn regression_49_predictions() {
    let predictions = compute_scorecard(50);
    assert_eq!(predictions.len(), 49, "Expected 49 predictions, got {}", predictions.len());

    // Hardcoded reference values (50-digit precision, captured 2026-02-23)
    // Updated: cross-rep vertex correction (952/893)^{3/20} applied to Σ_down
    let expected: &[(&str, f64)] = &[
        // #1: 1/α(0)
        ("1/α(0)", 1.370359991770789e2),
        // #2: sin²θ_W(M_Z)
        ("sin²θ_W(M_Z)", 2.312159268200373e-1),
        // #3: α_s(M_Z)
        ("α_s(M_Z)", 1.179413093345547e-1),
        // #4: sin²θ_W(GUT)
        ("sin²θ_W(GUT)", 3.750000000000000e-1),
        // #5: m_e
        ("m_e", 5.109740402604934e-1),
        // #6: m_μ
        ("m_μ", 1.056542637585243e2),
        // #7: m_τ
        ("m_τ", 1.776898345471456e3),
        // #8: m_u
        ("m_u", 2.207042894153784e0),
        // #9: m_c
        ("m_c", 1.264513385366838e3),
        // #10: m_t
        ("m_t", 1.725024484437401e5),
        // #11: m_d (vertex-corrected)
        ("m_d", 4.680704806077117e0),
        // #12: m_s (vertex-corrected)
        ("m_s", 9.343559248214324e1),
        // #13: m_b (vertex-corrected)
        ("m_b", 4.179632975423529e3),
        // #14: Σ_lep
        ("Σ_lep", 1.883063583270241e3),
        // #15: Σ_up
        ("Σ_up", 1.737691688720011e5),
        // #16: Σ_down (vertex-corrected: × (952/893)^{3/20})
        ("Σ_down", 4.277749272711750e3),
        // #17: Σ_ν
        ("Σ_ν", 5.856807214039605e1),
        // #18: V_ud (updated: down masses changed)
        ("V_ud", 9.745609132002185e-1),
        // #19: V_us
        ("V_us", 2.240933310545309e-1),
        // #20: V_ub
        ("V_ub", 3.633928898657079e-3),
        // #21: V_cd
        ("V_cd", 2.239645920971494e-1),
        // #22: V_cs
        ("V_cs", 9.737064811068042e-1),
        // #23: V_cb
        ("V_cb", 4.165993443780532e-2),
        // #24: V_td
        ("V_td", 8.419498150933135e-3),
        // #25: V_ts
        ("V_ts", 4.096178252333712e-2),
        // #26: V_tb
        ("V_tb", 9.991252396088276e-1),
        // #27: δ_CKM
        ("δ_CKM", 6.428571428571429e1),
        // #28: J (Jarlskog)
        ("J (Jarlskog)", 2.973240224236010e-5),
        // #29: sin²θ₁₂
        ("sin²θ₁₂", 3.110042339640731e-1),
        // #30: sin²θ₂₃
        ("sin²θ₂₃", 5.358983848622454e-1),
        // #31: sin²θ₁₃
        ("sin²θ₁₃", 2.232909936926022e-2),
        // #32: δ_PMNS
        ("δ_PMNS", 1.928571428571429e2),
        // #33: m_ν₁
        ("m_ν₁", 3.742568708979613e-1),
        // #34: m_ν₂
        ("m_ν₂", 8.699498019528253e0),
        // #35: m_ν₃
        ("m_ν₃", 4.949431724996984e1),
        // #36: Δm²₂₁
        ("Δm²₂₁", 7.554119758636166e-5),
        // #37: Δm²₃₁
        ("Δm²₃₁", 2.449547371835247e-3),
        // #38: m_H
        ("m_H", 1.251242618125252e2),
        // #39: λ_H
        ("λ_H", 1.315323374301730e-1),
        // #40: m_S
        ("m_S", 9.556523350522819e1),
        // #41: θ̄_QCD
        ("θ\u{304}_QCD", 0.000000000000000e0),
        // #42: |Φ_E8|
        ("|Φ_E8|", 2.400000000000000e2),
        // #43: Plaquettes
        ("Plaquettes", 2.240000000000000e3),
        // #44: Per root
        ("Per root", 2.800000000000000e1),
        // #45: Generations
        ("Generations", 3.000000000000000e0),
        // #46: d = 8
        ("d = 8", 8.000000000000000e0),
        // #47: y_t
        ("y_t", 9.908021368357685e-1),
        // #48: M_GUT
        ("M_GUT", 4.578337500000000e18),
        // #49: Confinement
        ("Confinement", 1.000000000000000e0),
    ];

    for (i, (exp_name, exp_val)) in expected.iter().enumerate() {
        let pred = &predictions[i];
        assert_eq!(
            pred.name, *exp_name,
            "Prediction #{} name mismatch: got {:?}, expected {:?}",
            i + 1, pred.name, exp_name
        );
        check(&format!("#{} {}", i + 1, exp_name), pred.predicted, *exp_val);
    }
}
