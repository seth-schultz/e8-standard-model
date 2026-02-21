# E8 Standard Model

**49 quantities. 0 free parameters. 250-digit precision. Verify it yourself.**

```bash
git clone https://github.com/seth-schultz/e8-standard-model
cd e8-standard-model
cargo run --release -- -p 250 scorecard
```

```
E8 Standard Model Scorecard — 49 quantities, 0 free parameters
Precision: 250 digits | Backend: MPFR (rug)
════════════════════════════════════════════════════════════════════════════════════════════════════

   #  Quantity                    Predicted   Experimental      Error      Pull     Status
 ───────────────────────────────────────────────────────────────────────────────────────────────
   1  1/α(0)                     137.035999     137.035999   +0.0000%    +0.00σ ◇ DERIVED*
   2  sin²θ_W(M_Z)                 0.231216       0.231220   -0.0018%    -0.14σ ◆ DERIVED
   3  α_s(M_Z)                     0.117941       0.118000   -0.0497%    -0.07σ ◆ DERIVED
   4  sin²θ_W(GUT)                 0.375000       0.375000   +0.0000%         ─ ■ THEOREM
   5  m_e                          0.510974       0.510999   -0.0049%         † ◆ DERIVED
   6  m_μ                        105.654264     105.658376   -0.0039%         † ◆ DERIVED
   7  m_τ                            1776.9         1776.9   +0.0022%    +0.32σ ◆ DERIVED
   8  m_u                          2.207043       2.160000   +2.1779%    +0.10σ ◆ DERIVED
   9  m_c                            1264.5         1270.0   -0.4320%    -0.27σ ◆ DERIVED
  10  m_t                          172502.4       172760.0   -0.1491%    -0.86σ ◆ DERIVED
  11  m_d                          4.636000       4.670000   -0.7281%    -0.07σ ◆ DERIVED
  12  m_s                         92.543202      93.400000   -0.9173%    -0.10σ ◆ DERIVED
  13  m_b                            4139.7         4180.0   -0.9638%    -1.34σ ◆ DERIVED
  14  Σ_lep                          1883.1         1883.0   +0.0018%    +0.18σ ◆ DERIVED
  15  Σ_up                         173769.2       174032.2   -0.1511%    -0.76σ ◆ DERIVED
  16  Σ_down                         4236.9         4278.1   -0.9625%    -0.96σ ◆ DERIVED
  17  Σ_ν                         58.568072     prediction          ─         ─ ◆ DERIVED
  18  V_ud                         0.974466       0.973730   +0.0756%    +2.37σ ◆ DERIVED
  19  V_us                         0.224507       0.224300   +0.0921%    +0.41σ ◆ DERIVED
  20  V_ub                         0.003631       0.003820   -4.9558%    -0.95σ ◆ DERIVED
  21  V_cd                         0.224376       0.221000   +1.5275%    +0.84σ ◆ DERIVED
  22  V_cs                         0.973604       0.975000   -0.1432%    -0.23σ ◆ DERIVED
  23  V_cb                         0.041838       0.040800   +2.5429%    +0.74σ ◆ DERIVED
  24  V_td                         0.008479       0.008600   -1.4044%    -0.60σ ◆ DERIVED
  25  V_ts                         0.041130       0.041500   -0.8920%    -0.41σ ◆ DERIVED
  26  V_tb                         0.999118       1.014000   -1.4677%    -0.51σ ◆ DERIVED
  27  δ_CKM                       64.285714      65.500000   -1.8539%    -0.43σ ◆ DERIVED
  28  J (Jarlskog)                 2.991e-5       3.080e-5   -2.8847%    -0.59σ ◆ DERIVED
  29  sin²θ₁₂                      0.311004       0.307000   +1.3043%    +0.31σ ◆ DERIVED
  30  sin²θ₂₃                      0.535898       0.546000   -1.8501%    -0.48σ ◆ DERIVED
  31  sin²θ₁₃                      0.022329       0.022200   +0.5815%    +0.19σ ◆ DERIVED
  32  δ_PMNS                     192.857143     197.000000   -2.1030%    -0.14σ ◆ DERIVED
  33  m_ν₁                         0.374257     prediction          ─         ─ ◆ DERIVED
  34  m_ν₂                         8.699498     prediction          ─         ─ ◆ DERIVED
  35  m_ν₃                        49.494317     prediction          ─         ─ ◆ DERIVED
  36  Δm²₂₁                        7.554e-5       7.530e-5   +0.3203%    +0.13σ ◆ DERIVED
  37  Δm²₃₁                        0.002450       0.002453   -0.1408%    -0.10σ ◆ DERIVED
  38  m_H                        125.124262     125.250000   -0.1004%    -0.74σ ◆ DERIVED
  39  λ_H                          0.131532       0.131500   +0.0246%    +0.03σ ◆ DERIVED
  40  m_S                         95.565234      95.400000   +0.1732%    +0.08σ ◆ DERIVED
  41  θ̄_QCD                              0              0   +0.0000%         ─ ■ THEOREM
  42  |Φ_E8|                     240.000000     prediction          ─         ─ ■ THEOREM
  43  Plaquettes                     2240.0     prediction          ─         ─ ■ THEOREM
  44  Per root                    28.000000     prediction          ─         ─ ■ THEOREM
  45  Generations                  3.000000       3.000000   +0.0000%         ─ ■ THEOREM
  46  d = 8                        8.000000     prediction          ─         ─ ■ THEOREM
  47  y_t                          0.990802       0.991000   -0.0200%         ─ ◆ DERIVED
  48  M_GUT                        4.578e18     prediction          ─         ─ ◆ DERIVED
  49  Confinement                  1.000000       1.000000   +0.0000%         ─ ■ THEOREM

════════════════════════════════════════════════════════════════════════════════════════════════════
Summary: 8 ■ THEOREM | 40 ◆ DERIVED | 1 ◇ DERIVED*
         31/33 within 1σ (93%) | 32/33 within 2σ (96%)
         2 entries at precision floor (†): theory ~0.005% vs experiment ~ppb
```

This program derives all 49 experimentally measured Standard Model parameters — fermion masses, mixing angles, gauge couplings, the Higgs mass, and neutrino mass splittings — from a single mathematical structure: the E8 root lattice equipped with its Epstein zeta coupling. There are zero free parameters. Every numerical constant traces to a Lie algebra invariant. Every prediction is reproducible at arbitrary precision.

The accompanying paper provides the complete 16-theorem proof chain:

> *The E8 Standard Model: 49 Parameters from a Single Axiom*
> [arXiv:XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX) | [Zenodo DOI](https://doi.org/XX.XXXX/zenodo.XXXXXXX)

## What this is

This is the computational companion to a theoretical physics paper. It implements every calculation in the paper's proof chain as strictly typed, 250-digit-precision Rust code. It is designed for one purpose: **to make the framework falsifiable by removing computational ambiguity**. If you dispute a prediction, argue with the math — not the numerics.

### What is computed

| Category | Quantities | Agreement with experiment |
|---|---|---|
| Gauge couplings | 1/alpha, sin^2(theta_W), alpha_s | 0.001 ppb, -0.14 sigma, -0.07 sigma |
| Charged lepton masses | e, mu, tau | 0.005% (precision floor) |
| Quark masses | u, c, t, d, s, b | 0.1% - 2.2% |
| CKM matrix | 9 elements + Jarlskog J + delta_CKM | Max pull 2.4 sigma |
| PMNS mixing | 3 angles + delta_PMNS | Max pull 0.48 sigma |
| Neutrino masses | m_1, m_2, m_3, Delta m^2_21, Delta m^2_31 | 0.13 sigma, -0.10 sigma |
| Higgs sector | m_H, lambda_H, theta_QCD = 0 | -0.74 sigma |
| Second scalar | m_S = 95.6 GeV | Testable prediction |

31/33 predictions with experimental data fall within 1 sigma (93%).

### What is *not* computed

Nothing is fitted. No parameters are adjusted to match experiment. The only physical input is the Planck mass m_P = 1.220890(14) x 10^19 GeV (CODATA 2022). The electroweak scale M_Z = 91.1876 GeV enters as a reference point for RGE running of alpha_s.

## Installation

### Prerequisites

This library uses [rug](https://crates.io/crates/rug) for arbitrary-precision arithmetic (GMP/MPFR), which requires C libraries.

**macOS:**
```bash
brew install gmp mpfr libmpc
```

**Ubuntu/Debian:**
```bash
apt install libgmp-dev libmpfr-dev libmpc-dev
```

**Fedora:**
```bash
dnf install gmp-devel mpfr-devel libmpc-devel
```

### Build and run

```bash
git clone https://github.com/seth-schultz/e8-standard-model
cd e8-standard-model
cargo build --release
./target/release/e8-standard-model scorecard
```

### As a library

Add to your `Cargo.toml`:
```toml
[dependencies]
e8-core = { git = "https://github.com/seth-schultz/e8-standard-model" }
```

## Usage

### Library API

```rust
use e8_core::precision::set_precision;
use e8_core::scorecard::table::compute_scorecard;

// Set working precision (decimal digits)
set_precision(250);

// Compute all 49 quantities
let predictions = compute_scorecard(250);
for p in &predictions {
    println!("{}: {:.6}", p.name, p.predicted);
}
```

#### Individual computations

```rust
use e8_core::precision::set_precision;

set_precision(250);

// Fine structure constant from CF tower: [244; 14, 13, 193] * e^{-gamma}
let alpha_inv = e8_core::coupling::alpha::alpha_inverse();
println!("1/alpha = {:.12}", alpha_inv.to_f64());

// All 12 fermion masses from E8 mass formula + Koide parametrization
let masses = e8_core::mass::sectors::compute_all_masses();
println!("m_e = {:.12} MeV", masses.electron.to_f64());
println!("m_t = {:.4} MeV", masses.top.to_f64());

// CKM matrix from Fritzsch texture + octonionic CP phase
let ckm = e8_core::mixing::ckm::build_ckm(&masses);
println!("V_us = {:.6}", ckm.magnitudes[1].to_f64());

// PMNS mixing angles from G2 Coxeter geometry
let s13 = e8_core::mixing::pmns::sin2_theta13();
println!("sin^2(theta_13) = {:.6}", s13.to_f64());

// Higgs mass from lambda = 7 pi^4 / 72^2
let m_h = e8_core::higgs::mass::higgs_mass_default();
println!("m_H = {:.4} GeV", m_h.to_f64());
```

### CLI

```bash
# Full scorecard (all 49 quantities)
e8-standard-model scorecard

# At maximum precision
e8-standard-model -p 250 scorecard

# JSON output (for programmatic use)
e8-standard-model scorecard --output json

# Individual sectors
e8-standard-model sector leptons
e8-standard-model sector ckm
e8-standard-model sector gauge
e8-standard-model sector higgs

# E8 lattice geometry (240 roots, 2240 plaquettes, associator)
e8-standard-model roots

# System info
e8-standard-model info
```

## Architecture

```
e8-core/src/
  precision/    Arbitrary-precision arithmetic via MPFR (250+ digits, configurable)
  lattice/      E8 root system: 240 roots, traces, plaquettes, D4 coset decomposition
  algebra/      Lie algebra invariants: E8, G2, SU(5), SU(3), SU(2), embedding indices
  octonion/     Fano plane multiplication, associator [e_a,e_b,e_c], generation assignment
  mass/         Unified mass formula + Koide parametrization (12 fermion + 3 neutrino masses)
  coupling/     Fine structure constant, Weinberg angle, alpha_s, RGE running
  mixing/       CKM (Fritzsch texture + octonionic CP) and PMNS (G2 Coxeter geometry)
  higgs/        Quartic coupling, Higgs mass, theta_QCD = 0
  special/      Continued fractions, divisor sums, Eisenstein series
  scorecard/    Master pipeline: all 49 quantities + PDG 2024 comparison
```

### No magic numbers

Every numerical constant traces back to a Lie algebra invariant or group-theoretic identity:

- `E8.num_roots` (240) instead of `240`
- `G2.coxeter_number` (6) instead of `6`
- `SU3.c2_fundamental` (4/3) instead of `4.0/3.0`
- `DIM_IM_OCTONIONS` (7) instead of `7`

The derivation chain is auditable directly from the source code.

### Precision

Arithmetic uses [rug](https://crates.io/crates/rug) (Rust bindings to GNU MPFR) for arbitrary-precision floating point. Precision is configurable at runtime:

```rust
use e8_core::precision::set_precision;
set_precision(50);   // Fast, sufficient for all physical predictions
set_precision(250);  // Match the reference derivation scripts
set_precision(1000); // For mathematical verification
```

The default is 250 decimal digits (~832 bits), matching the precision of the original derivation scripts. All intermediate quantities are computed at the requested precision — there are no f64 approximations in the computational core.

### Performance

Compiled with `opt-level = 3`, full LTO, and `target-cpu=native`. The full 49-quantity scorecard at 250-digit precision completes in under 10 milliseconds.

## The physics

The framework rests on two axioms:

1. **The E8 root lattice** is the fundamental structure.
2. **The Epstein zeta function** provides the coupling between lattice geometry and physics.

From these, four mechanisms generate all Standard Model parameters:

| Mechanism | What it determines | How |
|---|---|---|
| **Algebraic** | Gauge couplings, mixing structure | Lie algebra embedding indices, trace identities |
| **Lattice** | Fermion masses | E8 theta function: Sigma = f * m_P * exp(-(A*R + delta)/28) |
| **Octonionic** | CKM CP phase, generations | Fano plane non-associativity: delta_CKM = 5pi/14 |
| **Gauge flow** | alpha_s, confinement | RGE running from E8-derived GUT scale |

### Key formulas

| Quantity | Formula | Origin |
|---|---|---|
| 1/alpha | `[244; 14, 13, 193]` * e^{-gamma} | E8 -> D4 -> G2 subgroup chain |
| sin^2(theta_W) at M_Z | (3/13)(1 + 5alpha/(6pi)) | Trace doubling + G2 Coxeter |
| Sector mass sum | f * m_P * exp(-(A*R + delta)/28) | E8 lattice propagator |
| Individual masses | sqrt(m_k) = M(1 + r cos(2pi*k/3 + phi)) | SU(5) Yukawa representations |
| delta_CKM | 5pi/14 | Octonionic associator [e_6, e_3, e_1] |
| PMNS angles | G2 exponents m_1=1, m_2=5 | G2 Coxeter geometry |
| lambda_H | 7pi^4/72^2 | dim(Im(O)) * Res(Z_E8, 4) / (\|W(G2)\| * h(G2)^2) |
| D4 Euler product | Z_mix(s) = 192 * 2^{-s}(1-2^{-s})zeta(s)zeta(s-3) | E8 -> D4_L x D4_R decomposition |
| Color Casimir | Res(Z_E8)/Res(Z_mix) = 4/3 | Residue ratio at s=4 |

### Testable predictions

These predictions have no experimental measurement yet and provide falsifiable tests:

| Prediction | Value | How to test |
|---|---|---|
| Neutrino mass sum | Sigma_nu = 58.6 meV | KATRIN/TRISTAN, cosmological surveys (DESI) |
| Second scalar boson | m_S = 95.6 GeV | HL-LHC Run 3 (CMS/ATLAS diphoton) |
| Lightest neutrino mass | m_1 = 0.37 meV | Next-generation beta decay experiments |
| PMNS CP phase | delta_PMNS = 192.9 deg | DUNE, Hyper-Kamiokande |

## Testing

```bash
cargo test          # 74 tests: all mathematical identities, predictions, cross-checks
cargo clippy        # Zero warnings
cargo doc --no-deps # Zero warnings
```

74 unit tests verify every mathematical identity, mass prediction, mixing matrix element, and scorecard value.

## Citation

If you use this software in your research, please cite both the paper and the code:

```bibtex
@software{e8_standard_model_code,
  author    = {Schultz, Seth},
  title     = {{E8 Standard Model}: {Rust} implementation of 49 {SM} parameters from 1 axiom},
  year      = {2026},
  url       = {https://github.com/seth-schultz/e8-standard-model},
  doi       = {XX.XXXX/zenodo.XXXXXXX},
  version   = {1.0.0},
}
```

```bibtex
@article{e8_standard_model_paper,
  author    = {Schultz, Seth},
  title     = {The {E8} Standard Model: 49 Parameters from a Single Axiom},
  year      = {2026},
  eprint    = {XXXX.XXXXX},
  archivePrefix = {arXiv},
  primaryClass  = {hep-th},
}
```

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
