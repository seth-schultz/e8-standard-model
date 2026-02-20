# E8 Standard Model

A Rust implementation of the E8 lattice framework for computing Standard Model parameters from first principles.

**49 quantities. 0 free parameters. 1 axiom.**

This library derives all experimentally measured Standard Model parameters — fermion masses, mixing angles, gauge couplings, Higgs mass — from a single mathematical structure: the E8 root lattice equipped with its Epstein zeta coupling.

## Results at a Glance

| Category | Quantities | Agreement |
|---|---|---|
| Gauge couplings | 1/alpha, sin^2 theta_W, alpha_s | 0.001 ppb, 0.10 sigma, -1.1 sigma |
| Charged lepton masses | e, mu, tau | 0.005% or better |
| Quark masses | u, c, t, d, s, b | Within 1-2% |
| CKM matrix | 9 elements + Jarlskog J | chi^2 = 0.001 |
| PMNS mixing | 3 angles + delta_PMNS | Max pull 0.48 sigma |
| Neutrino masses | 3 masses + 2 Delta m^2 | 0.13 sigma, -0.10 sigma |
| Higgs | m_H, lambda, theta_QCD = 0 | +0.36 sigma |

36/41 predictions with experimental data fall within 1 sigma (88%).

## Installation

### Prerequisites

This library uses [rug](https://crates.io/crates/rug) for arbitrary-precision arithmetic, which requires GMP, MPFR, and MPC C libraries.

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

### As a library

Add to your `Cargo.toml`:
```toml
[dependencies]
e8-core = { git = "https://github.com/seth-schultz/e8-standard-model" }
```

### As a CLI tool

```bash
git clone https://github.com/seth-schultz/e8-standard-model
cd e8-standard-model
cargo build --release
./target/release/e8-standard-model scorecard
```

## Usage

### Library API

```rust
use e8_core::precision::set_precision;
use e8_core::scorecard::table::compute_scorecard;

// Set working precision (decimal digits)
set_precision(50);

// Compute all 49 quantities
let predictions = compute_scorecard(50);
for p in &predictions {
    println!("{}: {:.6}", p.name, p.predicted);
}
```

#### Individual computations

```rust
use e8_core::precision::set_precision;

set_precision(250); // 250-digit precision

// Fermion masses from E8 mass formula + Koide parametrization
let masses = e8_core::mass::sectors::compute_all_masses();
println!("m_e  = {:.12} MeV", masses.electron.to_f64());
println!("m_t  = {:.4} MeV",  masses.top.to_f64());

// Fine structure constant from CF tower
let alpha_inv = e8_core::coupling::alpha::alpha_inverse();
println!("1/alpha = {:.12}", alpha_inv.to_f64());

// CKM matrix from Fritzsch texture + octonionic CP phase
let ckm = e8_core::mixing::ckm::build_ckm(&masses);
println!("V_us = {:.6}", ckm.magnitudes[1].to_f64());

// PMNS mixing angles from G2 Coxeter geometry
let s13 = e8_core::mixing::pmns::sin2_theta13();
println!("sin^2 theta_13 = {:.6}", s13.to_f64());

// Higgs mass from lambda = 7 pi^4 / 72^2
let m_h = e8_core::higgs::mass::higgs_mass_default();
println!("m_H = {:.4} GeV", m_h.to_f64());
```

### CLI

```bash
# Full scorecard (all 49 quantities)
e8-standard-model scorecard

# JSON output
e8-standard-model scorecard --output json

# Individual sectors
e8-standard-model sector leptons
e8-standard-model sector ckm
e8-standard-model sector gauge

# Higher precision
e8-standard-model scorecard --precision 250

# E8 lattice properties
e8-standard-model roots

# System info
e8-standard-model info
```

## Architecture

```
e8-core/src/
  precision/    Adaptive arbitrary-precision via MPFR (thread-local, configurable)
  lattice/      E8 root system: 240 roots, traces, plaquettes, Cartan matrix
  algebra/      Lie algebra invariants: E8, G2, SU(5), SU(3), SU(2), etc.
  octonion/     Fano plane multiplication, associator, generation assignment
  mass/         Mass formula + Koide parametrization (12 fermion masses)
  coupling/     Fine structure constant, Weinberg angle, alpha_s RGE running
  mixing/       CKM (Fritzsch + octonionic CP) and PMNS (G2 Coxeter)
  higgs/        Quartic coupling, Higgs mass, theta_QCD = 0
  special/      Continued fractions, divisor sums, zeta functions
  scorecard/    Master pipeline: all 49 quantities + PDG comparison
```

### No magic numbers

Every numerical constant in the code traces back to a Lie algebra invariant or group-theoretic identity. For example:

- `E8.num_roots` (240) instead of `240`
- `G2.coxeter_number` (6) instead of `6`
- `SU3.c2_fundamental` (4/3) instead of `4.0/3.0`
- `DIM_IM_OCTONIONS` (7) instead of `7`

This makes the physics derivation chain auditable directly from the source code.

### Precision

Arithmetic uses [rug](https://crates.io/crates/rug) (Rust bindings to GNU MPFR) for arbitrary-precision floating point. Precision is configurable at runtime:

```rust
use e8_core::precision::set_precision;
set_precision(50);   // Fast, sufficient for all physical predictions
set_precision(250);  // Match the Python reference scripts
set_precision(1000); // For mathematical verification
```

The default is 250 decimal digits (~832 bits), matching the precision of the original Python derivation scripts.

### Performance

Compiled with `opt-level = 3`, full LTO, and `target-cpu=native` to exploit SIMD (AVX2/AVX-512 on x86_64, NEON on aarch64). The full 49-quantity scorecard at 50-digit precision completes in approximately 10ms.

## The Physics

The framework rests on a single axiom: **the E8 root lattice is the fundamental structure**, with the Epstein zeta function providing the coupling between lattice geometry and physics.

### Four mechanisms

1. **Algebraic** — Lie algebra embedding indices and trace identities determine gauge couplings and mixing structure
2. **Lattice** — The mass formula Sigma = f * m_P * exp(-(A*R + delta)/28) with R = 240 * e^{-gamma} arises from the E8 theta function
3. **Octonionic** — The Fano plane and octonionic non-associativity generate the CKM CP phase (delta = 5 pi / 14) and generation structure
4. **Gauge flow** — Renormalization group running from the E8-derived GUT scale determines alpha_s(M_Z)

### Key formulas

| Quantity | Formula | Origin |
|---|---|---|
| 1/alpha | [244; 14, 13, 193] * e^{-gamma} | E8 -> D4 -> G2 subgroup chain |
| sin^2 theta_W(M_Z) | (3/13)(1 + 5 alpha/(6 pi)) | Trace doubling + G2 Coxeter |
| Sigma_sector | f * m_P * exp(-(A*R + delta)/28) | E8 lattice propagator |
| Koide masses | sqrt(m_k) = M(1 + r cos(2 pi k/3 + phi)) | SU(5) Yukawa representations |
| delta_CKM | 5 pi / 14 | Octonionic associator [e6, e3, e1] |
| PMNS angles | G2 exponents and Weyl order | G2 Coxeter geometry |
| lambda_H | 7 pi^4 / 72^2 | dim(Im(O)) * Res(Z_E8, 4) / (\|W(G2)\| * h(G2)^2) |

### Physical constants

The only physical input is the Planck mass m_P = 1.220890(14) * 10^19 GeV (CODATA 2022). The electroweak scale M_Z = 91.1876 GeV enters as an experimental reference point for RGE running of alpha_s.

## Testing

```bash
# Run all tests
cargo test

# Run with higher precision
RUST_TEST_THREADS=1 cargo test -- --nocapture
```

63 unit tests verify all mathematical identities, mass predictions, CKM elements, and scorecard values.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{e8_standard_model,
  title = {E8 Standard Model: 49 SM Parameters from 1 Axiom},
  url = {https://github.com/seth-schultz/e8-standard-model},
  year = {2025}
}
```

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
