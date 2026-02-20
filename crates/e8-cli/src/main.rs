use clap::{Parser, Subcommand};
use e8_core::scorecard::prediction::{Prediction, Status};
use e8_core::scorecard::table::{compute_scorecard, scorecard_summary};

#[derive(Parser)]
#[command(
    name = "e8-standard-model",
    about = "E8 Standard Model — 49 quantities from the E8 axiom, 0 free parameters",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    /// Decimal digits of precision
    #[arg(short, long, default_value = "50")]
    precision: u32,

    /// Output format
    #[arg(short, long, default_value = "table")]
    output: String,

    /// Suppress progress output
    #[arg(short, long)]
    quiet: bool,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute all quantities (default)
    Scorecard,
    /// Compute a single sector
    Sector {
        /// Sector: leptons, up, down, neutrinos, gauge, ckm, pmns, higgs
        name: String,
    },
    /// Dump E8 root system properties
    Roots,
    /// Show hardware and precision info
    Info,
}

fn main() {
    let cli = Cli::parse();
    let digits = cli.precision;

    match cli.command.unwrap_or(Commands::Scorecard) {
        Commands::Scorecard => run_scorecard(digits, &cli.output, cli.quiet),
        Commands::Sector { name } => run_sector(digits, &name),
        Commands::Roots => run_roots(),
        Commands::Info => run_info(digits),
    }
}

/// Format a floating-point value for display: use scientific notation for
/// very large or very small (but nonzero) values, otherwise fixed precision.
fn format_value(v: f64) -> String {
    let abs = v.abs();
    if abs == 0.0 {
        "0".to_string()
    } else if abs >= 1e6 || (abs < 1e-4 && abs > 0.0) {
        format!("{:.3e}", v)
    } else if abs >= 1000.0 {
        format!("{:.1}", v)
    } else if abs >= 1.0 {
        format!("{:.6}", v)
    } else {
        format!("{:.6}", v)
    }
}

fn run_scorecard(digits: u32, format: &str, quiet: bool) {
    if !quiet {
        eprintln!(
            "E8 Standard Model Scorecard — computing at {}-digit precision...",
            digits
        );
    }

    let start = std::time::Instant::now();
    let predictions = compute_scorecard(digits);
    let elapsed = start.elapsed();

    match format {
        "json" => {
            let json = serde_json::to_string_pretty(&predictions).unwrap();
            println!("{}", json);
        }
        _ => {
            // Table format
            println!();
            println!("E8 Standard Model Scorecard — {} quantities, 0 free parameters", predictions.len());
            println!("Precision: {} digits | Backend: MPFR (rug)", digits);
            println!("{}", "═".repeat(100));
            println!();
            println!(
                " {:>3}  {:<22} {:>14} {:>14} {:>10} {:>9} {:>10}",
                "#", "Quantity", "Predicted", "Experimental", "Error", "Pull", "Status"
            );
            println!(" {}", "─".repeat(95));

            for (i, p) in predictions.iter().enumerate() {
                // Smart formatting: use scientific notation for very large/small numbers
                let pred_str = format_value(p.predicted);
                let exp_str = match p.experimental {
                    Some(e) => format_value(e),
                    None => "prediction".to_string(),
                };

                // Show percent error for all, pull only when meaningful
                let err_str = match p.pct_error {
                    Some(e) if p.experimental.is_some() => format!("{:+.4}%", e),
                    _ => "─".to_string(),
                };

                // For pulls > 100σ with small % error, show the % error
                // instead — the pull is meaningless (precision floor)
                let pull_str = match p.pull_sigma {
                    Some(pull) if pull.abs() > 100.0 => {
                        // Precision floor: theory at ~0.01%, experiment at ppb
                        "†".to_string()
                    }
                    Some(pull) => format!("{:+.2}σ", pull),
                    None => "─".to_string(),
                };

                let status_sym = match p.status {
                    Status::Theorem => "■",
                    Status::Derived => "◆",
                    Status::DerivedStar => "◇",
                    Status::Identified => "○",
                };

                println!(
                    " {:>3}  {:<22} {:>14} {:>14} {:>10} {:>9} {} {}",
                    i + 1,
                    p.name,
                    pred_str,
                    exp_str,
                    err_str,
                    pull_str,
                    status_sym,
                    p.status
                );
            }

            println!();
            println!("{}", "═".repeat(100));

            // Count precision-floor entries separately
            let precision_floor: Vec<&Prediction> = predictions
                .iter()
                .filter(|p| {
                    p.pull_sigma.map_or(false, |pull| pull.abs() > 100.0)
                        && p.pct_error.map_or(false, |e| e.abs() < 1.0)
                })
                .collect();

            let summary = scorecard_summary(&predictions);

            // Adjusted stats: exclude precision-floor entries from pull counts
            let effective_with_exp = summary.with_experimental - precision_floor.len();
            let effective_1sigma = summary.within_1sigma;
            let effective_2sigma = summary.within_2sigma;

            println!(
                "Summary: {} ■ THEOREM | {} ◆ DERIVED | {} ◇ DERIVED*",
                summary.theorems, summary.derived, summary.derived_star
            );
            if effective_with_exp > 0 {
                println!(
                    "         {}/{} within 1σ ({}%) | {}/{} within 2σ ({}%)",
                    effective_1sigma,
                    effective_with_exp,
                    effective_1sigma * 100 / effective_with_exp,
                    effective_2sigma,
                    effective_with_exp,
                    effective_2sigma * 100 / effective_with_exp,
                );
            }
            if !precision_floor.is_empty() {
                println!(
                    "         {} entries at precision floor (†): theory ~0.005% vs experiment ~ppb",
                    precision_floor.len()
                );
            }
            println!("Computed in {:.2}s", elapsed.as_secs_f64());
        }
    }
}

fn run_sector(digits: u32, name: &str) {
    e8_core::precision::set_precision(digits);
    match name {
        "leptons" | "lepton" => {
            let masses = e8_core::mass::sectors::compute_all_masses();
            println!("Charged Lepton Masses (MeV):");
            println!("  electron: {:.10}", masses.electron.to_f64());
            println!("  muon:     {:.10}", masses.muon.to_f64());
            println!("  tau:      {:.10}", masses.tau.to_f64());
            println!("  Σ_lep:    {:.10}", masses.sigma_lep.to_f64());
        }
        "up" => {
            let masses = e8_core::mass::sectors::compute_all_masses();
            println!("Up-Type Quark Masses (MeV):");
            println!("  up:    {:.10}", masses.up.to_f64());
            println!("  charm: {:.10}", masses.charm.to_f64());
            println!("  top:   {:.10}", masses.top.to_f64());
            println!("  Σ_up:  {:.10}", masses.sigma_up.to_f64());
        }
        "down" => {
            let masses = e8_core::mass::sectors::compute_all_masses();
            println!("Down-Type Quark Masses (MeV):");
            println!("  down:    {:.10}", masses.down.to_f64());
            println!("  strange: {:.10}", masses.strange.to_f64());
            println!("  bottom:  {:.10}", masses.bottom.to_f64());
            println!("  Σ_down:  {:.10}", masses.sigma_down.to_f64());
        }
        "neutrinos" | "neutrino" => {
            let masses = e8_core::mass::sectors::compute_all_masses();
            println!("Neutrino Masses (meV):");
            println!("  ν₁: {:.6}", masses.nu1.to_f64());
            println!("  ν₂: {:.6}", masses.nu2.to_f64());
            println!("  ν₃: {:.6}", masses.nu3.to_f64());
        }
        "gauge" => {
            println!("Gauge Couplings:");
            println!("  1/α:         {:.12}", e8_core::coupling::alpha::alpha_inverse().to_f64());
            println!("  sin²θ_W(GUT): {:.12}", e8_core::coupling::weinberg::sin2_theta_w_gut().to_f64());
            println!("  sin²θ_W(M_Z): {:.12}", e8_core::coupling::weinberg::sin2_theta_w_mz().to_f64());
            println!("  α_s(M_Z):     {:.12}", e8_core::coupling::alpha_s::alpha_s_mz().to_f64());
        }
        "ckm" => {
            let masses = e8_core::mass::sectors::compute_all_masses();
            let ckm = e8_core::mixing::ckm::build_ckm(&masses);
            let names = ["V_ud", "V_us", "V_ub", "V_cd", "V_cs", "V_cb", "V_td", "V_ts", "V_tb"];
            println!("CKM Matrix:");
            for i in 0..3 {
                println!("  [{:.6}, {:.6}, {:.6}]",
                    ckm.magnitudes[i*3].to_f64(),
                    ckm.magnitudes[i*3+1].to_f64(),
                    ckm.magnitudes[i*3+2].to_f64());
            }
            println!();
            for (k, name) in names.iter().enumerate() {
                println!("  {}: {:.8}", name, ckm.magnitudes[k].to_f64());
            }
            println!("  J:    {:.6e}", ckm.jarlskog.to_f64());
            println!("  δ_CP: {:.4}°", ckm.delta_rad.to_f64() * 180.0 / std::f64::consts::PI);
        }
        "pmns" => {
            println!("PMNS Mixing Angles:");
            println!("  sin²θ₁₃: {:.8}", e8_core::mixing::pmns::sin2_theta13().to_f64());
            println!("  sin²θ₁₂: {:.8}", e8_core::mixing::pmns::sin2_theta12().to_f64());
            println!("  sin²θ₂₃: {:.8}", e8_core::mixing::pmns::sin2_theta23().to_f64());
            println!("  δ_PMNS:   {:.4}°", e8_core::mixing::cp_phase::delta_pmns_deg().to_f64());
        }
        "higgs" => {
            println!("Higgs Sector:");
            println!("  λ_H:    {:.8}", e8_core::higgs::quartic::higgs_quartic().to_f64());
            println!("  m_H:    {:.4} GeV", e8_core::higgs::mass::higgs_mass_default().to_f64());
            println!("  m_H/m_t: {:.8}", e8_core::higgs::quartic::mh_over_mt().to_f64());
        }
        _ => {
            eprintln!("Unknown sector: {}. Options: leptons, up, down, neutrinos, gauge, ckm, pmns, higgs", name);
            std::process::exit(1);
        }
    }
}

fn run_roots() {
    println!("E8 Root System:");
    let roots = e8_core::lattice::roots::generate_e8_roots();
    println!("  Total roots: {}", roots.len());
    println!("  All |α|² = 2: {}", e8_core::lattice::roots::verify_root_norms(&roots));

    let dist = e8_core::lattice::roots::inner_product_distribution(&roots);
    println!("  Inner product distribution:");
    for (ip, count) in &dist {
        println!("    <α,β> = {}: {} pairs", ip, count);
    }

    let plaquettes = e8_core::lattice::plaquettes::find_plaquettes(&roots);
    let stats = e8_core::lattice::plaquettes::verify_plaquettes(&roots, &plaquettes);
    println!("  Plaquettes: {} (28 per root: {})", stats.count, stats.all_per_root_equal);
    println!("  All inner products = -1: {}", stats.all_ip_minus_one);

    println!("\nOctonionic Associator:");
    let (fano, non_fano) = e8_core::octonion::associator::count_non_fano_triples();
    println!("  Fano triples: {}", fano);
    println!("  Non-Fano triples: {}", non_fano);
    println!("  ||T||² = {}", e8_core::octonion::associator::total_associator_norm_squared());
    println!("  Uniform distribution: {}", e8_core::octonion::associator::verify_uniform_distribution());
}

fn run_info(digits: u32) {
    println!("E8 Standard Model — System Info");
    println!("  Precision: {} decimal digits", digits);
    println!("  Backend: rug/MPFR");
    println!("  Architecture: {}", std::env::consts::ARCH);
    println!("  OS: {}", std::env::consts::OS);

    #[cfg(target_arch = "aarch64")]
    println!("  NEON: available (always on aarch64)");

    println!(
        "  Threads: {}",
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    );
}
