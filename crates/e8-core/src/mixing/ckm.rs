//! CKM matrix from Fritzsch texture + octonionic CP phase.
//!
//! Down sector: real symmetric Fritzsch texture (10×5̄ = inner product)
//! Up sector: complex Hermitian Fritzsch texture (10×10 = octonionic cross product)
//!
//! Self-energy corrections:
//!   Down: D₁ = -m_u (U(1)), D₂ = -dim(su(3))×m_u (SU(3))
//!   Up:   D₁ = C₂(SU3,fund)×m_u, D₂ = C₂(SU3)×(|W(G₂)|+1)/h(G₂) × m_c

use rug::ops::NegAssign;
use rug::Float;

use crate::algebra::groups::identities::{DIM_IM_OCTONIONS, W_G2_PLUS_1};
use crate::algebra::groups::{G2, SU3};
use crate::mass::sectors::AllMasses;
use crate::precision::{pi, precision_bits, sqrt};

/// CKM matrix result: 3×3 complex matrix + derived quantities.
#[derive(Debug, Clone)]
pub struct CKMResult {
    /// |V_ij| magnitudes, row-major: [ud, us, ub, cd, cs, cb, td, ts, tb]
    pub magnitudes: [Float; 9],
    /// Jarlskog invariant J
    pub jarlskog: Float,
    /// CP phase δ in radians
    pub delta_rad: Float,
}

/// Build the full CKM matrix from quark masses (zero free parameters).
#[allow(clippy::needless_range_loop)]
pub fn build_ckm(masses: &AllMasses<Float>) -> CKMResult {
    let prec = precision_bits();

    let m_u = &masses.up;
    let m_c = &masses.charm;
    let m_t = &masses.top;
    let m_d = &masses.down;
    let m_s = &masses.strange;
    let m_b = &masses.bottom;

    // ════════════════════════════════════════════
    // DOWN SECTOR: Real symmetric Fritzsch texture
    // ════════════════════════════════════════════
    // Self-energy corrections:
    //   D₁ = -m_u (U(1) self-energy)
    //   D₂ = -dim(su(3)) × m_u (SU(3) self-energy)
    let d1_d = Float::with_val(prec, -m_u);
    let d2_d = Float::with_val(prec, -(SU3.dimension as i32)) * m_u;

    // Symmetric functions with sign convention: m_d, -m_s, m_b
    let s1_d = Float::with_val(prec, Float::with_val(prec, m_d - m_s) + m_b);
    let md_ms = Float::with_val(prec, m_d * m_s);
    let md_mb = Float::with_val(prec, m_d * m_b);
    let ms_mb = Float::with_val(prec, m_s * m_b);
    let neg_md_ms = Float::with_val(prec, -&md_ms);
    let s12_d = Float::with_val(prec, Float::with_val(prec, &neg_md_ms + &md_mb) - &ms_mb);
    let s123_d = Float::with_val(prec, -(Float::with_val(prec, &md_ms) * m_b));

    let d3_d = Float::with_val(prec, Float::with_val(prec, &s1_d - &d1_d) - &d2_d);

    // Solve for A², B²
    let q_d = Float::with_val(prec, Float::with_val(prec, &d1_d * &d2_d + &d1_d * &d3_d) + &d2_d * &d3_d);
    let ab_d = Float::with_val(prec, &q_d - &s12_d);
    let rhs_d = Float::with_val(prec, Float::with_val(prec, &d1_d * &d2_d) * &d3_d - &s123_d);
    let numerator = Float::with_val(prec, &rhs_d - Float::with_val(prec, &d1_d * &ab_d));
    let denominator = Float::with_val(prec, &d3_d - &d1_d);
    let a_sq_d = Float::with_val(prec, numerator / denominator);
    let b_sq_d = Float::with_val(prec, &ab_d - &a_sq_d);

    let a_d = sqrt(&a_sq_d);
    let b_d = sqrt(&b_sq_d);

    // Build 3×3 real symmetric Fritzsch matrix
    // [D1, A, 0; A, D2, B; 0, B, D3]
    let md_matrix = [
        [d1_d.clone(), a_d.clone(), Float::with_val(prec, 0)],
        [a_d, d2_d.clone(), b_d.clone()],
        [Float::with_val(prec, 0), b_d, d3_d.clone()],
    ];

    // Diagonalize real symmetric matrix → eigenvectors
    let (evals_d, evecs_d) = diag_real_symmetric_3x3(&md_matrix);

    // Physical ordering: eigenvalues are [-m_s, m_d, m_b] → reorder to [d, s, b]
    let order_d = physical_order_down(&evals_d);

    // ════════════════════════════════════════════
    // UP SECTOR: Complex Hermitian Fritzsch texture
    // ════════════════════════════════════════════
    // D₁ = C₂(SU3,fund) × m_u = (4/3) × m_u
    let (c2_num, c2_den) = SU3.c2_fundamental;
    let c2_fund = Float::with_val(prec, c2_num) / Float::with_val(prec, c2_den);
    let d1_u = Float::with_val(prec, &c2_fund * m_u);
    // D₂ = C₂(SU3,fund) × (|W(G₂)|+1) / h(G₂) × m_c = (4/3)×13/6 × m_c = (26/9) × m_c
    let d2_u_factor = Float::with_val(prec, &c2_fund * Float::with_val(prec, W_G2_PLUS_1))
        / Float::with_val(prec, G2.coxeter_number);
    let d2_u = Float::with_val(prec, d2_u_factor * m_c);

    let s1_u = Float::with_val(prec, Float::with_val(prec, m_u - m_c) + m_t);
    let d3_u = Float::with_val(prec, Float::with_val(prec, &s1_u - &d1_u) - &d2_u);

    // CP phase from octonionic associator [e₆,e₃,e₁]
    // arg(C) = π/dim(Im(O)) = π/7
    let arg_c = pi() / Float::with_val(prec, DIM_IM_OCTONIONS);

    // |C| = √(C₂(fund)) × √(m_u × m_t)
    let c_mag = sqrt(&(Float::with_val(prec, Float::with_val(prec, &c2_fund * m_u) * m_t)));

    // C as complex number
    let c_re = Float::with_val(prec, &c_mag * crate::precision::cos(&arg_c));
    let c_im = Float::with_val(prec, &c_mag * crate::precision::sin(&arg_c));

    let c_sq = Float::with_val(prec, &c_mag * &c_mag); // |C|²

    // Solve for A, B using eigenvalue constraints
    let mu_mc = Float::with_val(prec, m_u * m_c);
    let mu_mt = Float::with_val(prec, m_u * m_t);
    let mc_mt = Float::with_val(prec, m_c * m_t);
    let neg_mu_mc = Float::with_val(prec, -&mu_mc);
    let s12_u = Float::with_val(prec, Float::with_val(prec, &neg_mu_mc + &mu_mt) - &mc_mt);
    let sig_u_temp = Float::with_val(prec, Float::with_val(prec, &d1_u * &d2_u + &d1_u * &d3_u) + &d2_u * &d3_u);
    let sig_u = Float::with_val(prec, Float::with_val(prec, sig_u_temp - &s12_u) - &c_sq);

    let rhs_u = Float::with_val(prec, Float::with_val(prec, &d1_u * &d2_u) * &d3_u + Float::with_val(prec, m_u * m_c) * m_t);
    let numerator_u = Float::with_val(prec, &rhs_u - Float::with_val(prec, &d1_u * &sig_u));
    let denominator_u = Float::with_val(prec, &d3_u - &d1_u);
    let a_sq_u_guess = Float::with_val(prec, numerator_u / denominator_u);

    // Newton refinement for A
    let a_u = newton_solve_a(&d1_u, &d2_u, &d3_u, &sig_u, &c_re, &c_sq, m_u, m_c, m_t, &a_sq_u_guess);
    let b_sq_u = Float::with_val(prec, &sig_u - &a_u * &a_u);

    let b_u = sqrt(&b_sq_u);

    // Build complex Hermitian matrix:
    // [D1,    A,       C    ]
    // [A,     D2,      B    ]
    // [C*,    B,       D3   ]
    // Diagonalize via real reduction (M → Re+Im decomposition)
    let (evals_u, evecs_u) = diag_complex_hermitian_3x3(
        &d1_u, &d2_u, &d3_u, &a_u, &b_u, &c_re, &c_im,
    );

    let order_u = physical_order_up(&evals_u);

    // Phase matrix P† = diag(1, e^{4πi/dim(Im(O))}, e^{4πi/dim(Im(O))})
    // 4π/7 from Fano plane phase assignment: (5-4-0)×π/7 → 4π/7
    let phase_angle_temp = Float::with_val(prec, Float::with_val(prec, 4) * pi());
    let phase_angle = phase_angle_temp / Float::with_val(prec, DIM_IM_OCTONIONS);
    let p_re = [
        Float::with_val(prec, 1),
        crate::precision::cos(&phase_angle),
        crate::precision::cos(&phase_angle),
    ];
    let p_im = [
        Float::with_val(prec, 0),
        crate::precision::sin(&phase_angle),
        crate::precision::sin(&phase_angle),
    ];

    // CKM = Q_u† × P† × Q_d
    // V[i][j] = Σ_k conj(Q_u[k][order_u[i]]) × P†[k] × Q_d[k][order_d[j]]
    let mut v_re = [[Float::with_val(prec, 0), Float::with_val(prec, 0), Float::with_val(prec, 0)],
                    [Float::with_val(prec, 0), Float::with_val(prec, 0), Float::with_val(prec, 0)],
                    [Float::with_val(prec, 0), Float::with_val(prec, 0), Float::with_val(prec, 0)]];
    let mut v_im = [[Float::with_val(prec, 0), Float::with_val(prec, 0), Float::with_val(prec, 0)],
                    [Float::with_val(prec, 0), Float::with_val(prec, 0), Float::with_val(prec, 0)],
                    [Float::with_val(prec, 0), Float::with_val(prec, 0), Float::with_val(prec, 0)]];

    for i in 0..3 {
        for j in 0..3 {
            let mut sum_re = Float::with_val(prec, 0);
            let mut sum_im = Float::with_val(prec, 0);
            for k in 0..3 {
                // conj(Q_u[k][i]) = (q_u_re, -q_u_im)
                let qu_re = &evecs_u[k][order_u[i]].0;
                let qu_im_neg = Float::with_val(prec, -&evecs_u[k][order_u[i]].1);

                // P†[k] = (p_re[k], p_im[k])
                // Product: conj(Q_u) × P†
                let cp_re = Float::with_val(prec, Float::with_val(prec, qu_re * &p_re[k]) - &qu_im_neg * &p_im[k]);
                let cp_im = Float::with_val(prec, Float::with_val(prec, qu_re * &p_im[k]) + &qu_im_neg * &p_re[k]);

                // × Q_d[k][j] (real for down sector)
                let qd = &evecs_d[k][order_d[j]];

                sum_re += &cp_re * qd;
                sum_im += &cp_im * qd;
            }
            v_re[i][j] = sum_re;
            v_im[i][j] = sum_im;
        }
    }

    // Phase convention: make V_ud real positive
    let v_ud_abs = sqrt(&(Float::with_val(prec, Float::with_val(prec, &v_re[0][0] * &v_re[0][0]) + &v_im[0][0] * &v_im[0][0])));
    if v_ud_abs > Float::with_val(prec, 1e-30f64) {
        // phase factor = conj(V_ud)/|V_ud|
        let pf_re = Float::with_val(prec, &v_re[0][0] / &v_ud_abs);
        let neg_vim = Float::with_val(prec, -&v_im[0][0]);
        let pf_im = Float::with_val(prec, neg_vim / &v_ud_abs);
        for i in 0..3 {
            for j in 0..3 {
                let new_re = Float::with_val(prec, Float::with_val(prec, &v_re[i][j] * &pf_re) - &v_im[i][j] * &pf_im);
                let new_im = Float::with_val(prec, Float::with_val(prec, &v_re[i][j] * &pf_im) + &v_im[i][j] * &pf_re);
                v_re[i][j] = new_re;
                v_im[i][j] = new_im;
            }
        }
    }

    // Extract magnitudes
    let mut magnitudes = std::array::from_fn(|_| Float::with_val(prec, 0));
    for i in 0..3 {
        for j in 0..3 {
            magnitudes[i * 3 + j] =
                sqrt(&(Float::with_val(prec, Float::with_val(prec, &v_re[i][j] * &v_re[i][j]) + &v_im[i][j] * &v_im[i][j])));
        }
    }

    // Jarlskog invariant: J = Im(V_ud V_cs V*_us V*_cd)
    // V_ud × V_cs
    let prod1_re = Float::with_val(prec, Float::with_val(prec, &v_re[0][0] * &v_re[1][1]) - &v_im[0][0] * &v_im[1][1]);
    let prod1_im = Float::with_val(prec, Float::with_val(prec, &v_re[0][0] * &v_im[1][1]) + &v_im[0][0] * &v_re[1][1]);

    // conj(V_us) × conj(V_cd) = conj(V_us × V_cd)
    let prod2_re = Float::with_val(prec, Float::with_val(prec, &v_re[0][1] * &v_re[1][0]) - &v_im[0][1] * &v_im[1][0]);
    let prod2_im = Float::with_val(prec, Float::with_val(prec, &v_re[0][1] * &v_im[1][0]) + &v_im[0][1] * &v_re[1][0]);
    // conjugate:
    let prod2c_re = prod2_re;
    let mut prod2c_im = prod2_im;
    prod2c_im.neg_assign();

    // Full product
    let j_im = Float::with_val(prec, Float::with_val(prec, &prod1_re * &prod2c_im) + &prod1_im * &prod2c_re);
    let jarlskog = j_im.abs();

    CKMResult {
        magnitudes,
        jarlskog,
        delta_rad: super::cp_phase::delta_ckm_rad(),
    }
}

#[allow(clippy::needless_range_loop)]
/// Diagonalize a 3×3 real symmetric matrix via Jacobi rotations.
/// Returns (eigenvalues sorted ascending, eigenvector columns).
fn diag_real_symmetric_3x3(m: &[[Float; 3]; 3]) -> (Vec<Float>, Vec<Vec<Float>>) {
    let prec = precision_bits();
    let n = 3;

    // Copy matrix
    let mut a = vec![vec![Float::with_val(prec, 0); n]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = m[i][j].clone();
        }
    }

    // Eigenvector matrix (starts as identity)
    let mut v = vec![vec![Float::with_val(prec, 0); n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = Float::with_val(prec, 1);
    }

    // Jacobi iteration
    for _iter in 0..100 {
        // Find largest off-diagonal element
        let mut max_val = Float::with_val(prec, 0);
        let mut p = 0;
        let mut q = 1;
        for i in 0..n {
            for j in (i + 1)..n {
                let aij_abs = a[i][j].clone().abs();
                if aij_abs > max_val {
                    max_val = aij_abs;
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < Float::with_val(prec, 1e-200f64) {
            break;
        }

        // Compute rotation angle
        let diff = Float::with_val(prec, &a[q][q] - &a[p][p]);
        let theta = if diff.clone().abs() < Float::with_val(prec, 1e-200f64) {
            pi() / Float::with_val(prec, 4)
        } else {
            let tau = Float::with_val(prec, Float::with_val(prec, 2) * &a[p][q]) / &diff;
            Float::with_val(prec, tau.atan() / Float::with_val(prec, 2))
        };
        let c = crate::precision::cos(&theta);
        let s = crate::precision::sin(&theta);

        // Apply rotation
        let mut new_a = a.clone();
        let c_sq = Float::with_val(prec, &c * &c);
        let s_sq = Float::with_val(prec, &s * &s);
        let cs2 = Float::with_val(prec, Float::with_val(prec, 2) * Float::with_val(prec, &c * &s));
        new_a[p][p] = Float::with_val(
            prec,
            Float::with_val(prec, Float::with_val(prec, &c_sq * &a[p][p]) - Float::with_val(prec, &cs2 * &a[p][q])) + &s_sq * &a[q][q],
        );
        new_a[q][q] = Float::with_val(
            prec,
            Float::with_val(prec, Float::with_val(prec, &s_sq * &a[p][p]) + Float::with_val(prec, &cs2 * &a[p][q])) + &c_sq * &a[q][q],
        );
        new_a[p][q] = Float::with_val(prec, 0);
        new_a[q][p] = Float::with_val(prec, 0);

        for r in 0..n {
            if r != p && r != q {
                new_a[p][r] = Float::with_val(prec, &c * &a[p][r] - &s * &a[q][r]);
                new_a[r][p] = new_a[p][r].clone();
                new_a[q][r] = Float::with_val(prec, &s * &a[p][r] + &c * &a[q][r]);
                new_a[r][q] = new_a[q][r].clone();
            }
        }
        a = new_a;

        // Update eigenvectors
        let mut new_v = v.clone();
        for r in 0..n {
            new_v[r][p] = Float::with_val(prec, &c * &v[r][p] - &s * &v[r][q]);
            new_v[r][q] = Float::with_val(prec, &s * &v[r][p] + &c * &v[r][q]);
        }
        v = new_v;
    }

    // Eigenvalues are on diagonal
    let eigenvalues: Vec<Float> = (0..n).map(|i| a[i][i].clone()).collect();

    (eigenvalues, v)
}

/// Diagonalize a 3×3 complex Hermitian matrix.
/// Uses characteristic polynomial (cubic formula) + null space for eigenvectors.
/// Returns (real eigenvalues sorted ascending, complex eigenvectors as (re, im) pairs).
/// Eigenvector format: evecs[row][col] where col indexes eigenvalue.
fn diag_complex_hermitian_3x3(
    d1: &Float,
    d2: &Float,
    d3: &Float,
    a: &Float,  // real off-diagonal (0,1)
    b: &Float,  // real off-diagonal (1,2)
    c_re: &Float,
    c_im: &Float,
) -> (Vec<Float>, Vec<Vec<(Float, Float)>>) {
    let prec = precision_bits();

    // M = [[d1, a, c], [a, d2, b], [c*, b, d3]]
    // Characteristic polynomial: λ³ - Sλ² + Qλ - R = 0
    // S = tr(M)
    let s_trace = Float::with_val(prec, Float::with_val(prec, d1 + d2) + d3);

    // Q = sum of 2×2 principal minors
    let a_sq = Float::with_val(prec, a * a);
    let b_sq = Float::with_val(prec, b * b);
    let c_sq = Float::with_val(prec, Float::with_val(prec, c_re * c_re) + c_im * c_im);
    let q_cofac = Float::with_val(
        prec,
        Float::with_val(prec, Float::with_val(prec, d1 * d2 + d1 * d3) + d2 * d3)
            - Float::with_val(prec, Float::with_val(prec, &a_sq + &b_sq) + &c_sq),
    );

    // R = det(M) = d1(d2*d3-b²) - a(a*d3-b*Re(c)) + Re(c)(a*b-d2*Re(c)) + Im(c)(...)
    // det = d1*d2*d3 - d1*b² - a²*d3 + 2*a*b*Re(c) - d2*|c|²
    let det_r = Float::with_val(
        prec,
        Float::with_val(
            prec,
            Float::with_val(
                prec,
                Float::with_val(prec, Float::with_val(prec, d1 * d2) * d3)
                    - Float::with_val(prec, d1 * &b_sq),
            ) - Float::with_val(prec, &a_sq * d3),
        ) + Float::with_val(
            prec,
            Float::with_val(prec, Float::with_val(prec, 2) * Float::with_val(prec, a * b) * c_re)
                - Float::with_val(prec, d2 * &c_sq),
        ),
    );

    // Depressed cubic: t³ + pt + q = 0 where λ = t + S/3
    // p = Q - S²/3
    // q = -R + S*Q/3 - 2S³/27
    let s_over_3 = Float::with_val(prec, &s_trace / Float::with_val(prec, 3));
    let s_sq = Float::with_val(prec, &s_trace * &s_trace);
    let p_dep = Float::with_val(prec, &q_cofac - Float::with_val(prec, &s_sq / Float::with_val(prec, 3)));
    let sq_term = Float::with_val(prec, Float::with_val(prec, &s_trace * &q_cofac) / Float::with_val(prec, 3));
    let cu_term = Float::with_val(
        prec,
        Float::with_val(prec, Float::with_val(prec, Float::with_val(prec, 2) * &s_sq) * &s_trace) / Float::with_val(prec, 27),
    );
    let q_dep = Float::with_val(
        prec,
        Float::with_val(prec, Float::with_val(prec, -&det_r) + &sq_term) - &cu_term,
    );

    // Trigonometric solution: t_k = 2√(-p/3) × cos(φ/3 - 2πk/3)
    // where cos(φ) = (3q)/(2p) × √(-3/p)
    let neg_p = Float::with_val(prec, -&p_dep);
    let neg_p_over_3 = Float::with_val(prec, &neg_p / Float::with_val(prec, 3));
    let m = Float::with_val(prec, Float::with_val(prec, 2) * sqrt(&neg_p_over_3));

    // cos(φ) = (3q)/(2p) × √(-3/p)
    let three_q = Float::with_val(prec, Float::with_val(prec, 3) * &q_dep);
    let two_p = Float::with_val(prec, Float::with_val(prec, 2) * &p_dep);
    let ratio = Float::with_val(prec, &three_q / &two_p);
    let neg_3_over_p = Float::with_val(prec, Float::with_val(prec, -3) / &p_dep);
    let cos_phi_arg = Float::with_val(prec, &ratio * sqrt(&neg_3_over_p));

    // Clamp to [-1, 1] for numerical safety
    let cos_phi_clamped = if cos_phi_arg > Float::with_val(prec, 1) {
        Float::with_val(prec, 1)
    } else if cos_phi_arg < Float::with_val(prec, -1) {
        Float::with_val(prec, -1)
    } else {
        cos_phi_arg
    };
    let phi = Float::with_val(prec, cos_phi_clamped.acos());

    let two_pi_3 = Float::with_val(prec, Float::with_val(prec, Float::with_val(prec, 2) * pi()) / Float::with_val(prec, 3));

    let phi_over_3 = Float::with_val(prec, &phi / Float::with_val(prec, 3));
    let mut evals = Vec::with_capacity(3);
    for k in 0..3u32 {
        let shift = Float::with_val(prec, &two_pi_3 * k);
        let angle = Float::with_val(prec, &phi_over_3 - &shift);
        let t = Float::with_val(prec, &m * crate::precision::cos(&angle));
        evals.push(Float::with_val(prec, &t + &s_over_3));
    }

    // Sort eigenvalues ascending
    let mut idx: Vec<usize> = vec![0, 1, 2];
    idx.sort_by(|&i, &j| evals[i].partial_cmp(&evals[j]).unwrap_or(std::cmp::Ordering::Equal));
    let sorted_evals: Vec<Float> = idx.iter().map(|&i| evals[i].clone()).collect();

    // Compute eigenvectors from null space of (M - λI)
    // For eigenvalue λ, row 0 = (d1-λ, a, c), row 1 = (a, d2-λ, b)
    // Cross product: v = row0 × row1
    let mut result_evecs = vec![
        vec![(Float::with_val(prec, 0), Float::with_val(prec, 0)); 3]; 3
    ];

    for (col, lam) in sorted_evals.iter().enumerate() {
        let n00 = Float::with_val(prec, d1 - lam); // d1 - λ (real)
        let n11 = Float::with_val(prec, d2 - lam); // d2 - λ (real)
        let n22 = Float::with_val(prec, d3 - lam); // d3 - λ (real)

        // Row 0: (n00, a, c_re+ic_im)
        // Row 1: (a, n11, b)  (all real)
        // v = row0 × row1:
        // v₀ = row0[1]*row1[2] - row0[2]*row1[1]
        //    = a*b - (c_re+ic_im)*n11
        let v0_re = Float::with_val(prec, Float::with_val(prec, a * b) - c_re * &n11);
        let v0_im = Float::with_val(prec, -(Float::with_val(prec, c_im * &n11)));

        // v₁ = row0[2]*row1[0] - row0[0]*row1[2]
        //    = (c_re+ic_im)*a - n00*b
        let v1_re = Float::with_val(prec, Float::with_val(prec, c_re * a) - &n00 * b);
        let v1_im = Float::with_val(prec, c_im * a);

        // v₂ = row0[0]*row1[1] - row0[1]*row1[0]
        //    = n00*n11 - a² (real)
        let v2_re = Float::with_val(prec, Float::with_val(prec, &n00 * &n11) - &a_sq);
        let v2_im = Float::with_val(prec, 0);

        // Check if cross product is too small (rows nearly parallel); try other pair
        let norm_sq = Float::with_val(
            prec,
            Float::with_val(
                prec,
                Float::with_val(prec, &v0_re * &v0_re + &v0_im * &v0_im)
                    + Float::with_val(prec, &v1_re * &v1_re + &v1_im * &v1_im),
            ) + Float::with_val(prec, &v2_re * &v2_re),
        );

        let (f0_re, f0_im, f1_re, f1_im, f2_re, f2_im) =
            if norm_sq > Float::with_val(prec, 1e-100f64) {
                (v0_re, v0_im, v1_re, v1_im, v2_re, v2_im)
            } else {
                // Fallback: cross product of row 0 and row 2
                // Row 2: (c_re-ic_im, b, n22)  (note: conjugate of c)
                // v₀ = row0[1]*row2[2] - row0[2]*row2[1]
                //    = a*n22 - (c_re+ic_im)*b
                let w0_re = Float::with_val(prec, Float::with_val(prec, a * &n22) - c_re * b);
                let w0_im = Float::with_val(prec, -(Float::with_val(prec, c_im * b)));
                // v₁ = row0[2]*row2[0] - row0[0]*row2[2]
                //    = (c_re+ic_im)*(c_re-ic_im) - n00*n22
                //    = |c|² - n00*n22
                let w1_re = Float::with_val(prec, &c_sq - Float::with_val(prec, &n00 * &n22));
                let w1_im = Float::with_val(prec, 0);
                // v₂ = row0[0]*row2[1] - row0[1]*row2[0]
                //    = n00*b - a*(c_re-ic_im)
                let w2_re = Float::with_val(prec, Float::with_val(prec, &n00 * b) - a * c_re);
                let w2_im = Float::with_val(prec, a * c_im);
                (w0_re, w0_im, w1_re, w1_im, w2_re, w2_im)
            };

        // Normalize
        let norm_sq2 = Float::with_val(
            prec,
            Float::with_val(
                prec,
                Float::with_val(prec, &f0_re * &f0_re + &f0_im * &f0_im)
                    + Float::with_val(prec, &f1_re * &f1_re + &f1_im * &f1_im),
            ) + Float::with_val(prec, &f2_re * &f2_re + &f2_im * &f2_im),
        );
        let norm = sqrt(&norm_sq2);

        result_evecs[0][col] = (
            Float::with_val(prec, &f0_re / &norm),
            Float::with_val(prec, &f0_im / &norm),
        );
        result_evecs[1][col] = (
            Float::with_val(prec, &f1_re / &norm),
            Float::with_val(prec, &f1_im / &norm),
        );
        result_evecs[2][col] = (
            Float::with_val(prec, &f2_re / &norm),
            Float::with_val(prec, &f2_im / &norm),
        );
    }

    (sorted_evals, result_evecs)
}

/// Evaluate the determinant residual: det(M) - (-m_u m_c m_t) for given A value.
#[allow(clippy::too_many_arguments)]
fn det_residual(
    a_val: &Float,
    d1: &Float,
    d2: &Float,
    d3: &Float,
    sig: &Float,
    c_re: &Float,
    c_sq: &Float,
    target_det: &Float,
) -> Option<Float> {
    let prec = precision_bits();
    let b_sq = Float::with_val(prec, sig - a_val * a_val);
    if b_sq < Float::with_val(prec, 0) {
        return None;
    }
    let b_val = sqrt(&b_sq);

    let det = Float::with_val(
        prec,
        Float::with_val(
            prec,
            Float::with_val(
                prec,
                Float::with_val(prec, Float::with_val(prec, d1 * d2) * d3)
                    - Float::with_val(prec, d1 * &b_sq),
            ) - Float::with_val(prec, Float::with_val(prec, a_val * a_val) * d3),
        ) + Float::with_val(
            prec,
            Float::with_val(prec, Float::with_val(prec, 2) * Float::with_val(prec, a_val * &b_val)) * c_re
                - Float::with_val(prec, d2 * c_sq),
        ),
    );
    Some(Float::with_val(prec, &det - target_det))
}

/// Newton solve for A in the up-sector determinant constraint.
/// Matches Python's findroot: start from sqrt(max(0, Asq_guess)), iterate
/// with Newton's method using numerical derivative.
#[allow(clippy::too_many_arguments)]
fn newton_solve_a(
    d1: &Float,
    d2: &Float,
    d3: &Float,
    sig: &Float,
    c_re: &Float,
    c_sq: &Float,
    m_u: &Float,
    m_c: &Float,
    m_t: &Float,
    a_sq_guess: &Float,
) -> Float {
    let prec = precision_bits();
    let target_det = Float::with_val(prec, -(Float::with_val(prec, m_u * m_c) * m_t));

    // Match Python: start from sqrt(max(0, Asq_guess))
    let mut a_val = if *a_sq_guess > Float::with_val(prec, 0) {
        sqrt(a_sq_guess)
    } else {
        Float::with_val(prec, 0)
    };

    // Numerical derivative epsilon: ~half the working precision digits
    // For N-digit precision, optimal step for numerical derivative ≈ 10^{-N/2}
    let digits = prec * 3 / 10; // rough bits → digits
    let half_digits = digits / 2;
    // Build 10^{-half_digits}
    let mut eps_scale = Float::with_val(prec, 1);
    for _ in 0..half_digits {
        eps_scale /= Float::with_val(prec, 10);
    }

    for _iter in 0..200 {
        let residual = match det_residual(&a_val, d1, d2, d3, sig, c_re, c_sq, &target_det) {
            Some(r) => r,
            None => break,
        };

        if residual.clone().abs() < Float::with_val(prec, &eps_scale * &eps_scale) {
            break;
        }

        // Numerical derivative: use max(|a_val|, 1) * eps_scale as step
        let a_abs = Float::with_val(prec, a_val.clone().abs());
        let base = if a_abs > Float::with_val(prec, 1) {
            a_abs
        } else {
            Float::with_val(prec, 1)
        };
        let eps = Float::with_val(prec, &base * &eps_scale);
        let a_plus = Float::with_val(prec, &a_val + &eps);
        let residual_plus = match det_residual(&a_plus, d1, d2, d3, sig, c_re, c_sq, &target_det) {
            Some(r) => r,
            None => break,
        };

        let deriv = Float::with_val(prec, Float::with_val(prec, &residual_plus - &residual) / &eps);
        if deriv.clone().abs() < Float::with_val(prec, 1e-250f64) {
            break;
        }

        let step = Float::with_val(prec, &residual / &deriv);
        a_val -= &step;

        // Clamp to valid range: 0 ≤ A² < Sig
        let sig_sqrt = sqrt(sig);
        if a_val < Float::with_val(prec, 0) {
            a_val = Float::with_val(prec, 0);
        } else if a_val > sig_sqrt {
            a_val = Float::with_val(prec, &sig_sqrt * Float::with_val(prec, 0.99f64));
        }
    }

    a_val
}

/// Physical eigenvalue ordering for down sector:
/// Eigenvalues are approximately [-m_s, m_d, m_b] sorted ascending.
/// Physical order: [d, s, b] → indices [1, 0, 2] into sorted eigenvalues.
fn physical_order_down(evals: &[Float]) -> [usize; 3] {
    // Sort indices by eigenvalue
    let mut indices: Vec<usize> = (0..3).collect();
    indices.sort_by(|&a, &b| evals[a].partial_cmp(&evals[b]).unwrap_or(std::cmp::Ordering::Equal));
    // Sorted ascending: [most negative, middle, most positive]
    // Most negative ≈ -m_s → physical s
    // Middle ≈ m_d → physical d
    // Most positive ≈ m_b → physical b
    // Physical order: [d, s, b] = [indices[1], indices[0], indices[2]]
    [indices[1], indices[0], indices[2]]
}

/// Physical eigenvalue ordering for up sector (analogous).
fn physical_order_up(evals: &[Float]) -> [usize; 3] {
    let mut indices: Vec<usize> = (0..evals.len()).collect();
    indices.sort_by(|&a, &b| evals[a].partial_cmp(&evals[b]).unwrap_or(std::cmp::Ordering::Equal));
    [indices[1], indices[0], indices[2]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mass::sectors::compute_all_masses;
    use crate::precision::set_precision;

    #[test]
    fn test_ckm_magnitudes() {
        set_precision(50);
        let masses = compute_all_masses::<Float>();
        let ckm = build_ckm(&masses);

        // Python reference: V_ud=0.974466, V_us=0.224507, V_ub=0.003631
        let v_ud = ckm.magnitudes[0].to_f64();
        assert!(
            (v_ud - 0.9745).abs() < 0.005,
            "V_ud = {} (expected ~0.974)",
            v_ud
        );

        let v_us = ckm.magnitudes[1].to_f64();
        assert!(
            (v_us - 0.2245).abs() < 0.005,
            "V_us = {} (expected ~0.2245)",
            v_us
        );

        let v_ub = ckm.magnitudes[2].to_f64();
        assert!(
            (v_ub - 0.0036).abs() < 0.001,
            "V_ub = {} (expected ~0.0036)",
            v_ub
        );

        // Unitarity: |V_ud|² + |V_us|² + |V_ub|² ≈ 1
        let row1_sum = ckm.magnitudes[0].to_f64().powi(2)
            + ckm.magnitudes[1].to_f64().powi(2)
            + ckm.magnitudes[2].to_f64().powi(2);
        assert!(
            (row1_sum - 1.0).abs() < 0.001,
            "Row 1 unitarity: {}",
            row1_sum
        );

        // Jarlskog: J ≈ 3.0e-5
        let j = ckm.jarlskog.to_f64();
        assert!(
            (j - 3.0e-5).abs() < 1.0e-5,
            "J = {} (expected ~3.0e-5)",
            j
        );
    }
}
