//! CKM matrix from Fritzsch texture + octonionic CP phase.
//!
//! Down sector: real symmetric Fritzsch texture (10×5̄ = inner product)
//! Up sector: complex Hermitian Fritzsch texture (10×10 = octonionic cross product)
//!
//! Self-energy corrections:
//!   Down: D₁ = -m_u (U(1)), D₂ = -dim(su(3))×m_u (SU(3))
//!   Up:   D₁ = C₂(SU3,fund)×m_u, D₂ = C₂(SU3)×(|W(G₂)|+1)/h(G₂) × m_c

use crate::algebra::groups::identities::{DIM_IM_OCTONIONS, W_G2_PLUS_1};
use crate::algebra::groups::{G2, SU3};
use crate::mass::sectors::AllMasses;
use crate::precision::scalar::Scalar;

/// CKM matrix result: 3×3 complex matrix + derived quantities.
#[derive(Debug, Clone)]
pub struct CKMResult<S: Scalar> {
    /// |V_ij| magnitudes, row-major: [ud, us, ub, cd, cs, cb, td, ts, tb]
    pub magnitudes: [S; 9],
    /// Jarlskog invariant J
    pub jarlskog: S,
    /// CP phase δ in radians
    pub delta_rad: S,
}

/// Build the full CKM matrix from quark masses (zero free parameters).
#[allow(clippy::needless_range_loop)]
pub fn build_ckm<S: Scalar>(masses: &AllMasses<S>) -> CKMResult<S> {
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
    let d1_d = m_u.clone().neg();
    let d2_d = S::from_i64(-(SU3.dimension as i64)) * m_u.clone();

    // Symmetric functions with sign convention: m_d, -m_s, m_b
    let s1_d = m_d.clone() - m_s.clone() + m_b.clone();
    let md_ms = m_d.clone() * m_s.clone();
    let md_mb = m_d.clone() * m_b.clone();
    let ms_mb = m_s.clone() * m_b.clone();
    let neg_md_ms = md_ms.clone().neg();
    let s12_d = neg_md_ms + md_mb - ms_mb;
    let s123_d = (md_ms * m_b.clone()).neg();

    let d3_d = s1_d.clone() - d1_d.clone() - d2_d.clone();

    // Solve for A², B²
    let q_d = d1_d.clone() * d2_d.clone() + d1_d.clone() * d3_d.clone() + d2_d.clone() * d3_d.clone();
    let ab_d = q_d - s12_d;
    let rhs_d = d1_d.clone() * d2_d.clone() * d3_d.clone() - s123_d;
    let numerator = rhs_d - d1_d.clone() * ab_d.clone();
    let denominator = d3_d.clone() - d1_d.clone();
    let a_sq_d = numerator / denominator;
    let b_sq_d = ab_d - a_sq_d.clone();

    let a_d = a_sq_d.sqrt();
    let b_d = b_sq_d.sqrt();

    // Build 3×3 real symmetric Fritzsch matrix
    // [D1, A, 0; A, D2, B; 0, B, D3]
    let md_matrix = [
        [d1_d.clone(), a_d.clone(), S::zero()],
        [a_d, d2_d.clone(), b_d.clone()],
        [S::zero(), b_d, d3_d.clone()],
    ];

    // Diagonalize real symmetric matrix → eigenvectors
    let (evals_d, evecs_d) = diag_real_symmetric_3x3(&md_matrix);

    // Physical ordering: eigenvalues are [-m_s, m_d, m_b] → reorder to [d, s, b]
    let order_d = physical_order(&evals_d);

    // ════════════════════════════════════════════
    // UP SECTOR: Complex Hermitian Fritzsch texture
    // ════════════════════════════════════════════
    // D₁ = C₂(SU3,fund) × m_u = (4/3) × m_u
    let (c2_num, c2_den) = SU3.c2_fundamental;
    let c2_fund = S::from_u64(c2_num as u64) / S::from_u64(c2_den as u64);
    let d1_u = c2_fund.clone() * m_u.clone();
    // D₂ = C₂(SU3,fund) × (|W(G₂)|+1) / h(G₂) × m_c = (4/3)×13/6 × m_c = (26/9) × m_c
    let d2_u_factor = c2_fund.clone() * S::from_u64(W_G2_PLUS_1)
        / S::from_u64(G2.coxeter_number as u64);
    let d2_u = d2_u_factor * m_c.clone();

    let s1_u = m_u.clone() - m_c.clone() + m_t.clone();
    let d3_u = s1_u - d1_u.clone() - d2_u.clone();

    // CP phase from octonionic associator [e₆,e₃,e₁]
    // arg(C) = π/dim(Im(O)) = π/7
    let arg_c = S::pi() / S::from_u64(DIM_IM_OCTONIONS as u64);

    // |C| = √(C₂(fund)) × √(m_u × m_t)
    let c_mag = (c2_fund.clone() * m_u.clone() * m_t.clone()).sqrt();

    // C as complex number
    let c_re = c_mag.clone() * arg_c.cos();
    let c_im = c_mag.clone() * arg_c.sin();

    let c_sq = c_mag.clone() * c_mag; // |C|²

    // Solve for A, B using eigenvalue constraints
    let mu_mc = m_u.clone() * m_c.clone();
    let mu_mt = m_u.clone() * m_t.clone();
    let mc_mt = m_c.clone() * m_t.clone();
    let s12_u = mu_mt - mu_mc.clone() - mc_mt;
    let sig_u_temp = d1_u.clone() * d2_u.clone() + d1_u.clone() * d3_u.clone() + d2_u.clone() * d3_u.clone();
    let sig_u = sig_u_temp - s12_u - c_sq.clone();

    let rhs_u = d1_u.clone() * d2_u.clone() * d3_u.clone() + mu_mc * m_t.clone();
    let numerator_u = rhs_u - d1_u.clone() * sig_u.clone();
    let denominator_u = d3_u.clone() - d1_u.clone();
    let a_sq_u_guess = numerator_u / denominator_u;

    // Newton refinement for A
    let a_u = newton_solve_a(&d1_u, &d2_u, &d3_u, &sig_u, &c_re, &c_sq, m_u, m_c, m_t, &a_sq_u_guess);
    let b_sq_u = sig_u - a_u.clone() * a_u.clone();

    let b_u = b_sq_u.sqrt();

    // Build complex Hermitian matrix:
    // [D1,    A,       C    ]
    // [A,     D2,      B    ]
    // [C*,    B,       D3   ]
    // Diagonalize via real reduction (M → Re+Im decomposition)
    let (evals_u, evecs_u) = diag_complex_hermitian_3x3(
        &d1_u, &d2_u, &d3_u, &a_u, &b_u, &c_re, &c_im,
    );

    let order_u = physical_order(&evals_u);

    // Phase matrix P† = diag(1, e^{4πi/dim(Im(O))}, e^{4πi/dim(Im(O))})
    // 4π/7 from Fano plane phase assignment: (5-4-0)×π/7 → 4π/7
    let phase_angle = S::from_u64(4) * S::pi() / S::from_u64(DIM_IM_OCTONIONS as u64);
    let p_cos = phase_angle.cos();
    let p_sin = phase_angle.sin();
    let p_re = [S::one(), p_cos.clone(), p_cos];
    let p_im = [S::zero(), p_sin.clone(), p_sin];

    // CKM = Q_u† × P† × Q_d
    // V[i][j] = Σ_k conj(Q_u[k][order_u[i]]) × P†[k] × Q_d[k][order_d[j]]
    let mut v_re = [[S::zero(), S::zero(), S::zero()],
                    [S::zero(), S::zero(), S::zero()],
                    [S::zero(), S::zero(), S::zero()]];
    let mut v_im = [[S::zero(), S::zero(), S::zero()],
                    [S::zero(), S::zero(), S::zero()],
                    [S::zero(), S::zero(), S::zero()]];

    for i in 0..3 {
        for j in 0..3 {
            let mut sum_re = S::zero();
            let mut sum_im = S::zero();
            for k in 0..3 {
                // conj(Q_u[k][i]) = (q_u_re, -q_u_im)
                let qu_re = evecs_u[k][order_u[i]].0.clone();
                let qu_im_neg = evecs_u[k][order_u[i]].1.clone().neg();

                // P†[k] = (p_re[k], p_im[k])
                // Product: conj(Q_u) × P†
                let cp_re = qu_re.clone() * p_re[k].clone() - qu_im_neg.clone() * p_im[k].clone();
                let cp_im = qu_re * p_im[k].clone() + qu_im_neg * p_re[k].clone();

                // × Q_d[k][j] (real for down sector)
                let qd = evecs_d[k][order_d[j]].clone();

                sum_re = sum_re + cp_re * qd.clone();
                sum_im = sum_im + cp_im * qd;
            }
            v_re[i][j] = sum_re;
            v_im[i][j] = sum_im;
        }
    }

    // Phase convention: make V_ud real positive
    let v_ud_abs = (v_re[0][0].clone() * v_re[0][0].clone() + v_im[0][0].clone() * v_im[0][0].clone()).sqrt();
    if v_ud_abs > S::from_f64(1e-30) {
        // phase factor = conj(V_ud)/|V_ud|
        let pf_re = v_re[0][0].clone() / v_ud_abs.clone();
        let pf_im = v_im[0][0].clone().neg() / v_ud_abs;
        for i in 0..3 {
            for j in 0..3 {
                let new_re = v_re[i][j].clone() * pf_re.clone() - v_im[i][j].clone() * pf_im.clone();
                let new_im = v_re[i][j].clone() * pf_im.clone() + v_im[i][j].clone() * pf_re.clone();
                v_re[i][j] = new_re;
                v_im[i][j] = new_im;
            }
        }
    }

    // Extract magnitudes
    let mut magnitudes: [S; 9] = std::array::from_fn(|_| S::zero());
    for i in 0..3 {
        for j in 0..3 {
            magnitudes[i * 3 + j] =
                (v_re[i][j].clone() * v_re[i][j].clone() + v_im[i][j].clone() * v_im[i][j].clone()).sqrt();
        }
    }

    // Jarlskog invariant: J = Im(V_ud V_cs V*_us V*_cd)
    // V_ud × V_cs
    let prod1_re = v_re[0][0].clone() * v_re[1][1].clone() - v_im[0][0].clone() * v_im[1][1].clone();
    let prod1_im = v_re[0][0].clone() * v_im[1][1].clone() + v_im[0][0].clone() * v_re[1][1].clone();

    // conj(V_us) × conj(V_cd) = conj(V_us × V_cd)
    let prod2_re = v_re[0][1].clone() * v_re[1][0].clone() - v_im[0][1].clone() * v_im[1][0].clone();
    let prod2_im = v_re[0][1].clone() * v_im[1][0].clone() + v_im[0][1].clone() * v_re[1][0].clone();
    // conjugate:
    let prod2c_re = prod2_re;
    let prod2c_im = prod2_im.neg();

    // Full product
    let j_im = prod1_re * prod2c_im.clone() + prod1_im * prod2c_re;
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
fn diag_real_symmetric_3x3<S: Scalar>(m: &[[S; 3]; 3]) -> (Vec<S>, Vec<Vec<S>>) {
    let n = 3;

    // Copy matrix
    let mut a = vec![vec![S::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = m[i][j].clone();
        }
    }

    // Eigenvector matrix (starts as identity)
    let mut v = vec![vec![S::zero(); n]; n];
    for i in 0..n {
        v[i][i] = S::one();
    }

    // Jacobi iteration
    for _iter in 0..100 {
        // Find largest off-diagonal element
        let mut max_val = S::zero();
        let mut p = 0;
        let mut q = 1;
        for i in 0..n {
            for j in (i + 1)..n {
                let aij_abs = a[i][j].abs();
                if aij_abs > max_val {
                    max_val = aij_abs;
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < S::from_f64(1e-200) {
            break;
        }

        // Compute rotation angle
        let diff = a[q][q].clone() - a[p][p].clone();
        let theta = if diff.abs() < S::from_f64(1e-200) {
            S::pi() / S::from_u64(4)
        } else {
            let tau = S::from_u64(2) * a[p][q].clone() / diff;
            // atan via acos: atan(x) = acos(1/sqrt(1+x²)) with sign
            // Simpler: use the identity theta = atan(tau)/2
            // We compute atan as acos(1/sqrt(1+x²)) * sign(x)
            let tau_abs = tau.abs();
            let atan_val = (S::one() / (S::one() + tau_abs.clone() * tau_abs).sqrt()).acos();
            let signed = if tau < S::zero() { atan_val.neg() } else { atan_val };
            signed / S::from_u64(2)
        };
        let c = theta.cos();
        let s_rot = theta.sin();

        // Apply rotation
        let mut new_a = a.clone();
        let c_sq = c.clone() * c.clone();
        let s_sq = s_rot.clone() * s_rot.clone();
        let cs2 = S::from_u64(2) * c.clone() * s_rot.clone();
        new_a[p][p] = c_sq.clone() * a[p][p].clone() - cs2.clone() * a[p][q].clone() + s_sq.clone() * a[q][q].clone();
        new_a[q][q] = s_sq * a[p][p].clone() + cs2 * a[p][q].clone() + c_sq * a[q][q].clone();
        new_a[p][q] = S::zero();
        new_a[q][p] = S::zero();

        for r in 0..n {
            if r != p && r != q {
                new_a[p][r] = c.clone() * a[p][r].clone() - s_rot.clone() * a[q][r].clone();
                new_a[r][p] = new_a[p][r].clone();
                new_a[q][r] = s_rot.clone() * a[p][r].clone() + c.clone() * a[q][r].clone();
                new_a[r][q] = new_a[q][r].clone();
            }
        }
        a = new_a;

        // Update eigenvectors
        let mut new_v = v.clone();
        for r in 0..n {
            new_v[r][p] = c.clone() * v[r][p].clone() - s_rot.clone() * v[r][q].clone();
            new_v[r][q] = s_rot.clone() * v[r][p].clone() + c.clone() * v[r][q].clone();
        }
        v = new_v;
    }

    // Eigenvalues are on diagonal
    let eigenvalues: Vec<S> = (0..n).map(|i| a[i][i].clone()).collect();

    (eigenvalues, v)
}

/// Diagonalize a 3×3 complex Hermitian matrix.
/// Uses characteristic polynomial (cubic formula) + null space for eigenvectors.
/// Returns (real eigenvalues sorted ascending, complex eigenvectors as (re, im) pairs).
fn diag_complex_hermitian_3x3<S: Scalar>(
    d1: &S, d2: &S, d3: &S,
    a: &S,  // real off-diagonal (0,1)
    b: &S,  // real off-diagonal (1,2)
    c_re: &S, c_im: &S,
) -> (Vec<S>, Vec<Vec<(S, S)>>) {
    // Characteristic polynomial: λ³ - Sλ² + Qλ - R = 0
    let s_trace = d1.clone() + d2.clone() + d3.clone();

    let a_sq = a.clone() * a.clone();
    let b_sq = b.clone() * b.clone();
    let c_sq = c_re.clone() * c_re.clone() + c_im.clone() * c_im.clone();
    let q_cofac = d1.clone() * d2.clone() + d1.clone() * d3.clone() + d2.clone() * d3.clone()
        - a_sq.clone() - b_sq.clone() - c_sq.clone();

    // det = d1*d2*d3 - d1*b² - a²*d3 + 2*a*b*Re(c) - d2*|c|²
    let det_r = d1.clone() * d2.clone() * d3.clone()
        - d1.clone() * b_sq
        - a_sq.clone() * d3.clone()
        + S::from_u64(2) * a.clone() * b.clone() * c_re.clone()
        - d2.clone() * c_sq.clone();

    // Depressed cubic: t³ + pt + q = 0 where λ = t + S/3
    let s_over_3 = s_trace.clone() / S::from_u64(3);
    let s_sq = s_trace.clone() * s_trace.clone();
    let p_dep = q_cofac.clone() - s_sq.clone() / S::from_u64(3);
    let sq_term = s_trace.clone() * q_cofac / S::from_u64(3);
    let cu_term = S::from_u64(2) * s_sq * s_trace / S::from_u64(27);
    let q_dep = det_r.neg() + sq_term - cu_term;

    // Trigonometric solution
    let neg_p = p_dep.clone().neg();
    let neg_p_over_3 = neg_p / S::from_u64(3);
    let m = S::from_u64(2) * neg_p_over_3.clone().sqrt();

    let three_q = S::from_u64(3) * q_dep;
    let two_p = S::from_u64(2) * p_dep.clone();
    let ratio = three_q / two_p;
    let neg_3_over_p = S::from_i64(-3) / p_dep;
    let cos_phi_arg = ratio * neg_3_over_p.sqrt();

    // Clamp to [-1, 1] for numerical safety
    let cos_phi_clamped = if cos_phi_arg > S::one() {
        S::one()
    } else if cos_phi_arg < S::one().neg() {
        S::one().neg()
    } else {
        cos_phi_arg
    };
    let phi = cos_phi_clamped.acos();

    let two_pi_3 = S::from_u64(2) * S::pi() / S::from_u64(3);

    let phi_over_3 = phi / S::from_u64(3);
    let mut evals = Vec::with_capacity(3);
    for k in 0..3u64 {
        let shift = two_pi_3.clone() * S::from_u64(k);
        let angle = phi_over_3.clone() - shift;
        let t = m.clone() * angle.cos();
        evals.push(t + s_over_3.clone());
    }

    // Sort eigenvalues ascending
    let mut idx: Vec<usize> = vec![0, 1, 2];
    idx.sort_by(|&i, &j| evals[i].partial_cmp(&evals[j]).unwrap_or(std::cmp::Ordering::Equal));
    let sorted_evals: Vec<S> = idx.iter().map(|&i| evals[i].clone()).collect();

    // Compute eigenvectors from null space of (M - λI)
    let mut result_evecs = vec![
        vec![(S::zero(), S::zero()); 3]; 3
    ];

    for (col, lam) in sorted_evals.iter().enumerate() {
        let n00 = d1.clone() - lam.clone();
        let n11 = d2.clone() - lam.clone();
        let n22 = d3.clone() - lam.clone();

        // Row 0: (n00, a, c_re+ic_im)
        // Row 1: (a, n11, b)
        // v = row0 × row1
        let v0_re = a.clone() * b.clone() - c_re.clone() * n11.clone();
        let v0_im = (c_im.clone() * n11.clone()).neg();
        let v1_re = c_re.clone() * a.clone() - n00.clone() * b.clone();
        let v1_im = c_im.clone() * a.clone();
        let v2_re = n00.clone() * n11.clone() - a_sq.clone();
        let v2_im = S::zero();

        let norm_sq = v0_re.clone() * v0_re.clone() + v0_im.clone() * v0_im.clone()
            + v1_re.clone() * v1_re.clone() + v1_im.clone() * v1_im.clone()
            + v2_re.clone() * v2_re.clone();

        let (f0_re, f0_im, f1_re, f1_im, f2_re, f2_im) =
            if norm_sq > S::from_f64(1e-100) {
                (v0_re, v0_im, v1_re, v1_im, v2_re, v2_im)
            } else {
                // Fallback: cross product of row 0 and row 2
                let w0_re = a.clone() * n22.clone() - c_re.clone() * b.clone();
                let w0_im = (c_im.clone() * b.clone()).neg();
                let w1_re = c_sq.clone() - n00.clone() * n22.clone();
                let w1_im = S::zero();
                let w2_re = n00.clone() * b.clone() - a.clone() * c_re.clone();
                let w2_im = a.clone() * c_im.clone();
                (w0_re, w0_im, w1_re, w1_im, w2_re, w2_im)
            };

        // Normalize
        let norm_sq2 = f0_re.clone() * f0_re.clone() + f0_im.clone() * f0_im.clone()
            + f1_re.clone() * f1_re.clone() + f1_im.clone() * f1_im.clone()
            + f2_re.clone() * f2_re.clone() + f2_im.clone() * f2_im.clone();
        let norm = norm_sq2.sqrt();

        result_evecs[0][col] = (f0_re / norm.clone(), f0_im / norm.clone());
        result_evecs[1][col] = (f1_re / norm.clone(), f1_im / norm.clone());
        result_evecs[2][col] = (f2_re / norm.clone(), f2_im / norm);
    }

    (sorted_evals, result_evecs)
}

/// Evaluate the determinant residual: det(M) - (-m_u m_c m_t) for given A value.
#[allow(clippy::too_many_arguments)]
fn det_residual<S: Scalar>(
    a_val: &S, d1: &S, d2: &S, d3: &S,
    sig: &S, c_re: &S, c_sq: &S, target_det: &S,
) -> Option<S> {
    let b_sq = sig.clone() - a_val.clone() * a_val.clone();
    if b_sq < S::zero() {
        return None;
    }
    let b_val = b_sq.clone().sqrt();

    let det = d1.clone() * d2.clone() * d3.clone()
        - d1.clone() * b_sq
        - a_val.clone() * a_val.clone() * d3.clone()
        + S::from_u64(2) * a_val.clone() * b_val * c_re.clone()
        - d2.clone() * c_sq.clone();
    Some(det - target_det.clone())
}

/// Newton solve for A in the up-sector determinant constraint.
#[allow(clippy::too_many_arguments)]
fn newton_solve_a<S: Scalar>(
    d1: &S, d2: &S, d3: &S, sig: &S,
    c_re: &S, c_sq: &S,
    m_u: &S, m_c: &S, m_t: &S,
    a_sq_guess: &S,
) -> S {
    let target_det = (m_u.clone() * m_c.clone() * m_t.clone()).neg();

    // Start from sqrt(max(0, Asq_guess))
    let mut a_val = if a_sq_guess.clone() > S::zero() {
        a_sq_guess.sqrt()
    } else {
        S::zero()
    };

    // Numerical derivative epsilon: build 10^{-half_digits}
    // Use a conservative eps based on ~25 digits of precision
    let mut eps_scale = S::one();
    for _ in 0..25 {
        eps_scale = eps_scale / S::from_u64(10);
    }
    let eps_sq = eps_scale.clone() * eps_scale.clone();

    for _iter in 0..200 {
        let residual = match det_residual(&a_val, d1, d2, d3, sig, c_re, c_sq, &target_det) {
            Some(r) => r,
            None => break,
        };

        if residual.abs() < eps_sq.clone() {
            break;
        }

        // Numerical derivative
        let a_abs = a_val.abs();
        let base = if a_abs.clone() > S::one() { a_abs } else { S::one() };
        let eps = base * eps_scale.clone();
        let a_plus = a_val.clone() + eps.clone();
        let residual_plus = match det_residual(&a_plus, d1, d2, d3, sig, c_re, c_sq, &target_det) {
            Some(r) => r,
            None => break,
        };

        let deriv = (residual_plus - residual.clone()) / eps;
        if deriv.abs() < S::from_f64(1e-250) {
            break;
        }

        let step = residual / deriv;
        a_val = a_val - step;

        // Clamp to valid range
        let sig_sqrt = sig.sqrt();
        if a_val.clone() < S::zero() {
            a_val = S::zero();
        } else if a_val.clone() > sig_sqrt.clone() {
            a_val = sig_sqrt * S::from_f64(0.99);
        }
    }

    a_val
}

/// Physical eigenvalue ordering:
/// Eigenvalues are approximately [-m_s, m_d, m_b] sorted ascending.
/// Physical order: [d, s, b] → indices [1, 0, 2] into sorted eigenvalues.
fn physical_order<S: Scalar>(evals: &[S]) -> [usize; 3] {
    let mut indices: Vec<usize> = (0..3).collect();
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
        let masses = compute_all_masses::<rug::Float>();
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
