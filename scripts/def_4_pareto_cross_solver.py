"""Script 4: Curva de Pareto y Cross-Solver (Clarabel vs SCS).

Auditoria de robustez y validacion cruzada de solvers:
- x_optimo: solucion optima par para h=1e-3 (mu_cert = 0.692716).
- Barrido de lam_M en [1e-11, 1e-10, 1e-9, 1e-8, 1e-7].
- Resolucion dual idéntica con Clarabel y con SCS (eps=1e-12).
- Certificacion primal cruzada con mpmath a 50 digitos para ambos solvers.
- Tabla comparativa de mu_cert y lam_min(M) fisico.
"""

from __future__ import annotations

import math
import warnings
warnings.filterwarnings("ignore")

import cvxpy as cp
import mpmath
import numpy as np
import scipy.linalg

from kuramoto_templates.cases.s7 import compute_q_star_and_w_star
from kuramoto_templates.domain.conditions import (
    build_ac_matrix,
    build_ad_matrix,
    build_af_matrix,
    build_lw_matrix,
    get_pair_indices,
    p_matrix_to_vector,
    p_vector_to_matrix,
)
from kuramoto_templates.evaluators.scaled import (
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)

q_star, w_star = compute_q_star_and_w_star(dps=50)
a_star = float(1.0 / (7.0 + q_star))
x_W_base = float(1.0 - 6.0 * a_star)

h = 0.001
half_pi = math.pi / 2.0
pi = math.pi

phi = np.array([
    pi - w_star * h,      # 0
    -half_pi + h,         # 1
    -half_pi - h,         # 2
    0.0,                  # 3 (W fijo)
    half_pi + h,          # 4
    half_pi - h,          # 5
    -pi + w_star * h,     # 6
], dtype=float)

# Vector x_optimo obtenido en la optimizacion par exacta para h=1e-3
dx_0_6 = +1.112791e-05
dx_1_5 = +9.353764e-05
dx_2_4 = -1.348106e-04
dx_W = -2.0 * (dx_0_6 + dx_1_5 + dx_2_4)

x_opt = np.array([
    a_star + dx_0_6,
    a_star + dx_1_5,
    a_star + dx_2_4,
    x_W_base + dx_W,
    a_star + dx_2_4,
    a_star + dx_1_5,
    a_star + dx_0_6,
], dtype=float)


def compute_mpmath_M(xc: np.ndarray, phic: np.ndarray, P_mat: np.ndarray) -> float:
    """Calcula el autovalor minimo fisico lam_min(M) en mpmath a 50 digitos."""
    m_len = len(xc)
    p_v = p_matrix_to_vector(np.array(P_mat, dtype=float))
    L_w = build_lw_matrix(xc, phic, p_v)
    Q_np = scipy.linalg.null_space(np.ones((1, m_len), dtype=float))
    dim = m_len - 1
    mpmath.mp.dps = 50
    Q_mp = mpmath.mp.matrix(m_len, dim)
    for r in range(m_len):
        for c in range(dim):
            Q_mp[r, c] = mpmath.mp.mpf(str(Q_np[r, c]))
    for j in range(dim):
        v = Q_mp[:, j]
        for k in range(j):
            proj = mpmath.mp.fdot(Q_mp[:, k], v)
            for r in range(m_len):
                v[r] -= proj * Q_mp[r, k]
        norm_v = mpmath.sqrt(sum(v[r] ** 2 for r in range(m_len)))
        for r in range(m_len):
            Q_mp[r, j] = v[r] / norm_v
    L_mp = mpmath.mp.matrix(m_len, m_len)
    for r in range(m_len):
        for c in range(m_len):
            L_mp[r, c] = mpmath.mp.mpf(str(L_w[r, c]))
    M_mp = Q_mp.T * L_mp * Q_mp
    eigs = mpmath.mp.eigsy(M_mp, eigvals_only=True)
    return float(min(eigs))


def solve_with_solver(
    xc: np.ndarray,
    phic: np.ndarray,
    lam_M: float,
    solver_name: str,
) -> tuple[float | None, float | None, str]:
    """Resuelve la LMI reescalada con el solver especificado y certifica en mpmath."""
    m = len(xc)
    pairs = get_pair_indices(m)
    K = len(pairs)
    sheets, ref_sheet = compute_sheets(xc, phic)
    h_eff = float(compute_h_eff(xc, phic, sheets))

    T = build_scaled_basis_t(xc, phic, sheets, ref_sheet, h_eff)
    assert T is not None

    P_perp = np.eye(m) - np.ones((m, m)) / float(m)
    T_prime = P_perp @ T
    metric_M_rescaled = (T_prime.T @ T_prime) / (h_eff**2)

    A_F = build_af_matrix(xc, phic)
    A_c = build_ac_matrix(xc, phic)
    A_d = build_ad_matrix(xc)
    N = scipy.linalg.null_space(A_F)

    corr, _, _, _ = np.linalg.lstsq(A_F, A_F @ N, rcond=None)
    N = N - corr
    N, _ = np.linalg.qr(N)

    G_rescaled_basis = []
    for a, b in pairs:
        d_ab = np.zeros(m, dtype=float)
        d_ab[a] = 1.0
        d_ab[b] = -1.0
        Td = T.T @ d_ab
        weight = xc[a] * xc[b] * math.cos(phic[a] - phic[b])
        Gk = weight * np.outer(Td, Td) / (h_eff**2)
        G_rescaled_basis.append(0.5 * (Gk + Gk.T))

    z = cp.Variable(N.shape[1])
    t = cp.Variable()
    p = N @ z
    L_expr = sum(p[k] * G_rescaled_basis[k] for k in range(K))

    prob = cp.Problem(
        cp.Maximize(t),
        [
            xc + A_d @ p >= t,
            p >= 0.0,
            p <= 1.0,
            L_expr >> lam_M * metric_M_rescaled,
            xc + A_c @ p >= 1e-4,
        ],
    )

    if solver_name == "CLARABEL":
        prob.solve(
            solver=cp.CLARABEL,
            tol_feas=1e-10,
            tol_gap_abs=1e-10,
            tol_gap_rel=1e-10,
            verbose=False,
        )
    elif solver_name == "SCS":
        prob.solve(
            solver=cp.SCS,
            eps_abs=1e-10,
            eps_rel=1e-10,
            max_iters=25000,
            verbose=False,
        )
    else:
        raise ValueError(f"Solver desconocido: {solver_name}")

    if prob.status not in ("optimal", "optimal_inaccurate") or z.value is None:
        return None, None, prob.status

    p_clip = np.clip(N @ z.value, 0.0, 1.0)
    mu_cert = float(np.min(xc + A_d @ p_clip))
    P_mat = p_vector_to_matrix(p_clip, m)
    lam_mp = compute_mpmath_M(xc, phic, P_mat)
    return mu_cert, lam_mp, prob.status


def main() -> None:
    print("=" * 110)
    print("SCRIPT 4: CURVA DE PARETO Y CROSS-SOLVER (CLARABEL vs SCS)")
    print("Punto evaluado: x_optimo (deformacion par exacta, h=1e-3, delta_W=0.0)")
    print("Certificacion primal rigurosa en mpmath a 50 digitos sobre P recortada")
    print("=" * 110)

    lam_list = [1e-11, 1e-10, 1e-9, 1e-8, 1e-7]
    header = (
        f"{'lam_M':<10} | {'mu_cert (Clarabel)':<20} | {'lam_mp (Clarabel)':<20} | "
        f"{'mu_cert (SCS)':<18} | {'lam_mp (SCS)':<18} | {'Diferencia mu':<14}"
    )
    print(f"\n{header}")
    print("-" * len(header))

    for lam_test in lam_list:
        mu_cl, lam_cl, st_cl = solve_with_solver(x_opt, phi, lam_test, "CLARABEL")
        mu_scs, lam_scs, st_scs = solve_with_solver(x_opt, phi, lam_test, "SCS")

        d_mu = abs(mu_cl - mu_scs) if (mu_cl is not None and mu_scs is not None) else float("nan")

        str_mu_cl = f"{mu_cl:.8f}" if mu_cl is not None else "FAIL"
        str_lam_cl = f"{lam_cl:.2e}" if lam_cl is not None else "FAIL"
        str_mu_scs = f"{mu_scs:.8f}" if mu_scs is not None else "FAIL"
        str_lam_scs = f"{lam_scs:.2e}" if lam_scs is not None else "FAIL"
        str_d_mu = f"{d_mu:.2e}" if not math.isnan(d_mu) else "N/A"

        print(
            f"{lam_test:<10.0e} | {str_mu_cl:<20} | {str_lam_cl:<20} | "
            f"{str_mu_scs:<18} | {str_lam_scs:<18} | {str_d_mu:<14}"
        )

    print("-" * len(header))
    print("=" * 110)


if __name__ == "__main__":
    main()
