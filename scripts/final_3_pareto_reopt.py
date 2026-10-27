"""Script 3: El Frente de Pareto Robusto (Reoptimizacion Completa).

Auditoria de robustez en el frente de Pareto:
- Parametrizacion simetrica par (3 masas: dx_0_6, dx_1_5, dx_2_4) con delta_W = 0, h = 1e-3.
- Reoptimizacion completa desde cero para lam_M in [1e-10, 1e-9, 1e-8].
- Clarabel tol_feas=1e-10, certificacion primal mpmath 50 digitos.
- Tabla limpia: lam_M, mu_S7, mu_optimo_recalculado, Ganancia.
"""

from __future__ import annotations

import math
import warnings
warnings.filterwarnings("ignore")

import cvxpy as cp
import mpmath
import numpy as np
import scipy.linalg
from scipy.optimize import minimize

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
x_s7_manual = np.array([a_star, a_star, a_star, x_W_base, a_star, a_star, a_star], dtype=float)

h = 0.001
half_pi = math.pi / 2.0
pi = math.pi

phi = np.array([
    pi - w_star * h,      # 0
    -half_pi + h,         # 1
    -half_pi - h,         # 2
    0.0,                  # 3 (W fijo en 0)
    half_pi + h,          # 4
    half_pi - h,          # 5
    -pi + w_star * h,     # 6
], dtype=float)


def u_to_masses(u: np.ndarray | list[float]) -> np.ndarray:
    """Convierte u = [dx_0_6, dx_1_5, dx_2_4] [en 1e-4] a vector de masas completo."""
    dx_0_6 = float(u[0]) * 1e-4
    dx_1_5 = float(u[1]) * 1e-4
    dx_2_4 = float(u[2]) * 1e-4
    dx_W = -2.0 * (dx_0_6 + dx_1_5 + dx_2_4)
    return np.array([
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


def solve_exact(
    xc: np.ndarray,
    phic: np.ndarray,
    lam_M: float = 1e-10,
    tol: float = 1e-10,
    certify: bool = False,
) -> tuple[float | None, float | None, str]:
    """Resuelve la LMI reescalada y certifica primalmente en mpmath."""
    if np.min(xc) <= 1e-4:
        return None, None, "INFEASIBLE_MASS"

    m = len(xc)
    pairs = get_pair_indices(m)
    K = len(pairs)
    sheets, ref_sheet = compute_sheets(xc, phic)
    h_eff = float(compute_h_eff(xc, phic, sheets))
    if h_eff < 1e-6:
        return None, None, "SCALE_TOO_SMALL"

    T = build_scaled_basis_t(xc, phic, sheets, ref_sheet, h_eff)
    if T is None:
        return None, None, "NOT_SCALABLE"

    P_perp = np.eye(m) - np.ones((m, m)) / float(m)
    T_prime = P_perp @ T
    metric_M_rescaled = (T_prime.T @ T_prime) / (h_eff**2)

    A_F = build_af_matrix(xc, phic)
    A_c = build_ac_matrix(xc, phic)
    A_d = build_ad_matrix(xc)
    N = scipy.linalg.null_space(A_F)
    if N.shape[1] == 0:
        return None, None, "NO_NULLSPACE"

    res_inf = float(np.max(np.abs(A_F @ N)))
    if res_inf > 1e-14:
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
    try:
        prob.solve(
            solver=cp.CLARABEL,
            tol_feas=tol,
            tol_gap_abs=tol,
            tol_gap_rel=tol,
            verbose=False,
        )
    except Exception:
        return None, None, "SOLVER_EXCEPTION"

    if prob.status != "optimal":
        return None, None, prob.status

    p_clip = np.clip(N @ z.value, 0.0, 1.0)
    mu_cert = float(np.min(xc + A_d @ p_clip))
    if certify:
        P_mat = p_vector_to_matrix(p_clip, m)
        lam_mp = compute_mpmath_M(xc, phic, P_mat)
    else:
        lam_mp = 0.0
    return mu_cert, lam_mp, prob.status


def main() -> None:
    print("=" * 105)
    print("SCRIPT 3: EL FRENTE DE PARETO ROBUSTO (REOPTIMIZACION COMPLETA)")
    print("LMI REESCALADA (h_eff^2) | Clarabel tol_feas=1e-10 | Certificacion Primal mpmath 50d")
    print("Reoptimizacion de 3 masas simetricas para cada lam_M (h = 1e-3, delta_W = 0.0)")
    print("=" * 105)

    lam_targets = [1.0e-10, 1.0e-9, 1.0e-8]
    results = []

    for lam_val in lam_targets:
        mu_s7, lam_s7, st_s7 = solve_exact(x_s7_manual, phi, lam_M=lam_val, certify=True)

        def loss(u: np.ndarray) -> float:
            xc = u_to_masses(u)
            mu, _, st = solve_exact(xc, phi, lam_M=lam_val, certify=False)
            if st != "optimal" or mu is None:
                return 1e9
            return -mu

        u_start = np.array([+0.11, +0.94, -1.35])
        sim = [u_start]
        for d in range(3):
            pt = np.copy(u_start)
            pt[d] += 0.2
            sim.append(pt)
        sim = np.array(sim)

        res = minimize(
            loss,
            u_start,
            method="Nelder-Mead",
            options={
                "initial_simplex": sim,
                "maxiter": 1000,
                "xatol": 1e-5,
                "fatol": 1e-8,
            },
        )

        xc_opt = u_to_masses(res.x)
        mu_opt, lam_opt, st_opt = solve_exact(xc_opt, phi, lam_M=lam_val, certify=True)
        gain = (mu_opt - mu_s7) * 1e4

        results.append({
            "lam_M": lam_val,
            "mu_s7": mu_s7,
            "mu_opt": mu_opt,
            "gain": gain,
            "lam_opt": lam_opt,
            "u_opt": res.x,
        })

    print(f"\n{'lam_M':<10} | {'mu_cert S7':<16} | {'mu_opt Recalculado':<22} | {'Ganancia (x1e-4)':<18} | {'lam_min(M) Realizado'}")
    print("-" * 95)
    for r in results:
        print(
            f"{r['lam_M']:<10.0e} | {r['mu_s7']:<16.8f} | {r['mu_opt']:<22.8f} | "
            f"{r['gain']:+18.3f} | {r['lam_opt']:.2e}"
        )
    print("-" * 95)

    print("\nPARAMETROS OPTIMOS RECALCULADOS [1e-4]:")
    for r in results:
        u = r["u_opt"]
        print(f"  lam_M = {r['lam_M']:.0e} : dx_0_6={u[0]:+.4f}, dx_1_5={u[1]:+.4f}, dx_2_4={u[2]:+.4f}")
    print("=" * 105)


if __name__ == "__main__":
    main()
