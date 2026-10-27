"""Script 3: La Ley de Escala y la Familia Analitica.

Auditoria de escalabilidad topologica:
- Parametrizacion simetrica par (3 masas: dx_0_6, dx_1_5, dx_2_4) con delta_W = 0.
- Optimizacion Nelder-Mead para h = 5e-4, 1e-3 y 2e-3 a lam_M = 1e-10.
- Tabla limpia con: h, mu_cert_S7, mu_cert_optimo, Ganancia Neta, y (dx_0_6, dx_1_5, dx_2_4) / h.
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

half_pi = math.pi / 2.0
pi = math.pi


def build_phi(h_val: float) -> np.ndarray:
    """Fases exactas de S7(h) con delta_W = 0.0 estricto."""
    return np.array([
        pi - w_star * h_val,      # 0
        -half_pi + h_val,         # 1
        -half_pi - h_val,         # 2
        0.0,                      # 3 (W fijo)
        half_pi + h_val,          # 4
        half_pi - h_val,          # 5
        -pi + w_star * h_val,     # 6
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
    P_mat = p_vector_to_matrix(p_clip, m)
    lam_mp = compute_mpmath_M(xc, phic, P_mat)
    return mu_cert, lam_mp, prob.status


def main() -> None:
    print("=" * 105)
    print("SCRIPT 3: LA LEY DE ESCALA Y LA FAMILIA ANALITICA (SIMETRIA PAR R)")
    print("LMI REESCALADA (h_eff^2) | Clarabel tol_feas=1e-10 | Certificacion Primal mpmath 50d")
    print("Parametrizacion: dx_0_6, dx_1_5, dx_2_4 | delta_W = 0.0 estricto | lam_M = 1.0e-10")
    print("=" * 105)

    scales = [5.0e-4, 1.0e-3, 2.0e-3]
    results = []

    for h_val in scales:
        phi_val = build_phi(h_val)
        mu_s7, lam_s7, st_s7 = solve_exact(x_s7_manual, phi_val)

        def loss(u: np.ndarray) -> float:
            xc = u_to_masses(u)
            mu, _, st = solve_exact(xc, phi_val)
            if st != "optimal" or mu is None:
                return 1e9
            return -mu

        ratio = h_val / 0.001
        u_start = np.array([0.11 * ratio, 0.93 * ratio, -1.35 * ratio])
        sim = [u_start]
        for d in range(3):
            pt = np.copy(u_start)
            pt[d] += 0.2 * ratio
            sim.append(pt)
        sim = np.array(sim)

        res = minimize(
            loss,
            u_start,
            method="Nelder-Mead",
            options={
                "initial_simplex": sim,
                "maxiter": 250,
                "xatol": 1e-4,
                "fatol": 1e-7,
            },
        )
        mu_opt = -res.fun
        xc_opt = u_to_masses(res.x)
        mu_c, lam_c, st_c = solve_exact(xc_opt, phi_val)

        dx_opt = res.x * 1e-4
        dx_over_h = dx_opt / h_val
        gain = (mu_c - mu_s7) * 1e4

        results.append({
            "h": h_val,
            "mu_s7": mu_s7,
            "mu_opt": mu_c,
            "gain": gain,
            "dx_over_h": dx_over_h,
            "lam_min": lam_c,
            "status": st_c,
        })

    # Imprimir tabla limpia requerida
    print(f"\n{'Escala h':<10} | {'mu_cert S7':<14} | {'mu_cert Optimo':<16} | {'Ganancia Neta':<15} | {'(dx_0_6, dx_1_5, dx_2_4) / h'}")
    print("-" * 105)
    for r in results:
        v_str = f"[{r['dx_over_h'][0]:+8.4f}, {r['dx_over_h'][1]:+8.4f}, {r['dx_over_h'][2]:+8.4f}]"
        print(
            f"{r['h']:<10.1e} | {r['mu_s7']:<14.8f} | {r['mu_opt']:<16.8f} | "
            f"{r['gain']:+11.3f}e-4 | {v_str}"
        )
    print("-" * 105)

    print("\nVERIFICACION DE LINEALIDAD DE PRIMER ORDEN:")
    d05 = results[0]["dx_over_h"]
    d10 = results[1]["dx_over_h"]
    d20 = results[2]["dx_over_h"]
    print(f"  h = 5.0e-4 : dx/h = [{d05[0]:+8.4f}, {d05[1]:+8.4f}, {d05[2]:+8.4f}]")
    print(f"  h = 1.0e-3 : dx/h = [{d10[0]:+8.4f}, {d10[1]:+8.4f}, {d10[2]:+8.4f}]")
    print(f"  h = 2.0e-3 : dx/h = [{d20[0]:+8.4f}, {d20[1]:+8.4f}, {d20[2]:+8.4f}]")
    print("=" * 105)


if __name__ == "__main__":
    main()
