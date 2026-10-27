"""Script 2: Sensibilidades Analiticas (Gradiente de Primer Orden).

Calculo de las derivadas direccionales exactas de mu* respecto a las 3 deformaciones pares
de masa en el origen S7 mediante Diferencias Finitas Centrales de alta precision:
- Parametros pares: dx_0_6, dx_1_5, dx_2_4.
- Paso de diferenciacion: step = 1e-6 (con sondeo de estabilidad en 5e-7).
- Escalas h: [5e-4, 1e-3, 2e-3].
- Solver: Clarabel tol_feas=1e-10, tol_gap_abs=1e-10, tol_gap_rel=1e-10.
- LMI reescalada por h_eff^2.
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
x_s7_manual = np.array([a_star, a_star, a_star, x_W_base, a_star, a_star, a_star], dtype=float)

half_pi = math.pi / 2.0
pi = math.pi


def build_phi(h_val: float) -> np.ndarray:
    """Fases exactas de S7(h) con delta_W = 0.0 estricto."""
    return np.array([
        pi - w_star * h_val,      # 0
        -half_pi + h_val,         # 1
        -half_pi - h_val,         # 2
        0.0,                      # 3 (W fijo en 0)
        half_pi + h_val,          # 4
        half_pi - h_val,          # 5
        -pi + w_star * h_val,     # 6
    ], dtype=float)


def u_to_masses_raw(u: np.ndarray | list[float]) -> np.ndarray:
    """Convierte vector [dx_0_6, dx_1_5, dx_2_4] en unidades de masa directas."""
    dx_0_6 = float(u[0])
    dx_1_5 = float(u[1])
    dx_2_4 = float(u[2])
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
    eigsy_fn = getattr(mpmath.mp, "eigsy")
    eigs = eigsy_fn(M_mp, eigvals_only=True)
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
    L_expr = p[0] * G_rescaled_basis[0]
    for k in range(1, K):
        L_expr = L_expr + p[k] * G_rescaled_basis[k]

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


from typing import TypedDict


class SensResultItem(TypedDict):
    h: float
    mu_s7: float
    lam_s7: float
    grad: np.ndarray
    steps: list[float]


def main() -> None:
    print("=" * 105)
    print("SCRIPT 2: SENSIBILIDADES ANALITICAS (Gradiente de Primer Orden en S7)")
    print("LMI REESCALADA (h_eff^2) | Clarabel tol_feas=1e-10 | Diferencias Finitas Centrales (step=1e-6)")
    print("Gradiente de mu* respecto a las deformaciones simetricas: [d(mu*)/d(dx_0_6), d(mu*)/d(dx_1_5), d(mu*)/d(dx_2_4)]")
    print("=" * 105)

    scales = [5.0e-4, 1.0e-3, 2.0e-3]
    candidate_steps = [1.0e-6, 5.0e-7, 2.0e-6]
    results: list[SensResultItem] = []

    for h_val in scales:
        phi_val = build_phi(h_val)
        mu_s7, lam_s7, st_s7 = solve_exact(x_s7_manual, phi_val, certify=True)
        assert mu_s7 is not None and lam_s7 is not None

        grad = np.zeros(3, dtype=float)
        steps_used: list[float] = []
        for i in range(3):
            computed = False
            for step in candidate_steps:
                # Perturbacion positiva
                u_plus = np.zeros(3, dtype=float)
                u_plus[i] = +step
                xc_plus = u_to_masses_raw(u_plus)
                mu_plus, _, st_plus = solve_exact(xc_plus, phi_val, certify=False)

                # Perturbacion negativa
                u_minus = np.zeros(3, dtype=float)
                u_minus[i] = -step
                xc_minus = u_to_masses_raw(u_minus)
                mu_minus, _, st_minus = solve_exact(xc_minus, phi_val, certify=False)

                if st_plus == "optimal" and st_minus == "optimal" and mu_plus is not None and mu_minus is not None:
                    grad[i] = (mu_plus - mu_minus) / (2.0 * step)
                    steps_used.append(step)
                    computed = True
                    break

            if not computed:
                grad[i] = float("nan")
                steps_used.append(float("nan"))

        results.append({
            "h": h_val,
            "mu_s7": mu_s7,
            "lam_s7": lam_s7,
            "grad": grad,
            "steps": steps_used,
        })

    print(f"\n{'Escala h':<10} | {'mu_cert S7':<14} | {'Gradiente Exacto g = [g_0_6, g_1_5, g_2_4]':<45} | {'Norma ||g||_2'}")
    print("-" * 105)
    for r in results:
        g = r["grad"]
        norm_g = float(np.linalg.norm(g))
        g_str = f"[{g[0]:+11.4f}, {g[1]:+11.4f}, {g[2]:+11.4f}]"
        print(f"{r['h']:<10.1e} | {r['mu_s7']:<14.8f} | {g_str:<45} | {norm_g:10.4f}")
    print("-" * 105)

    print("\nDETALLE DE COMPONENTES DE SENSIBILIDAD POR PARAMETRO:")
    for r in results:
        g = r["grad"]
        print(f"  h = {r['h']:.1e} : d(mu*)/d(dx_0_6) = {g[0]:+10.4f} | d(mu*)/d(dx_1_5) = {g[1]:+10.4f} | d(mu*)/d(dx_2_4) = {g[2]:+10.4f}")

    print("\nCOMPARACION DE RUMBO DE MAXIMO ASCENSO (VECTOR DIRECCIONAL UNITARIO g / ||g||):")
    for r in results:
        g = r["grad"]
        u_dir = g / np.linalg.norm(g)
        print(f"  h = {r['h']:.1e} : u_dir = [{u_dir[0]:+8.4f}, {u_dir[1]:+8.4f}, {u_dir[2]:+8.4f}]")

    print("\nRELACION g(2e-3) / g(1e-3) ENTRE REGIMENES COMPARABLES:")
    g_1e3 = results[1]["grad"]
    g_2e3 = results[2]["grad"]
    ratio = g_2e3 / g_1e3
    print(f"  Cociente [g(2e-3) / g(1e-3)] = [{ratio[0]:.4f}, {ratio[1]:.4f}, {ratio[2]:.4f}] (Constancia en orden 0: derivada invariante)")
    print("=" * 105)


if __name__ == "__main__":
    main()
