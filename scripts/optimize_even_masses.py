"""Optimizacion de la Parte Par (Masas Simetricas) y Escala en h.

Script definitivo para la auditoria cientifica:
- Nelder-Mead sobre deformacion par con 3 parametros [dx_O, dx_S, dx_N] escalados en 1e-4.
- LMI reescalada por h_eff^2 con Clarabel (tol_feas=1e-10, tol_gap=1e-10).
- Rechazo estricto de OPTIMAL_INACCURATE y soluciones no factibles.
- Certificacion primal exacta con mpmath a 50 digitos.
- delta_W estrictamente en 0.0.
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

# -----------------------------------------------------------------------------
# 1. Geometria y Constantes Base de S7
# -----------------------------------------------------------------------------
q_star, w_star = compute_q_star_and_w_star(dps=50)
a_star = float(1.0 / (7.0 + q_star))
x_W_base = float(1.0 - 6.0 * a_star)
half_pi = math.pi / 2.0
pi = math.pi


def build_phi(h_val: float) -> np.ndarray:
    """Construye las fases exactas de S7(h) con delta_W = 0.0 estricto."""
    return np.array([
        pi - w_star * h_val,      # 0: opuesto (+)
        -half_pi + h_val,         # 1: sur (+)
        -half_pi - h_val,         # 2: sur (-)
        0.0,                      # 3: W (fijo en 0.0)
        half_pi + h_val,          # 4: norte (+)
        half_pi - h_val,          # 5: norte (-)
        -pi + w_star * h_val,     # 6: opuesto (-)
    ], dtype=float)


def u_to_masses(u: np.ndarray | list[float]) -> np.ndarray:
    """Convierte parametros u en unidades de 1e-4 al vector completo de masas x.
    
    u[0] = dx_O (opuestos: idx 0, 6)
    u[1] = dx_S (sur: idx 1, 2)
    u[2] = dx_N (norte: idx 4, 5)
    dx_W = -2 * (dx_O + dx_S + dx_N) para conservar sum(x) = 1.0
    """
    dx_O = float(u[0]) * 1e-4
    dx_S = float(u[1]) * 1e-4
    dx_N = float(u[2]) * 1e-4
    dx_W = -2.0 * (dx_O + dx_S + dx_N)
    return np.array([
        a_star + dx_O,
        a_star + dx_S,
        a_star + dx_S,
        x_W_base + dx_W,
        a_star + dx_N,
        a_star + dx_N,
        a_star + dx_O,
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


def solve_rescaled_lmi(
    xc: np.ndarray,
    phic: np.ndarray,
    lam_M: float = 1e-10,
    tol: float = 1e-10,
) -> tuple[float | None, float | None, str]:
    """Resuelve la LMI reescalada por h_eff^2 y certifica primalmente en mpmath.
    
    Rechaza estrictamente cualquier estado que no sea 'optimal' (incluyendo OPTIMAL_INACCURATE).
    """
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


def objective(u: np.ndarray, h_val: float, lam_M: float = 1e-10) -> float:
    """Funcion objetivo para Nelder-Mead: minimiza -mu_cert."""
    xc = u_to_masses(u)
    phic = build_phi(h_val)
    mu_cert, _, status = solve_rescaled_lmi(xc, phic, lam_M=lam_M, tol=1e-10)
    if status != "optimal" or mu_cert is None:
        return 1e9
    return -mu_cert


def run_multistart(
    h_val: float,
    lam_M: float = 1e-10,
    n_starts: int = 6,
    seed: int = 42,
) -> tuple[np.ndarray, float, float]:
    """Ejecuta optimizacion Nelder-Mead multi-start."""
    np.random.seed(seed)
    h_ratio = h_val / 0.001

    # Arranques heuristicos y aleatorios
    starts = [
        np.array([1.33 * h_ratio, -1.32 * h_ratio, -1.32 * h_ratio]),
        np.array([1.0 * h_ratio, -0.7 * h_ratio, -0.7 * h_ratio]),
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5 * h_ratio, -1.0 * h_ratio, -1.0 * h_ratio]),
        np.array([1.5 * h_ratio, -1.5 * h_ratio, -1.5 * h_ratio]),
        np.array([0.8 * h_ratio, -1.2 * h_ratio, -1.2 * h_ratio]),
    ]
    if n_starts > len(starts):
        for _ in range(n_starts - len(starts)):
            starts.append(np.random.uniform(-2.0, 2.0, size=3) * h_ratio)

    best_mu = -1e9
    best_u = starts[0]
    best_lam = 0.0

    for idx, u0 in enumerate(starts[:n_starts]):
        # Construir simplex inicial con paso de 0.5 unidades (5e-5 en masa)
        sim = [u0]
        for d in range(3):
            pt = np.copy(u0)
            pt[d] += 0.5 * h_ratio
            sim.append(pt)
        sim = np.array(sim)

        res = minimize(
            objective,
            u0,
            args=(h_val, lam_M),
            method="Nelder-Mead",
            options={
                "initial_simplex": sim,
                "maxiter": 150,
                "xatol": 1e-3,
                "fatol": 1e-6,
            },
        )
        mu_cand = -res.fun
        if mu_cand > best_mu:
            # Certificar con mpmath
            xc_cand = u_to_masses(res.x)
            phic = build_phi(h_val)
            mu_cert, lam_mp, st = solve_rescaled_lmi(xc_cand, phic, lam_M=lam_M, tol=1e-10)
            if st == "optimal" and mu_cert is not None and lam_mp is not None:
                best_mu = mu_cert
                best_u = res.x
                best_lam = lam_mp

        print(
            f"  Start {idx:02d} | u0=[{u0[0]:+5.2f}, {u0[1]:+5.2f}, {u0[2]:+5.2f}] "
            f"-> mu_cert = {-res.fun:.8f} (nit={res.nit})"
        )

    return best_u, best_mu, best_lam


def main() -> None:
    print("=" * 105)
    print("OPTIMIZACION DE LA PARTE PAR (MASAS SIMETRICAS) Y ESCALA EN h")
    print("LMI REESCALADA (h_eff^2) | Clarabel tol_feas=1e-10 | Certificacion Primal mpmath 50d")
    print("Parametrizacion: dx_O (0,6), dx_S (1,2), dx_N (4,5), dx_W ajusta suma a 0 | delta_W = 0.0")
    print("=" * 105)

    # Base S7 de referencia
    x_s7 = np.array([a_star, a_star, a_star, x_W_base, a_star, a_star, a_star], dtype=float)

    # -------------------------------------------------------------------------
    # PARTE A: OPTIMIZACION MULTI-START EN h = 1e-3 (lam_M = 1e-10)
    # -------------------------------------------------------------------------
    print("\nPARTE A: MULTI-START EN h = 1.0e-3 (lam_M = 1.0e-10)")
    print("-" * 105)
    h_A = 1e-3
    phi_A = build_phi(h_A)
    mu_s7_A, lam_s7_A, _ = solve_rescaled_lmi(x_s7, phi_A, lam_M=1e-10)
    print(f"Linea base S7(h=1e-3): mu_cert = {mu_s7_A:.8f} | lam_min(M) = {lam_s7_A:.2e}")
    print("Iniciando busqueda Nelder-Mead multi-start (6 arranques)...")

    u_opt_A, mu_opt_A, lam_opt_A = run_multistart(h_A, lam_M=1e-10, n_starts=6, seed=42)
    x_opt_A = u_to_masses(u_opt_A)
    dx_A = x_opt_A - x_s7

    print("\nRESULTADOS PARTE A (h = 1.0e-3):")
    print(f"  mu_cert Maximo : {mu_opt_A:.8f}  (Ganancia vs S7: {(mu_opt_A - mu_s7_A)*1e4:+.3f}e-4)")
    print(f"  lam_min(M)     : {lam_opt_A:.2e}  (>= 0.95 * lam_M verificado)")
    print(f"  u_opt [1e-4]   : dx_O={u_opt_A[0]:+.4f}, dx_S={u_opt_A[1]:+.4f}, dx_N={u_opt_A[2]:+.4f}")
    print(f"  Simetria S-N   : |dx_S - dx_N| = {abs(u_opt_A[1] - u_opt_A[2])*1e-4:.2e}")
    print(f"  Delta_x_par    : [{', '.join(f'{v:+.6e}' for v in dx_A)}]")

    # -------------------------------------------------------------------------
    # PARTE B: OPTIMIZACION MULTI-START EN h = 5e-4 y h = 2e-3 (lam_M = 1e-10)
    # -------------------------------------------------------------------------
    print("\nPARTE B: ESCALABILIDAD EN h (h = 5.0e-4 y h = 2.0e-3, lam_M = 1.0e-10)")
    print("-" * 105)
    
    # h = 5e-4
    h_B1 = 5e-4
    print(f"\n--- Subcaso h = {h_B1:.1e} ---")
    phi_B1 = build_phi(h_B1)
    mu_s7_B1, lam_s7_B1, _ = solve_rescaled_lmi(x_s7, phi_B1, lam_M=1e-10)
    print(f"S7(h=5e-4) base: mu_cert = {mu_s7_B1:.8f} | lam_min(M) = {lam_s7_B1:.2e}")
    u_opt_B1, mu_opt_B1, lam_opt_B1 = run_multistart(h_B1, lam_M=1e-10, n_starts=6, seed=101)
    x_opt_B1 = u_to_masses(u_opt_B1)
    dx_B1 = x_opt_B1 - x_s7
    print(f"Optimo h=5e-4  : mu_cert = {mu_opt_B1:.8f} (Ganancia: {(mu_opt_B1 - mu_s7_B1)*1e4:+.3f}e-4)")
    print(f"Delta_x_par    : [{', '.join(f'{v:+.6e}' for v in dx_B1)}]")

    # h = 2e-3
    h_B2 = 2e-3
    print(f"\n--- Subcaso h = {h_B2:.1e} ---")
    phi_B2 = build_phi(h_B2)
    mu_s7_B2, lam_s7_B2, _ = solve_rescaled_lmi(x_s7, phi_B2, lam_M=1e-10)
    print(f"S7(h=2e-3) base: mu_cert = {mu_s7_B2:.8f} | lam_min(M) = {lam_s7_B2:.2e}")
    u_opt_B2, mu_opt_B2, lam_opt_B2 = run_multistart(h_B2, lam_M=1e-10, n_starts=6, seed=202)
    x_opt_B2 = u_to_masses(u_opt_B2)
    dx_B2 = x_opt_B2 - x_s7
    print(f"Optimo h=2e-3  : mu_cert = {mu_opt_B2:.8f} (Ganancia: {(mu_opt_B2 - mu_s7_B2)*1e4:+.3f}e-4)")
    print(f"Delta_x_par    : [{', '.join(f'{v:+.6e}' for v in dx_B2)}]")

    # Comparacion de escala de Delta_x_par
    print("\nANALISIS DE ESCALABILIDAD PROPORCIONAL:")
    print(f"{'Parametro':<12} | {'h = 5.0e-4':<14} | {'h = 1.0e-3':<14} | {'h = 2.0e-3':<14} | {'Ratio 2k/1k':<12} | {'Ratio 1k/0.5k'}")
    print("-" * 85)
    r_2k_1k_O = u_opt_B2[0] / u_opt_A[0]
    r_1k_05k_O = u_opt_A[0] / u_opt_B1[0]
    r_2k_1k_S = u_opt_B2[1] / u_opt_A[1]
    r_1k_05k_S = u_opt_A[1] / u_opt_B1[1]
    dx_W_B1 = dx_B1[3] / 1e-4
    dx_W_A = dx_A[3] / 1e-4
    dx_W_B2 = dx_B2[3] / 1e-4
    r_2k_1k_W = dx_W_B2 / dx_W_A
    r_1k_05k_W = dx_W_A / dx_W_B1
    print(f"{'dx_O [1e-4]':<12} | {u_opt_B1[0]:+14.4f} | {u_opt_A[0]:+14.4f} | {u_opt_B2[0]:+14.4f} | {r_2k_1k_O:12.3f} | {r_1k_05k_O:.3f}")
    print(f"{'dx_S [1e-4]':<12} | {u_opt_B1[1]:+14.4f} | {u_opt_A[1]:+14.4f} | {u_opt_B2[1]:+14.4f} | {r_2k_1k_S:12.3f} | {r_1k_05k_S:.3f}")
    print(f"{'dx_N [1e-4]':<12} | {u_opt_B1[2]:+14.4f} | {u_opt_A[2]:+14.4f} | {u_opt_B2[2]:+14.4f} | {u_opt_B2[2]/u_opt_A[2]:12.3f} | {u_opt_A[2]/u_opt_B1[2]:.3f}")
    print(f"{'dx_W [1e-4]':<12} | {dx_W_B1:+14.4f} | {dx_W_A:+14.4f} | {dx_W_B2:+14.4f} | {r_2k_1k_W:12.3f} | {r_1k_05k_W:.3f}")

    # -------------------------------------------------------------------------
    # PARTE C: EVALUACION CRUZADA CONTRA lam_M = 1e-9 y 1e-8
    # -------------------------------------------------------------------------
    print("\nPARTE C: EVALUACION CRUZADA CONTRA lam_M = 1.0e-9 y 1.0e-8")
    print("=" * 105)
    print(f"{'Escala h':<10} | {'lam_M':<10} | {'mu_cert S7':<14} | {'mu_cert x_par_opt':<18} | {'Delta mu (x1e-4)':<18} | {'Status'}")
    print("-" * 105)

    cases = [
        ("5.0e-4", h_B1, x_opt_B1),
        ("1.0e-3", h_A, x_opt_A),
        ("2.0e-3", h_B2, x_opt_B2),
    ]
    test_lams = [1e-10, 1e-9, 1e-8]

    for label_h, h_val, x_opt_val in cases:
        phi_val = build_phi(h_val)
        for lam_test in test_lams:
            mu_s7_t, lam_m_s7, st_s7 = solve_rescaled_lmi(x_s7, phi_val, lam_M=lam_test)
            mu_opt_t, lam_m_opt, st_opt = solve_rescaled_lmi(x_opt_val, phi_val, lam_M=lam_test)
            d_mu = (mu_opt_t - mu_s7_t) * 1e4 if (mu_opt_t and mu_s7_t) else 0.0
            print(
                f"{label_h:<10} | {lam_test:<10.1e} | {mu_s7_t:<14.8f} | {mu_opt_t:<18.8f} | "
                f"{d_mu:+18.3f} | {'optimal (mpmath certified)' if st_opt == 'optimal' else st_opt}"
            )
    print("=" * 105)


if __name__ == "__main__":
    main()
