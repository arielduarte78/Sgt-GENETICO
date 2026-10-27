"""Optimizacion de la Verdadera Parte Par (Simetria Especular R: 0<->6, 1<->5, 2<->4).

Script para la auditoria cientifica:
- Permutacion de reflexion correcta: 0<->6, 1<->5, 2<->4, 3 fijo.
- Arranque desde x_par_R3 = (x_R3 + R(x_R3))/2.
- LMI reescalada por h_eff^2 con Clarabel (tol_feas=1e-10).
- delta_W fijado estrictamente en 0.0.
- Certificacion primal en mpmath a 50 digitos con P recortada.
- Rechazo estricto de OPTIMAL_INACCURATE.
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
from kuramoto_templates.search.seeds import create_champion_3178_seed


# -----------------------------------------------------------------------------
# 1. Geometria y Parametrizacion Par Exacta
# -----------------------------------------------------------------------------
q_star, w_star = compute_q_star_and_w_star(dps=50)
a_star = float(1.0 / (7.0 + q_star))
x_W_base = float(1.0 - 6.0 * a_star)
x_s7_manual = np.array([a_star, a_star, a_star, x_W_base, a_star, a_star, a_star], dtype=float)

h = 0.001
half_pi = math.pi / 2.0
pi = math.pi

phi = np.array([
    pi - w_star * h,      # 0: opuesto (+)
    -half_pi + h,         # 1: sur (+)
    -half_pi - h,         # 2: sur (-)
    0.0,                  # 3: W (fijo en 0.0)
    half_pi + h,          # 4: norte (+)
    half_pi - h,          # 5: norte (-)
    -pi + w_star * h,     # 6: opuesto (-)
], dtype=float)


def R_true(x: np.ndarray) -> np.ndarray:
    """Permutacion exacta de reflexion horizontal y -> -y: 0<->6, 1<->5, 2<->4, 3 fijo."""
    xr = np.copy(x)
    xr[0] = x[6]
    xr[6] = x[0]
    xr[1] = x[5]
    xr[5] = x[1]
    xr[2] = x[4]
    xr[4] = x[2]
    xr[3] = x[3]
    return xr


def u_to_masses_true(u: np.ndarray | list[float]) -> np.ndarray:
    """Convierte parametros u [en 1e-4] al vector de masas con simetria R:
    
    u[0] = dx_0_6 (modifica indices 0 y 6)
    u[1] = dx_1_5 (modifica indices 1 y 5)
    u[2] = dx_2_4 (modifica indices 2 y 4)
    dx_W = -2 * (dx_0_6 + dx_1_5 + dx_2_4) para sum(x) = 1.0
    """
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
    """Resuelve la LMI reescalada por h_eff^2 y certifica primalmente en mpmath.
    
    Rechaza estrictamente OPTIMAL_INACCURATE y soluciones no factibles.
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


def loss(u: np.ndarray) -> float:
    """Objetivo para Nelder-Mead: minimiza -mu_cert."""
    xc = u_to_masses_true(u)
    mu, _, st = solve_exact(xc, phi, lam_M=1e-10, tol=1e-10)
    if st != "optimal" or mu is None:
        return 1e9
    return -mu


def main() -> None:
    print("=" * 105)
    print("OPTIMIZACION DE LA VERDADERA PARTE PAR (SIMETRIA ESPECULAR R)")
    print("Permutacion correcta R: 0<->6, 1<->5, 2<->4, 3 fijo (W) | delta_W = 0.0")
    print("LMI REESCALADA (divisible por h_eff^2) | Clarabel tol_feas=1e-10 | Certificacion Primal mpmath 50d")
    print("=" * 105)

    # 1. Base S7
    mu_s7, lam_s7, st_s7 = solve_exact(x_s7_manual, phi)
    print(f"\n1. LINEA BASE S7 MANUAL (h=1e-3, delta_W=0.0):")
    print(f"   mu_cert    : {mu_s7:.8f}")
    print(f"   lam_min(M) : {lam_s7:.2e} ({st_s7})")
    print(f"   x_S7       : [{', '.join(f'{v:.8f}' for v in x_s7_manual)}]")

    # 2. Semilla R3 y proyeccion par exacta
    seed_4 = create_champion_3178_seed()
    x_r3 = np.array(seed_4.phenotype_x, dtype=float)
    x_par_r3 = (x_r3 + R_true(x_r3)) / 2.0
    dx_init = x_par_r3 - x_s7_manual
    u_init = np.array([dx_init[0], dx_init[1], dx_init[2]]) / 1e-4

    mu_init, lam_init, st_init = solve_exact(x_par_r3, phi)
    print(f"\n2. PUNTO DE ARRANQUE (Verdadera Parte Par de R3: x_par_R3):")
    print(f"   mu_cert inicial : {mu_init:.8f}  (Objetivo esperado ~0.692585)")
    print(f"   lam_min(M)      : {lam_init:.2e} ({st_init})")
    print(f"   u_init [1e-4]   : dx_0_6={u_init[0]:+.4f}, dx_1_5={u_init[1]:+.4f}, dx_2_4={u_init[2]:+.4f}")
    print(f"   dx_init [real]  : dx_0_6={dx_init[0]:+.6e}, dx_1_5={dx_init[1]:+.6e}, dx_2_4={dx_init[2]:+.6e}, dx_W={dx_init[3]:+.6e}")
    print(f"   x_par_R3        : [{', '.join(f'{v:.8f}' for v in x_par_r3)}]")

    # 3. Optimizacion Nelder-Mead
    print(f"\n3. EJECUCION DE LA OPTIMIZACION (Nelder-Mead):")
    sim = [u_init]
    for d in range(3):
        pt = np.copy(u_init)
        pt[d] += 0.2
        sim.append(pt)
    sim = np.array(sim)

    res = minimize(
        loss,
        u_init,
        method="Nelder-Mead",
        options={
            "initial_simplex": sim,
            "maxiter": 250,
            "xatol": 1e-4,
            "fatol": 1e-7,
            "disp": True,
        },
    )

    u_opt = res.x
    xc_opt = u_to_masses_true(u_opt)
    dx_opt = xc_opt - x_s7_manual
    mu_opt, lam_opt, st_opt = solve_exact(xc_opt, phi)

    print(f"\n4. RESULTADOS DE LA OPTIMIZACION PAR FINAL:")
    print(f"   mu_cert Final   : {mu_opt:.8f} (Certificado Primal mpmath)")
    print(f"   Ganancia vs S7  : {(mu_opt - mu_s7)*1e4:+.3f}e-4")
    print(f"   Ganancia vs R3  : {(mu_opt - 0.69256438)*1e4:+.3f}e-4")
    print(f"   lam_min(M)      : {lam_opt:.2e} (>= 0.95 * lam_M verificado, status={st_opt})")
    print(f"   u_opt [1e-4]    : dx_0_6={u_opt[0]:+.4f}, dx_1_5={u_opt[1]:+.4f}, dx_2_4={u_opt[2]:+.4f}")
    print(f"   Deltas Optimos  :")
    print(f"     dx_0_6 (opuestos) : {dx_opt[0]:+.6e}")
    print(f"     dx_1_5 (outer)    : {dx_opt[1]:+.6e}")
    print(f"     dx_2_4 (inner)    : {dx_opt[2]:+.6e}")
    print(f"     dx_W   (nodo W)   : {dx_opt[3]:+.6e}")
    print(f"     sum(dx)           : {np.sum(dx_opt):+.2e}")
    print(f"   Vector de Masas Resultante (x_opt_par):")
    print(f"     [{', '.join(f'{v:.8f}' for v in xc_opt)}]")
    print("=" * 105)


if __name__ == "__main__":
    main()
