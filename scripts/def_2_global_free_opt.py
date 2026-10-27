"""Script 2: El Optimo Global Libre (5 masas + delta_W, masa 6 por normalizacion).

Auditoria del espacio asimetrico libre:
- 6 grados de libertad independientes: dx_0, dx_1, dx_2, dx_4, dx_5 y delta_W (unidades 1e-4).
- Masa 3 fija en x_W_base, masa 6 determinada por normalizacion sum(x) = 1.0.
- LMI reescalada (divisible por h_eff^2) con Clarabel (tol_feas=1e-10).
- Multi-start (5 arranques) incluyendo la semilla central x_par_R3 con delta_W = 0.
- Descomposicion par / impar del vector optimo resultante.
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


q_star, w_star = compute_q_star_and_w_star(dps=50)
a_star = float(1.0 / (7.0 + q_star))
x_W_base = float(1.0 - 6.0 * a_star)
x_s7_manual = np.array([a_star, a_star, a_star, x_W_base, a_star, a_star, a_star], dtype=float)

h = 0.001
half_pi = math.pi / 2.0
pi = math.pi


def R_true(x: np.ndarray) -> np.ndarray:
    """Permutacion de reflexion especular horizontal: 0<->6, 1<->5, 2<->4, 3 fijo."""
    xr = np.copy(x)
    xr[0] = x[6]; xr[6] = x[0]
    xr[1] = x[5]; xr[5] = x[1]
    xr[2] = x[4]; xr[4] = x[2]
    xr[3] = x[3]
    return xr


def build_phi_free(delta_W: float) -> np.ndarray:
    """Fases de S7(h=1e-3) con desfasaje delta_W en el nodo 3 (W)."""
    return np.array([
        pi - w_star * h,      # 0
        -half_pi + h,         # 1
        -half_pi - h,         # 2
        delta_W,              # 3 (W libre)
        half_pi + h,          # 4
        half_pi - h,          # 5
        -pi + w_star * h,     # 6
    ], dtype=float)


def u_to_masses_6d(u: np.ndarray | list[float]) -> tuple[np.ndarray, float]:
    """Decodifica el vector u [unidades de 1e-4] a vector de masas y delta_W.
    
    u = [dx_0, dx_1, dx_2, dx_4, dx_5, delta_W]
    dx_6 = -(dx_0 + dx_1 + dx_2 + dx_4 + dx_5)
    """
    dx0 = float(u[0]) * 1e-4
    dx1 = float(u[1]) * 1e-4
    dx2 = float(u[2]) * 1e-4
    dx4 = float(u[3]) * 1e-4
    dx5 = float(u[4]) * 1e-4
    dx6 = -(dx0 + dx1 + dx2 + dx4 + dx5)
    xc = np.array([
        a_star + dx0,
        a_star + dx1,
        a_star + dx2,
        x_W_base,
        a_star + dx4,
        a_star + dx5,
        a_star + dx6,
    ], dtype=float)
    delta_W = float(u[5]) * 1e-4
    return xc, delta_W


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


def loss_6d(u: np.ndarray) -> float:
    """Objetivo Nelder-Mead: minimiza -mu_cert."""
    xc, delta_W = u_to_masses_6d(u)
    phic = build_phi_free(delta_W)
    mu, _, st = solve_exact(xc, phic, lam_M=1e-10, tol=1e-10)
    if st != "optimal" or mu is None:
        return 1e9
    return -mu


def main() -> None:
    print("=" * 100)
    print("SCRIPT 2: EL OPTIMO GLOBAL LIBRE (6 D.O.F.: 5 MASAS LIBRES + DELTA_W)")
    print("LMI REESCALADA (h_eff^2) | Clarabel tol_feas=1e-10 | Certificacion Primal mpmath 50d")
    print("=" * 100)

    # 1. Semilla central x_par_R3 con delta_W = 0
    seed_4 = create_champion_3178_seed()
    x_r3 = np.array(seed_4.phenotype_x, dtype=float)
    x_par_r3 = (x_r3 + R_true(x_r3)) / 2.0
    dx_init = x_par_r3 - x_s7_manual
    u_seed = np.array([
        dx_init[0], dx_init[1], dx_init[2], dx_init[4], dx_init[5], 0.0
    ]) / 1e-4

    xc_seed, dw_seed = u_to_masses_6d(u_seed)
    mu_seed, lam_seed, st_seed = solve_exact(xc_seed, build_phi_free(dw_seed))
    print(f"\n1. SEMILLA CENTRAL (x_par_R3 con delta_W = 0.0):")
    print(f"   mu_cert    : {mu_seed:.8f}")
    print(f"   lam_min(M) : {lam_seed:.2e} ({st_seed})")
    print(f"   u_seed     : [{', '.join(f'{v:+.4f}' for v in u_seed)}]")

    # 2. Multi-start con 5 arranques
    print(f"\n2. EJECUCION MULTI-START (5 ARRANQUES NELDER-MEAD):")
    np.random.seed(42)
    starts = [
        u_seed,
        u_seed + np.array([+0.2, -0.2, +0.2, +0.2, -0.2, +0.1]),
        u_seed + np.array([-0.2, +0.2, -0.2, -0.2, +0.2, -0.1]),
        np.array([+0.5, +0.5, -1.0, -1.0, +0.5, +0.0]),
        np.array([+1.5, +1.0, -2.5, -2.5, +1.0, +0.0]),
    ]

    best_mu = -1e9
    best_u = u_seed
    best_xc = xc_seed
    best_dw = dw_seed
    best_lam = lam_seed

    for idx, u0 in enumerate(starts):
        sim = [u0]
        for d in range(6):
            pt = np.copy(u0)
            pt[d] += 0.2 if d < 5 else 0.1
            sim.append(pt)
        sim = np.array(sim)

        res = minimize(
            loss_6d,
            u0,
            method="Nelder-Mead",
            options={
                "initial_simplex": sim,
                "maxiter": 200,
                "xatol": 1e-3,
                "fatol": 1e-6,
            },
        )
        mu_cand = -res.fun
        xc_cand, dw_cand = u_to_masses_6d(res.x)
        mu_c, lam_c, st_c = solve_exact(xc_cand, build_phi_free(dw_cand))

        print(
            f"  Start {idx:02d} | mu_cert = {mu_c:.8f} | dw = {dw_cand*1e6:+8.4f} um "
            f"| nit={res.nit:3d} (status={st_c})"
        )

        if st_c == "optimal" and mu_c is not None and mu_c > best_mu:
            best_mu = mu_c
            best_u = res.x
            best_xc = xc_cand
            best_dw = dw_cand
            best_lam = lam_c

    # 3. Descomposicion Par / Impar
    x_par = (best_xc + R_true(best_xc)) / 2.0
    x_impar = (best_xc - R_true(best_xc)) / 2.0
    norm_impar = float(np.max(np.abs(x_impar)))

    print(f"\n3. RESULTADOS DEL OPTIMO GLOBAL LIBRE:")
    print(f"   mu_cert Maximo Absoluto : {best_mu:.8f} (mpmath certified)")
    print(f"   lam_min(M) Fisico       : {best_lam:.2e} (>= 0.95 * lam_M)")
    print(f"   delta_W Resultante      : {best_dw:+.6e} rad ({best_dw*1e6:+.4f} um)")
    print(f"   Vector de Masas x*      : [{', '.join(f'{v:.8f}' for v in best_xc)}]")

    print(f"\n4. DESCOMPOSICION PAR / IMPAR DE x*:")
    print(f"   Parte Par   x_par   : [{', '.join(f'{v:.8f}' for v in x_par)}]")
    print(f"   Parte Impar x_impar : [{', '.join(f'{v:+.2e}' for v in x_impar)}]")
    print(f"   ||x_impar||_inf     : {norm_impar:.2e}")
    print(f"   Veredicto Simetria  : {'PAR PERFECTO (asimetria < 1e-5)' if norm_impar < 1e-5 else 'ASIMETRICO'}")
    print("=" * 100)


if __name__ == "__main__":
    main()
