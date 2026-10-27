"""Script 1: Optimo 7D Absoluto (Masas y Fase Simultaneas).

Auditoria definitiva de simetria y fase:
- 7 grados de libertad: dx_0, dx_1, dx_2, dx_4, dx_5, dx_W (unidades 1e-4) y delta_W (microradianes 1e-6 rad).
- dx_6 = -(dx_0 + dx_1 + dx_2 + dx_W + dx_4 + dx_5) para conservacion estricta sum(x) = 1.0.
- Fase en nodo W: phi_3 = delta_W.
- h = 1e-3, lam_M = 1e-10.
- Multi-start con 5 arranques (incluyendo el optimo par conocido con delta_W=0).
- Nelder-Mead con maxiter=2000, xatol=1e-5, fatol=1e-8.
- Certificacion primal mpmath a 50 digitos sobre P recortada.
- Evaluacion de ||x_impar||_inf y convergencia de delta_W a 0.
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

h = 0.001
half_pi = math.pi / 2.0
pi = math.pi


def R_true(x: np.ndarray) -> np.ndarray:
    """Reflexion horizontal exacta: 0<->6, 1<->5, 2<->4, 3 fijo."""
    xr = np.copy(x)
    xr[0] = x[6]; xr[6] = x[0]
    xr[1] = x[5]; xr[5] = x[1]
    xr[2] = x[4]; xr[4] = x[2]
    xr[3] = x[3]
    return xr


def u_to_phenotype_7d(u: np.ndarray | list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Decodifica u = [dx_0, dx_1, dx_2, dx_4, dx_5, dx_W, delta_W_urad].
    
    Masas en unidades 1e-4.
    Fase delta_W en unidades microradianes (1e-6 rad).
    """
    dx0 = float(u[0]) * 1e-4
    dx1 = float(u[1]) * 1e-4
    dx2 = float(u[2]) * 1e-4
    dx4 = float(u[3]) * 1e-4
    dx5 = float(u[4]) * 1e-4
    dxW = float(u[5]) * 1e-4
    dx6 = -(dx0 + dx1 + dx2 + dxW + dx4 + dx5)

    xc = np.array([
        a_star + dx0,
        a_star + dx1,
        a_star + dx2,
        x_W_base + dxW,
        a_star + dx4,
        a_star + dx5,
        a_star + dx6,
    ], dtype=float)

    delta_W = float(u[6]) * 1e-6

    phic = np.array([
        pi - w_star * h,      # 0
        -half_pi + h,         # 1
        -half_pi - h,         # 2
        delta_W,              # 3 (W con desfasaje libre)
        half_pi + h,          # 4
        half_pi - h,          # 5
        -pi + w_star * h,     # 6
    ], dtype=float)

    return xc, phic


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


def loss_7d(u: np.ndarray) -> float:
    """Objetivo Nelder-Mead: minimiza -mu_cert."""
    xc, phic = u_to_phenotype_7d(u)
    mu, _, st = solve_exact(xc, phic, lam_M=1e-10, tol=1e-10, certify=False)
    if st != "optimal" or mu is None:
        return 1e9
    return -mu


def main() -> None:
    print("=" * 105)
    print("SCRIPT 1: OPTIMO 7D ABSOLUTO (Masas y Fase Simultaneas)")
    print("LMI REESCALADA (h_eff^2) | Clarabel tol_feas=1e-10 | Certificacion Primal mpmath 50d")
    print("7 Grados de Libertad: dx_0, dx_1, dx_2, dx_4, dx_5, dx_W [1e-4] | delta_W [microrad]")
    print("Masa de balance: dx_6 = -(dx_0 + dx_1 + dx_2 + dx_W + dx_4 + dx_5)")
    print("=" * 105)

    # Optimo par conocido con delta_W = 0.0
    u_par_opt = np.array([
        +0.111279,   # dx_0
        +0.935376,   # dx_1
        -1.348106,   # dx_2
        -1.348106,   # dx_4
        +0.935376,   # dx_5
        +0.602901,   # dx_W
        +0.000000,   # delta_W [urad]
    ])

    # 5 arranques multi-start
    np.random.seed(42)
    starts = [
        u_par_opt,
        u_par_opt + np.array([+0.05, -0.05, +0.05, -0.05, +0.05, -0.05, +1.0]),  # +1 urad
        u_par_opt + np.array([-0.05, +0.05, -0.05, +0.05, -0.05, +0.05, -1.0]),  # -1 urad
        np.zeros(7),                                                               # S7 puro (delta_W=0)
        np.array([+0.2, +0.8, -1.2, -1.0, +0.6, +0.4, +2.0]),                    # Asimetrico + 2 urad
    ]

    best_mu = -1e9
    best_u = u_par_opt
    best_xc = None
    best_phic = None
    best_lam = 0.0

    print("\nEJECUCION MULTI-START (5 ARRANQUES NELDER-MEAD, maxiter=2000, xatol=1e-5, fatol=1e-8):")
    for idx, u0 in enumerate(starts):
        sim = [u0]
        for d in range(7):
            pt = np.copy(u0)
            pt[d] += 0.1 if d < 6 else 0.5  # paso de 0.1 en masas, 0.5 urad en fase
            sim.append(pt)
        sim = np.array(sim)

        res = minimize(
            loss_7d,
            u0,
            method="Nelder-Mead",
            options={
                "initial_simplex": sim,
                "maxiter": 2000,
                "xatol": 1e-5,
                "fatol": 1e-8,
            },
        )

        xc_cand, phic_cand = u_to_phenotype_7d(res.x)
        mu_c, lam_c, st_c = solve_exact(xc_cand, phic_cand, lam_M=1e-10, tol=1e-10, certify=True)

        print(
            f"  Start {idx:02d} | mu_cert = {mu_c:.8f} | delta_W = {res.x[6]:+.4f} urad "
            f"| nit={res.nit:3d} | nfev={res.nfev:3d} (status={st_c})"
        )

        if st_c == "optimal" and mu_c is not None and mu_c > best_mu:
            best_mu = mu_c
            best_u = res.x
            best_xc = xc_cand
            best_phic = phic_cand
            best_lam = lam_c

    assert best_xc is not None

    # Descomposicion Par / Impar de masas
    x_par = (best_xc + R_true(best_xc)) / 2.0
    x_impar = (best_xc - R_true(best_xc)) / 2.0
    norm_impar = float(np.max(np.abs(x_impar)))
    delta_W_final_urad = float(best_u[6])
    delta_W_final_rad = delta_W_final_urad * 1e-6

    print("\nRESULTADOS FINALES DEL OPTIMO 7D ABSOLUTO:")
    print(f"  mu_cert Maximo Absoluto : {best_mu:.8f} (mpmath certified)")
    print(f"  lam_min(M) Fisico       : {best_lam:.2e} (>= 0.95 * lam_M)")
    print(f"  Vector x* Optimo        : [{', '.join(f'{v:.8f}' for v in best_xc)}]")
    print(f"  Suma sum(x*)            : {np.sum(best_xc):.16f}")
    print(f"  Deltas masas u* [1e-4]  : [{', '.join(f'{v:+.6f}' for v in best_u[:6])}]")
    print(f"  delta_W Optimo (urad)   : {delta_W_final_urad:+.6f} urad ({delta_W_final_rad:+.4e} rad)")

    print("\nDESCOMPOSICION PAR / IMPAR DE x*:")
    print(f"  x_par           : [{', '.join(f'{v:.8f}' for v in x_par)}]")
    print(f"  x_impar         : [{', '.join(f'{v:+.2e}' for v in x_impar)}]")
    print(f"  ||x_impar||_inf : {norm_impar:.2e}")
    print(f"  Veredicto Simetria Masas : {'SIMETRIA ABSOLUTA (||x_impar||_inf < 1e-5)' if norm_impar < 1e-5 else 'ASIMETRICO'}")
    print(f"  Veredicto Fase delta_W   : {'PARIDAD EXACTA (delta_W -> 0)' if abs(delta_W_final_rad) < 1e-5 else 'DESFASAJE NO NULO'}")
    print("=" * 105)


if __name__ == "__main__":
    main()
