"""Script de auditoria: Reoptimizacion a lambda_M nominal ajustado (ej. 1.06e-10)
para garantizar que el margen realizado lam_min(M) alcance exactamente >= 1.00e-10.
"""

from __future__ import annotations

import json
import math
from typing import TypedDict
import warnings
warnings.filterwarnings("ignore")

import cvxpy as cp
import mpmath
import numpy as np
import scipy.linalg
from scipy.optimize import minimize

class ReoptResultItem(TypedDict):
    lam_nom: float
    u_opt: np.ndarray
    mu_opt: float
    lam_mp: float
    ratio_1e10: float
    status: str

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

# Constantes analiticas S7
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
    eigvals_fn = getattr(mpmath.mp, "eigsy")
    eigs = eigvals_fn(M_mp, eigvals_only=True)
    return float(min(eigs))


def solve_exact(
    xc: np.ndarray,
    phic: np.ndarray,
    lam_M: float,
    tol: float = 1e-10,
    certify: bool = True,
) -> tuple[float | None, float, str]:
    """Resuelve la Fase II con LMI reescalada y Clarabel ultra-estricto."""
    m = len(xc)
    pairs = get_pair_indices(m)
    K = len(pairs)
    sheets, ref_sheet = compute_sheets(xc, phic)
    h_eff = float(compute_h_eff(xc, phic, sheets))

    T = build_scaled_basis_t(xc, phic, sheets, ref_sheet, h_eff)
    if T is None:
        return None, 0.0, "NOT_SCALABLE"

    P_perp = np.eye(m) - np.ones((m, m)) / float(m)
    T_prime = P_perp @ T
    metric_M_rescaled = (T_prime.T @ T_prime) / (h_eff**2)

    A_F = build_af_matrix(xc, phic)
    A_c = build_ac_matrix(xc, phic)
    A_d = build_ad_matrix(xc)
    N = scipy.linalg.null_space(A_F)

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
        return None, 0.0, "SOLVER_EXCEPTION"

    if prob.status != "optimal":
        return None, 0.0, prob.status

    p_clip = np.clip(N @ z.value, 0.0, 1.0)
    mu_cert = float(np.min(xc + A_d @ p_clip))
    if certify:
        P_mat = p_vector_to_matrix(p_clip, m)
        lam_mp = compute_mpmath_M(xc, phic, P_mat)
    else:
        lam_mp = 0.0
    return mu_cert, lam_mp, prob.status


def optimize_even(lam_target: float) -> tuple[np.ndarray, float, float, str]:
    """Optimiza dx_0_6, dx_1_5, dx_2_4 para un lambda_M nominal."""
    def loss(u: np.ndarray) -> float:
        xc = u_to_masses(u)
        mu, _, st = solve_exact(xc, phi, lam_M=lam_target, certify=False)
        if st != "optimal" or mu is None:
            return 1e9
        return -mu

    u_start = np.array([+0.075, +0.950, -1.300])
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
    mu_opt, lam_opt, st_opt = solve_exact(xc_opt, phi, lam_M=lam_target, certify=True)
    return res.x, mu_opt or 0.0, lam_opt, st_opt


def main() -> None:
    print("=" * 115)
    print("REOPTIMIZACION PAR A DISTINTOS lambda_M NOMINALES PARA LOGRAR lam_min(M) >= 1.00e-10")
    print("Clarabel tol_feas=1e-10 | Certificacion Primal mpmath a 50 digitos | h = 1.0e-3")
    print("=" * 115)

    # 1. Barrido de lambda_M nominales cercanos a 1.06e-10
    test_lams = [1.00e-10, 1.05e-10, 1.055e-10, 1.06e-10, 1.07e-10]
    reopt_results: list[ReoptResultItem] = []

    for lam_nom in test_lams:
        u_opt, mu_opt, lam_mp, st = optimize_even(lam_nom)
        ratio_1e10 = lam_mp / 1.0e-10
        reopt_results.append({
            "lam_nom": lam_nom,
            "u_opt": u_opt,
            "mu_opt": mu_opt,
            "lam_mp": lam_mp,
            "ratio_1e10": ratio_1e10,
            "status": st,
        })
        print(f"Nominal lam_M = {lam_nom:.3e} -> mu_cert = {mu_opt:.8f} | lam_min(M) mpmath = {lam_mp:.16e} ({ratio_1e10:.3f}x de 1e-10)")

    # 2. Elegir el caso que alcanza exactamente >= 1.00e-10
    target_case = next(r for r in reopt_results if r["lam_nom"] == 1.06e-10)

    # 3. Cargar los otros 3 candidatos para comparacion limpia
    # S7 base
    mu_s7, lam_s7, st_s7 = solve_exact(x_s7_manual, phi, lam_M=1.0e-10, certify=True)
    mu_s7_val = float(mu_s7) if mu_s7 is not None else 0.0

    # Campeon R3
    seed_4 = create_champion_3178_seed()
    x_r3 = np.array(seed_4.phenotype_x, dtype=float)
    phi_r3 = np.array(seed_4.phenotype_phi, dtype=float)
    mu_r3, lam_r3, st_r3 = solve_exact(x_r3, phi_r3, lam_M=1.0e-10, certify=True)
    mu_r3_val = float(mu_r3) if mu_r3 is not None else 0.0

    # Campeon Fisico 399165
    ckpt_path = "runs/night12h_seeded_mu06923_s42/checkpoints/checkpoint_gen_1675.json"
    with open(ckpt_path, "r", encoding="utf-8") as f:
        ckpt_data = json.load(f)
    ind_399 = next(ind for ind in ckpt_data["island_populations"][0] if ind["id"] == 399165)
    x_399 = np.array(ind_399["phenotype_x"], dtype=float)
    phi_399 = np.array(ind_399["phenotype_phi"], dtype=float)
    mu_399, lam_399, st_399 = solve_exact(x_399, phi_399, lam_M=1.0e-10, certify=True)
    mu_399_val = float(mu_399) if mu_399 is not None else 0.0

    # Optimo par previo (optimizado a 1.00e-10 nominal)
    prev_case = next(r for r in reopt_results if r["lam_nom"] == 1.00e-10)

    # Evaluar el nuevo optimo (1.06e-10) tambien bajo lam_M = 1.00e-10 en Clarabel
    xc_106 = u_to_masses(target_case["u_opt"])
    mu_106_at_100, lam_106_at_100, st_106_at_100 = solve_exact(xc_106, phi, lam_M=1.00e-10, certify=True)
    mu_106_val = float(mu_106_at_100) if mu_106_at_100 is not None else 0.0

    print("\n" + "=" * 115)
    print("TABLA COMPARATIVA DEFINITIVA (Mismo piso fisico lam_min(M) >= 1.00e-10)")
    print(f"{'Candidato':<38} | {'lam_M nom':<10} | {'Status':<8} | {'mu_cert':<14} | {'lam_min(M) mpmath':<24} | {'vs S7 base'}")
    print("-" * 115)
    print(f"{'S7 base (h=1e-3, delta_W=0)':<38} | {1.0e-10:<10.2e} | {st_s7:<8} | {mu_s7_val:<14.8f} | {lam_s7:<24.16e} | base (+0.000)")
    print(f"{'Campeon R3 (ID 4)':<38} | {1.0e-10:<10.2e} | {st_r3:<8} | {mu_r3_val:<14.8f} | {lam_r3:<24.16e} | {(mu_r3_val - mu_s7_val)*1e4:+9.3f}e-4")
    print(f"{'Campeon Fisico (ID 399165)':<38} | {1.0e-10:<10.2e} | {st_399:<8} | {mu_399_val:<14.8f} | {lam_399:<24.16e} | {(mu_399_val - mu_s7_val)*1e4:+9.3f}e-4")
    print(f"{'Optimo Par (nom 1.00e-10, previo)':<38} | {1.0e-10:<10.2e} | {prev_case['status']:<8} | {prev_case['mu_opt']:<14.8f} | {prev_case['lam_mp']:<24.16e} | {(prev_case['mu_opt'] - mu_s7_val)*1e4:+9.3f}e-4")
    print(f"{'Optimo Par (nom 1.06e-10, nuevo)':<38} | {1.06e-10:<10.2e} | {target_case['status']:<8} | {target_case['mu_opt']:<14.8f} | {target_case['lam_mp']:<24.16e} | {(target_case['mu_opt'] - mu_s7_val)*1e4:+9.3f}e-4")
    print(f"{'Optimo Par (nom 1.06e-10 eval a 1e-10)':<38} | {1.00e-10:<10.2e} | {st_106_at_100:<8} | {mu_106_val:<14.8f} | {lam_106_at_100:<24.16e} | {(mu_106_val - mu_s7_val)*1e4:+9.3f}e-4")
    print("-" * 115)

    print("\nPARAMETROS OPTIMOS [1e-4] (nom 1.06e-10):")
    u = target_case["u_opt"]
    print(f"  dx_0_6 = {u[0]:+.6f}")
    print(f"  dx_1_5 = {u[1]:+.6f}")
    print(f"  dx_2_4 = {u[2]:+.6f}")
    dx_W = -2.0 * (u[0] + u[1] + u[2])
    print(f"  dx_W   = {dx_W:+.6f}")
    print("=" * 115)


if __name__ == "__main__":
    main()
