"""Comprehensive technical audit and cheap experiments requested by user:

1. Permutation and gauge-aligned Delta x and Delta phi of R1 and R3 against S7(1e-3).
2. Exact check of A_F * p in 50 digits and strict check of 0 <= P <= 1 bounds.
3. Exact LMI with metric tensor T'^T T':
       T^T L_w(p) T >= lambda_M (T'^T T')
   and comparison table against simplified LMI.
4. Block eigenvalues of G_hat (eta block vs zeta block) at the optimum.
5. Nelder-Mead envelope optimization of S7 over (h, t, kappa) at each lambda_M.
6. Evaluation with seed lock at 2e-3.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import mpmath
import numpy as np
import scipy.linalg  # type: ignore # pyright: ignore[reportMissingTypeStubs]
import scipy.optimize  # type: ignore # pyright: ignore[reportMissingTypeStubs]

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.domain.conditions import (
    build_ac_matrix,
    build_ad_matrix,
    build_af_matrix,
    build_lw_matrix,
    get_pair_indices,
    p_matrix_to_vector,
)
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)
from kuramoto_templates.search.genome import wrap_minus_pi_to_pi

# Load report
with open("results/overnight_12h_verified_report.json", "r", encoding="utf-8") as f:
    report = json.load(f)

r1 = report["runs"][0]
r3 = report["runs"][2]

# ==============================================================================
# 1. PERMUTATION & GAUGE-ALIGNED DELTA X AND DELTA PHI
# ==============================================================================
print("=" * 80)
print("1. VOLCADO DE DELTA x Y DELTA phi DE R1 Y R3 CONTRA S7(1e-3)")
print("=" * 80)

s7_ref = get_s7_appendix_a_template(h=0.001, t=0.0)
xs7 = np.array(s7_ref.x)
phis7 = np.array(s7_ref.phi)
cluster_labels = ["Am", "Ap", "Bm", "Bp", "Cm", "Cp", "W"]


def align_to_s7(
    x_cand: tuple[float, ...] | list[float] | np.ndarray,
    phi_cand: tuple[float, ...] | list[float] | np.ndarray,
    xs7: np.ndarray,
    phis7: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    m = len(x_cand)
    x_c = np.array(x_cand, dtype=float)
    phi_c = np.array(phi_cand, dtype=float)
    best_dist = float("inf")

    # Test all phase shifts that match cluster W (largest mass)
    w_cand_idx = int(np.argmax(x_c))
    alpha = wrap_minus_pi_to_pi(phis7[6] - phi_c[w_cand_idx])
    phi_rot = np.array([wrap_minus_pi_to_pi(p + alpha) for p in phi_c])

    # Cost matrix for assignment
    cost = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            cost[i, j] = 10.0 * (x_c[i] - xs7[j])**2 + (wrap_minus_pi_to_pi(phi_rot[i] - phis7[j]))**2
    row_ind, col_ind = scipy.optimize.linear_sum_assignment(cost)

    x_aligned = np.zeros(m)
    phi_aligned = np.zeros(m)
    for r, c in zip(row_ind, col_ind):
        x_aligned[c] = x_c[r]
        phi_aligned[c] = phi_rot[r]
    return x_aligned, phi_aligned

x1_al, phi1_al = align_to_s7(r1["phenotype_x"], r1["phenotype_phi"], xs7, phis7)
x3_al, phi3_al = align_to_s7(r3["phenotype_x"], r3["phenotype_phi"], xs7, phis7)

dx1 = x1_al - xs7
dphi1 = np.array([wrap_minus_pi_to_pi(p1 - ps) for p1, ps in zip(phi1_al, phis7)])
dx3 = x3_al - xs7
dphi3 = np.array([wrap_minus_pi_to_pi(p3 - ps) for p3, ps in zip(phi3_al, phis7)])

print("\n--- Campeón R1 (ID 125678) alineado contra S7(1e-3) ---")
print(f"{'Cluster':<8} | {'x_S7':<12} | {'x_R1':<12} | {'Delta x':<14} | {'phi_S7':<12} | {'phi_R1':<12} | {'Delta phi (rad)':<16}")
print("-" * 96)
for i in range(7):
    print(f"{cluster_labels[i]:<8} | {xs7[i]:<12.6f} | {x1_al[i]:<12.6f} | {dx1[i]:<+14.4e} | {phis7[i]:<12.6f} | {phi1_al[i]:<12.6f} | {dphi1[i]:<+16.4e}")

print(f"\nNorma L2: ||Delta x|| = {np.linalg.norm(dx1):.4e}, ||Delta phi|| = {np.linalg.norm(dphi1):.4e}")

print("\n--- Campeón R3 (ID 4) alineado contra S7(1e-3) ---")
print(f"{'Cluster':<8} | {'x_S7':<12} | {'x_R3':<12} | {'Delta x':<14} | {'phi_S7':<12} | {'phi_R3':<12} | {'Delta phi (rad)':<16}")
print("-" * 96)
for i in range(7):
    print(f"{cluster_labels[i]:<8} | {xs7[i]:<12.6f} | {x3_al[i]:<12.6f} | {dx3[i]:<+14.4e} | {phis7[i]:<12.6f} | {phi3_al[i]:<12.6f} | {dphi3[i]:<+16.4e}")
print(f"\nNorma L2: ||Delta x|| = {np.linalg.norm(dx3):.4e}, ||Delta phi|| = {np.linalg.norm(dphi3):.4e}")


# ==============================================================================
# 2. EXACT MP-MATH 50-DIGIT CHECK OF A_F * p AND STRICT 0 <= P <= 1 BOUNDS
# ==============================================================================
print("\n" + "=" * 80)
print("2. VERIFICACION ESTRICTA: A_F * p EN 50 DIGITOS Y COTAS 0 <= P <= 1")
print("=" * 80)

evaluator = ScaledSDPEvaluator(
    EvalConfig(solver_name="Clarabel", solver_tol=1e-9, h_eff_min=1e-4),
    h_eff_min=1e-4,
)

cands_to_check = [
    ("Campeon R1 (ID 125678)", tuple(r1["phenotype_x"]), tuple(r1["phenotype_phi"])),
    ("Campeon R2 (ID 2)", tuple(report["runs"][1]["phenotype_x"]), tuple(report["runs"][1]["phenotype_phi"])),
    ("Campeon R3 (ID 4)", tuple(r3["phenotype_x"]), tuple(r3["phenotype_phi"])),
    ("S7(1e-3)", s7_ref.x, s7_ref.phi),
]

mpmath.mp.dps = 50

for cname, xc, phic in cands_to_check:
    res = evaluator.evaluate_phase1(xc, phic)
    P_mat = np.array(res.P, dtype=float)
    min_p = float(np.min(P_mat))
    max_p = float(np.max(P_mat))
    exceed_upper = max(0.0, max_p - 1.0)
    exceed_lower = max(0.0, -min_p)

    # 50-digit mpmath calculation of A_F * p
    m = len(xc)
    pairs = get_pair_indices(m)
    K = len(pairs)
    p_vec = p_matrix_to_vector(P_mat)

    # mpmath matrices
    F_mp = [mpmath.mp.mpf(0) for _ in range(m)]
    for a in range(m):
        for k_idx, (i, j) in enumerate(pairs):
            if a == i:
                diff = mpmath.mp.mpf(str(phic[i])) - mpmath.mp.mpf(str(phic[j]))
                term = mpmath.mp.mpf(str(xc[j])) * mpmath.sin(diff) * mpmath.mp.mpf(str(p_vec[k_idx]))
                F_mp[a] += term
            elif a == j:
                diff = mpmath.mp.mpf(str(phic[j])) - mpmath.mp.mpf(str(phic[i]))
                term = mpmath.mp.mpf(str(xc[i])) * mpmath.sin(diff) * mpmath.mp.mpf(str(p_vec[k_idx]))
                F_mp[a] += term

    max_F_mp = float(max(abs(f_val) for f_val in F_mp))

    print(f"\n--- {cname} ---")
    print(f"  min(P): {min_p:.15f} | Violacion inferior (<0): {exceed_lower:.2e}")
    print(f"  max(P): {max_p:.15f} | Violacion superior (>1): {exceed_upper:.2e}")
    print(f"  max|A_F * p| en mpmath (50 digitos): {max_F_mp:.2e}")
    print(f"  Cotas estrictas en [0, 1] dentro de 1e-9: {exceed_upper <= 1e-9 and exceed_lower <= 1e-9}")


# ==============================================================================
# 3. EXACT GENERALIZED LMI WITH T'^T T' (METRIC TENSOR)
# ==============================================================================
print("\n" + "=" * 80)
print("3. LMI EXACTA CON METRICA T'^T T' (RAYLEIGH QUOTIENT EXACTO)")
print("=" * 80)

def solve_phase2_exact_lmi(
    x: tuple[float, ...] | list[float] | np.ndarray,
    phi: tuple[float, ...] | list[float] | np.ndarray,
    lam_M: float,
) -> tuple[float, float]:
    x_arr = np.array(x, dtype=float)
    phi_arr = np.array(phi, dtype=float)
    m = len(x_arr)
    pairs = get_pair_indices(m)
    K = len(pairs)
    sheets, ref_sheet = compute_sheets(x_arr, phi_arr)
    heff = float(compute_h_eff(x_arr, phi_arr, sheets))
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, heff)
    if T is None:
        return float("nan"), float("nan")

    # T' = P_perp T = (I - 1/m 1 1^T) T
    P_perp = np.eye(m) - np.ones((m, m)) / float(m)
    T_prime = P_perp @ T
    metric_M = T_prime.T @ T_prime  # (m-1) x (m-1)

    A_F = build_af_matrix(x_arr, phi_arr)
    A_c = build_ac_matrix(x_arr, phi_arr)
    A_d = build_ad_matrix(x_arr)

    N = scipy.linalg.null_space(A_F)
    res_inf = float(np.max(np.abs(A_F @ N)))
    if res_inf > 1e-14:
        corr, _, _, _ = np.linalg.lstsq(A_F, A_F @ N, rcond=None)
        N = N - corr
        N, _ = np.linalg.qr(N)

    d = N.shape[1]
    z = cp.Variable(d)
    t_var = cp.Variable()
    p_z = N @ z

    # Build T^T L_w T basis (without division by heff^2)
    G_raw_basis = []
    for a, b in pairs:
        d_ab = np.zeros(m, dtype=float)
        d_ab[a] = 1.0
        d_ab[b] = -1.0
        Td = T.T @ d_ab
        weight = x_arr[a] * x_arr[b] * math.cos(phi_arr[a] - phi_arr[b])
        Gk = weight * np.outer(Td, Td)
        Gk = 0.5 * (Gk + Gk.T)
        G_raw_basis.append(Gk)

    L_T_expr = sum(p_z[k] * G_raw_basis[k] for k in range(K))

    # EXACT LMI: T^T L_w(p) T >= lam_M * (T'^T T')
    constraints = [
        x_arr + A_d @ p_z >= t_var,
        p_z >= 0.0,
        p_z <= 1.0,
        L_T_expr >> lam_M * metric_M,
        x_arr + A_c @ p_z >= 1e-4,
    ]

    prob = cp.Problem(cp.Maximize(t_var), constraints)
    prob.solve(
        solver=cp.CLARABEL,
        tol_gap_abs=1e-9,
        tol_gap_rel=1e-9,
        tol_feas=1e-9,
        tol_infeas_abs=1e-9,
        tol_infeas_rel=1e-9,
        verbose=False,
    )
    mu_val = float(t_var.value) if t_var.value is not None else float("nan")

    # Realized lambda_min(M)
    if p_z.value is not None:
        p_opt = np.clip(N @ z.value, 0.0, 1.0)
        L_w = build_lw_matrix(x_arr, phi_arr, p_opt)
        Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
        M = Q.T @ L_w @ Q
        real_lam = float(np.min(np.linalg.eigvalsh(M)))
    else:
        real_lam = float("nan")

    return mu_val, real_lam

print("\nComparación: LMI Simplificada vs LMI Exacta con Métrica T'^T T'")
lam_test_targets = [1e-11, 1e-10, 1e-9]

print(f"{'Candidato':<26} | {'lambda_M':<10} | {'mu* (LMI Aprox)':<16} | {'mu* (LMI Exacta)':<16} | {'lambda_M Realizado':<18}")
print("-" * 95)
for cname, xc, phic in cands_to_check:
    heff = float(compute_h_eff(xc, phic, compute_sheets(xc, phic)[0]))
    for lm in lam_test_targets:
        res_aprox = evaluator.evaluate_phase2(xc, phic, eps=lm / (heff**2))
        mu_aprox = res_aprox.mu
        mu_exact, real_lam = solve_phase2_exact_lmi(xc, phic, lm)
        print(f"{cname:<26} | {lm:<10.0e} | {mu_aprox:<16.8f} | {mu_exact:<16.8f} | {real_lam:<18.4e}")


# ==============================================================================
# 4. AUTOVALORES DE G_HAT POR BLOQUE (ETA Y ZETA)
# ==============================================================================
print("\n" + "=" * 80)
print("4. AUTOVALORES DE G_HAT POR BLOQUE (eta: macroscopico, zeta: intra-hoja)")
print("=" * 80)

for cname, xc, phic in cands_to_check:
    res = evaluator.evaluate_phase2(xc, phic, eps=1e-10 / (heff**2))
    p_opt = p_matrix_to_vector(np.array(res.P, dtype=float))
    x_arr = np.array(xc, dtype=float)
    phi_arr = np.array(phic, dtype=float)
    m = len(x_arr)
    sheets, ref_sheet = compute_sheets(x_arr, phi_arr)
    heff = float(compute_h_eff(x_arr, phi_arr, sheets))
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, heff)
    if T is None:
        continue

    # Number of eta columns vs zeta columns
    n_eta = len(sheets) - 1  # 4 sheets - 1 = 3
    n_zeta = (m - 1) - n_eta  # 6 - 3 = 3

    L_w = build_lw_matrix(x_arr, phi_arr, p_opt)
    G_hat = (T.T @ L_w @ T) / (heff**2)

    G_eta = G_hat[:n_eta, :n_eta]
    G_zeta = G_hat[n_eta:, n_eta:]
    G_cross = G_hat[:n_eta, n_eta:]

    eigs_eta = np.linalg.eigvalsh(G_eta)
    eigs_zeta = np.linalg.eigvalsh(G_zeta)
    eigs_full = np.linalg.eigvalsh(G_hat)

    print(f"\n--- {cname} (G_hat dim: 6x6) ---")
    print(f"  Bloque eta (inter-hoja 3x3) autovalores : {eigs_eta}")
    print(f"  Bloque zeta (intra-hoja 3x3) autovalores: {eigs_zeta}")
    print(f"  Acoplamiento inter-intra ||G_cross||_F  : {np.linalg.norm(G_cross):.4e}")
    print(f"  Autovalores globales de G_hat          : {eigs_full}")


# ==============================================================================
# 5. NELDER-MEAD ENVELOPE OVER (h, t, kappa) FOR S7
# ==============================================================================
print("\n" + "=" * 80)
print("5. ENVOLVENTE DE S7 OPTIMIZANDO (h, t, kappa) CON NELDER-MEAD")
print("=" * 80)

def evaluate_s7_mu(
    params: list[float] | np.ndarray,
    target_lam_M: float,
) -> float:
    h_v, t_v, kap_v = params
    if h_v < 1e-4 or h_v > 5e-3 or t_v < -0.01 or t_v > 0.05 or kap_v < 0.5 or kap_v > 5.0:
        return 1e6
    try:
        tmpl = get_s7_appendix_a_template(h=h_v, t=t_v, kappa=kap_v)
        mu_val, _ = solve_phase2_exact_lmi(tmpl.x, tmpl.phi, target_lam_M)
        if math.isnan(mu_val):
            return 1e6
        return -mu_val
    except Exception:
        return 1e6

print(f"{'Target lambda_M':<16} | {'mu* S7 Defecto':<16} | {'mu* S7 Envolvente':<18} | {'(h*, t*, kappa*) Optimos':<30} | {'Campeon R1':<14}")
print("-" * 100)

for lm in [1e-11, 1e-10, 1e-9]:
    # Defecto S7
    tmpl_def = get_s7_appendix_a_template(h=0.002, t=0.0, kappa=1.9)
    mu_def, _ = solve_phase2_exact_lmi(tmpl_def.x, tmpl_def.phi, lm)

    # Nelder-Mead optimization
    init_guess = [0.001, 0.0, 1.9]
    res_nm = scipy.optimize.minimize(
        evaluate_s7_mu,
        init_guess,
        args=(lm,),
        method="Nelder-Mead",
        options={"maxiter": 60, "xatol": 1e-4, "fatol": 1e-5},
    )
    mu_env = -res_nm.fun
    h_opt, t_opt, k_opt = res_nm.x

    # R1
    mu_r1, _ = solve_phase2_exact_lmi(r1["phenotype_x"], r1["phenotype_phi"], lm)
    print(f"{lm:<16.0e} | {mu_def:<16.8f} | {mu_env:<18.8f} | h={h_opt:.4e}, t={t_opt:.2e}, k={k_opt:.2f} | {mu_r1:<14.8f}")
