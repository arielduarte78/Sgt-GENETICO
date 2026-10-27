"""Comprehensive Audit Suite addressing user experiments in exact priority order:

1. Primal certification of R1, R3, S7(1e-3), and Polished Point:
   - mu_cert (min degree) vs realized lambda_min(M) in mpmath 50d.
2. Slack diagnosis:
   - Clarabel solver status, residuals, and lambda_min of slack matrix:
     S_slack = T^T L_w(p) T - lambda_M (T'^T T').
3. Multi-start local polish (20 random starts) with:
   - Joint optimization: (Delta x, Delta phi)
   - Ablation 1: Delta x only (Delta phi = 0)
   - Ablation 2: Delta phi only (Delta x = 0)
4. Sweep of mu*(lambda_M) with realized lambda_min(M) on the table.
"""

from __future__ import annotations

import json
import math
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
    p_vector_to_matrix,
)
from kuramoto_templates.evaluators.scaled import (
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)

# Load reports
with open("results/overnight_12h_verified_report.json", "r", encoding="utf-8") as f:
    report = json.load(f)

r1 = report["runs"][0]
r3 = report["runs"][2]


def fast_solve_phase2(
    x_arr: np.ndarray,
    phi_arr: np.ndarray,
    lam_M: float,
) -> tuple[float, np.ndarray | None, np.ndarray | None, np.ndarray | None, str]:
    """Fast Phase II evaluation using Clarabel (pure NumPy/CVXPY, no mpmath)."""
    m = len(x_arr)
    pairs = get_pair_indices(m)
    K = len(pairs)
    sheets, ref_sheet = compute_sheets(x_arr, phi_arr)
    heff = float(compute_h_eff(x_arr, phi_arr, sheets))
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, heff)
    if T is None:
        return float("nan"), None, None, None, "NO_LEAF_BASIS"

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

    constraints = [
        x_arr + A_d @ p_z >= t_var,
        p_z >= 0.0,
        p_z <= 1.0,
        L_T_expr >> lam_M * metric_M,
        x_arr + A_c @ p_z >= 1e-4,
    ]

    prob = cp.Problem(cp.Maximize(t_var), constraints)
    try:
        prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-9, tol_gap_rel=1e-9, verbose=False)
        if prob.status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            p_val = np.clip(p_z.value, 0.0, 1.0)
            return float(prob.value), p_val, T, metric_M, prob.status
        return float("nan"), None, T, metric_M, prob.status
    except Exception as e:
        return float("nan"), None, None, None, f"EXCEPTION: {e}"


def certify_primal_50d(
    x_arr: np.ndarray,
    phi_arr: np.ndarray,
    p_vec: np.ndarray,
    lam_M: float,
    T: np.ndarray,
    metric_M: np.ndarray,
) -> dict[str, Any]:
    """Rigorous primal certification in mpmath 50 digits."""
    m = len(x_arr)
    A_F = build_af_matrix(x_arr, phi_arr)
    A_c = build_ac_matrix(x_arr, phi_arr)
    A_d = build_ad_matrix(x_arr)

    d_vec = x_arr + A_d @ p_vec
    mu_cert = float(np.min(d_vec))
    F_vec = A_F @ p_vec
    max_F = float(np.max(np.abs(F_vec)))
    c_vec = x_arr + A_c @ p_vec
    min_c = float(np.min(c_vec))

    # Slack matrix
    L_w_num = build_lw_matrix(x_arr, phi_arr, p_vec)
    T_L_T = T.T @ L_w_num @ T
    S_slack = T_L_T - lam_M * metric_M
    S_slack = 0.5 * (S_slack + S_slack.T)
    eigs_slack = np.linalg.eigvalsh(S_slack)

    # mpmath 50 digits on Q^T L_w Q
    mpmath.mp.dps = 50
    Q_np = scipy.linalg.null_space(np.ones((1, m), dtype=float))
    dim = m - 1
    Q_mp = mpmath.mp.matrix(m, dim)
    for r in range(m):
        for c in range(dim):
            Q_mp[r, c] = mpmath.mp.mpf(str(Q_np[r, c]))

    for j in range(dim):
        v = Q_mp[:, j]
        for k in range(j):
            proj = mpmath.mp.fdot(Q_mp[:, k], v)
            for r in range(m):
                v[r] -= proj * Q_mp[r, k]
        norm_v = mpmath.sqrt(sum(v[r] ** 2 for r in range(m)))
        for r in range(m):
            Q_mp[r, j] = v[r] / norm_v

    L_mp = mpmath.mp.matrix(m, m)
    for r in range(m):
        for c in range(m):
            L_mp[r, c] = mpmath.mp.mpf(str(L_w_num[r, c]))

    M_mp = Q_mp.T * L_mp * Q_mp
    for r in range(dim):
        for c in range(dim):
            val = (M_mp[r, c] + M_mp[c, r]) / 2
            M_mp[r, c] = val
            M_mp[c, r] = val
    eigs_M_mp = mpmath.mp.eigsy(M_mp, eigvals_only=True)
    lam_min_realized = float(min(eigs_M_mp))

    return {
        "mu_cert": mu_cert,
        "max_F": max_F,
        "min_c": min_c,
        "min_slack_eig": float(np.min(eigs_slack)),
        "lam_min_realized": lam_min_realized,
        "p_min": float(np.min(p_vec)),
        "p_max": float(np.max(p_vec)),
    }


# ==============================================================================
# SECTION 3: MULTI-START LOCAL POLISH (20 ARRANQUES) Y ABLACIONES
# ==============================================================================
print("=" * 105)
print("3. MULTI-START DEL PULIDO LOCAL (20 ARRANQUES) Y ABLACIONES (lambda_M = 1e-10)")
print("=" * 105)

s7_base = get_s7_appendix_a_template(h=0.001, t=0.0, kappa=1.9)
x_base = np.array(s7_base.x, dtype=float)
phi_base = np.array(s7_base.phi, dtype=float)

rng = np.random.default_rng(42)
N_STARTS = 20
target_lam = 1e-10

def polish_obj_joint(params: np.ndarray) -> float:
    dx = params[:6]
    dphi = params[6:]
    x_new = np.zeros(7)
    x_new[:6] = x_base[:6] + dx
    x_new[6] = 1.0 - np.sum(x_new[:6])
    if np.any(x_new <= 0.01):
        return 1e6
    phi_new = np.zeros(7)
    phi_new[:6] = phi_base[:6] + dphi
    phi_new[6] = phi_base[6]
    mu_val, _, _, _, _ = fast_solve_phase2(x_new, phi_new, target_lam)
    if math.isnan(mu_val):
        return 1e6
    return -mu_val

print(f"\n--- A. Optimizacion Conjunta (Delta x, Delta phi) - {N_STARTS} arranques ---")
joint_results = []
best_joint_p = None
best_joint_mu = -1.0

for s_idx in range(N_STARTS):
    p0 = np.zeros(12) if s_idx == 0 else np.concatenate([rng.normal(0.0, 1e-4, size=6), rng.normal(0.0, 2e-6, size=6)])
    res = scipy.optimize.minimize(polish_obj_joint, p0, method="Nelder-Mead", options={"maxiter": 150, "xatol": 1e-4, "fatol": 1e-5})
    mu_val = -res.fun
    joint_results.append(mu_val)
    if mu_val > best_joint_mu:
        best_joint_mu = mu_val
        best_joint_p = res.x

print(f"  Resultados Joint (20 starts): Min mu*={min(joint_results):.8f} | Media mu*={np.mean(joint_results):.8f} | Max mu*={max(joint_results):.8f}")

# Extract best joint point
assert best_joint_p is not None
dx_best = best_joint_p[:6]
dphi_best = best_joint_p[6:]
x_pol = np.zeros(7)
x_pol[:6] = x_base[:6] + dx_best
x_pol[6] = 1.0 - np.sum(x_pol[:6])
phi_pol = np.zeros(7)
phi_pol[:6] = phi_base[:6] + dphi_best
phi_pol[6] = phi_base[6]

# Ablation 1: Delta x only
def polish_obj_dx(dx: np.ndarray) -> float:
    x_new = np.zeros(7)
    x_new[:6] = x_base[:6] + dx
    x_new[6] = 1.0 - np.sum(x_new[:6])
    if np.any(x_new <= 0.01):
        return 1e6
    mu_val, _, _, _, _ = fast_solve_phase2(x_new, phi_base, target_lam)
    if math.isnan(mu_val):
        return 1e6
    return -mu_val

print(f"\n--- B. Ablacion 1: Solo Delta x (Delta phi = 0 fijo) - {N_STARTS} arranques ---")
dx_results = []
for s_idx in range(N_STARTS):
    p0 = np.zeros(6) if s_idx == 0 else rng.normal(0.0, 1e-4, size=6)
    res = scipy.optimize.minimize(polish_obj_dx, p0, method="Nelder-Mead", options={"maxiter": 100, "xatol": 1e-4, "fatol": 1e-5})
    dx_results.append(-res.fun)

print(f"  Resultados Solo Delta x: Min mu*={min(dx_results):.8f} | Media mu*={np.mean(dx_results):.8f} | Max mu*={max(dx_results):.8f}")

# Ablation 2: Delta phi only
def polish_obj_dphi(dphi: np.ndarray) -> float:
    phi_new = np.zeros(7)
    phi_new[:6] = phi_base[:6] + dphi
    phi_new[6] = phi_base[6]
    mu_val, _, _, _, _ = fast_solve_phase2(x_base, phi_new, target_lam)
    if math.isnan(mu_val):
        return 1e6
    return -mu_val

print(f"\n--- C. Ablacion 2: Solo Delta phi (Delta x = 0 fijo) - {N_STARTS} arranques ---")
dphi_results = []
for s_idx in range(N_STARTS):
    p0 = np.zeros(6) if s_idx == 0 else rng.normal(0.0, 2e-6, size=6)
    res = scipy.optimize.minimize(polish_obj_dphi, p0, method="Nelder-Mead", options={"maxiter": 100, "xatol": 1e-4, "fatol": 1e-5})
    dphi_results.append(-res.fun)

print(f"  Resultados Solo Delta phi: Min mu*={min(dphi_results):.8f} | Media mu*={np.mean(dphi_results):.8f} | Max mu*={max(dphi_results):.8f}")


# ==============================================================================
# SECTION 1 & 2: PRIMAL CERTIFICATION & SLACK DIAGNOSTIC
# ==============================================================================
print("\n" + "=" * 105)
print("1 & 2. CERTIFICACION PRIMAL Y DIAGNOSTICO DE HOLGURA A lambda_M = 1e-10")
print("=" * 105)

points_to_cert = [
    ("Campeon R1 (ID 125678)", np.array(r1["phenotype_x"]), np.array(r1["phenotype_phi"])),
    ("Campeon R3 (ID 4)", np.array(r3["phenotype_x"]), np.array(r3["phenotype_phi"])),
    ("S7(h=1e-3, t=0)", x_base, phi_base),
    ("Punto Pulido (Multi-start)", x_pol, phi_pol),
]

print(f"{'Candidato':<28} | {'Status':<10} | {'mu_cert':<12} | {'lam_real (mpmath 50d)':<22} | {'min eig(Holgura)':<18} | {'max|F|':<10}")
print("-" * 110)

diag_records = {}
for cname, xc, phic in points_to_cert:
    mu_val, p_val, T, metric_M, status = fast_solve_phase2(xc, phic, 1e-10)
    assert p_val is not None and T is not None and metric_M is not None
    cert = certify_primal_50d(xc, phic, p_val, 1e-10, T, metric_M)
    diag_records[cname] = {**cert, "status": status}
    print(f"{cname:<28} | {status:<10} | {cert['mu_cert']:<12.8f} | {cert['lam_min_realized']:<22.12e} | {cert['min_slack_eig']:<18.12e} | {cert['max_F']:<10.2e}")


# ==============================================================================
# SECTION 4: BARRIDO DE MU* CON LAMBDA_REALIZADO EN EL EJE
# ==============================================================================
print("\n" + "=" * 105)
print("4. BARRIDO DE mu*(lambda_M) CON lambda_M REALIZADO EN EL EJE")
print("=" * 105)

targets = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
print(f"{'Target lam':<12} | {'Candidato':<28} | {'mu_cert':<14} | {'lambda_M Realizado (mpmath 50d)':<32}")
print("-" * 92)

realized_table = []
for lam in targets:
    for cname, xc, phic in points_to_cert:
        mu_val, p_val, T, metric_M, status = fast_solve_phase2(xc, phic, lam)
        if p_val is not None and T is not None and metric_M is not None:
            cert = certify_primal_50d(xc, phic, p_val, lam, T, metric_M)
            print(f"{lam:<12.0e} | {cname:<28} | {cert['mu_cert']:<14.8f} | {cert['lam_min_realized']:<32.16e}")
            realized_table.append({
                "target_lam": lam,
                "candidate": cname,
                "mu_cert": cert["mu_cert"],
                "lam_realized": cert["lam_min_realized"],
                "status": status,
            })
        else:
            print(f"{lam:<12.0e} | {cname:<28} | {'INFEASIBLE':<14} | {'N/A':<32}")

# Save all results
results = {
    "diagnostics_1e10": diag_records,
    "multi_start_stats": {
        "joint_min": float(min(joint_results)),
        "joint_mean": float(np.mean(joint_results)),
        "joint_max": float(max(joint_results)),
        "dx_min": float(min(dx_results)),
        "dx_mean": float(np.mean(dx_results)),
        "dx_max": float(max(dx_results)),
        "dphi_min": float(min(dphi_results)),
        "dphi_mean": float(np.mean(dphi_results)),
        "dphi_max": float(max(dphi_results)),
    },
    "realized_sweep": realized_table,
    "best_polished_point": {
        "phenotype_x": list(x_pol),
        "phenotype_phi": list(phi_pol),
        "norm_dx": float(np.linalg.norm(dx_best)),
        "norm_dphi": float(np.linalg.norm(dphi_best)),
    },
}

with open("results/comprehensive_audit_results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)

print("\nResultados completos guardados en results/comprehensive_audit_results.json")
