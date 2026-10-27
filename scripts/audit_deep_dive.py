"""Deep Dive Audit Script:

1. Realized lambda_min(M) with exact LMI calculated in mpmath to 50 digits.
2. Full curve mu*(lambda_M) from 1e-12 to 1e-8 with exact LMI for:
   - R1 (ID 125678)
   - R3 (ID 4)
   - S7(h=5e-4)
   - S7(h=1e-3)
   - S7(h=2e-3)
3. Multi-start of S7 envelope over (h, t, kappa).
4. Direct local polish of (x, phi) starting from S7 using SLSQP / Nelder-Mead on Phase II mu*.
   Comparison against R1.
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
)
from kuramoto_templates.evaluators.scaled import (
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


def solve_phase2_exact_lmi(
    x: tuple[float, ...] | list[float] | np.ndarray,
    phi: tuple[float, ...] | list[float] | np.ndarray,
    lam_M: float,
) -> tuple[float, np.ndarray | None]:
    """Solve Phase II maximizing mu with exact metric constraint T^T L_w T >= lam_M (T'^T T')."""
    x_arr = np.array(x, dtype=float)
    phi_arr = np.array(phi, dtype=float)
    m = len(x_arr)
    pairs = get_pair_indices(m)
    K = len(pairs)
    sheets, ref_sheet = compute_sheets(x_arr, phi_arr)
    heff = float(compute_h_eff(x_arr, phi_arr, sheets))
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, heff)
    if T is None:
        return float("nan"), None

    P_perp = np.eye(m) - np.ones((m, m)) / float(m)
    T_prime = P_perp @ T
    metric_M = T_prime.T @ T_prime

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
            return float(prob.value), p_val
    except Exception:
        pass
    return float("nan"), None


def compute_mpmath_50d_lam_min(
    x: tuple[float, ...] | list[float] | np.ndarray,
    phi: tuple[float, ...] | list[float] | np.ndarray,
    p_vec: np.ndarray,
) -> float:
    """Compute exact lam_min of Q^T L_w Q in mpmath with 50 decimal digits."""
    mpmath.mp.dps = 50
    x_arr = np.array(x, dtype=float)
    phi_arr = np.array(phi, dtype=float)
    m = len(x_arr)
    L_w = build_lw_matrix(x_arr, phi_arr, p_vec)

    # Q orthonormal basis for 1^perp
    Q_np = scipy.linalg.null_space(np.ones((1, m), dtype=float))
    # Convert Q and L_w to mpmath
    dim = m - 1
    Q_mp = mpmath.mp.matrix(m, dim)
    for r in range(m):
        for c in range(dim):
            Q_mp[r, c] = mpmath.mp.mpf(str(Q_np[r, c]))

    # Reorthonormalize Q in 50 digits via Gram-Schmidt
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
            L_mp[r, c] = mpmath.mp.mpf(str(L_w[r, c]))

    # M_mp = Q^T L_w Q
    M_mp = Q_mp.T * L_mp * Q_mp
    # Symmetrize
    for r in range(dim):
        for c in range(dim):
            val = (M_mp[r, c] + M_mp[c, r]) / 2
            M_mp[r, c] = val
            M_mp[c, r] = val

    eigs = mpmath.mp.eigsy(M_mp, eigvals_only=True)
    return float(min(eigs))


# ==============================================================================
# SECTION 1: REALIZED LAMBDA_MIN(M) EN MPMATH 50 DIGITOS
# ==============================================================================
print("=" * 90)
print("1. RECALCULO DE LAMBDA_MIN(M) REALIZADO CON LMI EXACTA EN MPMATH (50 DIGITOS)")
print("=" * 90)

targets_lam = [1e-11, 1e-10, 1e-9]
candidates = [
    ("Campeon R1 (ID 125678)", r1["phenotype_x"], r1["phenotype_phi"]),
    ("Campeon R3 (ID 4)", r3["phenotype_x"], r3["phenotype_phi"]),
    ("S7(h=1e-3)", get_s7_appendix_a_template(h=0.001, t=0.0).x, get_s7_appendix_a_template(h=0.001, t=0.0).phi),
]

print(f"{'Candidato':<26} | {'Target lam_M':<14} | {'mu* Obtenido':<14} | {'lam_min(M) mpmath 50d':<24}")
print("-" * 85)

for cname, xc, phic in candidates:
    for lam_target in targets_lam:
        mu_val, p_sol = solve_phase2_exact_lmi(xc, phic, lam_target)
        if p_sol is not None:
            lam_real_mp = compute_mpmath_50d_lam_min(xc, phic, p_sol)
            print(f"{cname:<26} | {lam_target:<14.0e} | {mu_val:<14.8f} | {lam_real_mp:<24.16e}")
        else:
            print(f"{cname:<26} | {lam_target:<14.0e} | {'INFEASIBLE':<14} | {'N/A':<24}")


# ==============================================================================
# SECTION 2: CURVA COMPLETA MU*(LAMBDA_M) DE 1E-12 A 1E-8
# ==============================================================================
print("\n" + "=" * 90)
print("2. CURVA COMPLETA mu*(lambda_M) DE 1e-12 A 1e-8 (LMI EXACTA)")
print("=" * 90)

lam_grid = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
configs_to_sweep = [
    ("Campeon R1", r1["phenotype_x"], r1["phenotype_phi"]),
    ("Campeon R3", r3["phenotype_x"], r3["phenotype_phi"]),
    ("S7(h=5e-4)", get_s7_appendix_a_template(h=0.0005, t=0.0).x, get_s7_appendix_a_template(h=0.0005, t=0.0).phi),
    ("S7(h=1e-3)", get_s7_appendix_a_template(h=0.0010, t=0.0).x, get_s7_appendix_a_template(h=0.0010, t=0.0).phi),
    ("S7(h=2e-3)", get_s7_appendix_a_template(h=0.0020, t=0.0).x, get_s7_appendix_a_template(h=0.0020, t=0.0).phi),
]

header = f"{'lambda_M':<12} | " + " | ".join(f"{cfg[0]:<14}" for cfg in configs_to_sweep)
print(header)
print("-" * len(header))

curve_data: dict[str, list[float]] = {cfg[0]: [] for cfg in configs_to_sweep}

for lam in lam_grid:
    row_str = f"{lam:<12.0e} | "
    for cname, xc, phic in configs_to_sweep:
        mu_val, _ = solve_phase2_exact_lmi(xc, phic, lam)
        curve_data[cname].append(mu_val)
        row_str += f"{mu_val:<14.8f} | "
    print(row_str)


# ==============================================================================
# SECTION 3: MULTI-START DEL ENVOLVENTE S7(h, t, kappa)
# ==============================================================================
print("\n" + "=" * 90)
print("3. MULTI-START DEL ENVOLVENTE S7(h, t, kappa) (10 SEMILLAS)")
print("=" * 90)

target_test_lam = 1e-10

def s7_obj(params: list[float] | np.ndarray) -> float:
    h_v, t_v, kap_v = params
    if not (1e-4 <= h_v <= 5e-3 and -0.01 <= t_v <= 0.05 and 0.5 <= kap_v <= 4.0):
        return 1e6
    try:
        tmpl = get_s7_appendix_a_template(h=h_v, t=t_v, kappa=kap_v)
        mu_val, _ = solve_phase2_exact_lmi(tmpl.x, tmpl.phi, target_test_lam)
        if math.isnan(mu_val):
            return 1e6
        return -mu_val
    except Exception:
        return 1e6

rng = np.random.default_rng(2026)
multistart_seeds = [
    [0.001, 0.0, 1.9],
    [0.0005, 0.0, 1.9],
    [0.002, 0.0, 1.9],
    [0.0008, 0.01, 2.5],
    [0.0006, -0.005, 2.2],
]
for _ in range(5):
    multistart_seeds.append([
        float(rng.uniform(0.0003, 0.0025)),
        float(rng.uniform(-0.005, 0.02)),
        float(rng.uniform(1.2, 3.2)),
    ])

best_s7_mu = -1.0
best_s7_params = None

print(f"{'Seed':<5} | {'h_init':<10} | {'t_init':<10} | {'k_init':<8} | {'mu* Final':<14} | {'(h*, t*, kappa*) Optimos'}")
print("-" * 85)

for s_idx, s_init in enumerate(multistart_seeds):
    res_opt = scipy.optimize.minimize(
        s7_obj,
        s_init,
        method="Nelder-Mead",
        options={"maxiter": 80, "xatol": 1e-4, "fatol": 1e-5},
    )
    mu_final = -res_opt.fun
    h_f, t_f, k_f = res_opt.x
    print(f"{s_idx+1:<5} | {s_init[0]:<10.4e} | {s_init[1]:<10.2e} | {s_init[2]:<8.2f} | {mu_final:<14.8f} | h={h_f:.4e}, t={t_f:.2e}, k={k_f:.2f}")
    if mu_final > best_s7_mu:
        best_s7_mu = mu_final
        best_s7_params = (h_f, t_f, k_f)

print(f"\nMejor mu* de Envolvente S7 Multi-start a lam_M=1e-10: {best_s7_mu:.8f}")
print(f"Parametros: h={best_s7_params[0]:.4e}, t={best_s7_params[1]:.2e}, kappa={best_s7_params[2]:.2f}")


# ==============================================================================
# SECTION 4: PULIDO LOCAL DIRECTO DE (x, phi) DESDE S7 CON OPTIMIZACION CONTINUA
# ==============================================================================
print("\n" + "=" * 90)
print("4. PULIDO LOCAL CONTINUO DE (x, phi) DESDE S7 CON NELDER-MEAD / SLSQP")
print("=" * 90)

# Start from S7(h=1e-3)
s7_base = get_s7_appendix_a_template(h=0.001, t=0.0, kappa=1.9)
x_base = np.array(s7_base.x, dtype=float)
phi_base = np.array(s7_base.phi, dtype=float)
m = 7

# Parameterization:
# delta_x for first 6 clusters (cluster W determined by 1 - sum(x[:6]))
# delta_phi for first 6 clusters (gauge phi[6] fixed at 3*pi/2)
dim_polish = 12

eval_counter = 0

def polish_obj(params: np.ndarray) -> float:
    global eval_counter
    eval_counter += 1
    dx = params[:6]
    dphi = params[6:]

    x_new = np.zeros(m)
    x_new[:6] = x_base[:6] + dx
    x_new[6] = 1.0 - np.sum(x_new[:6])

    if np.any(x_new <= 0.01):
        return 1e6

    phi_new = np.zeros(m)
    phi_new[:6] = phi_base[:6] + dphi
    phi_new[6] = phi_base[6]

    try:
        mu_val, _ = solve_phase2_exact_lmi(x_new, phi_new, target_test_lam)
        if math.isnan(mu_val):
            return 1e6
        return -mu_val
    except Exception:
        return 1e6

init_params = np.zeros(dim_polish)
print(f"Optimizando {dim_polish} variables continuas (dx, dphi) partiendo de S7(1e-3)...")

res_polish = scipy.optimize.minimize(
    polish_obj,
    init_params,
    method="Nelder-Mead",
    options={"maxiter": 250, "xatol": 1e-4, "fatol": 1e-5},
)

mu_polished = -res_polish.fun
dx_opt = res_polish.x[:6]
dphi_opt = res_polish.x[6:]

print(f"Evaluaciones realizadas: {eval_counter}")
print(f"mu* Inicial (S7 base): {-polish_obj(init_params):.8f}")
print(f"mu* Pulido Local     : {mu_polished:.8f}")
print(f"mu* Campeon R1       : {solve_phase2_exact_lmi(r1['phenotype_x'], r1['phenotype_phi'], target_test_lam)[0]:.8f}")
print(f"Norma ||dx|| opt     : {np.linalg.norm(dx_opt):.4e}")
print(f"Norma ||dphi|| opt   : {np.linalg.norm(dphi_opt):.4e}")

# Save all results to json
results = {
    "curve_data": curve_data,
    "lam_grid": lam_grid,
    "best_s7_envelope": {
        "mu_star": best_s7_mu,
        "params": best_s7_params,
    },
    "local_polish": {
        "mu_star": mu_polished,
        "evaluations": eval_counter,
        "norm_dx": float(np.linalg.norm(dx_opt)),
        "norm_dphi": float(np.linalg.norm(dphi_opt)),
    },
}

with open("results/audit_deep_dive.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)

print("\nResultados guardados exitosamente en results/audit_deep_dive.json")
