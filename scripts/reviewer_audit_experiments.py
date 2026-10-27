"""Reviewer Audit Experiments Script:

1. Verbose Clarabel runs of S7(1e-3) at lam_M = 1e-10:
   - Unscaled LMI vs Rescaled LMI (with h_eff^2)
2. Phase I comparison across key evolutionary milestones:
   - ID 4 (Seed), ID 125678 (HoF champ), ID 399165 (Island dominant)
   - lam_hat, raw lam, and physical lam_min(M) in mpmath 50d.
3. Sweep of ell_W from -6.0 to -4.0 on genome 274309:
   - lam_hat, h_eff, physical lam_min(M) (mpmath 50d), and Phase II mu*.
4. Phase II of distinct deduplicated phenotypes in Top-20:
   - 399165, 125678, 274309, 80548, 15223, 17863.
"""

from __future__ import annotations

import json
import math
from typing import Any

import cvxpy as cp
import mpmath
import numpy as np
import scipy.linalg  # type: ignore # pyright: ignore[reportMissingTypeStubs]

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
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)
from kuramoto_templates.search.genome import Genome, to_phenotype

# ==============================================================================
# EXP 1: VERBOSE CLARABEL RUNS (UNSCALED VS RESCALED LMI ON S7)
# ==============================================================================
print("=" * 90)
print("EXPERIMENTO 1: CLARABEL VERBOSE=TRUE EN S7(1e-3) (SIN ESCALAR VS REESCALADA)")
print("=" * 90)

s7 = get_s7_appendix_a_template(h=0.001, t=0.0, kappa=1.9)
x_s7 = np.array(s7.x, dtype=float)
phi_s7 = np.array(s7.phi, dtype=float)
m = 7
pairs = get_pair_indices(m)
K = len(pairs)
sheets, ref_sheet = compute_sheets(x_s7, phi_s7)
heff = float(compute_h_eff(x_s7, phi_s7, sheets))
T = build_scaled_basis_t(x_s7, phi_s7, sheets, ref_sheet, heff)
assert T is not None

P_perp = np.eye(m) - np.ones((m, m)) / float(m)
T_prime = P_perp @ T
metric_M_unscaled = T_prime.T @ T_prime
metric_M_rescaled = metric_M_unscaled / (heff**2)

A_F = build_af_matrix(x_s7, phi_s7)
A_c = build_ac_matrix(x_s7, phi_s7)
A_d = build_ad_matrix(x_s7)
N = scipy.linalg.null_space(A_F)

# Basis for unscaled LMI: T^T L_w T
G_unscaled_basis = []
for a, b in pairs:
    d_ab = np.zeros(m, dtype=float)
    d_ab[a] = 1.0; d_ab[b] = -1.0
    Td = T.T @ d_ab
    weight = x_s7[a] * x_s7[b] * math.cos(phi_s7[a] - phi_s7[b])
    Gk = weight * np.outer(Td, Td)
    G_unscaled_basis.append(0.5 * (Gk + Gk.T))

# Basis for rescaled LMI: (T^T L_w T) / heff^2
G_rescaled_basis = [Gk / (heff**2) for Gk in G_unscaled_basis]

lam_target = 1e-10

print("\n>>> A. LMI SIN ESCALAR: T^T L_w T >= lam_M * (T'^T T') <<<")
z_u = cp.Variable(N.shape[1])
t_u = cp.Variable()
p_u = N @ z_u
L_u = sum(p_u[k] * G_unscaled_basis[k] for k in range(K))
prob_u = cp.Problem(cp.Maximize(t_u), [
    x_s7 + A_d @ p_u >= t_u,
    p_u >= 0.0, p_u <= 1.0,
    L_u >> lam_target * metric_M_unscaled,
    x_s7 + A_c @ p_u >= 1e-4,
])
prob_u.solve(solver=cp.CLARABEL, tol_gap_abs=1e-9, tol_gap_rel=1e-9, verbose=True)

print("\n>>> B. LMI REESCALADA: G_hat >= lam_M * M_hat (con division por heff^2) <<<")
z_r = cp.Variable(N.shape[1])
t_r = cp.Variable()
p_r = N @ z_r
L_r = sum(p_r[k] * G_rescaled_basis[k] for k in range(K))
prob_r = cp.Problem(cp.Maximize(t_r), [
    x_s7 + A_d @ p_r >= t_r,
    p_r >= 0.0, p_r <= 1.0,
    L_r >> lam_target * metric_M_rescaled,
    x_s7 + A_c @ p_r >= 1e-4,
])
prob_r.solve(solver=cp.CLARABEL, tol_gap_abs=1e-9, tol_gap_rel=1e-9, verbose=True)


# ==============================================================================
# EXP 2: COMPARACION FISICA DE FASE I (ID 4, 125678, 399165)
# ==============================================================================
print("\n" + "=" * 90)
print("EXPERIMENTO 2: COMPARACION DEL MARGEN FISICO lam_min(M) DE FASE I")
print("=" * 90)

ckpt_path = "runs/night12h_seeded_mu06923_s42/checkpoints/checkpoint_gen_1675.json"
with open(ckpt_path, "r", encoding="utf-8") as f:
    ckpt_data = json.load(f)

# Find 125678 and 399165
ind_125 = next(ind for ind in ckpt_data["hall_of_fame_global"] if ind["id"] == 125678)
ind_399 = next(ind for ind in ckpt_data["island_populations"][0] if ind["id"] == 399165)

# Seed 4
from kuramoto_templates.search.seeds import create_champion_3178_seed
seed_4 = create_champion_3178_seed()
x_4, phi_4 = seed_4.phenotype_x, seed_4.phenotype_phi

compare_candidates = [
    ("ID 4 (Seed champ)", x_4, phi_4, seed_4.genome),
    ("ID 125678 (HoF champ)", tuple(ind_125["phenotype_x"]), tuple(ind_125["phenotype_phi"]), None),
    ("ID 399165 (Island dominant)", tuple(ind_399["phenotype_x"]), tuple(ind_399["phenotype_phi"]), None),
]

evaluator = ScaledSDPEvaluator(EvalConfig(solver_tol=1e-9, lam_floor_hat=1e-5))

print(f"{'Candidato':<28} | {'h_eff':<12} | {'lam_hat (Fase 1)':<18} | {'lam_min(M) mpmath 50d':<26} | {'Status'}")
print("-" * 95)

def compute_mpmath_M(xc, phic, P_mat):
    m_len = len(xc)
    x_c = np.array(xc, dtype=float)
    phi_c = np.array(phic, dtype=float)
    p_v = p_matrix_to_vector(np.array(P_mat, dtype=float))
    L_w = build_lw_matrix(x_c, phi_c, p_v)
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

for cname, xc, phic, gen in compare_candidates:
    sheets, _ = compute_sheets(np.array(xc), np.array(phic))
    h_val = float(compute_h_eff(np.array(xc), np.array(phic), sheets))
    res = evaluator.evaluate_phase1(xc, phic)
    assert res.P is not None
    lam_mp = compute_mpmath_M(xc, phic, res.P)
    print(f"{cname:<28} | {h_val:<12.6e} | {res.lam_normalized:<18.10e} | {lam_mp:<26.16e} | {res.status}")


# ==============================================================================
# EXP 3: BARRIDO DE ell_W DE -6.0 A -4.0 EN GENOMA 274309
# ==============================================================================
print("\n" + "=" * 90)
print("EXPERIMENTO 3: BARRIDO DE ell_W DE -6.0 A -4.0 EN GENOMA 274309")
print("=" * 90)

ind_274 = next(ind for ind in ckpt_data["hall_of_fame_global"] if ind["id"] == 274309)
base_g = ind_274["genome"]

ell_grid = np.linspace(-6.0, -4.0, 11)

print(f"{'ell_W':<8} | {'phi_W (rad)':<14} | {'h_eff':<12} | {'lam_hat (Fase 1)':<18} | {'lam_min(M) mpmath':<22} | {'mu* (Fase 2)'}")
print("-" * 95)

for ell_val in ell_grid:
    new_ell = list(base_g["ell"])
    new_ell[3] = float(ell_val)
    g_sweep = Genome(
        m=base_g["m"],
        K_ref=base_g["K_ref"],
        k=tuple(base_g["k"]),
        sgn=tuple(base_g["sgn"]),
        ell=tuple(new_ell),
        u=tuple(base_g["u"]),
    )
    x_sw, phi_sw = to_phenotype(g_sweep)
    sheets_sw, _ = compute_sheets(np.array(x_sw), np.array(phi_sw))
    heff_sw = float(compute_h_eff(np.array(x_sw), np.array(phi_sw), sheets_sw))

    # Phase I
    res_p1 = evaluator.evaluate_phase1(x_sw, phi_sw)
    if res_p1.P is not None:
        lam_M_val = compute_mpmath_M(x_sw, phi_sw, res_p1.P)
    else:
        lam_M_val = float("nan")

    # Phase II with exact LMI at 1e-10
    res_p2 = evaluator.evaluate_phase2(x_sw, phi_sw, eps=1e-10 / (heff_sw**2))
    mu_val = res_p2.mu if res_p2.mu is not None else float("nan")
    phi_W_rad = phi_sw[3]

    print(f"{ell_val:<8.2f} | {phi_W_rad:<14.4e} | {heff_sw:<12.6e} | {res_p1.lam_normalized:<18.10e} | {lam_M_val:<22.14e} | {mu_val:<12.8f}")


# ==============================================================================
# EXP 4: FASE II DE LOS FENOTIPOS DISTINTOS (DEDUPLICADOS) DEL TOP-20
# ==============================================================================
print("\n" + "=" * 90)
print("EXPERIMENTO 4: FASE II DE LOS FENOTIPOS DISTINTOS EN TOP-20 (INCLUYENDO 399165)")
print("=" * 90)

# Deduplicate Top-20 by phenotype
distinct_ids = [399165, 125678, 274309, 80548, 15223, 16123, 17863]
all_inds = {ind["id"]: ind for ind in ckpt_data["island_populations"][0] + ckpt_data["hall_of_fame_global"]}

print(f"{'ID':<8} | {'f (Fase 1)':<18} | {'mu* (Fase 2)':<14} | {'lam_min(M) mpmath':<24} | {'Status Fase 2'}")
print("-" * 80)

for did in distinct_ids:
    if did not in all_inds:
        continue
    ind = all_inds[did]
    xc = tuple(ind["phenotype_x"])
    phic = tuple(ind["phenotype_phi"])
    sheets_d, _ = compute_sheets(np.array(xc), np.array(phic))
    heff_d = float(compute_h_eff(np.array(xc), np.array(phic), sheets_d))

    res_p2 = evaluator.evaluate_phase2(xc, phic, eps=1e-10 / (heff_d**2))
    mu_val = res_p2.mu if res_p2.mu is not None else float("nan")
    if res_p2.P is not None:
        lam_M_val = compute_mpmath_M(xc, phic, res_p2.P)
    else:
        lam_M_val = float("nan")

    print(f"{did:<8} | {ind['f']:<18.10e} | {mu_val:<14.8f} | {lam_M_val:<24.16e} | {res_p2.status}")
