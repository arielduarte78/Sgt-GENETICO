import json
import math
import cvxpy as cp
import mpmath
import numpy as np
import scipy.linalg

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
from kuramoto_templates.search.seeds import create_champion_3178_seed


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


# ==============================================================================
# 1. BARRIDO DE ell_W DE -6 A -9 Y phi_W = 0 EXACTO SOBRE 274309
# ==============================================================================
print("=" * 90)
print("EXPERIMENTO: BARRIDO DE ell_W DE -6.0 A -9.0 Y phi_W = 0 EXACTO (GENOMA 274309)")
print("=" * 90)

ckpt_path = "runs/night12h_seeded_mu06923_s42/checkpoints/checkpoint_gen_1675.json"
with open(ckpt_path, "r", encoding="utf-8") as f:
    ckpt_data = json.load(f)

ind_274 = next(ind for ind in ckpt_data["hall_of_fame_global"] if ind["id"] == 274309)
base_g = ind_274["genome"]

evaluator = ScaledSDPEvaluator(EvalConfig(solver_tol=1e-9, lam_floor_hat=1e-5))

sweep_points = [
    ("-6.0", -6.0, 1),
    ("-6.5", -6.5, 1),
    ("-7.0", -7.0, 1),
    ("-7.5", -7.5, 1),
    ("-8.0", -8.0, 1),
    ("-8.5", -8.5, 1),
    ("-9.0", -9.0, 1),
    ("phi_W=0", -6.0, 0),  # sgn = 0 -> exact 0
]

print(f"{'Caso':<12} | {'phi_W (rad)':<14} | {'h_eff':<12} | {'lam_hat (Fase 1)':<18} | {'lam_min(M) mpmath':<22} | {'mu* (Fase 2)'}")
print("-" * 95)

for label, ell_val, sgn_val in sweep_points:
    new_ell = list(base_g["ell"])
    new_sgn = list(base_g["sgn"])
    new_ell[3] = float(ell_val)
    new_sgn[3] = int(sgn_val)
    # Decode directly using decode_x and decode_phi to explore ell in [-9.0, -6.0]
    from kuramoto_templates.search.genome import decode_phi, decode_x
    x_sw = decode_x(tuple(base_g["u"]))
    phi_sw = decode_phi(base_g["K_ref"], tuple(base_g["k"]), tuple(new_sgn), tuple(new_ell))
    sheets_sw, _ = compute_sheets(np.array(x_sw), np.array(phi_sw))
    heff_sw = float(compute_h_eff(np.array(x_sw), np.array(phi_sw), sheets_sw))

    # Phase 1
    res_p1 = evaluator.evaluate_phase1(x_sw, phi_sw)
    lam_mp1 = compute_mpmath_M(x_sw, phi_sw, res_p1.P) if res_p1.P is not None else float("nan")

    # Phase 2 (exact rescaled LMI with lam_M = 1e-10)
    res_p2 = evaluator.evaluate_phase2(x_sw, phi_sw, eps=1e-10 / (heff_sw**2))
    mu_val = res_p2.mu if res_p2.mu is not None else float("nan")

    print(f"{label:<12} | {phi_sw[3]:<14.4e} | {heff_sw:<12.6e} | {res_p1.lam_normalized:<18.10e} | {lam_mp1:<22.14e} | {mu_val:<12.8f}")


# ==============================================================================
# 2. LAS TRES VARIANTES DE SDP (Q CRUDO, T CRUDO, T/h^2) EN S7(1e-3)
# ==============================================================================
print("\n" + "=" * 90)
print("EXPERIMENTO: TRES VARIANTES DE SDP (Q CRUDO, T CRUDO, T/h^2) A lam_M = 1e-10")
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
P_perp = np.eye(m) - np.ones((m, m)) / float(m)
T_prime = P_perp @ T
metric_M_unscaled = T_prime.T @ T_prime
metric_M_rescaled = metric_M_unscaled / (heff**2)

A_F = build_af_matrix(x_s7, phi_s7)
A_c = build_ac_matrix(x_s7, phi_s7)
A_d = build_ad_matrix(x_s7)
N = scipy.linalg.null_space(A_F)

# Basis Q crudo
Q_np = scipy.linalg.null_space(np.ones((1, m), dtype=float))
G_Q_basis = []
for a, b in pairs:
    d_ab = np.zeros(m, dtype=float)
    d_ab[a] = 1.0; d_ab[b] = -1.0
    Qd = Q_np.T @ d_ab
    w = x_s7[a] * x_s7[b] * math.cos(phi_s7[a] - phi_s7[b])
    Gk = w * np.outer(Qd, Qd)
    G_Q_basis.append(0.5 * (Gk + Gk.T))

# Basis T crudo
G_T_basis = []
for a, b in pairs:
    d_ab = np.zeros(m, dtype=float)
    d_ab[a] = 1.0; d_ab[b] = -1.0
    Td = T.T @ d_ab
    w = x_s7[a] * x_s7[b] * math.cos(phi_s7[a] - phi_s7[b])
    Gk = w * np.outer(Td, Td)
    G_T_basis.append(0.5 * (Gk + Gk.T))

# Basis T/h^2
G_Th2_basis = [Gk / (heff**2) for Gk in G_T_basis]

lam_target = 1e-10

sdp_variants = [
    ("1. Q crudo (Q^T L_w Q >= lam*I)", G_Q_basis, lam_target * np.eye(m-1)),
    ("2. T crudo (T^T L_w T >= lam*T'^T T')", G_T_basis, lam_target * metric_M_unscaled),
    ("3. T/h^2 reescalado (G_hat >= lam*M_hat)", G_Th2_basis, lam_target * metric_M_rescaled),
]

print(f"{'Variante SDP':<40} | {'Status':<12} | {'Iters':<6} | {'mu*':<14} | {'lam_min(M) mpmath':<22}")
print("-" * 100)

for vname, basis, bound_mat in sdp_variants:
    z_var = cp.Variable(N.shape[1])
    t_var = cp.Variable()
    p_vec = N @ z_var
    L_expr = sum(p_vec[k] * basis[k] for k in range(K))
    prob = cp.Problem(cp.Maximize(t_var), [
        x_s7 + A_d @ p_vec >= t_var,
        p_vec >= 0.0, p_vec <= 1.0,
        L_expr >> bound_mat,
        x_s7 + A_c @ p_vec >= 1e-4,
    ])
    prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-9, tol_gap_rel=1e-9, tol_feas=1e-8, verbose=False)
    p_opt = np.clip(N @ z_var.value, 0.0, 1.0)
    P_opt = p_vector_to_matrix(p_opt, m)
    lam_mp = compute_mpmath_M(x_s7, phi_s7, P_opt)
    iters = prob.solver_stats.num_iters if prob.solver_stats else "N/A"
    print(f"{vname:<40} | {prob.status:<12} | {iters:<6} | {t_var.value:<14.8f} | {lam_mp:<22.14e}")


# ==============================================================================
# 3. FRENTE mu(lam_realizado) PARA S7, R3 Y R1* (ID 399165)
# ==============================================================================
print("\n" + "=" * 90)
print("EXPERIMENTO: FRENTE mu(lam_realizado) PARA S7(1e-3), R3 (Semilla 4) y R1* (399165)")
print("=" * 90)

ind_399 = next(ind for ind in ckpt_data["island_populations"][0] if ind["id"] == 399165)
x_399 = tuple(ind_399["phenotype_x"])
phi_399 = tuple(ind_399["phenotype_phi"])

seed_4 = create_champion_3178_seed()
x_4, phi_4 = seed_4.phenotype_x, seed_4.phenotype_phi

cands = [
    ("S7 (h=1e-3)", tuple(x_s7), tuple(phi_s7)),
    ("R3 (Semilla 4)", tuple(x_4), tuple(phi_4)),
    ("R1* (ID 399165)", tuple(x_399), tuple(phi_399)),
]

lam_sweep = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8]

print(f"{'Candidato':<18} | {'lam_M target':<14} | {'mu* Fase 2':<12} | {'lam_min(M) mpmath':<24} | {'Status'}")
print("-" * 80)

for cname, xc, phic in cands:
    sheets_c, _ = compute_sheets(np.array(xc), np.array(phic))
    heff_c = float(compute_h_eff(np.array(xc), np.array(phic), sheets_c))
    for lam_tgt in lam_sweep:
        res_p2 = evaluator.evaluate_phase2(xc, phic, eps=lam_tgt / (heff_c**2))
        if res_p2.P is not None and res_p2.status in ("OPTIMAL", "SOLVED"):
            lam_mp = compute_mpmath_M(xc, phic, res_p2.P)
            mu_val = res_p2.mu
        else:
            lam_mp = float("nan")
            mu_val = float("nan")
        print(f"{cname:<18} | {lam_tgt:<14.1e} | {mu_val:<12.8f} | {lam_mp:<24.16e} | {res_p2.status}")
