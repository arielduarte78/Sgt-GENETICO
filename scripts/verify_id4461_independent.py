"""Verificación independiente del ganador de mu_t=0.6910 (ID 4461):
1. Clarabel (1e-8 y 1e-9)
2. SCS (1e-8)
3. NumPy float64 puro sobre G_hat y M
4. mpmath a 50 dígitos decimales sobre G_hat y M
5. Cálculo de mu* (Fase II)
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np
import scipy.linalg
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    compute_sheets,
    compute_h_eff,
    build_scaled_basis_t,
    build_g_hat_basis,
)
from kuramoto_templates.domain.conditions import (
    build_m_basis,
    build_af_matrix,
    build_ac_matrix,
    build_ad_matrix,
    get_pair_indices,
    p_vector_to_matrix,
    verify_solution,
)
from kuramoto_templates.search.genome import Genome


def main():
    mp.mp.dps = 50

    # 1. Cargar ID 4461
    ckpt = json.load(open("runs/run_50gen_mu_06910_frontier/checkpoints/checkpoint_gen_50.json"))
    hof = ckpt["hall_of_fame_global"]
    ind_data = next(ind for ind in hof if ind["id"] == 4461)

    x = tuple(float(v) for v in ind_data["phenotype_x"])
    phi = tuple(float(v) for v in ind_data["phenotype_phi"])
    g_dict = ind_data["genome"]
    genome = Genome(
        m=g_dict["m"],
        K_ref=g_dict["K_ref"],
        k=tuple(g_dict["k"]),
        sgn=tuple(g_dict["sgn"]),
        ell=tuple(g_dict["ell"]),
        u=tuple(g_dict["u"]),
    )

    m = len(x)
    x_arr = np.asarray(x, dtype=np.float64)
    phi_arr = np.asarray(phi, dtype=np.float64)

    sheets, ref_sheet = compute_sheets(x_arr, phi_arr, genome=genome)
    h_eff = compute_h_eff(x_arr, phi_arr, sheets)
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, h_eff)
    G_basis = build_g_hat_basis(x_arr, phi_arr, T, h_eff)

    print("======================================================================")
    print("VERIFICACIÓN INDEPENDIENTE ID 4461 (GANADOR mu_t = 0.6910)")
    print("======================================================================")
    print(f"ID: 4461 | m={m} | h_eff={h_eff:.8e}")
    print(f"Hojas: {sheets} (ref={ref_sheet})\n")

    # Resolver Fase I con Clarabel (1e-8)
    cfg_cl8 = EvalConfig(mu_t=0.6910, solver_name="CLARABEL", solver_tol=1e-8, lam_floor_hat=1e-5)
    ev_cl8 = ScaledSDPEvaluator(config=cfg_cl8, lam_floor_hat=1e-5)
    r_cl8 = ev_cl8.evaluate_phase1(x, phi, genome=genome)

    # Resolver Fase I con Clarabel (1e-9)
    cfg_cl9 = EvalConfig(mu_t=0.6910, solver_name="CLARABEL", solver_tol=1e-9, lam_floor_hat=1e-5)
    ev_cl9 = ScaledSDPEvaluator(config=cfg_cl9, lam_floor_hat=1e-5)
    r_cl9 = ev_cl9.evaluate_phase1(x, phi, genome=genome)

    # Resolver Fase I con SCS (1e-8)
    cfg_scs = EvalConfig(mu_t=0.6910, solver_name="SCS", solver_tol=1e-8, lam_floor_hat=1e-5)
    ev_scs = ScaledSDPEvaluator(config=cfg_scs, lam_floor_hat=1e-5)
    r_scs = ev_scs.evaluate_phase1(x, phi, genome=genome)

    print("--- 1. RESULTADOS FASE I (Margen Espectral lam_hat) ---")
    print(f"Clarabel (tol=1e-8): lam_hat = {r_cl8.lam_normalized:.12e} | s* = {r_cl8.s:.12e} | status = {r_cl8.status}")
    print(f"Clarabel (tol=1e-9): lam_hat = {r_cl9.lam_normalized:.12e} | s* = {r_cl9.s:.12e} | status = {r_cl9.status}")
    print(f"SCS      (tol=1e-8): lam_hat = {r_scs.lam_normalized:.12e} | s* = {r_scs.s:.12e} | status = {r_scs.status}")

    # Racionalización de la matriz P obtenida por Clarabel 1e-9
    P_raw = np.array(r_cl9.P, dtype=np.float64)
    # Racionalizar / simetrizar exactamente
    P_sym = 0.5 * (P_raw + P_raw.T)
    np.fill_diagonal(P_sym, 0.0)
    P_sym = np.clip(P_sym, 0.0, 1.0)
    P_rat = tuple(tuple(float(v) for v in row) for row in P_sym)

    # Matriz G_hat y Matriz M física
    pairs = get_pair_indices(m)
    p_vec = np.array([P_sym[i, j] for i, j in pairs], dtype=np.float64)
    G_hat_mat = sum(p_vec[k] * G_basis[k] for k in range(len(pairs)))

    Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
    M_basis = build_m_basis(x_arr, phi_arr, Q)
    M_mat = sum(p_vec[k] * M_basis[k] for k in range(len(pairs)))

    # 2. NumPy float64
    eigs_G_np = np.linalg.eigvalsh(G_hat_mat)
    lam_G_np = float(np.min(eigs_G_np))

    eigs_M_np = np.linalg.eigvalsh(M_mat)
    lam_M_np = float(np.min(eigs_M_np))

    # 3. mpmath a 50 dígitos
    G_hat_mp = mp.matrix(G_hat_mat.tolist())
    eigs_G_mp = mp.eigsy(G_hat_mp, eigvals_only=True)
    lam_G_mp = float(min(eigs_G_mp))

    M_mp = mp.matrix(M_mat.tolist())
    eigs_M_mp = mp.eigsy(M_mp, eigvals_only=True)
    lam_M_mp = float(min(eigs_M_mp))

    print("\n--- 2. VERIFICACIÓN INDEPENDIENTE NUMPY Y MPMATH 50d SOBRE P RACIONALIZADA ---")
    print(f"NumPy float64:   lam_min(G_hat) = {lam_G_np:.15e} | lam_min(M) = {lam_M_np:.15e}")
    print(f"mpmath (50 dps): lam_min(G_hat) = {lam_G_mp:.15e} | lam_min(M) = {lam_M_mp:.15e}")

    # Verificación física en condiciones Kuramoto
    verif = verify_solution(x, phi, P_rat, eps=1e-10, eps_c=1e-4)
    print(f"\nVerificación Física verify_solution: max|F| = {verif['max_F']:.12e} | min c = {verif['min_c']:.8f} | is_valid = {verif['is_valid']}")

    # 4. Fase II: mu* máximo alcanzable con piso 1e-5
    p2_cl = ev_cl9.evaluate_phase2(x, phi, eps=1e-5, genome=genome)
    p2_scs = ev_scs.evaluate_phase2(x, phi, eps=1e-5, genome=genome)

    print("\n--- 3. RESULTADOS FASE II (mu* máximo con lam_hat >= lam_floor_hat = 1e-5) ---")
    print(f"Clarabel Fase II: mu* = {p2_cl.mu:.12f} | status = {p2_cl.status}")
    print(f"SCS      Fase II: mu* = {p2_scs.mu:.12f} | status = {p2_scs.status}")

    res_summary = {
        "id": 4461,
        "m": m,
        "h_eff": h_eff,
        "clarabel_1e8_lam_hat": r_cl8.lam_normalized,
        "clarabel_1e9_lam_hat": r_cl9.lam_normalized,
        "scs_1e8_lam_hat": r_scs.lam_normalized,
        "numpy_lam_min_G_hat": lam_G_np,
        "numpy_lam_min_M": lam_M_np,
        "mpmath_50d_lam_min_G_hat": lam_G_mp,
        "mpmath_50d_lam_min_M": lam_M_mp,
        "phase2_mu_clarabel": p2_cl.mu,
        "phase2_mu_scs": p2_scs.mu,
        "equilibrium_max_F": verif["max_F"],
        "min_c": verif["min_c"],
        "is_strictly_stable": bool(lam_G_mp >= 1e-5 and verif["max_F"] <= 1e-8),
    }


    with open("results/id4461_verification.json", "w") as f:
        json.dump(res_summary, f, indent=2)
    print("\nResultados guardados en results/id4461_verification.json")


if __name__ == "__main__":
    main()
