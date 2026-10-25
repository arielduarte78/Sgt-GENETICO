"""Verificación rápida de los 6 campeones de las corridas con semillas:
- mu_t = 0.6915 (seeds: 42, 123, 999)
- mu_t = 0.6920 (seeds: 42, 123, 999)

TAREAS:
1. Extrae (x, phi) exactos del mejor individuo de cada corrida.
2. Verificación estricta:
   - verify_solution (NumPy)
   - lam_min(M) y lam_min(G_hat) con mpmath a 50 dígitos
   - Resolución con SCS
3. Curva de Fase II (mu* vs piso):
   - evaluate_phase2 en pisos 1e-5, 1e-4 y 1e-3 (Clarabel y SCS)
4. Distancia gauge-alineada respecto a S7(h=0.002, t=0).
5. Criterio de Falso Positivo (Bloqueante):
   - Alarma gigante si SCS y Clarabel difieren en mu* > 1e-5 o si lam_min(M) < 0 en mpmath.
6. Guarda en results/mu_vs_floor.json.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np
import scipy.linalg
from scipy.optimize import linear_sum_assignment
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
    get_pair_indices,
    verify_solution,
)
from kuramoto_templates.search.genome import Genome, wrap_minus_pi_to_pi
from kuramoto_templates.cases.s7 import get_s7_appendix_a_template


def compute_gauge_aligned_s7_distance(x, phi, x_s7, phi_s7) -> float:
    """Calcula la distancia cuadrática mínima en (x, phi) con alineación circular global y permutación óptima."""
    x_arr = np.asarray(x, dtype=np.float64)
    phi_arr = np.asarray(phi, dtype=np.float64)
    best_dist = float("inf")
    for i in range(7):
        alpha = wrap_minus_pi_to_pi(phi_s7[6] - phi_arr[i])
        phi_rot = np.array([wrap_minus_pi_to_pi(p + alpha) for p in phi_arr])
        cost = np.zeros((7, 7))
        for r in range(7):
            for c in range(7):
                cost[r, c] = (x_arr[r] - x_s7[c]) ** 2 + (wrap_minus_pi_to_pi(phi_rot[r] - phi_s7[c])) ** 2
        row_ind, col_ind = linear_sum_assignment(cost)
        total_sq = sum(cost[r, c] for r, c in zip(row_ind, col_ind))
        dist = float(np.sqrt(total_sq))
        if dist < best_dist:
            best_dist = dist
    return best_dist


def main():
    mp.mp.dps = 50

    s7_ref = get_s7_appendix_a_template(h=0.002, t=0.0)
    x_s7 = np.array(s7_ref.x, dtype=np.float64)
    phi_s7 = np.array([wrap_minus_pi_to_pi(p) for p in s7_ref.phi], dtype=np.float64)

    runs = [
        ("with_seeds_mu_0.6915_seed_42", 0.6915, 42),
        ("with_seeds_mu_0.6915_seed_123", 0.6915, 123),
        ("with_seeds_mu_0.6915_seed_999", 0.6915, 999),
        ("with_seeds_mu_0.6920_seed_42", 0.6920, 42),
        ("with_seeds_mu_0.6920_seed_123", 0.6920, 123),
        ("with_seeds_mu_0.6920_seed_999", 0.6920, 999),
    ]

    floors = [1e-5, 1e-4, 1e-3]
    all_results = []
    global_alarms = []

    print("=" * 95)
    print("VERIFICACIÓN INDEPENDIENTE DE LOS 6 CAMPEONES CON SEMILLAS (PRE-CORRIDA NOCTURNA)")
    print("=" * 95)

    for run_name, mu_t, seed in runs:
        ckpt_path = Path("runs") / run_name / "checkpoints" / "checkpoint_gen_50.json"
        if not ckpt_path.exists():
            print(f"ERROR: No se encontró {ckpt_path}")
            continue

        with open(ckpt_path) as f:
            ckpt_data = json.load(f)

        hof = ckpt_data.get("hall_of_fame_global", [])
        if not hof:
            print(f"ERROR: Hall of fame vacío en {run_name}")
            continue

        champ = hof[0]
        cid = champ["id"]
        x = tuple(float(v) for v in champ["phenotype_x"])
        phi = tuple(float(v) for v in champ["phenotype_phi"])
        g_dict = champ["genome"]
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
        Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
        M_basis = build_m_basis(x_arr, phi_arr, Q)
        pairs = get_pair_indices(m)

        # Distancia gauge-aligned a S7
        d_s7 = compute_gauge_aligned_s7_distance(x_arr, phi_arr, x_s7, phi_s7)

        # 1. Resolver Fase II con piso base 1e-5 para verificación estricta de P
        ev_cl_base = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_name="CLARABEL", solver_tol=1e-9, lam_floor_hat=1e-5))
        r2_base = ev_cl_base.evaluate_phase2(x, phi, eps=1e-5, genome=genome)
        P_sol = r2_base.P

        # a) verify_solution en NumPy puro
        verif = verify_solution(x, phi, P_sol, eps=1e-10, eps_c=1e-4)

        # b) mpmath a 50 dígitos sobre G_hat y M
        p_vec = np.array([P_sol[i][j] for i, j in pairs], dtype=np.float64)
        G_mat = sum(p_vec[k] * G_basis[k] for k in range(len(pairs)))
        M_mat = sum(p_vec[k] * M_basis[k] for k in range(len(pairs)))

        eigs_G_mp = mp.eigsy(mp.matrix(G_mat.tolist()), eigvals_only=True)
        eigs_M_mp = mp.eigsy(mp.matrix(M_mat.tolist()), eigvals_only=True)
        lam_G_mp = float(min(eigs_G_mp))
        lam_M_mp = float(min(eigs_M_mp))

        # c) Curva de Fase II (mu* vs floor) para pisos 1e-5, 1e-4 y 1e-3
        curve_data = {}
        champ_alarms = []

        if lam_M_mp < 0.0:
            champ_alarms.append(f"lam_min(M) en mpmath es NEGATIVO ({lam_M_mp:.3e})")

        for floor in floors:
            # Clarabel
            ev_cl = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_name="CLARABEL", solver_tol=1e-9, lam_floor_hat=floor))
            r_cl = ev_cl.evaluate_phase2(x, phi, eps=floor, genome=genome)
            mu_cl = float(r_cl.mu) if r_cl.mu is not None else None

            # SCS
            ev_scs = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_name="SCS", solver_tol=1e-8, lam_floor_hat=floor))
            r_scs = ev_scs.evaluate_phase2(x, phi, eps=floor, genome=genome)
            mu_scs = float(r_scs.mu) if r_scs.mu is not None else None

            diff_mu = abs(mu_cl - mu_scs) if (mu_cl is not None and mu_scs is not None) else None

            if diff_mu is not None and diff_mu > 1e-5:
                champ_alarms.append(
                    f"Piso {floor:.0e}: Clarabel mu*={mu_cl:.8f} vs SCS mu*={mu_scs:.8f} (diff={diff_mu:.2e} > 1e-5, status_scs={r_scs.status})"
                )

            curve_data[str(floor)] = {
                "floor": floor,
                "clarabel_mu": mu_cl,
                "clarabel_status": r_cl.status,
                "clarabel_lam_hat": floor,
                "scs_mu": mu_scs,
                "scs_status": r_scs.status,
                "diff_mu": diff_mu,
            }

        is_alarm = len(champ_alarms) > 0
        if is_alarm:
            global_alarms.append((run_name, cid, champ_alarms))

        champ_record = {
            "run_id": run_name,
            "mu_t": mu_t,
            "seed": seed,
            "champion_id": cid,
            "m": m,
            "h_eff": h_eff,
            "d_s7_gauge_aligned": d_s7,
            "verify_solution": {
                "is_valid": verif["is_valid"],
                "max_F": verif["max_F"],
                "min_c": verif["min_c"],
            },
            "mpmath_50d": {
                "lam_min_G_hat": lam_G_mp,
                "lam_min_M": lam_M_mp,
            },
            "phase2_curve": curve_data,
            "is_alarm": is_alarm,
            "alarm_messages": champ_alarms,
        }
        all_results.append(champ_record)

        print(f"\n>>> [{run_name}] Ganador ID {cid} (m={m}, h_eff={h_eff:.4e}, d_S7={d_s7:.6f}):")
        print(f"    Física (NumPy):  is_valid={verif['is_valid']} | max|F|={verif['max_F']:.2e} | min c_a={verif['min_c']:.6f}")
        print(f"    mpmath (50 dps): lam_min(M)={lam_M_mp:.6e} | lam_min(G_hat)={lam_G_mp:.6e}")
        print("    Curva Fase II (mu* vs piso):")
        for f_key in ["1e-05", "0.0001", "0.001"]:
            c_info = curve_data[str(float(f_key))]
            d_str = f"{c_info['diff_mu']:.2e}" if c_info['diff_mu'] is not None else "N/A"
            print(f"      Piso {float(f_key):<7.0e}: Clarabel mu* = {c_info['clarabel_mu']:.8f} | SCS mu* = {c_info['scs_mu']:.8f} | diff = {d_str} ({c_info['scs_status']})")

        if is_alarm:
            print("    [!] ALERTAS DEL CANDIDATO:")
            for a_msg in champ_alarms:
                print(f"       -> {a_msg}")

    # Guardar resultados
    out_file = Path("results/mu_vs_floor.json")
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2)

    # 2. Criterio de Falso Positivo (Bloqueante)
    print("\n" + "=" * 95)
    print("CRITERIO DE FALSO POSITIVO (BLOQUEANTE)")
    print("=" * 95)
    if global_alarms:
        print("\n" + "!" * 80)
        print("!!! ALARMA CRITICA DE FALSO POSITIVO DETECTADA !!!")
        print("!" * 80)
        for r_name, cid, msgs in global_alarms:
            print(f"\n[Corrida {r_name} | ID {cid}]:")
            for m_txt in msgs:
                print(f"  [ALERTA] {m_txt}")
        print("\n" + "!" * 80 + "\n")
    else:
        print("[OK] NINGUNA ALARMA: Todos los candidatos tienen lam_min(M) > 0 y SCS concuerda con Clarabel dentro de 1e-5.")

    # Resumen tabular solicitado
    print("\n" + "=" * 95)
    print("TABLA COMPARATIVA FINAL: 6 CAMPEONES Y SU CURVA DE PISOS")
    print("=" * 95)
    print(f"{'Corrida':<30} | {'ID':<5} | {'mu*(1e-5)':<12} | {'mu*(1e-4)':<12} | {'mu*(1e-3)':<12} | {'lam_min(M)':<12} | {'d_S7':<10}")
    print("-" * 95)
    max_mu_1e4 = 0.0
    best_champ_1e4 = None
    for res in all_results:
        p2 = res["phase2_curve"]
        m_1e5 = p2["1e-05"]["clarabel_mu"]
        m_1e4 = p2["0.0001"]["clarabel_mu"]
        m_1e3 = p2["0.001"]["clarabel_mu"]
        m_mp = res["mpmath_50d"]["lam_min_M"]
        d_val = res["d_s7_gauge_aligned"]
        if m_1e4 > max_mu_1e4:
            max_mu_1e4 = m_1e4
            best_champ_1e4 = (res["run_id"], res["champion_id"])
        print(f"{res['run_id']:<30} | {res['champion_id']:<5d} | {m_1e5:<12.8f} | {m_1e4:<12.8f} | {m_1e3:<12.8f} | {m_mp:<12.3e} | {d_val:<10.6f}")
    print("=" * 95)
    print(f"\n>>> VALOR MÁXIMO DE mu* LOGRADO CON EL PISO MÁS CONSERVADOR (1e-4):")
    print(f"    mu* = {max_mu_1e4:.8f} (Corrida: {best_champ_1e4[0]}, Campeón ID {best_champ_1e4[1]})")
    print("=" * 95)


if __name__ == "__main__":
    main()
