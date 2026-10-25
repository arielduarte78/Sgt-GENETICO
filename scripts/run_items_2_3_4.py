"""Script integral para ejecutar los Items 1, 2, 3 y 4:

1. Verificación de robustez de ID 4461:
   - P exacta de Phase I y Phase II
   - P redondeada a 1e-10, 1e-8, 1e-6
   - lam_min(M) y lam_hat en cada caso
   - Distancia en (x, phi) directa y gauge-aligned a S7 directa.

2. 3 semillas x 50 gens CON semilla S7 a mu_t=0.6915 y 0.6920 (use_seeds=True):
   - Reporta mu* maximo y si supera 0.69132.

3. 3 semillas x 50 gens SIN semillas con m restringido en [5, 8] a mu_t=0.6915 y 0.6920:
   - Reporta mu* maximo por seed.

4. Criterio cuantitativo de estructura S7:
   - Distancia gauge-aligned d_{S7}(x, phi) a la semilla S7 (h=0.002, t=0).
   - Recuento con umbrales d <= 0.01, 0.05, 0.10, 0.25 en todas las corridas.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time
import numpy as np
import scipy.linalg
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    compute_sheets,
    compute_h_eff,
    build_scaled_basis_t,
    build_g_hat_basis,
)
from kuramoto_templates.domain.conditions import build_m_basis, get_pair_indices, verify_solution
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine
from kuramoto_templates.search.genome import Genome, wrap_minus_pi_to_pi
from kuramoto_templates.cases.s7 import get_s7_appendix_a_template


# ==============================================================================
# ITEM 4: Función cuantitativa de distancia a S7
# ==============================================================================
def get_s7_reference():
    s7_ref = get_s7_appendix_a_template(h=0.002, t=0.0)
    x_s7 = np.array(s7_ref.x, dtype=np.float64)
    phi_s7 = np.array([wrap_minus_pi_to_pi(p) for p in s7_ref.phi], dtype=np.float64)
    return x_s7, phi_s7


def compute_gauge_aligned_s7_distance(x, phi, x_s7, phi_s7) -> float:
    """Calcula la distancia mínima cuadrática en (x, phi) bajo rotación global U(1)

    y permutación de clusters (Hungarian algorithm):
    d = min_{alpha, sigma} sqrt(||x_sigma - x_s7||_2^2 + ||dist_S1(phi_sigma + alpha, phi_s7)||_2^2)
    """
    if len(x) != 7:
        return float("inf")

    x_arr = np.asarray(x, dtype=np.float64)
    phi_arr = np.asarray(phi, dtype=np.float64)

    best_dist = float("inf")
    # Probar alinear cada cluster i al cluster de la hoja solitaria (idx 6 en s7)
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


# ==============================================================================
# ITEM 1: Robustez de ID 4461
# ==============================================================================
def run_item_1():
    print("=" * 80)
    print("ITEM 1: VERIFICACIÓN DE ROBUSTEZ ID 4461")
    print("=" * 80)
    ckpt = json.load(open("runs/run_50gen_mu_06910_frontier/checkpoints/checkpoint_gen_50.json"))
    ind = next(i for i in ckpt["hall_of_fame_global"] if i["id"] == 4461)
    x = tuple(ind["phenotype_x"])
    phi = tuple(ind["phenotype_phi"])
    g = ind["genome"]
    genome = Genome(m=g["m"], K_ref=g["K_ref"], k=tuple(g["k"]), sgn=tuple(g["sgn"]), ell=tuple(g["ell"]), u=tuple(g["u"]))

    ev = ScaledSDPEvaluator(EvalConfig(mu_t=0.6910, solver_name="CLARABEL", solver_tol=1e-9, lam_floor_hat=1e-5))
    r1 = ev.evaluate_phase1(x, phi, genome=genome)
    r2 = ev.evaluate_phase2(x, phi, eps=1e-5, genome=genome)

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

    def analyze_P(P_raw, name):
        results = {}
        for round_level, round_name in [(None, "exact"), (10, "1e-10"), (8, "1e-8"), (6, "1e-6")]:
            if round_level is None:
                P_mat = np.array(P_raw, dtype=np.float64)
            else:
                P_mat = np.round(np.array(P_raw, dtype=np.float64), round_level)

            p_vec = np.array([P_mat[i, j] for i, j in pairs], dtype=np.float64)
            G_mat = sum(p_vec[k] * G_basis[k] for k in range(len(pairs)))
            M_mat = sum(p_vec[k] * M_basis[k] for k in range(len(pairs)))
            lam_G = float(np.min(np.linalg.eigvalsh(G_mat)))
            lam_M = float(np.min(np.linalg.eigvalsh(M_mat)))
            results[round_name] = {
                "lam_min_M": lam_M,
                "lam_hat": lam_G,
            }
            print(f"  {name} [redondeo={round_name:<6}]: lam_min(M) = {lam_M:16.12e} | lam_hat = {lam_G:16.12e}")
        return results

    print("\n--- A. Matriz P de Fase I (Clarabel 1e-9) ---")
    p1_results = analyze_P(r1.P, "Phase I P")

    print("\n--- B. Matriz P de Fase II (optimizando mu* con lam_hat >= 1e-5) ---")
    p2_results = analyze_P(r2.P, "Phase II P")

    # Distancias a S7
    x_s7, phi_s7 = get_s7_reference()
    # Directa
    d_x_dir = float(np.linalg.norm(x_arr - x_s7))
    d_phi_dir = float(np.linalg.norm([wrap_minus_pi_to_pi(p1 - p2) for p1, p2 in zip(phi_arr, phi_s7)]))
    # Gauge-aligned
    d_gauge = compute_gauge_aligned_s7_distance(x_arr, phi_arr, x_s7, phi_s7)

    print("\n--- C. Distancia a S7 directa (Appendix A, h=0.002, t=0) ---")
    print(f"  Distancia directa:       ||dx||_2 = {d_x_dir:.6f} | ||dphi||_2 = {d_phi_dir:.6f}")
    print(f"  Distancia gauge-aligned: d_S7(x, phi) = {d_gauge:.8f}")

    item1_data = {
        "p1_results": p1_results,
        "p2_results": p2_results,
        "distance_direct": {"dx": d_x_dir, "dphi": d_phi_dir},
        "distance_gauge_aligned": d_gauge,
    }
    with open("results/id4461_robustness.json", "w") as f:
        json.dump(item1_data, f, indent=2)
    return item1_data


# ==============================================================================
# Helper de ejecución de GA
# ==============================================================================
def run_ga_experiment(run_id: str, cfg: GAConfig, x_s7, phi_s7) -> dict:
    out_dir = Path(f"runs/{run_id}")
    print(f"\n>>> Lanzando corrida: {run_id} (mu_t={cfg.mu_t}, seed={cfg.seed}, use_seeds={cfg.use_seeds}, m in [{cfg.min_m}, {cfg.max_m}])")
    t0 = time.perf_counter()
    engine = GAEngine(config=cfg, output_dir=out_dir, allow_dirty=True)
    engine.run()
    elapsed = time.perf_counter() - t0

    ckpt_path = out_dir / "checkpoints" / "checkpoint_gen_50.json"
    status_path = out_dir / "status.json"

    with open(ckpt_path) as f:
        ckpt_data = json.load(f)
    with open(status_path) as f:
        status_data = json.load(f)

    hof_global = ckpt_data.get("hall_of_fame_global", [])
    cache = ckpt_data.get("cache", {})

    ev_scaled = ScaledSDPEvaluator(EvalConfig(mu_t=cfg.mu_t, solver_tol=1e-8, lam_floor_hat=1e-5))

    max_mu_star = None
    best_ind = None

    for ind in hof_global:
        lam_hat = ind.get("lam_normalized", ind.get("lam_star", 0.0))
        if lam_hat is not None and lam_hat >= 1e-5:
            x = tuple(ind["phenotype_x"])
            phi = tuple(ind["phenotype_phi"])
            g = ind["genome"]
            g_obj = Genome(
                m=g["m"],
                K_ref=g["K_ref"],
                k=tuple(g["k"]),
                sgn=tuple(g["sgn"]),
                ell=tuple(g["ell"]),
                u=tuple(g["u"]),
            )
            r_p2 = ev_scaled.evaluate_phase2(x, phi, eps=1e-5, genome=g_obj)
            mu_val = r_p2.mu
            if mu_val is not None:
                if max_mu_star is None or mu_val > max_mu_star:
                    max_mu_star = mu_val
                    best_ind = ind

    # Métricas cuantitativas S7 en cache
    distances_m7 = []
    for key, entry in cache.items():
        parsed = json.loads(key)
        if len(parsed) == 21:  # m = 7
            x_cand = parsed[0:7]
            phi_cand = parsed[7:14]
            d = compute_gauge_aligned_s7_distance(x_cand, phi_cand, x_s7, phi_s7)
            distances_m7.append(d)

    d_arr = np.array(distances_m7) if distances_m7 else np.array([])
    min_dist_s7 = float(np.min(d_arr)) if len(d_arr) > 0 else float("inf")
    cnt_001 = int(np.sum(d_arr <= 0.01))
    cnt_005 = int(np.sum(d_arr <= 0.05))
    cnt_010 = int(np.sum(d_arr <= 0.10))
    cnt_025 = int(np.sum(d_arr <= 0.25))

    winner = hof_global[0] if hof_global else None
    res = {
        "run_id": run_id,
        "mu_t": cfg.mu_t,
        "seed": cfg.seed,
        "use_seeds": cfg.use_seeds,
        "min_m": cfg.min_m,
        "max_m": cfg.max_m,
        "elapsed_seconds": elapsed,
        "total_evaluations": status_data.get("total_evaluations", len(cache)),
        "cache_hits": status_data.get("cache_hits", 0),
        "hof_length": len(hof_global),
        "hof_ge_floor_count": sum(1 for ind in hof_global if (ind.get("lam_normalized") or 0.0) >= 1e-5),
        "max_mu_star": max_mu_star,
        "beats_4461": bool(max_mu_star is not None and max_mu_star > 0.69132),
        "s7_distance_metrics": {
            "m7_evaluations": len(d_arr),
            "min_distance_to_s7": min_dist_s7,
            "count_le_0.01": cnt_001,
            "count_le_0.05": cnt_005,
            "count_le_0.10": cnt_010,
            "count_le_0.25": cnt_025,
        },
        "winner": {
            "id": winner.get("id") if winner else None,
            "m": winner.get("genome", {}).get("m") if winner else None,
            "f": winner.get("f") if winner else None,
            "s_star": winner.get("s_star") if winner else None,
            "lam_hat": winner.get("lam_normalized") if winner else None,
            "mu_star": max_mu_star,
        } if winner else None,
    }

    print(f"  Finalizada en {elapsed:.1f} s | Max mu*: {max_mu_star} | Supera 0.69132: {res['beats_4461']}")
    print(f"  Distancia mínima a S7: {min_dist_s7:.6f} | <=0.01: {cnt_001}, <=0.05: {cnt_005}, <=0.10: {cnt_010}")
    return res


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    x_s7, phi_s7 = get_s7_reference()

    # 1. Item 1
    run_item_1()

    # 2. Item 2: CON semilla S7 (use_seeds=True) a mu_t = 0.6915 y 0.6920
    print("\n" + "=" * 80)
    print("ITEM 2: 3 SEMILLAS X 50 GENS CON SEMILLA S7 (USE_SEEDS=TRUE)")
    print("=" * 80)
    seeds = [42, 123, 999]
    item2_results = []
    for mu_t in [0.6915, 0.6920]:
        for seed in seeds:
            run_id = f"with_seeds_mu_{mu_t:.4f}_seed_{seed}"
            cfg = GAConfig(
                mu_t=mu_t,
                num_islands=4,
                island_pop_size=20,
                num_generations=50,
                checkpoint_interval=10,
                use_seeds=True,
                seed=seed,
                evaluator="scaled",
                lam_floor_hat=1e-5,
            )
            r = run_ga_experiment(run_id, cfg, x_s7, phi_s7)
            item2_results.append(r)

    with open("results/item2_with_seeds_summary.json", "w") as f:
        json.dump(item2_results, f, indent=2)

    # 3. Item 3: SIN semillas, m en [5, 8] a mu_t = 0.6915 y 0.6920
    print("\n" + "=" * 80)
    print("ITEM 3: 3 SEMILLAS X 50 GENS SIN SEMILLAS CON m EN [5, 8] (USE_SEEDS=FALSE)")
    print("=" * 80)
    item3_results = []
    for mu_t in [0.6915, 0.6920]:
        for seed in seeds:
            run_id = f"no_seeds_m5to8_mu_{mu_t:.4f}_seed_{seed}"
            cfg = GAConfig(
                mu_t=mu_t,
                num_islands=4,
                island_pop_size=20,
                num_generations=50,
                checkpoint_interval=10,
                use_seeds=False,
                seed=seed,
                evaluator="scaled",
                lam_floor_hat=1e-5,
                min_m=5,
                max_m=8,
            )
            r = run_ga_experiment(run_id, cfg, x_s7, phi_s7)
            item3_results.append(r)

    with open("results/item3_no_seeds_m5to8_summary.json", "w") as f:
        json.dump(item3_results, f, indent=2)

    # Resumen final consolidado
    print("\n" + "=" * 80)
    print("TABLA RESUMEN ITEM 2 (CON SEMILLA S7)")
    print("=" * 80)
    print(f"{'mu_t':<8} | {'seed':<6} | {'Ganador m':<10} | {'lam_hat':<15} | {'mu* (Fase II)':<15} | {'Supera 0.69132'}")
    print("-" * 80)
    for r in item2_results:
        w = r.get("winner", {})
        m_val = w.get("m") if w else "-"
        lam_val = f"{w.get('lam_hat', 0.0):.6e}" if w and w.get('lam_hat') is not None else "-"
        mu_val = f"{r.get('max_mu_star'):.8f}" if r.get('max_mu_star') is not None else "None"
        beats = "SÍ" if r["beats_4461"] else "NO"
        print(f"{r['mu_t']:<8.4f} | {r['seed']:<6d} | {str(m_val):<10} | {lam_val:<15} | {mu_val:<15} | {beats}")

    print("\n" + "=" * 80)
    print("TABLA RESUMEN ITEM 3 (SIN SEMILLAS, m EN [5, 8])")
    print("=" * 80)
    print(f"{'mu_t':<8} | {'seed':<6} | {'Ganador m':<10} | {'lam_hat':<15} | {'mu* (Fase II)':<15} | {'d_min(S7)':<12} | {'<=0.05'}")
    print("-" * 80)
    for r in item3_results:
        w = r.get("winner", {})
        m_val = w.get("m") if w else "-"
        lam_val = f"{w.get('lam_hat', 0.0):.6e}" if w and w.get('lam_hat') is not None else "-"
        mu_val = f"{r.get('max_mu_star'):.8f}" if r.get('max_mu_star') is not None else "None"
        d_min = f"{r['s7_distance_metrics']['min_distance_to_s7']:.4f}"
        cnt_05 = r['s7_distance_metrics']['count_le_0.05']
        print(f"{r['mu_t']:<8.4f} | {r['seed']:<6d} | {str(m_val):<10} | {lam_val:<15} | {mu_val:<15} | {d_min:<12} | {cnt_05}")
    print("=" * 80)


if __name__ == "__main__":
    main()
