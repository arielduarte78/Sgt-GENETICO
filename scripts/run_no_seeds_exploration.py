"""Experimento Item 3:
Exploración sin semillas (use_seeds=False) con 3 semillas estocásticas distintas:
- mu_t = 0.6915 (seeds: 42, 123, 999)
- mu_t = 0.6920 (seeds: 42, 123, 999)
Control separado: Semilla S7 escalable.

Reporta mu* máximo con lam_hat >= piso por semilla y cuántas veces reencuentra m=7 con estructura S7.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    compute_sheets,
    compute_h_eff,
)
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine
from kuramoto_templates.search.seeds import create_s7_direct_seed
from kuramoto_templates.search.genome import Genome


def is_s7_structure(ind_dict: dict) -> bool:
    """Verifica si un individuo posee estructura topológica S7:

    - m = 7 clusters
    - Agrupamiento en 4 hojas con distribución de tamaños {2, 2, 2, 1}
    - Escala h_eff en [1e-4, 0.05]
    - Margen espectral lam_hat >= lam_floor_hat (1e-5)
    """
    g = ind_dict.get("genome", {})
    m = g.get("m")
    if m != 7:
        return False

    lam_norm = ind_dict.get("lam_normalized", ind_dict.get("lam_star", -1.0))
    if lam_norm is None or lam_norm < 1e-5:
        return False

    x = tuple(ind_dict["phenotype_x"])
    phi = tuple(ind_dict["phenotype_phi"])
    genome_obj = Genome(
        m=g["m"],
        K_ref=g["K_ref"],
        k=tuple(g["k"]),
        sgn=tuple(g["sgn"]),
        ell=tuple(g["ell"]),
        u=tuple(g["u"]),
    )

    sheets, _ = compute_sheets(x, phi, genome=genome_obj)
    sheet_sizes = sorted([len(idxs) for idxs in sheets.values()])
    if sheet_sizes != [1, 2, 2, 2]:
        return False

    h_eff = compute_h_eff(x, phi, sheets)
    if not (1e-4 <= h_eff <= 0.05):
        return False

    return True


def run_single(mu_t: float, seed: int):
    run_id = f"no_seeds_mu_{mu_t:.4f}_seed_{seed}"
    out_dir = Path(f"runs/{run_id}")
    print(f"\n>>> Lanzando corrida: {run_id} (mu_t={mu_t}, seed={seed}, 50 gens, use_seeds=False)")

    cfg = GAConfig(
        mu_t=mu_t,
        num_islands=4,
        island_pop_size=20,  # pop = 80
        num_generations=50,
        checkpoint_interval=10,
        use_seeds=False,
        seed=seed,
        evaluator="scaled",
        lam_floor_hat=1e-5,
    )

    t0 = time.perf_counter()
    engine = GAEngine(config=cfg, output_dir=out_dir, allow_dirty=True)
    engine.run()
    elapsed = time.perf_counter() - t0

    # Analizar resultados
    ckpt_path = out_dir / "checkpoints" / "checkpoint_gen_50.json"
    status_path = out_dir / "status.json"

    with open(ckpt_path) as f:
        ckpt_data = json.load(f)
    with open(status_path) as f:
        status_data = json.load(f)

    hof_global = ckpt_data.get("hall_of_fame_global", [])
    cache = ckpt_data.get("cache", {})

    # Evaluar Fase II en el mejor del HoF si es estable
    ev_scaled = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, lam_floor_hat=1e-5))

    max_mu_star = None
    best_stable_ind = None

    # Buscar mejor mu* entre todos los estables del HoF
    for ind in hof_global:
        lam_hat = ind.get("lam_normalized", ind.get("lam_star", 0.0))
        if lam_hat is not None and lam_hat >= 1e-5:
            # Evaluar Fase II
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
                    best_stable_ind = ind

    # Contar cuántas veces se reencuentra m=7 con estructura S7
    # (Tanto en el HoF como en todo el historial evaluado en la cache)
    s7_rediscovered_count_hof = sum(1 for ind in hof_global if is_s7_structure(ind))

    # Conteo en cache
    s7_rediscovered_count_cache = 0
    for key, entry in cache.items():
        lam_hat = entry.get("lam_normalized", entry.get("lam_star", -1.0))
        if lam_hat >= 1e-5:
            parsed = json.loads(key)
            if len(parsed) == 21:  # m = 7 (3 * 7 = 21: x, phi, k)
                k_list = parsed[14:21]
                counts = sorted([k_list.count(c) for c in set(k_list)])
                if counts == [1, 2, 2, 2]:
                    s7_rediscovered_count_cache += 1

    winner = hof_global[0] if hof_global else None
    result = {
        "run_id": run_id,
        "mu_t": mu_t,
        "seed": seed,
        "elapsed_seconds": elapsed,
        "total_evaluations": status_data.get("total_evaluations", len(cache)),
        "cache_hits": status_data.get("cache_hits", 0),
        "hof_length": len(hof_global),
        "hof_ge_floor_count": sum(1 for ind in hof_global if (ind.get("lam_normalized") or 0.0) >= 1e-5),
        "max_mu_star": max_mu_star,
        "s7_rediscovered_hof": s7_rediscovered_count_hof,
        "s7_rediscovered_cache": s7_rediscovered_count_cache,
        "winner": {
            "id": winner.get("id") if winner else None,
            "m": winner.get("genome", {}).get("m") if winner else None,
            "f": winner.get("f") if winner else None,
            "s_star": winner.get("s_star") if winner else None,
            "lam_hat": winner.get("lam_normalized") if winner else None,
            "mu_star": max_mu_star,
        } if winner else None,
    }

    print(f"  Corrida {run_id} finalizada en {elapsed:.1f} s")
    print(f"  HoF ge piso (1e-5): {result['hof_ge_floor_count']}/{len(hof_global)}")
    print(f"  Max mu* (Fase II): {max_mu_star}")
    print(f"  Reencuentros estructura S7: {s7_rediscovered_count_hof} en HoF, {s7_rediscovered_count_cache} en Cache")
    return result


def main():
    print("=" * 75)
    print("CONTROL SEPARADO S7 ESCALABLE")
    print("=" * 75)
    s7 = create_s7_direct_seed(h=0.002, t=0.0)
    for mu in [0.6915, 0.6920]:
        ev = ScaledSDPEvaluator(EvalConfig(mu_t=mu, solver_tol=1e-8, lam_floor_hat=1e-5))
        r1 = ev.evaluate_phase1(s7.phenotype_x, s7.phenotype_phi, genome=s7.genome)
        r2 = ev.evaluate_phase2(s7.phenotype_x, s7.phenotype_phi, eps=1e-5, genome=s7.genome)
        print(f"S7 Control @ mu_t = {mu:.4f}: f={r1.f:.8e}, s*={r1.s:.8e}, lam_hat={r1.lam_normalized:.8e}, mu*(Fase II)={r2.mu:.8f}")

    print("\n" + "=" * 75)
    print("EXPLORACIÓN SIN SEMILLAS (USE_SEEDS = FALSE): 3 SEMILLAS X 2 NIVELES DE MU_T")
    print("=" * 75)

    seeds = [42, 123, 999]
    mu_targets = [0.6915, 0.6920]

    all_results = []
    for mu_t in mu_targets:
        for seed in seeds:
            res = run_single(mu_t=mu_t, seed=seed)
            all_results.append(res)

    out_file = Path("results/no_seeds_exploration_summary.json")
    with open(out_file, "w") as f:
        json.dump(all_results, f, indent=2)

    print("\n" + "=" * 75)
    print("TABLA RESUMEN DE RESULTADOS SIN SEMILLAS (ITEM 3)")
    print("=" * 75)
    print(f"{'mu_t':<8} | {'seed':<6} | {'Ganador m':<10} | {'lam_hat':<15} | {'mu* (Fase II)':<15} | {'Reencuentra S7'}")
    print("-" * 75)
    for r in all_results:
        w = r.get("winner", {})
        m_val = w.get("m") if w else "-"
        lam_val = f"{w.get('lam_hat', 0.0):.6e}" if w and w.get('lam_hat') is not None else "-"
        mu_val = f"{r.get('max_mu_star'):.8f}" if r.get('max_mu_star') is not None else "None (< piso)"
        s7_found = f"{r['s7_rediscovered_hof']} en HoF ({r['s7_rediscovered_cache']} en cache)"
        print(f"{r['mu_t']:<8.4f} | {r['seed']:<6d} | {str(m_val):<10} | {lam_val:<15} | {mu_val:<15} | {s7_found}")
    print("=" * 75)


if __name__ == "__main__":
    main()
