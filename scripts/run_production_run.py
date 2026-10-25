"""Execute the production run specified in Punto 5:
mu_t=0.6929077, 8 islands, mu=64, 200 generations, with seeds.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

# Ensure src is in python path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine


def main():
    mu_t = 0.6929077
    out_dir = Path("runs/prod_mu_0.6929077_8isl_mu64_gen200")
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = GAConfig(
        mu_t=mu_t,
        num_islands=8,
        island_pop_size=64,
        num_generations=200,
        checkpoint_interval=10,
        use_seeds=True,
        seed=42,
    )
    cfg_path = out_dir / "config.yaml"
    cfg.to_yaml(cfg_path)

    print("=================================================================")
    print(f"Starting Production GA Run:")
    print(f"  mu_t:            {mu_t}")
    print(f"  num_islands:     8")
    print(f"  island_pop_size: 64")
    print(f"  generations:     200")
    print(f"  seeds enabled:   True")
    print(f"  output dir:      {out_dir}")
    print("=================================================================", flush=True)

    t0 = time.perf_counter()
    engine = GAEngine(config=cfg, output_dir=out_dir, allow_dirty=True)

    # Wrap generation progress reporting if possible or run directly
    res = engine.run()
    t1 = time.perf_counter()
    total_elapsed = t1 - t0

    print("\n=================================================================")
    print(f"Production Run Completed in {total_elapsed:.2f} seconds ({total_elapsed/60:.2f} min)")
    print(f"Total Unique Evaluations (cache size): {len(engine.cache)}")
    print(f"Total Cache Hits:                      {engine.cache_hits}")
    print(f"Global Best Fitness:                   {res['best_f']:.12e} (m={res['best_m']})")
    print(f"Global Best lam*:                      {res['best_lam']:.12e}")
    print("=================================================================\n", flush=True)

    # Extract evolution of best lam* per generation from run.jsonl
    log_file = out_dir / "run.jsonl"
    gen_evolution = {}
    if log_file.exists():
        with open(log_file, "r", encoding="utf-8") as f:
            for line in f:
                rec = json.loads(line)
                if rec.get("type") == "generation":
                    g = rec["generation"]
                    lam = rec["best_lam_raw"]
                    f_val = rec["best_f"]
                    # Track global best across islands for this generation
                    if g not in gen_evolution or f_val > gen_evolution[g]["best_f"]:
                        gen_evolution[g] = {
                            "generation": g,
                            "best_f": f_val,
                            "best_lam": lam,
                            "best_m": rec["best_m"],
                        }

    # Extract best per m from hall_of_fame_per_m
    best_per_m = {}
    for m in range(2, 11):
        if m in engine.hall_of_fame_per_m:
            ind = engine.hall_of_fame_per_m[m]
            val_reported = ind.lam_star if ind.s_star <= cfg.tol else -ind.s_star
            metric_label = "lam*" if ind.s_star <= cfg.tol else "-s*"
            best_per_m[m] = {
                "reported_value": val_reported,
                "metric_label": metric_label,
                "s_star": ind.s_star,
                "lam_star": ind.lam_star,
                "lam_normalized": ind.lam_normalized,
                "f": ind.f,
                "id": ind.id,
                "operator": ind.operator_name,
            }
        else:
            best_per_m[m] = None

    summary = {
        "mu_t": mu_t,
        "num_islands": 8,
        "island_pop_size": 64,
        "num_generations": 200,
        "total_elapsed_seconds": total_elapsed,
        "new_evaluations": len(engine.cache),
        "cache_hits": engine.cache_hits,
        "best_global_f": res["best_f"],
        "best_global_m": res["best_m"],
        "best_global_lam": res["best_lam"],
        "best_per_m": best_per_m,
        "generation_evolution": [gen_evolution[g] for g in sorted(gen_evolution.keys())],
    }

    res_json = Path("results/prod_run_0.6929077.json")
    with open(res_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved complete summary to {res_json}", flush=True)


if __name__ == "__main__":
    main()
