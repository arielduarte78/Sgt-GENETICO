"""Re-run Smoke Test 3b at mu_t = 0.69291 with scale penalty and diverse scalable seeds."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.evaluators.scaled import compute_h_eff, compute_sheets
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine
from kuramoto_templates.search.genome import Genome


def main():
    print("=" * 70)
    print("RE-EJECUCIÓN SMOKE TEST 3b: Frontera mu_t = 0.69291")
    print("=" * 70, flush=True)

    cfg_path = Path("configs/smoke_frontier_scaled.yaml")
    config = GAConfig.from_yaml(cfg_path)
    output_dir = Path("runs/smoke_frontier_scaled_rerun")
    output_dir.mkdir(parents=True, exist_ok=True)

    engine = GAEngine(config=config, output_dir=output_dir, allow_dirty=True)

    t0 = time.perf_counter()
    results = engine.run()
    elapsed = time.perf_counter() - t0

    # Analyze run.jsonl
    runlog_path = output_dir / "run.jsonl"
    evals = []
    false_positives = []
    with runlog_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("type") == "eval":
                evals.append(rec)
            elif rec.get("type") == "false_positive":
                false_positives.append(rec)

    total_evals = len(evals)
    status_counts: dict[str, int] = {}
    for ev in evals:
        st = str(ev.get("solver_status", "UNKNOWN"))
        status_counts[st] = status_counts.get(st, 0) + 1

    not_scalable_count = status_counts.get("NOT_SCALABLE", 0)
    not_scalable_rate = (not_scalable_count / total_evals) if total_evals > 0 else 0.0

    # Best individual in Hall of Fame
    checkpoints = sorted(output_dir.glob("checkpoints/checkpoint_gen_*.json"))
    final_cp = checkpoints[-1] if checkpoints else None
    hof_global = []
    best_h_eff = None
    best_ind_dict = None

    if final_cp and final_cp.exists():
        with final_cp.open("r", encoding="utf-8") as f:
            cp_data = json.load(f)
            hof_global = cp_data.get("hall_of_fame_global", [])
            if hof_global:
                best_ind_dict = hof_global[0]
                g_d = best_ind_dict["genome"]
                g_obj = Genome(
                    m=int(g_d["m"]),
                    K_ref=int(g_d["K_ref"]),
                    k=tuple(int(v) for v in g_d["k"]),
                    sgn=tuple(int(v) for v in g_d["sgn"]),
                    ell=tuple(float(v) for v in g_d["ell"]),
                    u=tuple(float(v) for v in g_d["u"]),
                )
                sheets, _ = compute_sheets(
                    best_ind_dict["phenotype_x"],
                    best_ind_dict["phenotype_phi"],
                    genome=g_obj,
                )
                best_h_eff = compute_h_eff(
                    best_ind_dict["phenotype_x"],
                    best_ind_dict["phenotype_phi"],
                    sheets,
                )

    print(f"\n--- Resultados Re-ejecución Smoke Test 3b ---")
    print(f"Tiempo total: {elapsed:.2f} s")
    print(f"Total evaluaciones logged: {total_evals}")
    print(f"Evaluaciones nuevas (cache size): {len(engine.cache)}")
    print(f"Cache hits: {engine.cache_hits}")
    print(f"Velocidad: {total_evals / elapsed:.2f} evals/s (log total), {len(engine.cache) / elapsed:.2f} nuevas/s")
    print(f"Desglose de estados solver: {status_counts}")
    print(f"Tasa NOT_SCALABLE: {not_scalable_count}/{total_evals} ({not_scalable_rate*100:.2f}%)")
    print(f"Falsos positivos detectados: {len(false_positives)}")
    print(f"Mejor f global: {results['best_f']:.12e}")
    print(f"Mejor m: {results['best_m']}")
    print(f"Mejor h_eff: {best_h_eff:.8f}" if best_h_eff is not None else "Mejor h_eff: N/A")
    print(f"Mejor slack s*: {results['best_s']:.12e}")
    print(f"Mejor lam_hat: {results['best_lam']:.12e}")
    print(f"Mejor mu_star (Fase II): {results.get('best_mu_star')}")
    if best_ind_dict:
        print(f"Detalle top HoF: ID={best_ind_dict.get('id')}, op={best_ind_dict.get('operator_name')}")

    out_summary = {
        "elapsed_seconds": elapsed,
        "total_evals": total_evals,
        "unique_evals": len(engine.cache),
        "cache_hits": engine.cache_hits,
        "evals_per_second": total_evals / elapsed,
        "status_breakdown": status_counts,
        "not_scalable_count": not_scalable_count,
        "not_scalable_rate": not_scalable_rate,
        "false_positives_detected": len(false_positives),
        "best_f": results["best_f"],
        "best_m": results["best_m"],
        "best_h_eff": best_h_eff,
        "best_s_star": results["best_s"],
        "best_lam_hat": results["best_lam"],
        "best_mu_star": results.get("best_mu_star"),
        "best_individual": best_ind_dict,
    }

    out_file = Path("results/smoke_test_3b_recomputed.json")
    with out_file.open("w", encoding="utf-8") as f:
        json.dump(out_summary, f, indent=2)
    print(f"\nResumen de 3b guardado en {out_file}")


if __name__ == "__main__":
    main()
