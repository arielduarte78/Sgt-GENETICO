"""Run and analyze 50-generation smoke tests:
3a: Control run with S7 at mu_t = 0.690
3b: Frontier exploration at mu_t = 0.69291
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

# Ensure src is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine


def analyze_run(output_dir: Path) -> dict:
    runlog_path = output_dir / "run.jsonl"
    status_path = output_dir / "status.json"

    evals = []
    false_positives = []
    gen_summaries = []

    if runlog_path.exists():
        with runlog_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                record = json.loads(line)
                rec_type = record.get("type")
                if rec_type == "eval":
                    evals.append(record)
                elif rec_type == "false_positive":
                    false_positives.append(record)
                elif rec_type == "generation":
                    gen_summaries.append(record)

    with status_path.open("r", encoding="utf-8") as f:
        status_data = json.load(f)

    # Checkpoint analysis for Hall of Fame
    checkpoints = sorted(output_dir.glob("checkpoints/checkpoint_gen_*.json"))
    final_cp = checkpoints[-1] if checkpoints else None
    hof_global = []
    hof_per_m = {}
    if final_cp and final_cp.exists():
        with final_cp.open("r", encoding="utf-8") as f:
            cp_data = json.load(f)
            hof_global = cp_data.get("hall_of_fame_global", [])
            hof_per_m = cp_data.get("hall_of_fame_per_m", {})

    # Solver status breakdown
    status_counts = {}
    for ev in evals:
        st = ev.get("solver_status", "UNKNOWN")
        status_counts[st] = status_counts.get(st, 0) + 1

    total_evals = len(evals)
    not_scalable_count = status_counts.get("NOT_SCALABLE", 0)
    not_scalable_rate = (not_scalable_count / total_evals) if total_evals > 0 else 0.0

    return {
        "status_data": status_data,
        "total_logged_evals": total_evals,
        "status_counts": status_counts,
        "not_scalable_count": not_scalable_count,
        "not_scalable_rate": not_scalable_rate,
        "false_positive_count": len(false_positives),
        "false_positives": false_positives,
        "hall_of_fame_global": hof_global,
        "hall_of_fame_per_m": hof_per_m,
    }


def run_test_3a():
    print("=" * 70)
    print("RUNNING SMOKE TEST 3a: Control con semilla S7 en mu_t = 0.690")
    print("=" * 70, flush=True)

    cfg_path = Path("configs/smoke_control_scaled.yaml")
    config = GAConfig.from_yaml(cfg_path)
    output_dir = Path("runs/smoke_control_scaled")
    output_dir.mkdir(parents=True, exist_ok=True)

    engine = GAEngine(config=config, output_dir=output_dir, allow_dirty=True)

    t0 = time.perf_counter()
    results = engine.run()
    elapsed = time.perf_counter() - t0

    analysis = analyze_run(output_dir)

    print(f"\n--- Resultados Smoke Test 3a ---")
    print(f"Tiempo total: {elapsed:.2f} s")
    print(f"Total evaluaciones nuevas (cache size): {len(engine.cache)}")
    print(f"Total cache hits: {engine.cache_hits}")
    print(f"Velocidad: {analysis['total_logged_evals'] / elapsed:.2f} evals/s (log total), {len(engine.cache) / elapsed:.2f} nuevas/s")
    print(f"Falsos positivos rechazados: {analysis['false_positive_count']}")
    print(f"Mejor f global: {results['best_f']:.12e} (m={results['best_m']})")
    print(f"Mejor lam_hat: {results['best_lam']:.12e}")
    print(f"Mejor mu_star (Fase II): {results.get('best_mu_star')}")

    # Check S7 in HoF
    hof_global = analysis["hall_of_fame_global"]
    if hof_global:
        best_hof = hof_global[0]
        print(f"Individuo top en Hall of Fame: ID={best_hof.get('id')}, op={best_hof.get('operator_name')}, f={best_hof.get('f')}, lam_hat={best_hof.get('lam_star')}, mu_star={best_hof.get('mu_star')}")

    return {
        "elapsed_sec": elapsed,
        "results": results,
        "analysis": analysis,
    }


def run_test_3b():
    print("\n" + "=" * 70)
    print("RUNNING SMOKE TEST 3b: Exploracion en la frontera mu_t = 0.69291")
    print("=" * 70, flush=True)

    cfg_path = Path("configs/smoke_frontier_scaled.yaml")
    config = GAConfig.from_yaml(cfg_path)
    output_dir = Path("runs/smoke_frontier_scaled")
    output_dir.mkdir(parents=True, exist_ok=True)

    engine = GAEngine(config=config, output_dir=output_dir, allow_dirty=True)

    t0 = time.perf_counter()
    results = engine.run()
    elapsed = time.perf_counter() - t0

    analysis = analyze_run(output_dir)

    print(f"\n--- Resultados Smoke Test 3b ---")
    print(f"Tiempo total: {elapsed:.2f} s")
    print(f"Total evaluaciones nuevas (cache size): {len(engine.cache)}")
    print(f"Total cache hits: {engine.cache_hits}")
    print(f"Velocidad: {analysis['total_logged_evals'] / elapsed:.2f} evals/s (log total), {len(engine.cache) / elapsed:.2f} nuevas/s")
    print(f"Desglose de estados solver: {analysis['status_counts']}")
    print(f"Tasa NOT_SCALABLE: {analysis['not_scalable_count']}/{analysis['total_logged_evals']} ({analysis['not_scalable_rate']*100:.2f}%)")
    print(f"Falsos positivos detectados: {analysis['false_positive_count']}")
    print(f"Mejor f global: {results['best_f']:.12e} (m={results['best_m']}, slack s*={results['best_s']:.12e}, lam_hat={results['best_lam']:.12e})")
    print(f"Mejor mu_star (Fase II): {results.get('best_mu_star')}")

    hof_global = analysis["hall_of_fame_global"]
    if hof_global:
        best_hof = hof_global[0]
        print(f"Individuo top en Hall of Fame: ID={best_hof.get('id')}, op={best_hof.get('operator_name')}, f={best_hof.get('f')}, s*={best_hof.get('s_star')}, lam_hat={best_hof.get('lam_star')}, mu_star={best_hof.get('mu_star')}")

    return {
        "elapsed_sec": elapsed,
        "results": results,
        "analysis": analysis,
    }


def main():
    res_3a = run_test_3a()
    res_3b = run_test_3b()

    summary = {
        "run_3a_control": {
            "elapsed_seconds": res_3a["elapsed_sec"],
            "best_f": res_3a["results"]["best_f"],
            "best_m": res_3a["results"]["best_m"],
            "best_lam_hat": res_3a["results"]["best_lam"],
            "best_s_star": res_3a["results"]["best_s"],
            "best_mu_star": res_3a["results"].get("best_mu_star"),
            "unique_evals": res_3a["analysis"]["status_data"]["total_evaluations"],
            "total_evals": res_3a["analysis"]["total_logged_evals"],
            "cache_hits": res_3a["analysis"]["status_data"]["cache_hits"],
            "evals_per_second": res_3a["analysis"]["total_logged_evals"] / res_3a["elapsed_sec"],
            "false_positives_rejected": res_3a["analysis"]["false_positive_count"],
            "s7_survived_50_gens": bool(res_3a["results"]["best_f"] > 0.0),
        },
        "run_3b_frontier": {
            "elapsed_seconds": res_3b["elapsed_sec"],
            "best_f": res_3b["results"]["best_f"],
            "best_m": res_3b["results"]["best_m"],
            "best_lam_hat": res_3b["results"]["best_lam"],
            "best_s_star": res_3b["results"]["best_s"],
            "best_mu_star": res_3b["results"].get("best_mu_star"),
            "total_evals": res_3b["analysis"]["total_logged_evals"],
            "cache_hits": res_3b["analysis"]["status_data"]["cache_hits"],
            "unique_evals": res_3b["analysis"]["status_data"]["total_evaluations"],
            "evals_per_second": res_3b["analysis"]["total_logged_evals"] / res_3b["elapsed_sec"],
            "status_breakdown": res_3b["analysis"]["status_counts"],
            "not_scalable_count": res_3b["analysis"]["not_scalable_count"],
            "not_scalable_rate": res_3b["analysis"]["not_scalable_rate"],
            "false_positives_detected": res_3b["analysis"]["false_positive_count"],
        },
    }

    out_file = Path("results/smoke_tests_50gens_summary.json")
    with out_file.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResumen completo guardado en {out_file}")


if __name__ == "__main__":
    main()
