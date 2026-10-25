"""Ejecuta los dos experimentos de 50 generaciones requeridos en el Item 4:
1. mu_t = 0.690 (control) con nueva fitness ordenada.
2. mu_t = 0.6910 con nueva fitness ordenada.

Evalúa si el GA conserva/mejora S7 y reporta el mejor mu alcanzable con lam_hat >= piso (Fase II).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import ScaledSDPEvaluator
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine
from kuramoto_templates.search.seeds import create_s7_direct_seed


def run_experiment(mu_t: float, run_id: str) -> None:
    print("=" * 75)
    print(f"INICIANDO EXPERIMENTO 50 GENS: {run_id} (mu_t = {mu_t:.4f})")
    print("=" * 75)

    out_dir = Path(f"runs/{run_id}")
    cfg = GAConfig(
        mu_t=mu_t,
        num_islands=4,
        island_pop_size=25,  # pop total = 100
        num_generations=50,
        checkpoint_interval=10,
        use_seeds=True,
        seed=2026,
        evaluator="scaled",
        lam_floor_hat=1e-5,
    )

    t0 = time.perf_counter()
    engine = GAEngine(config=cfg, output_dir=out_dir, allow_dirty=True)
    engine.run()
    elapsed = time.perf_counter() - t0

    # Analizar resultados
    status_path = out_dir / "status.json"
    ckpt_path = out_dir / "checkpoints" / "checkpoint_gen_50.json"
    log_path = out_dir / "evaluations.jsonl"

    with open(status_path) as f:
        status_data = json.load(f)

    with open(ckpt_path) as f:
        ckpt_data = json.load(f)

    hof_global = ckpt_data.get("hall_of_fame_global", [])

    # Contar evaluaciones
    total_evals = status_data.get("total_evaluations", 0)
    cache_hits = status_data.get("cache_hits", 0)

    # Identificar ganador
    winner = hof_global[0] if hof_global else None

    # Evaluar Fase II en el ganador y en S7 inicial
    ev_scaled = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, lam_floor_hat=1e-5))

    # Baseline S7 inicial
    s7_seed = create_s7_direct_seed(h=0.002, t=0.0)
    s7_phase1 = ev_scaled.evaluate_phase1(
        s7_seed.phenotype_x, s7_seed.phenotype_phi, genome=s7_seed.genome
    )
    s7_phase2 = ev_scaled.evaluate_phase2(
        s7_seed.phenotype_x, s7_seed.phenotype_phi, eps=1e-5, genome=s7_seed.genome
    )

    winner_phase2_mu = None
    if winner:
        w_x = tuple(winner["phenotype_x"])
        w_phi = tuple(winner["phenotype_phi"])
        # Genome dict to Genome dataclass if needed
        from kuramoto_templates.search.genome import Genome
        w_g_dict = winner["genome"]
        w_genome = Genome(
            m=w_g_dict["m"],
            K_ref=w_g_dict["K_ref"],
            k=tuple(w_g_dict["k"]),
            sgn=tuple(w_g_dict["sgn"]),
            ell=tuple(w_g_dict["ell"]),
            u=tuple(w_g_dict["u"]),
        )
        r_p2 = ev_scaled.evaluate_phase2(w_x, w_phi, eps=1e-5, genome=w_genome)
        winner_phase2_mu = r_p2.mu

    # Contar falsos positivos en el log
    false_positives = 0
    if log_path.exists():
        with open(log_path) as f:
            for line in f:
                if "falso_positivo" in line:
                    false_positives += 1

    # Contar cuántos en HoF tienen lam_hat >= lam_floor_hat
    lam_floor_hat = 1e-5
    hof_ge_floor = sum(
        1 for ind in hof_global
        if (ind.get("lam_normalized") is not None and ind.get("lam_normalized") >= lam_floor_hat)
        or (ind.get("lam_star") is not None and ind.get("lam_star") >= lam_floor_hat)
    )

    report = {
        "run_id": run_id,
        "mu_t": mu_t,
        "elapsed_seconds": elapsed,
        "total_evaluations": total_evals,
        "cache_hits": cache_hits,
        "false_positives_detected": false_positives,
        "hof_global_count": len(hof_global),
        "hof_ge_floor_count": hof_ge_floor,
        "s7_baseline": {
            "f": s7_phase1.f,
            "lam_hat": s7_phase1.lam_normalized,
            "s_star": s7_phase1.s,
            "mu_phase2": s7_phase2.mu,
            "status": s7_phase1.status,
        },
        "winner": {
            "id": winner.get("id") if winner else None,
            "m": winner.get("genome", {}).get("m") if winner else None,
            "f": winner.get("f") if winner else None,
            "s_star": winner.get("s_star") if winner else None,
            "lam_hat": winner.get("lam_normalized") if winner else None,
            "mu_phase2": winner_phase2_mu,
            "operator": winner.get("operator_name") if winner else None,
        } if winner else None,
        "hof_top5": [
            {
                "id": ind.get("id"),
                "m": ind.get("genome", {}).get("m"),
                "f": ind.get("f"),
                "s_star": ind.get("s_star"),
                "lam_hat": ind.get("lam_normalized"),
                "operator": ind.get("operator_name"),
            }
            for ind in hof_global[:5]
        ],
    }

    # Guardar reporte
    out_report = Path(f"results/{run_id}_summary.json")
    with open(out_report, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\n--- Resumen Corrida {run_id} ---")
    print(f"Tiempo: {elapsed:.2f} s | Total Evals: {total_evals} | Cache Hits: {cache_hits}")
    print(f"Falsos positivos: {false_positives}")
    print(f"HoF con lam_hat >= piso ({lam_floor_hat}): {hof_ge_floor} de {len(hof_global)}")
    print(f"S7 Baseline Phase II (mu*): {s7_phase2.mu:.8f}")
    if winner:
        w_m = winner.get("genome", {}).get("m")
        w_f = winner.get("f")
        w_s = winner.get("s_star")
        w_lam = winner.get("lam_normalized")
        print(
            f"Ganador: ID {winner.get('id')}, m={w_m}, "
            f"f={w_f:.8e}, s*={w_s:.6e}, lam_hat={w_lam:.6e}"
        )
        p2_str = f"{winner_phase2_mu:.8f}" if winner_phase2_mu else "None"
        print(f"Ganador Phase II (mu* max con lam_hat >= 1e-5): {p2_str}")
        conserva_s7 = winner.get("f") >= s7_phase1.f
        print(
            f"¿GA conserva/mejora S7?: {conserva_s7} "
            f"(f_ganador={w_f:.8e} vs f_S7={s7_phase1.f:.8e})"
        )
    print("=" * 75 + "\n")


def main() -> None:
    # 1. Corrida control a mu_t = 0.690
    run_experiment(mu_t=0.690, run_id="run_50gen_mu_0690_control")

    # 2. Corrida frontera a mu_t = 0.6910
    run_experiment(mu_t=0.6910, run_id="run_50gen_mu_06910_frontier")

    # Cargar y mostrar comparativa final
    p0 = Path("results/run_50gen_mu_0690_control_summary.json")
    p1 = Path("results/run_50gen_mu_06910_frontier_summary.json")
    if p0.exists() and p1.exists():
        with open(p0) as f0:
            rep_0690 = json.load(f0)
        with open(p1) as f1:
            rep_06910 = json.load(f1)

        print("\n" + "=" * 75)
        print("RESULTADOS FINALES COMPARATIVOS (ITEM 4)")
        print("=" * 75)
        w0 = rep_0690.get("winner") or {}
        w1 = rep_06910.get("winner") or {}
        print(f"mu_t = 0.6900: Ganador m={w0.get('m')}, f={w0.get('f')}, mu*={w0.get('mu_phase2')}")
        print(f"mu_t = 0.6910: Ganador m={w1.get('m')}, f={w1.get('f')}, mu*={w1.get('mu_phase2')}")




if __name__ == "__main__":
    main()
