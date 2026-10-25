"""Execute investigation runs:
1. 50 generations run to count false positives caught by SCS + NumPy HoF gate.
2. 200 generations control run at mu_t=0.690 without seeds to check if any verified lam* >= 1e-7 is found.
3. Precise S7 evaluation at h=0.002 and h=0.02 under mu_t=0.690.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.fitness import compute_continuous_fitness
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine
from kuramoto_templates.search.seeds import create_s7_direct_seed


def main():
    print("=" * 70)
    print("1. EVALUACION DE S7 EN mu_t=0.690 CON h=0.002 Y h=0.02")
    print("=" * 70)
    evaluator = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-8, lam_floor=1e-7))

    for h_val in [0.002, 0.02]:
        s7_seed = create_s7_direct_seed(h=h_val, t=0.0)
        res = evaluator.evaluate_phase1(s7_seed.phenotype_x, s7_seed.phenotype_phi)
        f_val = compute_continuous_fitness(lam=res.lam, s=res.s, tol=1e-8, lam_floor=1e-7)
        print(f"S7 (h={h_val:5.3f}):")
        print(f"  lam*:      {res.lam:16.8e}")
        print(f"  s*:        {res.s:16.8e}")
        print(f"  status:    {res.status}")
        print(f"  fitness f: {f_val:16.8e} (f > 0? {f_val > 0})")

    print("\n" + "=" * 70)
    print("2. CORRIDA DE 50 GENERACIONES PARA CONTAR FALSOS POSITIVOS")
    print("=" * 70)
    out_dir_50 = Path("runs/test_50gen_false_positives")
    cfg_50 = GAConfig(
        mu_t=0.6929077,
        num_islands=4,
        island_pop_size=16,
        num_generations=50,
        checkpoint_interval=10,
        use_seeds=True,
        seed=42,
        lam_floor=1e-7,
    )
    t0 = time.perf_counter()
    engine_50 = GAEngine(config=cfg_50, output_dir=out_dir_50, allow_dirty=True)
    engine_50.run()
    t_50 = time.perf_counter() - t0

    # Count false positives in log
    fp_records = []
    log_50 = out_dir_50 / "run.jsonl"
    with open(log_50, "r", encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("type") == "falso_positivo":
                fp_records.append(rec)

    print(f"Corrida 50 gens completada en {t_50:.2f}s.")
    print(f"Total evaluaciones: {len(engine_50.cache)}, Hits: {engine_50.cache_hits}")
    print(f"Falsos positivos detectados y rechazados del HoF: {len(fp_records)}")
    for fp in fp_records[:5]:
        print(f"  ID {fp['eval_id']} (m={fp['m']}): lam_scs={fp['lam_scs']:12.4e}, lam_min_M={fp['lam_min_M']:12.4e}")

    print("\n" + "=" * 70)
    print("3. RE-CORRIDA CONTROL mu_t=0.690 SIN SEMILLAS (200 GENERACIONES)")
    print("=" * 70)
    out_dir_ctrl = Path("runs/ctrl_mu_0.690_no_seeds_200gen")
    cfg_ctrl = GAConfig(
        mu_t=0.690,
        num_islands=4,
        island_pop_size=16,
        num_generations=200,
        checkpoint_interval=20,
        use_seeds=False,
        seed=42,
        lam_floor=1e-7,
    )
    t0_ctrl = time.perf_counter()
    engine_ctrl = GAEngine(config=cfg_ctrl, output_dir=out_dir_ctrl, no_seeds=True, allow_dirty=True)
    res_ctrl = engine_ctrl.run()
    t_ctrl = time.perf_counter() - t0_ctrl

    print(f"Control 200 gens completado en {t_ctrl:.2f}s ({t_ctrl/60:.2f} min).")
    print(f"Evaluaciones unicas: {len(engine_ctrl.cache)}, Hits: {engine_ctrl.cache_hits}")
    print(f"Mejor fitness global: {res_ctrl['best_f']:.12e} (m={res_ctrl['best_m']})")
    print(f"Mejor lam* global:    {res_ctrl['best_lam']:.12e}")
    print(f"Mejor s* global:      {res_ctrl['best_s']:.12e}")

    # Inspect Hall of Fame
    print("\nMejores por m en el control sin semillas:")
    verified_positive = []
    for m in range(2, 11):
        if m in engine_ctrl.hall_of_fame_per_m:
            ind = engine_ctrl.hall_of_fame_per_m[m]
            print(f"  m={m:2d}: f={ind.f:14.6e}, s*={ind.s_star:12.4e}, lam*={ind.lam_star:14.6e}, op={ind.operator_name}")
            if ind.f > 0:
                verified_positive.append((m, ind))
        else:
            print(f"  m={m:2d}: Ninguno")

    print(f"\nIndividuos con lam* >= 1e-7 verificados (f > 0): {len(verified_positive)}")
    for m, ind in verified_positive:
        print(f"  m={m}: ID={ind.id}, lam*={ind.lam_star:.12e}, f={ind.f:.12e}")


if __name__ == "__main__":
    main()
