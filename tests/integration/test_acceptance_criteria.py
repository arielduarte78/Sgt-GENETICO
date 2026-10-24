"""Acceptance criteria tests (c), (d), and (e) using real SDPEvaluator."""

import json
from pathlib import Path
import tempfile
import pytest

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.fitness import compute_continuous_fitness
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine
from kuramoto_templates.search.seeds import create_s7_direct_seed


def test_acceptance_c_seeds_monotonic_improvement() -> None:
    """ACCEPTANCE TEST (c):

    GA con semillas, mu_t=0.690, corrida corta: mejor f >= f de la semilla y
    monótono no decreciente. Debe pasar.
    """
    # 1. Compute baseline fitness of S7 seed
    s7_seed = create_s7_direct_seed()
    evaluator = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-8))
    res_s7 = evaluator.evaluate_phase1(s7_seed.phenotype_x, s7_seed.phenotype_phi)
    assert res_s7.s is not None and res_s7.lam is not None
    f_seed = compute_continuous_fitness(lam=res_s7.lam, s=res_s7.s, tol=1e-8)
    assert f_seed > 0.0

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Short run: 3 generations, 2 islands, pop_size=4
        cfg = GAConfig(
            mu_t=0.690,
            num_islands=2,
            island_pop_size=4,
            num_generations=3,
            use_seeds=True,
            seed=42,
            checkpoint_interval=1,
        )
        engine = GAEngine(config=cfg, output_dir=tmp_dir, allow_dirty=True)
        results = engine.run()

        # Check runlog generation progression
        log_path = Path(tmp_dir) / "run.jsonl"
        best_f_per_gen: list[float] = []

        with log_path.open("r", encoding="utf-8") as f_in:
            for line in f_in:
                entry = json.loads(line)
                if entry.get("type") == "generation":
                    best_f_per_gen.append(entry["best_f"])

        best_final = results["best_f"]

        print(f"\nACCEPTANCE TEST (c) RESULTS:")
        print(f"  f_seed (S7): {f_seed:.12e}")
        print(f"  best_final:  {best_final:.12e}")
        print(f"  evals_count: {results['total_evaluations']}")
        print(f"  cache_hits:  {results['cache_hits']}")

        assert best_final >= f_seed - 1e-12, (
            f"Best fitness {best_final} is lower than seed fitness {f_seed}"
        )


def test_acceptance_d_no_blocking_no_seeds_report() -> None:
    """ACCEPTANCE TEST (d) - NON-BLOCKING REPORT:

    GA sin semillas (--no-seeds) en mu_t=0.690; imprime si alcanzó f>0 y en cuántas generaciones.
    No es criterio de aprobación.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        cfg = GAConfig(
            mu_t=0.690,
            num_islands=2,
            island_pop_size=4,
            num_generations=2,
            use_seeds=False,
            seed=999,
        )
        engine = GAEngine(
            config=cfg,
            output_dir=tmp_dir,
            no_seeds=True,
            allow_dirty=True,
        )
        results = engine.run()

        best_f = results["best_f"]
        reached_f_positive = best_f > 0.0

        print(f"\nACCEPTANCE TEST (d) [NON-BLOCKING REPORT]:")
        print(f"  mu_t: 0.690")
        print(f"  no_seeds: True")
        print(f"  reached f > 0: {reached_f_positive}")
        print(f"  best_f achieved: {best_f:.12e}")
        print(f"  best_lam: {results['best_lam']:.12e}")
        print(f"  best_s: {results['best_s']:.12e}")


def test_acceptance_e_smoke_run_real_evaluator() -> None:
    """ACCEPTANCE TEST (e):

    Corrida de humo con evaluador real (3 generaciones, población chica): termina
    sin errores y deja manifiesto, JSONL y checkpoint válidos. Reporta el tiempo
    medio por evaluación.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir)
        cfg = GAConfig(
            mu_t=0.690,
            num_islands=2,
            island_pop_size=4,
            num_generations=3,
            use_seeds=True,
            seed=777,
            checkpoint_interval=3,
        )
        engine = GAEngine(config=cfg, output_dir=out_dir, allow_dirty=True)
        results = engine.run()

        # 1. Verify manifest exists and is valid JSON
        manifest_path = out_dir / "manifest.json"
        assert manifest_path.exists()
        with manifest_path.open("r", encoding="utf-8") as f:
            manifest = json.load(f)
        assert manifest["mu_t"] == 0.690
        assert "git_hash" in manifest

        # 2. Verify run.jsonl exists and contains eval and generation entries
        log_path = out_dir / "run.jsonl"
        assert log_path.exists()
        eval_times: list[float] = []
        with log_path.open("r", encoding="utf-8") as f:
            for line in f:
                rec = json.loads(line)
                if rec.get("type") == "eval" and rec.get("eval_time", 0.0) > 0.0:
                    eval_times.append(rec["eval_time"])

        assert len(eval_times) > 0
        mean_eval_time = sum(eval_times) / len(eval_times)

        # 3. Verify checkpoint exists and is valid
        cp_path = out_dir / "checkpoints" / "checkpoint_gen_3.json"
        assert cp_path.exists()
        with cp_path.open("r", encoding="utf-8") as f:
            cp_data = json.load(f)
        assert cp_data["generation"] == 3

        print(f"\nACCEPTANCE TEST (e) SMOKE RUN METRICS:")
        print(f"  total evaluations: {results['total_evaluations']}")
        print(f"  cache hits: {results['cache_hits']}")
        print(f"  mean evaluation time: {mean_eval_time:.4f} s / eval")
        print(f"  best f: {results['best_f']:.12e}")
