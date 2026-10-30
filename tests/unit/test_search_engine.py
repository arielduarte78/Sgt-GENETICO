"""Unit and deterministic reproduction tests for GA search engine."""

import tempfile
from pathlib import Path

import pytest

from kuramoto_templates.evaluators.results import EvalResult
from kuramoto_templates.search.checkpoint import load_checkpoint
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine


def mock_analytical_evaluator(
    x: tuple[float, ...], phi: tuple[float, ...]
) -> EvalResult:
    """Mock evaluator with known optimum for acceptance test (a).

    Maximizes target function when phases are close to pi/2, with s* = 0.
    """
    m = len(x)
    target_phi = 0.5
    dev = sum((p - target_phi) ** 2 for p in phi) / float(m)
    lam = float(max(1.0 - dev, 0.0))
    s = 0.0  # strictly feasible
    return EvalResult(
        f=lam,
        mu=0.7,
        Phi=0.0,
        lam=lam,
        s=s,
        lam_normalized=lam,
        status="OPTIMAL",
    )


def test_acceptance_a_mock_evaluator_finds_optimum() -> None:
    """ACCEPTANCE TEST (a):

    Evaluador falso (función analítica barata): el GA halla el óptimo conocido con seed fija.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        cfg = GAConfig(
            mu_t=0.690,
            num_islands=2,
            island_pop_size=10,
            num_generations=15,
            use_seeds=False,
            seed=42,
            checkpoint_interval=5,
        )
        engine = GAEngine(
            config=cfg,
            output_dir=tmp_dir,
            evaluator_fn=mock_analytical_evaluator,
            allow_dirty=True,
        )
        results = engine.run()

        # Engine should find fitness very close to the analytical optimum (1.0)
        best_f = results["best_f"]
        assert best_f >= 0.95, f"Mock evaluator optimum not found, best_f = {best_f}"


def test_engine_bit_for_bit_checkpoint_resume_identity() -> None:
    """MANDATORY TEST:

    20 generaciones seguidas == 10 + resume + 10, comparando archivos bit a bit.
    """
    with tempfile.TemporaryDirectory() as base_tmp:
        dir_full = Path(base_tmp) / "full_20"
        dir_split = Path(base_tmp) / "split_10_10"

        cfg_full = GAConfig(
            mu_t=0.690,
            num_islands=2,
            island_pop_size=8,
            num_generations=20,
            use_seeds=False,
            seed=12345,
            checkpoint_interval=10,
        )

        # 1. Run 20 generations continuously
        engine_full = GAEngine(
            config=cfg_full,
            output_dir=dir_full,
            evaluator_fn=mock_analytical_evaluator,
            run_id="run_reproducibility",
            allow_dirty=True,
        )
        engine_full.run()

        # 2. Run 10 generations, then resume for 10 more
        cfg_part1 = GAConfig(
            mu_t=0.690,
            num_islands=2,
            island_pop_size=8,
            num_generations=10,
            use_seeds=False,
            seed=12345,
            checkpoint_interval=10,
        )
        engine_part1 = GAEngine(
            config=cfg_part1,
            output_dir=dir_split,
            evaluator_fn=mock_analytical_evaluator,
            run_id="run_reproducibility",
            allow_dirty=True,
        )
        engine_part1.run()

        # Resume from checkpoint at generation 10 to generation 20
        cp_path = dir_split / "checkpoints" / "checkpoint_gen_10.json"
        engine_resume = GAEngine.resume(
            checkpoint_path=cp_path,
            total_generations=20,
            evaluator_fn=mock_analytical_evaluator,
        )
        engine_resume.run()

        # 3. Compare final checkpoint state bit a bit
        cp_full = dir_full / "checkpoints" / "checkpoint_gen_20.json"
        cp_resumed = dir_split / "checkpoints" / "checkpoint_gen_20.json"

        assert cp_full.exists()
        assert cp_resumed.exists()

        state_full = load_checkpoint(cp_full)
        state_resumed = load_checkpoint(cp_resumed)

        # Generations, cache hits, next_id
        assert state_full.generation == state_resumed.generation == 20
        assert state_full.cache_hits == state_resumed.cache_hits
        assert state_full.next_individual_id == state_resumed.next_individual_id

        # Verify island populations have exact identical genetic contents
        for pop_f, pop_r in zip(
            state_full.island_populations, state_resumed.island_populations, strict=False
        ):
            assert len(pop_f) == len(pop_r)
            for ind_f, ind_r in zip(pop_f, pop_r, strict=False):
                assert ind_f["id"] == ind_r["id"]
                assert ind_f["genome"] == ind_r["genome"]
                assert ind_f["phenotype_x"] == ind_r["phenotype_x"]
                assert ind_f["phenotype_phi"] == ind_r["phenotype_phi"]
                assert ind_f["f"] == pytest.approx(ind_r["f"])
                assert ind_f["s_star"] == pytest.approx(ind_r["s_star"])
                assert ind_f["lam_star"] == pytest.approx(ind_r["lam_star"])
                assert ind_f["parent_ids"] == ind_r["parent_ids"]
                assert ind_f["operator_name"] == ind_r["operator_name"]

        # Verify Hall of Fame bit-for-bit
        assert len(state_full.hall_of_fame_global) == len(state_resumed.hall_of_fame_global)
        for hof_f, hof_r in zip(
            state_full.hall_of_fame_global, state_resumed.hall_of_fame_global, strict=False
        ):
            assert hof_f["id"] == hof_r["id"]
            assert hof_f["genome"] == hof_r["genome"]
            assert hof_f["f"] == pytest.approx(hof_r["f"])


def test_hall_of_fame_reverification_rejects_false_positive() -> None:
    """Valida que un individuo con f>0 pero falso positivo en SCS/NumPy sea rechazado del HoF."""
    import json

    from kuramoto_templates.search.genome import Genome
    from kuramoto_templates.search.selection import Individual
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        cfg = GAConfig(
            mu_t=0.6929077,
            num_islands=1,
            island_pop_size=4,
            num_generations=1,
            lam_floor=1e-7,
        )
        engine = GAEngine(config=cfg, output_dir=tmp_dir, allow_dirty=True)

        # Usar coordenadas de ID 39260 (que sabemos da falso positivo en SCS y NumPy)
        x = (0.2494426094558034, 0.25018825224215374, 0.25053408179500397, 0.2498350565070389)
        phi = (-3.141591616161418, -1.5707953267948964, 1.000000000139778e-06, 1.5707973267948967)

        # Crear un individuo artificial con f > 0 pretendiendo entrar al HoF
        fake_ind = Individual(
            id=999,
            genome=Genome(
                m=4,
                K_ref=4,
                k=(0, 0, 0, 0),
                sgn=(1, 1, 1, 1),
                ell=(-2.0, -2.0, -2.0, -2.0),
                u=(0.0, 0.0, 0.0, 0.0),
            ),
            phenotype_x=x,
            phenotype_phi=phi,
            s_star=0.0,
            lam_star=2e-7,
            lam_normalized=2e-7,
            f=2e-7,  # f > 0
            f_sel=2e-7,
            parent_ids=(),
            operator_name="test_operator",
            generation=1,
        )

        engine._update_hall_of_fame(fake_ind)

        # Debe ser rechazado del HoF global y por m
        assert len(engine.hall_of_fame_global) == 0 or all(
            ind.id != 999 for ind in engine.hall_of_fame_global
        )
        assert 4 not in engine.hall_of_fame_per_m or engine.hall_of_fame_per_m[4].id != 999

        # Debe haberse registrado en run.jsonl como falso_positivo
        log_path = Path(tmp_dir) / "run.jsonl"
        assert log_path.exists()
        found_fp = False
        with log_path.open("r", encoding="utf-8") as f:
            for line in f:
                rec = json.loads(line)
                if rec.get("type") == "falso_positivo" and rec.get("eval_id") == 999:
                    found_fp = True
                    assert "lam_scs" in rec
                    assert "lam_min_M" in rec
                    assert rec["lam_min_M"] < cfg.lam_floor
        assert found_fp, "No se encontró registro de falso_positivo en run.jsonl"


def test_checkpoint_pruning_bounds_disk_usage(tmp_path: Path) -> None:
    """Verify that checkpoint pruning retains only the most recent N checkpoints."""
    cfg = GAConfig(
        mu_t=0.690,
        num_islands=2,
        island_pop_size=4,
        num_generations=6,
        use_seeds=False,
        seed=42,
        checkpoint_interval=2,
        checkpoint_keep_count=2,
    )
    engine = GAEngine(
        config=cfg,
        output_dir=tmp_path,
        evaluator_fn=mock_analytical_evaluator,
        allow_dirty=True,
    )
    engine.run()

    ckpt_dir = tmp_path / "checkpoints"
    ckpts = sorted(ckpt_dir.glob("checkpoint_gen_*.json"))
    # Only 2 checkpoints should remain: gen 4 and gen 6
    assert len(ckpts) == 2
    assert ckpts[0].name == "checkpoint_gen_4.json"
    assert ckpts[1].name == "checkpoint_gen_6.json"
    assert not (ckpt_dir / "checkpoint_gen_2.json").exists()


