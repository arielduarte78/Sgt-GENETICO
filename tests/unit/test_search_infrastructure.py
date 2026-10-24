"""Unit tests for runlog, manifest generation, and atomic checkpointing."""

import tempfile
from pathlib import Path

from kuramoto_templates.search.checkpoint import (
    CheckpointState,
    load_checkpoint,
    save_checkpoint_atomic,
)
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.manifest import create_run_manifest
from kuramoto_templates.search.runlog import RunLogger, reconstruct_lineage


def test_runlog_logging_and_lineage_reconstruction() -> None:
    """RunLogger must write JSONL entries and allow full ancestor lineage reconstruction."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        log_path = Path(tmp_dir) / "run.jsonl"
        logger = RunLogger(log_path, run_id="test_run_01")

        # Initial individual (ID 1, no parents)
        logger.log_evaluation(
            eval_id=1,
            parent_ids=(),
            operator="init",
            f=-0.05,
            s_star=0.05,
            lam_raw=0.0,
            lam_norm=0.0,
            solver_status="OPTIMAL",
            eval_time=0.01,
            m=4,
            island=0,
            generation=0,
        )

        # Child individual (ID 2, parent 1)
        logger.log_evaluation(
            eval_id=2,
            parent_ids=(1,),
            operator="gauss_u",
            f=0.02,
            s_star=0.0,
            lam_raw=0.02,
            lam_norm=0.05,
            solver_status="OPTIMAL",
            eval_time=0.015,
            m=4,
            island=0,
            generation=1,
        )

        # Grandchild individual (ID 3, parent 2)
        logger.log_evaluation(
            eval_id=3,
            parent_ids=(2,),
            operator="phase_dev",
            f=0.03,
            s_star=0.0,
            lam_raw=0.03,
            lam_norm=0.08,
            solver_status="OPTIMAL",
            eval_time=0.012,
            m=4,
            island=0,
            generation=2,
        )

        logger.log_generation_summary(
            island=0,
            generation=2,
            best_f=0.03,
            best_m=4,
            best_lam_raw=0.03,
            mean_f=0.01,
            cache_hits=5,
            evals_count=3,
        )

        lineage = reconstruct_lineage(log_path, target_eval_id=3)
        assert len(lineage) == 3
        assert lineage[0]["eval_id"] == 1
        assert lineage[1]["eval_id"] == 2
        assert lineage[2]["eval_id"] == 3


def test_manifest_creation() -> None:
    """Manifest records git commit, config SHA, versions, mu_t, ablation flags, and seeds."""
    cfg = GAConfig(mu_t=0.690, num_islands=2, island_pop_size=5, num_generations=2)
    manifest = create_run_manifest(
        config=cfg,
        run_id="run_test",
        allow_dirty=True,
        no_seeds=False,
        no_split_merge=False,
    )

    assert "git_hash" in manifest
    assert "config_sha256" in manifest
    assert manifest["mu_t"] == 0.690
    assert "versions" in manifest
    assert "numpy" in manifest["versions"]
    assert "clarabel" in manifest["versions"]
    assert "cvxpy" in manifest["versions"]
    assert "spawned_seeds" in manifest
    assert manifest["ablation_flags"]["no_seeds"] is False


def test_atomic_checkpoint_save_and_load() -> None:
    """Atomic checkpoint uses temporary file + replace and loads bit-for-bit identical state."""
    state = CheckpointState(
        run_id="run_cp_test",
        generation=5,
        island_populations=[[{"id": 1, "f": 0.5}], [{"id": 2, "f": 0.6}]],
        hall_of_fame_global=[{"id": 2, "f": 0.6}],
        hall_of_fame_per_m={4: {"id": 2, "f": 0.6}},
        cache={"dummy_key": {"lam_raw": 0.5}},
        cache_hits=12,
        next_individual_id=10,
    )

    with tempfile.TemporaryDirectory() as tmp_dir:
        cp_path = Path(tmp_dir) / "checkpoint_gen_5.json"
        save_checkpoint_atomic(state, cp_path)

        assert cp_path.exists()
        # Ensure temporary file does not linger
        tmp_pattern = list(Path(tmp_dir).glob("*.tmp"))
        assert len(tmp_pattern) == 0

        loaded_state = load_checkpoint(cp_path)
        assert loaded_state.run_id == state.run_id
        assert loaded_state.generation == state.generation
        assert loaded_state.cache_hits == state.cache_hits
        assert loaded_state.island_populations == state.island_populations
