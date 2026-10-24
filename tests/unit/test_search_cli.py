"""Unit tests for search CLI commands: run, resume, status, and tail."""

import json
import tempfile
from pathlib import Path

import pytest

from kuramoto_templates.search.cli import main
from kuramoto_templates.search.config import GAConfig


def test_cli_run_status_and_tail(capsys: pytest.CaptureFixture[str]) -> None:
    """CLI ga run executes search, produces status.json, and allows status/tail inspection."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir) / "run_out"
        cfg_path = Path(tmp_dir) / "test_config.yaml"

        cfg = GAConfig(
            mu_t=0.690,
            num_islands=1,
            island_pop_size=4,
            num_generations=2,
            use_seeds=True,
            checkpoint_interval=1,
        )
        cfg.to_yaml(cfg_path)

        # 1. Test 'ga run'
        exit_code = main([
            "run",
            str(cfg_path),
            "--output-dir",
            str(out_dir),
            "--allow-dirty",
        ])
        assert exit_code == 0

        # Check status.json was created
        status_file = out_dir / "status.json"
        assert status_file.exists()
        with status_file.open("r", encoding="utf-8") as f:
            status_data = json.load(f)
        assert status_data["status"] == "completed"
        assert status_data["generation"] == 2

        # 2. Test 'ga status'
        capsys.readouterr()
        status_exit = main(["status", "--output-dir", str(out_dir)])
        assert status_exit == 0
        captured = capsys.readouterr()
        assert "completed" in captured.out

        # 3. Test 'ga tail'
        tail_exit = main(["tail", "--output-dir", str(out_dir), "-n", "3"])
        assert tail_exit == 0
        captured_tail = capsys.readouterr()
        assert "generation" in captured_tail.out or "eval" in captured_tail.out
