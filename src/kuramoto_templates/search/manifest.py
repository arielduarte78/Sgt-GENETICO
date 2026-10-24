"""Run manifest generation capturing git provenance, environment versions, and seeds."""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import clarabel  # type: ignore[import-untyped]
import cvxpy
import numpy as np
import scipy  # type: ignore[import-untyped]

from kuramoto_templates.search.config import GAConfig


def _get_git_info(allow_dirty: bool = False) -> tuple[str, bool]:
    """Retrieve git commit hash and verify working tree cleanliness."""
    try:
        status_res = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True,
            text=True,
            check=True,
        )
        is_dirty = bool(status_res.stdout.strip())
        if is_dirty and not allow_dirty:
            raise RuntimeError(
                "Git working tree is dirty. Commit changes or explicitly pass --allow-dirty."
            )

        hash_res = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        git_hash = hash_res.stdout.strip()
        return git_hash, is_dirty
    except (subprocess.SubprocessError, FileNotFoundError):
        return "UNKNOWN_OR_UNTRACKED", True


def create_run_manifest(
    config: GAConfig,
    run_id: str,
    allow_dirty: bool = False,
    no_seeds: bool = False,
    no_split_merge: bool = False,
    spawned_seeds: list[int] | None = None,
) -> dict[str, Any]:
    """Construct full audit manifest for the search run."""
    git_hash, is_dirty = _get_git_info(allow_dirty=allow_dirty)

    if spawned_seeds is None:
        # Spawn seeds deterministically using SeedSequence
        seq = np.random.SeedSequence(config.seed)
        total_spawns = config.num_islands * config.num_generations
        children = seq.spawn(total_spawns)
        spawned_seeds = [int(c.generate_state(1)[0]) for c in children]

    manifest: dict[str, Any] = {
        "run_id": run_id,
        "date_utc": datetime.now(UTC).isoformat(),
        "git_hash": git_hash,
        "is_dirty": is_dirty,
        "config": asdict(config),
        "config_sha256": config.compute_sha256(),
        "mu_t": config.mu_t,
        "num_workers": config.num_workers,
        "ablation_flags": {
            "no_seeds": no_seeds,
            "no_split_merge": no_split_merge,
        },
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": getattr(scipy, "__version__", "unknown"),
            "cvxpy": getattr(cvxpy, "__version__", "unknown"),
            "clarabel": getattr(clarabel, "__version__", "unknown"),
        },
        "spawned_seeds": spawned_seeds,
    }
    return manifest


def save_manifest(manifest: dict[str, Any], filepath: Path | str) -> None:
    """Save manifest to JSON file."""
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
