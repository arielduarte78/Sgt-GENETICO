"""Atomic checkpointing for GA search state with bit-for-bit resume capability."""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class CheckpointState:
    """Serializable snapshot of the full genetic algorithm state."""

    run_id: str
    generation: int
    island_populations: list[list[dict[str, Any]]]
    hall_of_fame_global: list[dict[str, Any]]
    hall_of_fame_per_m: dict[int, dict[str, Any]]
    cache: dict[str, dict[str, Any]]
    cache_hits: int = 0
    next_individual_id: int = 0
    rng_states: list[dict[str, Any]] = field(default_factory=list)


def save_checkpoint_atomic(state: CheckpointState, filepath: Path | str) -> None:
    """Save checkpoint atomically via temporary file and atomic replace."""
    target_path = Path(filepath)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = target_path.with_suffix(".tmp")

    data = asdict(state)
    # Serialize to temp file (compact JSON to save disk space and I/O time)
    with temp_path.open("w", encoding="utf-8") as f:
        json.dump(data, f, separators=(",", ":"))

    # Atomic rename / replace
    os.replace(temp_path, target_path)


def load_checkpoint(filepath: Path | str) -> CheckpointState:
    """Load checkpoint state from JSON file."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Checkpoint file does not exist: {filepath}")

    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    # Convert keys of hall_of_fame_per_m back to int
    if "hall_of_fame_per_m" in data:
        data["hall_of_fame_per_m"] = {
            int(k): v for k, v in data["hall_of_fame_per_m"].items()
        }

    return CheckpointState(**data)
