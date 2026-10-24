"""JSONL logging for per-evaluation events, generation summaries, and lineage reconstruction."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class RunLogger:
    """Thread-safe / process-safe appending logger writing structured JSONL."""

    def __init__(self, filepath: Path | str, run_id: str) -> None:
        self.filepath = Path(filepath)
        self.run_id = run_id
        self.filepath.parent.mkdir(parents=True, exist_ok=True)

    def log_evaluation(
        self,
        eval_id: int,
        parent_ids: tuple[int, ...],
        operator: str,
        f: float,
        s_star: float,
        lam_raw: float,
        lam_norm: float,
        solver_status: str,
        eval_time: float,
        m: int,
        island: int,
        generation: int,
    ) -> None:
        """Log a single evaluation record to JSONL."""
        entry: dict[str, Any] = {
            "type": "eval",
            "run_id": self.run_id,
            "eval_id": eval_id,
            "parent_ids": list(parent_ids),
            "operator": operator,
            "f": f,
            "s_star": s_star,
            "lam_raw": lam_raw,
            "lam_norm": lam_norm,
            "solver_status": solver_status,
            "eval_time": eval_time,
            "m": m,
            "island": island,
            "generation": generation,
        }
        with self.filepath.open("a", encoding="utf-8") as f_out:
            f_out.write(json.dumps(entry) + "\n")

    def log_generation_summary(
        self,
        island: int,
        generation: int,
        best_f: float,
        best_m: int,
        best_lam_raw: float,
        mean_f: float,
        cache_hits: int,
        evals_count: int,
    ) -> None:
        """Log an end-of-generation summary for an island."""
        entry: dict[str, Any] = {
            "type": "generation",
            "run_id": self.run_id,
            "island": island,
            "generation": generation,
            "best_f": best_f,
            "best_m": best_m,
            "best_lam_raw": best_lam_raw,
            "mean_f": mean_f,
            "cache_hits": cache_hits,
            "evals_count": evals_count,
        }
        with self.filepath.open("a", encoding="utf-8") as f_out:
            f_out.write(json.dumps(entry) + "\n")


def reconstruct_lineage(
    log_path: Path | str, target_eval_id: int
) -> list[dict[str, Any]]:
    """Reconstruct complete ancestor lineage leading to a specific individual."""
    path = Path(log_path)
    if not path.exists():
        raise FileNotFoundError(f"Runlog not found at {log_path}")

    eval_by_id: dict[int, dict[str, Any]] = {}

    with path.open("r", encoding="utf-8") as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            if entry.get("type") == "eval":
                eval_by_id[int(entry["eval_id"])] = entry

    if target_eval_id not in eval_by_id:
        raise KeyError(f"Individual ID {target_eval_id} not found in log.")

    # Traverse backwards using BFS / recursion
    visited_ids: set[int] = set()
    queue = [target_eval_id]

    while queue:
        curr_id = queue.pop(0)
        if curr_id in visited_ids:
            continue
        visited_ids.add(curr_id)
        parent_ids = eval_by_id[curr_id].get("parent_ids", [])
        for pid in parent_ids:
            if pid in eval_by_id and pid not in visited_ids:
                queue.append(pid)

    # Return sorted chronologically by eval_id
    lineage = [eval_by_id[eid] for eid in sorted(visited_ids)]
    return lineage
