"""Command-line interface (CLI) for genetic algorithm search engine."""

from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path
from typing import Any

from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine


def _write_status(
    output_dir: Path,
    status: str,
    gen: int = 0,
    total_gens: int = 0,
    best_f: float = -1.0,
    best_m: int = 0,
    best_lam: float = 0.0,
    cache_hits: int = 0,
    evals: int = 0,
) -> None:
    """Update status.json atomically."""
    status_path = output_dir / "status.json"
    temp_path = output_dir / "status.tmp"

    data: dict[str, Any] = {
        "status": status,
        "pid": os.getpid(),
        "generation": gen,
        "total_generations": total_gens,
        "best_f": best_f,
        "best_m": best_m,
        "best_lam": best_lam,
        "cache_hits": cache_hits,
        "total_evaluations": evals,
        "updated_at": time.time(),
    }
    with temp_path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    os.replace(temp_path, status_path)


def handle_run(args: argparse.Namespace) -> int:
    """Execute search run."""
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: configuration file {config_path} not found.", file=sys.stderr)
        return 1

    config = GAConfig.from_yaml(config_path)
    output_dir = Path(args.output_dir or f"runs/{config_path.stem}")
    output_dir.mkdir(parents=True, exist_ok=True)

    engine = GAEngine(
        config=config,
        output_dir=output_dir,
        allow_dirty=args.allow_dirty,
        no_seeds=args.no_seeds,
        no_split_merge=args.no_split_merge,
    )

    # Signal handler for clean checkpointing on SIGINT or SIGTERM
    def sig_handler(signum: int, frame: Any) -> None:
        print(
            "\nTermination signal received. Saving final emergency checkpoint...",
            file=sys.stderr,
        )
        engine.save_checkpoint(engine.current_generation)
        engine.update_status_file(status="interrupted")
        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    if hasattr(signal, "SIGTERM"):
        signal.signal(signal.SIGTERM, sig_handler)
    if hasattr(signal, "SIGBREAK"):
        signal.signal(signal.SIGBREAK, sig_handler)

    engine.update_status_file(status="running")

    try:
        results = engine.run()
        engine.update_status_file(status="completed")
        print(f"Run completed successfully. Results saved in {output_dir}")
        msg = (
            f"Best fitness f: {results['best_f']:.12e} "
            f"(m={results['best_m']}, lam*={results['best_lam']:.12e})"
        )
        print(msg)
        return 0
    except Exception as exc:
        _write_status(output_dir, status="error")
        print(f"Fatal error during execution: {exc}", file=sys.stderr)
        return 1


def handle_resume(args: argparse.Namespace) -> int:
    """Resume search from an existing checkpoint."""
    cp_path = Path(args.checkpoint)
    if not cp_path.exists():
        print(f"Error: checkpoint {cp_path} does not exist.", file=sys.stderr)
        return 1

    engine = GAEngine.resume(cp_path, total_generations=args.generations)
    results = engine.run()
    print(f"Resumed run completed. Best fitness: {results['best_f']:.12e}")
    return 0


def handle_status(args: argparse.Namespace) -> int:
    """Display status of a run."""
    output_dir = Path(args.output_dir)
    status_file = output_dir / "status.json"
    if not status_file.exists():
        print(f"No status.json found in {output_dir}")
        return 1

    with status_file.open("r", encoding="utf-8") as f:
        data = json.load(f)

    print("=" * 40)
    print(f"Status: {data.get('status')}")
    print(f"Progress: Gen {data.get('generation')}/{data.get('total_generations')}")
    print(f"Best f: {data.get('best_f')}")
    print(f"Best lam*: {data.get('best_lam')}")
    print(f"Cache hits: {data.get('cache_hits')}")
    print("=" * 40)
    return 0


def handle_tail(args: argparse.Namespace) -> int:
    """Display latest lines of run.jsonl log."""
    output_dir = Path(args.output_dir)
    log_path = output_dir / "run.jsonl"
    if not log_path.exists():
        print(f"No run.jsonl log found in {output_dir}")
        return 1

    with log_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    n = args.lines
    for line in lines[-n:]:
        print(line.strip())
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser."""
    parser = argparse.ArgumentParser(description="Kuramoto GA Search CLI")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # ga run
    p_run = subparsers.add_parser("run", help="Start a new GA search run")
    p_run.add_argument("config", help="Path to GA YAML configuration file")
    p_run.add_argument("--output-dir", "-o", default=None, help="Directory to save run results")
    p_run.add_argument(
        "--allow-dirty", action="store_true", help="Allow uncommitted git working tree"
    )
    p_run.add_argument(
        "--no-seeds", action="store_true", help="Ablation: do not inject analytical seeds"
    )
    p_run.add_argument(
        "--no-split-merge", action="store_true", help="Ablation: disable split/merge"
    )
    p_run.set_defaults(func=handle_run)

    # ga resume
    p_resume = subparsers.add_parser("resume", help="Resume from checkpoint")
    p_resume.add_argument("checkpoint", help="Path to checkpoint JSON file")
    p_resume.add_argument(
        "--generations", "-g", type=int, default=None, help="New total generations target"
    )
    p_resume.set_defaults(func=handle_resume)

    # ga status
    p_status = subparsers.add_parser("status", help="Inspect status of a run")
    p_status.add_argument("--output-dir", "-o", default="runs/latest", help="Run directory")
    p_status.set_defaults(func=handle_status)

    # ga tail
    p_tail = subparsers.add_parser("tail", help="Tail run JSONL log")
    p_tail.add_argument("--output-dir", "-o", default="runs/latest", help="Run directory")
    p_tail.add_argument("-n", "--lines", type=int, default=10, help="Number of lines to tail")
    p_tail.set_defaults(func=handle_tail)

    return parser


def main(args: list[str] | None = None) -> int:
    """Entry point for CLI execution."""
    parser = build_parser()
    parsed = parser.parse_args(args)
    return int(parsed.func(parsed))


if __name__ == "__main__":
    sys.exit(main())
