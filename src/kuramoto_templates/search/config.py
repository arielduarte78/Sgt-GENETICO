"""Configuration specification, YAML serialization, and SHA-256 for GA search."""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass, field
from pathlib import Path

import yaml


@dataclass(frozen=True)
class GAConfig:
    """Hyperparameters and execution configuration for Genetic Algorithm search.

    Critical requirement: mu_t is configured explicitly and never hardcoded.
    """

    mu_t: float
    fitness_noise: float = 1e-9
    tol: float = 1e-8
    lam_floor: float = 1e-7
    lam_floor_hat: float = 1e-5
    num_islands: int = 4
    island_pop_size: int = 20
    num_generations: int = 20
    migration_interval: int = 5
    migration_count: int = 2
    tournament_k: int = 3
    use_seeds: bool = True
    allow_split_merge: bool = True
    seed: int = 42
    num_workers: int = 1
    checkpoint_interval: int = 5
    evaluator: str = "plain"
    min_m: int = 2
    max_m: int = 10
    operator_probs: dict[str, float] = field(
        default_factory=lambda: {
            "gauss_u": 0.25,
            "phase_dev": 0.25,
            "anchor_jump": 0.20,
            "split": 0.075,
            "merge": 0.075,
            "crossover": 0.15,
        }
    )

    def __post_init__(self) -> None:
        """Validate parameter bounds."""
        if self.evaluator not in ("plain", "scaled"):
            raise ValueError(f"evaluator must be 'plain' or 'scaled', got '{self.evaluator}'.")
        if not (self.fitness_noise >= 0.0):
            raise ValueError(f"fitness_noise must be non-negative, got {self.fitness_noise}.")
        if not (self.tol > 0.0):
            raise ValueError(f"tol must be strictly positive, got {self.tol}.")
        if not (self.lam_floor > 0.0):
            raise ValueError(f"lam_floor must be strictly positive, got {self.lam_floor}.")
        if not (self.lam_floor_hat > 0.0):
            raise ValueError(f"lam_floor_hat must be strictly positive, got {self.lam_floor_hat}.")
        if self.num_islands < 1:
            raise ValueError(f"num_islands must be >= 1, got {self.num_islands}.")
        if self.island_pop_size < 2:
            raise ValueError(f"island_pop_size must be >= 2, got {self.island_pop_size}.")
        if self.num_generations < 1:
            raise ValueError(f"num_generations must be >= 1, got {self.num_generations}.")
        if self.tournament_k < 1:
            raise ValueError(f"tournament_k must be >= 1, got {self.tournament_k}.")
        if self.min_m < 2:
            raise ValueError(f"min_m must be >= 2, got {self.min_m}.")
        if self.max_m < self.min_m:
            raise ValueError(f"max_m must be >= min_m, got {self.max_m} < {self.min_m}.")

    def to_yaml(self, path: Path | str) -> None:
        """Serialize configuration to a YAML file."""
        target_path = Path(path)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        data = asdict(self)
        with target_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(data, f, sort_keys=True, indent=2)

    @classmethod
    def from_yaml(cls, path: Path | str) -> GAConfig:
        """Load and instantiate configuration from a YAML file."""
        target_path = Path(path)
        with target_path.open("r", encoding="utf-8") as f:
            raw = yaml.safe_load(f)
        if not isinstance(raw, dict):
            raise ValueError(f"Invalid YAML content in {path}, expected dict.")
        return cls(**raw)

    def compute_sha256(self) -> str:
        """Compute deterministic SHA-256 hex digest of the canonical YAML representation."""
        data = asdict(self)
        canonical_str = yaml.safe_dump(data, sort_keys=True, indent=2)
        return hashlib.sha256(canonical_str.encode("utf-8")).hexdigest()
