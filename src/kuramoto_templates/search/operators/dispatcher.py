"""Operator dispatcher selecting and applying operators via single rng.choice call."""

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from kuramoto_templates.search.genome import Genome
from kuramoto_templates.search.operators.anchor_jump import anchor_jump
from kuramoto_templates.search.operators.crossover import crossover
from kuramoto_templates.search.operators.gauss_u import gauss_u
from kuramoto_templates.search.operators.merge import merge
from kuramoto_templates.search.operators.phase_dev import phase_dev
from kuramoto_templates.search.operators.split import split

OperatorFunc = Callable[..., Genome]


@dataclass
class OperatorRegistry:
    """Registry of mutation and crossover operators with associated selection weights."""

    probabilities: dict[str, float] = field(
        default_factory=lambda: {
            "gauss_u": 0.25,
            "phase_dev": 0.25,
            "anchor_jump": 0.20,
            "split": 0.075,
            "merge": 0.075,
            "crossover": 0.15,
        }
    )
    allow_split_merge: bool = True

    def get_normalized_probabilities(self) -> tuple[list[str], np.ndarray]:
        """Return operator names and strictly normalized selection probabilities."""
        names = []
        weights = []

        for name, prob in self.probabilities.items():
            if not self.allow_split_merge and name in ("split", "merge"):
                continue
            if prob > 0.0:
                names.append(name)
                weights.append(prob)

        weights_arr = np.array(weights, dtype=np.float64)
        total = float(np.sum(weights_arr))
        if total <= 0.0:
            raise ValueError("Total operator selection probability must be positive.")
        normalized = weights_arr / total
        return names, normalized


def select_and_apply_operator(
    parent1: Genome,
    parent2: Genome,
    registry: OperatorRegistry,
    rng: np.random.Generator,
) -> tuple[Genome, str]:
    """Select exactly one operator using a single call to rng.choice and apply it.

    Returns the mutated/crossover child and the name of the chosen operator.
    """
    names, probs = registry.get_normalized_probabilities()

    # CRITICAL: exactly ONE call to rng.choice with normalized probabilities
    chosen_op = str(rng.choice(names, p=probs))

    if chosen_op == "gauss_u":
        child = gauss_u(parent1, rng)
    elif chosen_op == "phase_dev":
        child = phase_dev(parent1, rng)
    elif chosen_op == "anchor_jump":
        child = anchor_jump(parent1, rng)
    elif chosen_op == "split":
        child = split(parent1, rng)
    elif chosen_op == "merge":
        child = merge(parent1, rng)
    elif chosen_op == "crossover":
        child = crossover(parent1, parent2, rng)
    else:
        raise ValueError(f"Unknown operator: {chosen_op}")

    return child, chosen_op
