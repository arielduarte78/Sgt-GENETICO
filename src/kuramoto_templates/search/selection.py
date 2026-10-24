"""Selection mechanisms, fitness noise tie-breaking, and deterministic crowding."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kuramoto_templates.search.genome import Genome, wrap_minus_pi_to_pi


@dataclass(frozen=True)
class Individual:
    """Evaluated candidate individual within the population."""

    id: int
    genome: Genome
    phenotype_x: tuple[float, ...]
    phenotype_phi: tuple[float, ...]
    s_star: float
    lam_star: float
    lam_normalized: float
    f: float
    f_sel: float
    parent_ids: tuple[int, ...]
    operator_name: str
    generation: int = 0
    solver_status: str = "OPTIMAL"
    eval_time: float = 0.0
    mu_star: float | None = None


def compute_f_sel(s_star: float, lam_norm: float, tol: float = 1e-8) -> float:
    """Compute selection fitness f_sel.

    f_sel = -s* if s* > tol else lam_norm.
    """
    if s_star > tol:
        return -float(s_star)
    return float(lam_norm)


def compare_individuals(
    ind1: Individual,
    ind2: Individual,
    fitness_noise: float = 1e-9,
) -> int:
    """Compare two individuals considering fitness noise and tie-breaking rules.

    Returns:
        > 0 if ind1 is better than ind2.
        < 0 if ind2 is better than ind1.
        0 if they are identical.

    Rules:
        - If |f_sel_1 - f_sel_2| < fitness_noise:
            - Smaller number of clusters m wins (parsimony).
            - If m is equal, smaller unique id wins (strict determinism).
        - Otherwise, higher f_sel wins.
    """
    diff = ind1.f_sel - ind2.f_sel
    if abs(diff) < fitness_noise:
        # Tie-break 1: smaller m wins
        if ind1.genome.m != ind2.genome.m:
            return 1 if ind1.genome.m < ind2.genome.m else -1
        # Tie-break 2: smaller id wins
        if ind1.id != ind2.id:
            return 1 if ind1.id < ind2.id else -1
        return 0

    return 1 if diff > 0.0 else -1


def tournament_selection(
    population: list[Individual],
    k: int,
    fitness_noise: float,
    rng: np.random.Generator,
) -> Individual:
    """Perform k-tournament selection returning the winning individual."""
    pop_size = len(population)
    if pop_size == 0:
        raise ValueError("Cannot perform selection on empty population.")

    k_eff = min(k, pop_size)
    chosen_indices = rng.choice(pop_size, size=k_eff, replace=False)

    best = population[int(chosen_indices[0])]
    for idx in chosen_indices[1:]:
        candidate = population[int(idx)]
        if compare_individuals(candidate, best, fitness_noise=fitness_noise) > 0:
            best = candidate

    return best


def phenotype_distance(ind1: Individual, ind2: Individual) -> float:
    """Compute circular phenotype distance between two individuals.

    Measures phase and mass discrepancies.
    """
    m1, m2 = ind1.genome.m, ind2.genome.m
    phi1, phi2 = ind1.phenotype_phi, ind2.phenotype_phi
    x1, x2 = ind1.phenotype_x, ind2.phenotype_x

    # If dimension matches, direct pairwise distance
    if m1 == m2:
        phase_dist = sum(
            abs(wrap_minus_pi_to_pi(p1 - p2)) for p1, p2 in zip(phi1, phi2, strict=True)
        )
        mass_dist = sum(abs(v1 - v2) for v1, v2 in zip(x1, x2, strict=True))
        return float(phase_dist + mass_dist)

    # If dimension differs, circular distance to nearest cluster
    d_total = abs(m1 - m2) * 1.0
    for p1 in phi1:
        d_min = min(abs(wrap_minus_pi_to_pi(p1 - p2)) for p2 in phi2)
        d_total += d_min
    return float(d_total)


def deterministic_crowding_replace(
    parent1: Individual,
    parent2: Individual,
    child1: Individual,
    child2: Individual,
    fitness_noise: float = 1e-9,
) -> tuple[Individual, Individual]:
    """Deterministic crowding replacement between parent pair and children pair."""
    d11 = phenotype_distance(parent1, child1)
    d22 = phenotype_distance(parent2, child2)
    d12 = phenotype_distance(parent1, child2)
    d21 = phenotype_distance(parent2, child1)

    if (d11 + d22) <= (d12 + d21):
        # child1 competes with parent1, child2 with parent2
        next1 = child1 if compare_individuals(child1, parent1, fitness_noise) > 0 else parent1
        next2 = child2 if compare_individuals(child2, parent2, fitness_noise) > 0 else parent2
    else:
        # child2 competes with parent1, child1 with parent2
        next1 = child2 if compare_individuals(child2, parent1, fitness_noise) > 0 else parent1
        next2 = child1 if compare_individuals(child1, parent2, fitness_noise) > 0 else parent2

    return next1, next2
