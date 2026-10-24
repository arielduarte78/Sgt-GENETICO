"""Anchor jump mutation operator shifting cluster reference anchors k."""

from __future__ import annotations

import numpy as np

from kuramoto_templates.search.genome import Genome


def anchor_jump(
    genome: Genome,
    rng: np.random.Generator,
    jump_type_prob: float = 0.5,
) -> Genome:
    """Mutate reference anchors k by discrete shifts modulo K_ref.

    Pure function: returns a new Genome instance without mutating inputs.
    """
    new_k = list(genome.k)
    idx = int(rng.integers(0, genome.m))

    if rng.uniform(0.0, 1.0) < jump_type_prob:
        step = int(rng.choice([-1, 1]))
        new_k[idx] = (new_k[idx] + step) % genome.K_ref
    else:
        new_k[idx] = int(rng.integers(0, genome.K_ref))

    return Genome(
        m=genome.m,
        K_ref=genome.K_ref,
        k=tuple(new_k),
        sgn=genome.sgn,
        ell=genome.ell,
        u=genome.u,
    )
