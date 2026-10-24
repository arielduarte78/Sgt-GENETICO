"""Phase deviation mutation operator modifying ell and sgn."""

from __future__ import annotations

import numpy as np

from kuramoto_templates.search.genome import Genome


def phase_dev(
    genome: Genome,
    rng: np.random.Generator,
    sigma_ell: float = 0.5,
    flip_sgn_prob: float = 0.2,
) -> Genome:
    """Mutate phase deviations ell and signs sgn.

    Pure function: returns a new Genome instance without mutating inputs.
    """
    new_ell = list(genome.ell)
    new_sgn = list(genome.sgn)

    num_to_mutate = int(rng.integers(1, genome.m + 1))
    indices = rng.choice(genome.m, size=num_to_mutate, replace=False)

    for idx in indices:
        delta = float(rng.normal(loc=0.0, scale=sigma_ell))
        new_ell[idx] = float(np.clip(new_ell[idx] + delta, -6.0, -0.5))

        if rng.uniform(0.0, 1.0) < flip_sgn_prob:
            new_sgn[idx] = -new_sgn[idx] if new_sgn[idx] != 0 else 1

    return Genome(
        m=genome.m,
        K_ref=genome.K_ref,
        k=genome.k,
        sgn=tuple(new_sgn),
        ell=tuple(new_ell),
        u=genome.u,
    )
