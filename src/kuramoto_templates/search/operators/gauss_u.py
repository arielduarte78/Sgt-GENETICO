"""Log-multiscale Gaussian mutation on mass weights u."""

from __future__ import annotations

import numpy as np

from kuramoto_templates.search.genome import Genome


def gauss_u(
    genome: Genome,
    rng: np.random.Generator,
    log_scale_min: float = -3.0,
    log_scale_max: float = 0.0,
) -> Genome:
    """Apply log-multiscale Gaussian mutation to mass logits u.

    Pure function: returns a new Genome instance without mutating inputs.
    The scale sigma is drawn log-uniformly from [10^log_scale_min, 10^log_scale_max].
    """
    log_sigma = float(rng.uniform(log_scale_min, log_scale_max))
    sigma = float(10.0**log_sigma)

    delta = rng.normal(loc=0.0, scale=sigma, size=genome.m)
    new_u = tuple(float(u_val + d) for u_val, d in zip(genome.u, delta, strict=False))

    return Genome(
        m=genome.m,
        K_ref=genome.K_ref,
        k=genome.k,
        sgn=genome.sgn,
        ell=genome.ell,
        u=new_u,
    )
