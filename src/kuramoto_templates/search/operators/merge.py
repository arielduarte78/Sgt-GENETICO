"""Merge mutation operator combining two clusters into one preserving mass and phase."""

from __future__ import annotations

import math

import numpy as np

from kuramoto_templates.search.genome import (
    Genome,
    to_phenotype,
    wrap_minus_pi_to_pi,
)
from kuramoto_templates.search.operators.split import _encode_phase_to_genes


def merge(
    genome: Genome,
    rng: np.random.Generator,
    min_m: int = 2,
) -> Genome:
    """Merge two neighboring clusters into one.

    Preserves combined mass and circular barycenter of the merged pair.
    Pure function: returns a new Genome without mutating inputs.
    If genome.m <= min_m, returns an identical copy.
    """
    if genome.m <= min_m:
        return Genome(
            m=genome.m,
            K_ref=genome.K_ref,
            k=genome.k,
            sgn=genome.sgn,
            ell=genome.ell,
            u=genome.u,
        )

    x, phi = to_phenotype(genome)

    # Find the pair of clusters with minimum circular distance
    m = genome.m
    min_dist = float("inf")
    best_pair = (0, 1)

    for i in range(m):
        for j in range(i + 1, m):
            dist = abs(wrap_minus_pi_to_pi(phi[i] - phi[j]))
            if dist < min_dist:
                min_dist = dist
                best_pair = (i, j)

    idx1, idx2 = best_pair
    x1, x2 = x[idx1], x[idx2]
    phi1, phi2 = phi[idx1], phi[idx2]

    # Combined mass
    x_c = x1 + x2
    # Combined circular barycenter: Arg(x1*exp(i*phi1) + x2*exp(i*phi2))
    c_x = x1 * math.cos(phi1) + x2 * math.cos(phi2)
    c_y = x1 * math.sin(phi1) + x2 * math.sin(phi2)
    phi_c = wrap_minus_pi_to_pi(math.atan2(c_y, c_x))

    # Construct new phenotype with m-1 elements
    new_m = m - 1
    hi_idx, lo_idx = max(idx1, idx2), min(idx1, idx2)

    new_x = np.delete(np.array(x), [lo_idx, hi_idx])
    new_x = np.append(new_x, x_c)
    new_x = new_x / np.sum(new_x)

    new_phi = np.delete(np.array(phi), [lo_idx, hi_idx])
    new_phi = np.append(new_phi, phi_c)

    # Convert new_x to logits u
    xmin = 1e-3
    p = np.maximum(new_x - xmin, 1e-12) / (1.0 - new_m * xmin)
    new_u_arr = np.log(p)
    new_u_arr = new_u_arr - np.mean(new_u_arr)

    new_k: list[int] = []
    new_sgn: list[int] = []
    new_ell: list[float] = []

    for i in range(new_m):
        k_val, sgn_val, ell_val = _encode_phase_to_genes(float(new_phi[i]), genome.K_ref)
        new_k.append(k_val)
        new_sgn.append(sgn_val)
        new_ell.append(ell_val)

    return Genome(
        m=new_m,
        K_ref=genome.K_ref,
        k=tuple(new_k),
        sgn=tuple(new_sgn),
        ell=tuple(new_ell),
        u=tuple(float(v) for v in new_u_arr),
    )
