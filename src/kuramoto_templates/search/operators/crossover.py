"""Crossover operator pairing clusters between two parents by circular distance."""

from __future__ import annotations

import numpy as np

from kuramoto_templates.search.genome import (
    Genome,
    to_phenotype,
    wrap_minus_pi_to_pi,
)


def crossover(
    g1: Genome,
    g2: Genome,
    rng: np.random.Generator,
) -> Genome:
    """Pair clusters between parents by circular distance and recombine genes.

    Pure function: returns a new Genome instance without mutating inputs.
    """
    if rng.uniform(0.0, 1.0) < 0.5:
        p_main, p_aux = g1, g2
    else:
        p_main, p_aux = g2, g1

    m = p_main.m
    _, phi_main = to_phenotype(p_main)
    _, phi_aux = to_phenotype(p_aux)

    # For each cluster in p_main, find closest cluster in p_aux by circular distance
    m_aux = p_aux.m
    paired_indices: list[int] = []
    for i in range(m):
        best_j = 0
        min_d = float("inf")
        for j in range(m_aux):
            dist = abs(wrap_minus_pi_to_pi(phi_main[i] - phi_aux[j]))
            if dist < min_d:
                min_d = dist
                best_j = j
        paired_indices.append(best_j)

    new_k = list(p_main.k)
    new_sgn = list(p_main.sgn)
    new_ell = list(p_main.ell)
    new_u = list(p_main.u)

    for i in range(m):
        j = paired_indices[i]

        # Recombine mass logits u
        if rng.uniform(0.0, 1.0) < 0.5:
            new_u[i] = p_aux.u[j]

        # Recombine phase deviations / anchors
        if rng.uniform(0.0, 1.0) < 0.5:
            if p_main.K_ref == p_aux.K_ref:
                new_k[i] = p_aux.k[j]
                new_sgn[i] = p_aux.sgn[j]
                new_ell[i] = p_aux.ell[j]
            else:
                aux_phi = phi_aux[j]
                grid_spacing = 2.0 * np.pi / float(p_main.K_ref)
                k_val = int(round(aux_phi / grid_spacing)) % p_main.K_ref
                anchor_phi = wrap_minus_pi_to_pi(k_val * grid_spacing)
                diff = wrap_minus_pi_to_pi(aux_phi - anchor_phi)
                new_k[i] = k_val
                new_sgn[i] = 1 if diff >= 0 else -1
                new_ell[i] = float(np.clip(np.log10(max(abs(diff), 1e-6)), -6.0, -0.5))

    return Genome(
        m=m,
        K_ref=p_main.K_ref,
        k=tuple(new_k),
        sgn=tuple(new_sgn),
        ell=tuple(new_ell),
        u=tuple(new_u),
    )
