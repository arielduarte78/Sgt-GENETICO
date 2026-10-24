"""Split mutation operator dividing a cluster into two preserving mass and phase."""

from __future__ import annotations

import math

import numpy as np

from kuramoto_templates.search.genome import (
    Genome,
    to_phenotype,
    wrap_minus_pi_to_pi,
)


def _encode_phase_to_genes(
    phi: float, K_ref: int
) -> tuple[int, int, float]:
    """Find closest anchor k, sign sgn, and log-deviation ell for target phase phi."""
    grid_spacing = 2.0 * math.pi / float(K_ref)
    k = int(round(phi / grid_spacing)) % K_ref
    anchor_phi = wrap_minus_pi_to_pi(k * grid_spacing)

    diff = wrap_minus_pi_to_pi(phi - anchor_phi)
    sgn = 1 if diff >= 0 else -1
    abs_diff = max(abs(diff), 1e-6)
    ell = float(np.clip(math.log10(abs_diff), -6.0, -0.5))
    return k, sgn, ell


def split(
    genome: Genome,
    rng: np.random.Generator,
    delta_phi: float = 0.005,
) -> Genome:
    """Split a cluster into two subclusters.

    Preserves total mass and circular barycenter of the split pair.
    Pure function: returns a new Genome without mutating inputs.
    If genome.m >= 10, returns an identical copy.
    """
    if genome.m >= 10:
        return Genome(
            m=genome.m,
            K_ref=genome.K_ref,
            k=genome.k,
            sgn=genome.sgn,
            ell=genome.ell,
            u=genome.u,
        )

    x, phi = to_phenotype(genome)

    # Pick the cluster with largest mass to split
    idx = int(np.argmax(x))

    x_orig = x[idx]
    x_half = x_orig / 2.0
    phi_orig = phi[idx]

    # Symmetric split in angle preserves exact circular barycenter
    # since x_half * exp(i*(phi+d)) + x_half * exp(i*(phi-d)) = x_orig * cos(d) * exp(i*phi)
    d_phi = min(delta_phi, 0.05)
    phi_1 = wrap_minus_pi_to_pi(phi_orig + d_phi)
    phi_2 = wrap_minus_pi_to_pi(phi_orig - d_phi)

    # Build new x array with m+1 elements
    new_m = genome.m + 1
    new_x = np.delete(np.array(x), idx)
    new_x = np.insert(new_x, idx, [x_half, x_half])
    new_x = new_x / np.sum(new_x)

    # Build new phi array with m+1 elements
    new_phi = np.delete(np.array(phi), idx)
    new_phi = np.insert(new_phi, idx, [phi_1, phi_2])

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
