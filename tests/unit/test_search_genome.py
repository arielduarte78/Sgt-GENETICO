"""Pruebas unitarias de la estructura del genoma, decodificación y canonicalización (Gauge)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from kuramoto_templates.search.genome import (
    Genome,
    canonicalize,
    decode_phi,
    decode_x,
    to_phenotype,
)


def test_genome_creation_and_bounds() -> None:
    """Valida la instanciación e invariantes del genoma."""
    m = 4
    K_ref = 8
    k = (0, 2, 4, 6)
    sgn = (1, -1, 1, -1)
    ell = (-2.0, -3.0, -2.5, -4.0)
    u = (0.0, 1.0, -0.5, 0.2)

    g = Genome(m=m, K_ref=K_ref, k=k, sgn=sgn, ell=ell, u=u)
    assert g.m == 4
    assert g.K_ref == 8
    assert len(g.k) == 4
    assert len(g.sgn) == 4
    assert len(g.ell) == 4
    assert len(g.u) == 4


def test_genome_validation_errors() -> None:
    """Valida que parámetros inválidos lancen ValueError."""
    # m fuera de [2, 10]
    with pytest.raises(ValueError, match="m debe estar en"):
        Genome(m=1, K_ref=4, k=(0,), sgn=(1,), ell=(-2.0,), u=(0.0,))

    with pytest.raises(ValueError, match="m debe estar en"):
        Genome(
            m=11,
            K_ref=4,
            k=tuple(range(11)),
            sgn=(1,) * 11,
            ell=(-2.0,) * 11,
            u=(0.0,) * 11,
        )

    # K_ref inválido
    with pytest.raises(ValueError, match="K_ref debe pertenecer a"):
        Genome(m=3, K_ref=6, k=(0, 1, 2), sgn=(1, 1, 1), ell=(-2.0, -2.0, -2.0), u=(0.0, 0.0, 0.0))

    # Longitudes inconsistentes
    with pytest.raises(ValueError, match="Longitud"):
        Genome(m=3, K_ref=4, k=(0, 1), sgn=(1, 1, 1), ell=(-2.0, -2.0, -2.0), u=(0.0, 0.0, 0.0))


def test_decode_x_simplex_and_minimum_size() -> None:
    """Valida que decode_x genere puntos en el símplex con x_a >= xmin = 1e-3."""
    xmin = 1e-3
    u = (1.5, -2.0, 0.0, 3.2, -1.1)
    m = len(u)
    x = decode_x(u, xmin=xmin)

    assert len(x) == m
    assert math.isclose(sum(x), 1.0, abs_tol=1e-12)
    for val in x:
        assert val >= xmin - 1e-14
        assert val > 0.0


def test_decode_phi_wrapped_in_minus_pi_to_pi() -> None:
    """Valida que decode_phi envuelva correctamente en (-pi, pi]."""
    K_ref = 4
    k = (0, 1, 2, 3)
    sgn = (1, 1, -1, 0)
    ell = (-2.0, -1.0, -3.0, -2.0)

    phi = decode_phi(K_ref=K_ref, k=k, sgn=sgn, ell=ell)
    assert len(phi) == 4
    for p in phi:
        assert -math.pi < p <= math.pi


def test_to_phenotype_returns_x_and_phi() -> None:
    """Valida to_phenotype."""
    g = Genome(
        m=3,
        K_ref=4,
        k=(0, 1, 2),
        sgn=(1, -1, 1),
        ell=(-2.0, -2.0, -2.0),
        u=(0.0, 0.0, 0.0),
    )
    x, phi = to_phenotype(g)
    assert len(x) == 3
    assert len(phi) == 3
    assert math.isclose(sum(x), 1.0, abs_tol=1e-12)


def test_canonicalize_gauge_invariance_and_idempotence() -> None:
    """Valida la canonicalización: orden por fase, mayor cluster anclado, e idempotencia."""
    g = Genome(
        m=3,
        K_ref=4,
        k=(2, 0, 1),
        sgn=(1, 1, -1),
        ell=(-1.5, -2.0, -3.0),
        u=(1.0, 3.0, 0.5),  # cluster 1 tiene mayor masa
    )

    c1 = canonicalize(g)
    c2 = canonicalize(c1)

    # Idempotencia estricta: canonicalize(canonicalize(g)) == canonicalize(g)
    assert c1.m == c2.m
    assert c1.K_ref == c2.K_ref
    assert c1.k == c2.k
    assert c1.sgn == c2.sgn
    assert np.allclose(c1.ell, c2.ell, atol=1e-12)
    assert np.allclose(c1.u, c2.u, atol=1e-12)

    # Las fases decodificadas deben estar ordenadas
    _, phi_c = to_phenotype(c1)
    for a in range(len(phi_c) - 1):
        assert phi_c[a] <= phi_c[a + 1] + 1e-12
