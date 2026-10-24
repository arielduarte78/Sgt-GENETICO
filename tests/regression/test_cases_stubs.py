"""Pruebas de regresión para casos de prueba base y consistencia de stubs."""

from __future__ import annotations

import math

import pytest

from kuramoto_templates.cases.c4 import get_c4_template
from kuramoto_templates.cases.canale import get_canale_template
from kuramoto_templates.cases.s7 import (
    MU_ZERO_CITADO,
    Q_STAR_CITADO,
    W_STAR_CITADO,
    get_s7_template,
)
from kuramoto_templates.domain.template import Template


def test_c4_state_validity() -> None:
    """Verifica que el estado analítico C4 genera una Template válida."""
    c4 = get_c4_template()
    assert isinstance(c4, Template)
    assert c4.m == 4
    assert len(c4.x) == 4
    assert len(c4.phi) == 4
    assert c4.P is not None
    assert len(c4.P) == 4

    # Verificar tamaños simétricos de 1/4
    for xa in c4.x:
        assert math.isclose(xa, 0.25)

    # Grados analíticos: d_a = 1/4 + 1/4 + 1/4 = 3/4
    # Para cada fila de P, hay 3 unos y 1 cero
    for row in c4.P:
        assert sum(row) == 3.0


def test_canale_is_explicit_stub() -> None:
    """Verifica que el caso de Canale sea un stub honesto con TODO explícito."""
    with pytest.raises(NotImplementedError, match="TODO: Implementar get_canale_template"):
        get_canale_template()


def test_s7_template_structure() -> None:
    """Verifica que la familia S7 construya una Template válida con m=7 clusters."""
    t = get_s7_template(h=0.002, t=0.0)
    assert t.m == 7
    assert len(t.x) == 7
    assert len(t.phi) == 7
    assert t.P is not None
    assert len(t.P) == 7



def test_s7_cited_constants_algebraic_consistency() -> None:
    """Verifica la consistencia algebraica interna de las constantes citadas de O'Keeffe (2026).

    Evalúa:
      1. Polinomio irreducible de grado 8 (Ec. 3):
         P8(q) = 9q^8 + 72q^7 + 288q^6 + 448q^5 + 470q^4 + 88q^3 - 208q^2 - 128q - 15 = 0
      2. Razón de apertura w* (Ec. 11):
         w^2 = (q^2 + 6q + 1) / (q + 1)^2
      3. Cota asintótica mu_0 (Ec. 2 y 10):
         mu_0 = (4 + 2*q) / (7 + q)
    """
    q = Q_STAR_CITADO
    poly_eval = (
        9 * (q**8)
        + 72 * (q**7)
        + 288 * (q**6)
        + 448 * (q**5)
        + 470 * (q**4)
        + 88 * (q**3)
        - 208 * (q**2)
        - 128 * q
        - 15
    )
    assert abs(poly_eval) < 1e-9, (
        f"El valor citado q*={q} no anula el polinomio de grado 8: {poly_eval}"
    )
    expected_w_squared = (q**2 + 6 * q + 1) / ((q + 1) ** 2)
    assert math.isclose(W_STAR_CITADO**2, expected_w_squared, rel_tol=1e-9)

    expected_mu0 = (4 + 2 * q) / (7 + q)
    assert math.isclose(MU_ZERO_CITADO, expected_mu0, rel_tol=1e-9)
