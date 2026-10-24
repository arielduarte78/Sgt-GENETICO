"""Pruebas unitarias para las estructuras y validaciones del dominio matemático puro."""

from __future__ import annotations

import math
from dataclasses import FrozenInstanceError

import pytest

from kuramoto_templates.domain.conditions import (
    c_a,
    connectivity,
    degrees,
    equilibrium_residual,
    projected_laplacian,
)
from kuramoto_templates.domain.template import Template


def test_template_creation_and_immutability() -> None:
    """Verifica que Template es inmutable y almacena adecuadamente los campos."""
    t = Template(
        x=(0.5, 0.5),
        phi=(0.0, 1.57),
        P=((1.0, 0.5), (0.5, 1.0)),
    )
    assert t.m == 2
    assert t.x == (0.5, 0.5)

    with pytest.raises(FrozenInstanceError):
        # Intentar modificar un atributo en una dataclass frozen debe fallar
        t.x = (0.6, 0.4)  # type: ignore[misc]


def test_template_validation_m_less_than_two() -> None:
    """Verifica que m < 2 sea rechazado."""
    with pytest.raises(ValueError, match="al menos m=2"):
        Template(x=(1.0,), phi=(0.0,))


def test_template_validation_dimension_mismatch() -> None:
    """Verifica el rechazo si len(x) != len(phi)."""
    with pytest.raises(ValueError, match="Discrepancia dimensional"):
        Template(x=(0.5, 0.5), phi=(0.0,))


def test_template_validation_sizes_not_summing_to_one() -> None:
    """Verifica el rechazo si sum(x) != 1.0."""
    with pytest.raises(ValueError, match="suma de tamaños"):
        Template(x=(0.6, 0.6), phi=(0.0, 1.0))


def test_template_validation_strictly_positive_sizes() -> None:
    """Verifica el rechazo si algún x_a <= 0."""
    with pytest.raises(ValueError, match="x_a > 0"):
        Template(x=(1.0, 0.0), phi=(0.0, 1.0))


def test_template_validation_asymmetric_p() -> None:
    """Verifica que matrices P asimétricas sean rechazadas."""
    with pytest.raises(ValueError, match="Asimetría"):
        Template(
            x=(0.5, 0.5),
            phi=(0.0, 1.0),
            P=((1.0, 0.2), (0.8, 1.0)),
        )


def test_template_validation_p_bounds() -> None:
    """Verifica que valores P_ab fuera de [0, 1] sean rechazados."""
    with pytest.raises(ValueError, match="fuera del intervalo"):
        Template(
            x=(0.5, 0.5),
            phi=(0.0, 1.0),
            P=((1.0, 1.5), (1.5, 1.0)),
        )


def test_condition_signatures_raise_without_p() -> None:
    """Verifica que las firmas de condiciones requieran P."""
    t_no_p = Template(x=(0.5, 0.5), phi=(0.0, 1.0))

    with pytest.raises(ValueError, match="sin una matriz de densidades P"):
        degrees(t_no_p)

    with pytest.raises(ValueError, match="sin matriz P"):
        connectivity(t_no_p)

    with pytest.raises(ValueError, match="sin matriz P"):
        equilibrium_residual(t_no_p)

    with pytest.raises(ValueError, match="sin matriz P"):
        c_a(t_no_p)

    with pytest.raises(ValueError, match="sin matriz P"):
        projected_laplacian(t_no_p)


def test_condition_computations_with_p() -> None:
    """Verifica que las firmas de condiciones computen valores finitos y coherentes."""
    t = Template(x=(0.5, 0.5), phi=(0.0, 1.0), P=((1.0, 0.5), (0.5, 1.0)))

    deg = degrees(t)
    assert len(deg) == 2
    assert math.isclose(deg[0], 0.75)
    assert math.isclose(deg[1], 0.75)

    conn = connectivity(t)
    assert math.isclose(conn, 0.75)

    res = equilibrium_residual(t)
    assert len(res) == 2
    assert math.isclose(res[0] + res[1], 0.0, abs_tol=1e-12)

    c_vals = c_a(t)
    assert len(c_vals) == 2
    assert c_vals[0] > 0.0

    lap = projected_laplacian(t)
    assert len(lap) == 1
    assert len(lap[0]) == 1
    assert lap[0][0] > 0.0

