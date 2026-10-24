"""Pruebas unitarias de la función de fitness continua (Item 2).

Definición:
    si s* > tol:  f = -s*     (distancia a la conectividad objetivo, <= 0)
    si s* <= tol: f = +lam*   (margen espectral, >= 0)
Sin pesos 1e10 ni log10 de penalización.
"""

from __future__ import annotations

import math

import pytest

from kuramoto_templates.evaluators.fitness import compute_continuous_fitness


def test_continuous_fitness_when_slack_exceeds_tolerance() -> None:
    """Valida que si s* > tol, f = -s* (penalización proporcional a la distancia)."""
    tol = 1e-6
    s_val = 0.05
    lam_val = 0.02

    f = compute_continuous_fitness(lam=lam_val, s=s_val, tol=tol)
    assert f == -s_val
    assert f < 0.0


def test_continuous_fitness_when_slack_is_zero_or_below_tolerance() -> None:
    """Valida que si s* <= tol, f = +lam* (recompensa igual al margen espectral)."""
    tol = 1e-6
    lam_val = 0.035

    # s = 0 exacto
    f_zero = compute_continuous_fitness(lam=lam_val, s=0.0, tol=tol)
    assert math.isclose(f_zero, lam_val, abs_tol=1e-12)

    # s dentro de tolerancia
    f_tol = compute_continuous_fitness(lam=lam_val, s=1e-7, tol=tol)
    assert math.isclose(f_tol, lam_val, abs_tol=1e-12)


def test_continuous_fitness_exact_boundary() -> None:
    """Valida el comportamiento en la frontera exacta s == tol."""
    tol = 1e-6
    lam_val = 0.015

    f = compute_continuous_fitness(lam=lam_val, s=tol, tol=tol)
    assert math.isclose(f, lam_val, abs_tol=1e-12)


def test_continuous_fitness_rejects_negative_tolerance() -> None:
    """Valida que tol negativa o cero lance ValueError."""
    with pytest.raises(ValueError, match="tol debe ser estrictamente positiva"):
        compute_continuous_fitness(lam=0.1, s=0.0, tol=-1e-6)
