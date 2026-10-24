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


def test_continuous_fitness_with_lam_floor_below_floor() -> None:
    """Valida que si s* <= tol pero lam* < lam_floor, f = lam* - lam_floor < 0."""
    tol = 1e-8
    lam_floor = 1e-7
    # Caso como ID 39260: s* dentro de tolerancia, lam* < 1e-7
    s_val = 3.61e-9
    lam_val = 1.715719955466e-08

    f = compute_continuous_fitness(lam=lam_val, s=s_val, tol=tol, lam_floor=lam_floor)
    expected_f = lam_val - lam_floor
    assert f < 0.0
    assert math.isclose(f, expected_f, abs_tol=1e-12)


def test_continuous_fitness_with_lam_floor_above_or_equal_floor() -> None:
    """Valida que si s* <= tol y lam* >= lam_floor, f = lam* >= lam_floor > 0."""
    tol = 1e-8
    lam_floor = 1e-7
    s_val = 1e-9
    lam_val = 2.5e-7

    f = compute_continuous_fitness(lam=lam_val, s=s_val, tol=tol, lam_floor=lam_floor)
    assert f == lam_val
    assert f >= lam_floor


def test_continuous_fitness_rejects_non_positive_lam_floor() -> None:
    """Valida que lam_floor <= 0 lance ValueError."""
    with pytest.raises(ValueError, match="lam_floor debe ser estrictamente positiva"):
        compute_continuous_fitness(lam=0.1, s=0.0, tol=1e-6, lam_floor=0.0)


def test_evalconfig_has_lam_floor_default_and_validation() -> None:
    """Valida que EvalConfig tenga lam_floor=1e-7 por defecto y valide positividad."""
    from kuramoto_templates.evaluators.config import EvalConfig

    cfg = EvalConfig()
    assert hasattr(cfg, "lam_floor")
    assert cfg.lam_floor == 1e-7

    with pytest.raises(ValueError, match="lam_floor debe ser estrictamente positivo"):
        EvalConfig(lam_floor=-1e-7)


def test_id_39260_float64_phenotype_yields_negative_fitness() -> None:
    """TEST REQUERIDO: ID 39260 (x y phi float64 del checkpoint) debe dar f < 0 con el piso."""
    from kuramoto_templates.evaluators.config import EvalConfig
    from kuramoto_templates.evaluators.sdp import SDPEvaluator

    # Coordenadas exactas float64 de ID 39260
    x = (0.2494426094558034, 0.25018825224215374, 0.25053408179500397, 0.2498350565070389)
    phi = (-3.141591616161418, -1.5707953267948964, 1.000000000139778e-06, 1.5707973267948967)

    cfg = EvalConfig(mu_t=0.6929077, solver_tol=1e-8, lam_floor=1e-7)
    evaluator = SDPEvaluator(config=cfg)

    res = evaluator.evaluate_phase1(x, phi)
    assert res.f is not None
    # Con lam_floor=1e-7, aunque Clarabel reporte lam*~1.7e-8 y s*~3.6e-9, f = lam* - 1e-7 < 0
    assert res.f < 0.0, f"ID 39260 debió dar f < 0, pero dio f = {res.f}"

