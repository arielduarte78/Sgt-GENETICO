"""Pruebas de integración sobre la interfaz Evaluator y el contrato de EvalResult."""

from __future__ import annotations

import pytest

from kuramoto_templates.evaluators.base import Evaluator
from kuramoto_templates.evaluators.results import EvalResult
from kuramoto_templates.search.genome import Genome


def test_evaluator_cannot_be_instantiated_directly() -> None:
    """Verifica que Evaluator es una clase abstracta pura."""
    with pytest.raises(TypeError):
        Evaluator()  # type: ignore[abstract]


def test_evaluator_subclass_contract() -> None:
    """Verifica que una implementación mock del evaluador respete el contrato de tipos."""

    class MockEvaluator(Evaluator):
        def evaluate(
            self,
            x: tuple[float, ...],
            phi: tuple[float, ...],
            genome: Genome | None = None,
        ) -> EvalResult:
            # Simulación de resultado conforme a EvalResult
            return EvalResult(
                f=0.6875,
                mu=0.6875,
                Phi=0.0,
                P=((1.0, 0.5), (0.5, 1.0)),
                status="OPTIMAL",
                solver_info={"solve_time_sec": 0.012, "iterations": 8},
            )

    evaluator = MockEvaluator()
    res = evaluator.evaluate(x=(0.5, 0.5), phi=(0.0, 1.57))

    assert isinstance(res, EvalResult)
    assert res.f == 0.6875
    assert res.mu == 0.6875
    assert res.Phi == 0.0
    assert res.status == "OPTIMAL"
    assert res.solver_info["iterations"] == 8
