"""Pruebas unitarias de las propiedades matemáticas y cotas de la Fase I (Suite B).

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección IV (convex SDP para P dados x y phi).
    - Plan de Diseño del Algoritmo Genético, Sección 4.1 (Fase I elástica).
"""

from __future__ import annotations

import math

import numpy as np

from kuramoto_templates.cases.c4 import get_c4_template
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def test_phase1_elastic_bounds_on_random_candidates() -> None:
    """Valida que la métrica elástica Phi de la Fase I esté acotada estrictamente en [0, 1].

    Qué valida:
        Para 50 individuos aleatorios (x, phi) generados con seed fija, la función objetivo
        Phi(x, phi) calculada por el programa SDP elástico de Fase I satisface:
            0.0 <= Phi <= 1.0 + 1e-8.

    De dónde viene el valor esperado:
        Propiedad matemática formal de la formulación elástica:
            Phi = min { tau / (eps * theta_tol) + delta / eps }
        El punto p = 0 (grafo completamente vacío) siempre es factible eligiendo:
          - tau = 0 (pues F(0) = 0 idénticamente),
          - delta = eps (con lo cual M(0) = 0 >= (eps - delta) * I = 0 * I),
          - c_a(0) = x_a >= eps_c (válido para todo cluster no degenerado x_a >= 1e-4).
        El valor del objetivo en dicho punto es 0 + eps / eps = 1. Como el solver minimiza
        sobre todo p in [0, 1]^K, el óptimo no puede superar 1.0. Dado que tau >= 0 y delta >= 0,
        el objetivo no puede ser negativo, garantizando 0 <= Phi <= 1.
    """
    rng = np.random.default_rng(seed=2026)
    config = EvalConfig(eps=1e-6, eps_c=1e-4, theta_tol=1e-4, solver_tol=1e-10)
    evaluator = SDPEvaluator(config=config)

    for i in range(50):
        m = int(rng.integers(3, 7))
        # Generar x en el simplex con x_a >= 2 * eps_c
        raw_x = rng.exponential(scale=1.0, size=m)
        x_arr = raw_x / np.sum(raw_x)
        # Asegurar cota mínima estricta para evitar clusters degenerados
        x_arr = np.maximum(x_arr, 2.0 * config.eps_c)
        x_arr = x_arr / np.sum(x_arr)
        x = tuple(float(v) for v in x_arr)

        phi = tuple(float(v) for v in rng.uniform(-math.pi, math.pi, size=m))

        res = evaluator.evaluate_phase1(x, phi)

        assert res.Phi >= -1e-8, f"Individuo {i}: Phi={res.Phi} violó la cota inferior 0."
        assert res.Phi <= 1.0 + 1e-6, f"Individuo {i}: Phi={res.Phi} violó la cota superior 1."


def test_c4_phase1_rejection_of_stability() -> None:
    """Valida que la Fase I rechace certificar el estado C4 como estrictamente estable (Phi > 0).

    Qué valida:
        Para el estado de cuatro clusters C4, con margen espectral eps = 1e-6 > 0,
        el evaluador devuelve Phi > 0 (típicamente Phi >= 0.999), demostrando que NO certifica
        C4 como estrictamente estable debido a su degeneración espectral.

    De dónde viene el valor esperado:
        Teorema de degeneración lineal de C4 (Kassabov, Strogatz & Townsend 2021, Chaos 31, 073135):
        En C4, todo enlace no nulo conecta clusters a distancia pi/2, donde cos(pi/2) = 0.
        Por tanto L_w(p) == 0 y M(p) == 0 para cualquier matriz de densidades admisible con
        pares opuestos desconectados. Como M(p) no puede tener autovalores >= eps, la restricción
        LMI M(p) >> (eps - delta) I exige forzosamente delta >= eps, produciendo delta / eps = 1,
        y consecuentemente Phi >= 1.0 a precisión del solver.
    """
    c4 = get_c4_template()
    config = EvalConfig(eps=1e-6, eps_c=1e-4, theta_tol=1e-4, solver_tol=1e-10)
    evaluator = SDPEvaluator(config=config)

    res = evaluator.evaluate_phase1(c4.x, c4.phi)

    assert res.Phi > 0.5, (
        f"El evaluador no detectó la degeneración de C4: Phi={res.Phi} esperado ~1.0."
    )
    assert res.Phi <= 1.0 + 1e-6, f"Phi de C4 ({res.Phi}) excede la cota superior unitaria."
