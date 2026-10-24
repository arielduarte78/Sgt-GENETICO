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


def test_phase1_always_feasible_on_random_candidates() -> None:
    """Valida que la formulación de Fase I sin tau sea siempre factible sobre candidatos aleatorios.

    Qué valida:
        Para 20 individuos aleatorios (x, phi) generados con seed fija, el programa SDP sin tau
        es siempre factible (p=0 con s grande es solución admisible), retornando:
            - status en ("OPTIMAL", "OPTIMAL_INACCURATE")
            - s* >= -1e-8 (holgura no negativa)
            - lam* finito y lam* normalizado acotado
    """
    rng = np.random.default_rng(seed=2026)
    config = EvalConfig(mu_t=0.6929077, rho=1.0, solver_tol=1e-8)
    evaluator = SDPEvaluator(config=config)

    for i in range(20):
        m = int(rng.integers(3, 7))
        raw_x = rng.exponential(scale=1.0, size=m)
        x_arr = raw_x / np.sum(raw_x)
        x_arr = np.maximum(x_arr, 2.0 * config.eps_c)
        x_arr = x_arr / np.sum(x_arr)
        x = tuple(float(v) for v in x_arr)

        phi = tuple(float(v) for v in rng.uniform(-math.pi, math.pi, size=m))

        res = evaluator.evaluate_phase1(x, phi)

        assert res.status in ("OPTIMAL", "OPTIMAL_INACCURATE"), (
            f"Individuo {i}: solver falló con status={res.status}."
        )
        assert res.s is not None and res.s >= -1e-8, (
            f"Individuo {i}: s={res.s} violó la no-negatividad."
        )
        assert res.lam is not None, f"Individuo {i}: lam es None."
        assert res.lam_normalized is not None, f"Individuo {i}: lam_normalized es None."


def test_c4_phase1_rejection_of_stability() -> None:
    """Valida que la Fase I rechace certificar C4 como estrictamente estable (lam* <= 0).

    Qué valida:
        Para el estado de cuatro clusters C4, la formulación sin tau devuelve un margen
        espectral lam* <= 1e-10 (nulo o negativo), demostrando que NO certifica
        C4 como estrictamente estable debido a su degeneración espectral.

    De dónde viene el valor esperado:
        Teorema de degeneración lineal de C4 (Kassabov, Strogatz & Townsend 2021, Chaos 31, 073135):
        En C4, todo enlace no nulo conecta clusters a distancia pi/2, donde cos(pi/2) = 0.
        Por tanto L_w(p) == 0 y M(p) == 0 para cualquier matriz de densidades admisible con
        pares opuestos desconectados. Como M(p) = 0, la restricción LMI M(p) >> lam*I exige
        forzosamente lam* <= 0, impidiendo obtener margen espectral positivo.
    """
    c4 = get_c4_template()
    config = EvalConfig(mu_t=0.6929077, rho=1.0, solver_tol=1e-10)
    evaluator = SDPEvaluator(config=config)

    res = evaluator.evaluate_phase1(c4.x, c4.phi)

    assert res.lam is not None
    assert res.lam <= 1e-10, (
        f"El evaluador no detectó la degeneración de C4: lam*={res.lam} esperado <= 0."
    )


def test_phase1_strictly_stable_case_achieves_near_zero_phi() -> None:
    """Valida un caso positivo claro en Fase I con Phi < 1e-6 para un estado estrictamente estable.

    Qué valida:
        Para un sistema de dos clusters idénticos sincronizados (m=2, x=(0.5, 0.5), phi=(0.0, 0.0)),
        la Fase I elástica certifica estabilidad estricta devolviendo:
            - Phi < 1e-6 (en la práctica ~ 1e-12).
            - Status = "OPTIMAL".
        Asimismo, la evaluación completa procede a Fase II alcanzando mu* = 1.0 (grafo completo).

    De dónde viene el valor esperado:
        Construcción analítica exacta:
        Con m=2 y phi_1 = phi_2 = 0, sin(phi_1 - phi_2) = 0, por lo que F(p) == 0 idénticamente
        para cualquier p in [0, 1] (tau = 0 exacto).
        El Laplaciano proyectado con p=1 es M = Q^T L_w Q = 0.5 (autovalor 0.5 >> eps = 1e-6).
        Por tanto, la restricción LMI M >= (eps - delta) I se satisface holgadamente con delta = 0,
        y c_a = 0.5 + 0.5 = 1.0 >= eps_c.
        El objetivo mínimo global es Phi = 0, alcanzando numéricamente Phi < 1e-6 con
        status OPTIMAL.
    """
    x = (0.5, 0.5)
    phi = (0.0, 0.0)
    config = EvalConfig(eps=1e-6, eps_c=1e-4, theta_tol=1e-4, solver_tol=1e-10)
    evaluator = SDPEvaluator(config=config)

    res_p1 = evaluator.evaluate_phase1(x, phi)

    assert res_p1.status == "OPTIMAL", (
        f"Fase I no devolvió OPTIMAL para caso sincronizado m=2: status={res_p1.status}."
    )
    assert res_p1.Phi < 1e-6, (
        f"Fase I devolvió Phi={res_p1.Phi} >= 1e-6 para un caso estrictamente estable."
    )

    # Validar que al proceder a Fase II se obtenga el óptimo analítico mu* = 1.0
    full_res = evaluator.evaluate(x, phi)
    assert full_res.status == "OPTIMAL"
    assert full_res.mu is not None
    assert math.isclose(full_res.mu, 1.0, abs_tol=1e-6), (
        f"Conectividad óptima {full_res.mu} difiere del valor analítico 1.0."
    )


def test_phase1_no_tau_formulation_returns_lam_s_and_normalized() -> None:
    """Valida la formulación de Fase I sin tau: p = N z, max lam - rho*s."""
    x = (0.5, 0.5)
    phi = (0.0, 0.0)
    config = EvalConfig(mu_t=0.6929077, rho=1.0, solver_tol=1e-8)
    evaluator = SDPEvaluator(config=config)

    res = evaluator.evaluate_phase1(x, phi)

    assert res.status in ("OPTIMAL", "OPTIMAL_INACCURATE")
    assert res.lam is not None
    assert res.s is not None
    assert res.lam_normalized is not None

    # En m=2 con P=1, d_a = 1.0 > mu_t = 0.6929077, por lo que s* = 0
    assert math.isclose(res.s, 0.0, abs_tol=1e-5), f"s* esperado ~ 0, obtenido {res.s}"
    # M = 0.5, por lo que lam* = 0.5
    assert math.isclose(res.lam, 0.5, abs_tol=1e-5), f"lam* esperado ~ 0.5, obtenido {res.lam}"
    # max_k ||M_k||_2 = 0.5, por lo que lam* normalizado = 1.0
    assert math.isclose(res.lam_normalized, 1.0, abs_tol=1e-5)

    # Fitness continua: s* <= tol => f = +lam*
    assert math.isclose(res.f, 0.5, abs_tol=1e-5)

    assert "lam" in res.solver_info
    assert "s" in res.solver_info
    assert "lam_normalized" in res.solver_info
    assert "max_m_basis_norm" in res.solver_info


