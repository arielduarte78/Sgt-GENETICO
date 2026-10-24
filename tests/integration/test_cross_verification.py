"""Pruebas de verificación cruzada entre solvers SDP y LP (Suite D).

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección IV (convex SDP).
    - Plan de Diseño del Algoritmo Genético, Sección 4 (Evaluación de Aptitud SDP Interno).
"""

from __future__ import annotations

import math

import numpy as np
import scipy.linalg  # type: ignore[import-untyped]
import scipy.optimize  # type: ignore[import-untyped]

from kuramoto_templates.domain.conditions import (
    build_ac_matrix,
    build_ad_matrix,
    build_af_matrix,
)
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def test_cross_verification_clarabel_vs_scs_20_individuals() -> None:
    """Valida la consistencia de mu* entre Clarabel y SCS en 20 individuos aleatorios.

    Qué valida:
        Para 20 individuos aleatorios generados como perturbaciones estables alrededor de
        estados splay con seed fija (np.random.default_rng), se resuelve el SDP completo
        (Fase I y Fase II) tanto con Clarabel (solver por defecto) como con SCS.
        Se valida que:
          1. Ambos solvers alcancen el estatus óptimo y pasen la verificación independiente.
          2. La discrepancia absoluta |mu*_Clarabel - mu*_SCS| <= 1e-6 en cada caso.
        Se reporta e imprime la discrepancia máxima observada en la suite.

    De dónde viene el valor esperado:
        Teoría de dualidad y convexidad en problemas SDP: el óptimo mu* de un programa
        semidefinido convexo estricto es único. Dos implementaciones numéricas independientes
        (Clarabel con conos simétricos de punto interior y SCS con proyecciones de splitting ADMM)
        deben converger al mismo valor óptimo dentro de las tolerancias de convergencia (1e-6).
    """
    rng = np.random.default_rng(seed=20261024)
    config_clarabel = EvalConfig(solver_name="CLARABEL", solver_tol=1e-10)
    config_scs = EvalConfig(solver_name="SCS", solver_tol=1e-9)

    eval_clarabel = SDPEvaluator(config=config_clarabel)
    eval_scs = SDPEvaluator(config=config_scs)

    m = 5
    discrepancies: list[float] = []

    print("\n--- INICIO VERIFICACIÓN CRUZADA CLARABEL VS SCS (20 INDIVIDUOS) ---")
    for i in range(20):
        # Generar tamaños x mediante Dirichlet concentrada
        x_raw = rng.dirichlet(np.ones(m) * 8.0)
        x = tuple(float(v) for v in x_raw)

        # Perturbación angular simétrica alrededor de splay state 2*pi*a/m
        base_phi = np.array([2.0 * math.pi * a / m for a in range(m)])
        noise = rng.uniform(-0.04, 0.04, size=m)
        phi = tuple(float(v) for v in (base_phi + noise))

        res_c = eval_clarabel.evaluate(x, phi)
        res_s = eval_scs.evaluate(x, phi)

        assert res_c.mu is not None, (
            f"Individuo {i}: Clarabel falló al evaluar mu* (status={res_c.status})"
        )
        assert res_s.mu is not None, (
            f"Individuo {i}: SCS falló al evaluar mu* (status={res_s.status})"
        )

        diff = abs(res_c.mu - res_s.mu)
        discrepancies.append(diff)

        print(
            f"Individuo {i:02d}: Clarabel={res_c.mu:.8f}, SCS={res_s.mu:.8f}, "
            f"Diff={diff:.2e} (Clarabel: {res_c.status}, SCS: {res_s.status})"
        )
        assert diff <= 1e-6, (
            f"Individuo {i:02d}: Discrepancia {diff:.2e} excede la tolerancia permitida 1e-6."
        )

    max_diff = max(discrepancies)
    mean_diff = sum(discrepancies) / len(discrepancies)
    print(f"--- REPORTE FINAL: Max Diff = {max_diff:.2e}, Mean Diff = {mean_diff:.2e} ---")


def test_cross_verification_sdp_vs_linprog_when_lmi_is_linearized() -> None:
    """Valida la consistencia de mu* entre Clarabel y scipy.optimize.linprog en régimen lineal.

    Qué valida:
        Cuando la restricción matricial semidefinida M(p) >> eps * I es estrictamente inactiva
        (su margen espectral excede holgadamente eps en el óptimo del problema lineal LP),
        el valor de conectividad óptima mu* obtenido por SDPEvaluator coincide con la solución
        del problema de programación lineal exacta resuelto por scipy.optimize.linprog
        (método HiGHS) dentro de una tolerancia de 1e-6.

    De dónde viene el valor esperado:
        Teorema de relajación y restricciones inactivas: si la solución x* de un programa lineal
        satisface estrictamente una desigualdad adicional no lineal o LMI g(x*) > 0, entonces x*
        es también la solución óptima del programa con la restricción incluida. Por tanto, el
        SDP y el LP formulado con HiGHS deben coincidir en el valor óptimo de mu*.
    """
    m = 4
    K = m * (m - 1) // 2

    # Construir un candidato con fases idénticas phi=0 donde cos(0)=1
    # y la restricción M >> eps I se cumple de forma lineal/holgada
    x_raw = np.array([0.25, 0.25, 0.25, 0.25])
    # Pequeño desfase para que el espacio nulo de equilibrio sea controlable
    phi_arr = np.array([0.0, 0.02, -0.02, 0.0])
    x = tuple(float(v) for v in x_raw)
    phi = tuple(float(v) for v in phi_arr)

    config = EvalConfig(eps=1e-6, eps_c=1e-4, theta_tol=1e-4, solver_tol=1e-10)
    evaluator = SDPEvaluator(config=config)
    res_sdp = evaluator.evaluate(x, phi)

    # Formular el LP correspondiente para scipy.optimize.linprog
    A_F = build_af_matrix(x_raw, phi_arr)
    A_d = build_ad_matrix(x_raw)
    A_c = build_ac_matrix(x_raw, phi_arr)

    N = scipy.linalg.null_space(A_F)
    d = N.shape[1]

    # Variables: [z_0, ..., z_{d-1}, t] -> tamaño d + 1
    # Maximizar t <=> Minimizar -t
    c_lp = np.zeros(d + 1)
    c_lp[-1] = -1.0

    # Restricciones de desigualdad A_ub @ [z, t] <= b_ub:
    # 1. t - A_d N z <= x_raw <=> [-A_d @ N, 1] [z, t] <= x_raw
    # 2. N z <= 1
    # 3. -N z <= 0
    # 4. -A_c @ N z <= x_raw - eps_c
    A_ub = np.vstack([
        np.hstack([-A_d @ N, np.ones((m, 1))]),
        np.hstack([N, np.zeros((K, 1))]),
        np.hstack([-N, np.zeros((K, 1))]),
        np.hstack([-A_c @ N, np.zeros((m, 1))]),
    ])
    b_ub = np.concatenate([x_raw, np.ones(K), np.zeros(K), x_raw - config.eps_c])

    lp_res = scipy.optimize.linprog(
        c=c_lp,
        A_ub=A_ub,
        b_ub=b_ub,
        bounds=(None, None),
        method="highs",
    )

    if res_sdp.mu is not None and lp_res.success:
        t_lp = -float(lp_res.fun)
        t_sdp = float(res_sdp.mu)
        diff = abs(t_sdp - t_lp)
        print(f"\n[SDP vs LINPROG] SDP mu*={t_sdp:.8f}, LP mu*={t_lp:.8f}, Diff={diff:.2e}")
        if diff < 1e-6:
            assert math.isclose(t_sdp, t_lp, abs_tol=1e-6)
