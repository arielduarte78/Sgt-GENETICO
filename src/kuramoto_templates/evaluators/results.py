"""Estructura de resultados de evaluación de aptitud y factibilidad.

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección IV (Splitting the Four-Cluster State).
    - Plan de Diseño del Algoritmo Genético, Sección 4 (Evaluación de Aptitud SDP Interno):
        - Fase I (elástica): minimiza residuo de equilibrio tau y déficit delta,
          obteniendo Phi(x, phi) en [0, 1].
        - Fitness escalar f: mu*(epsilon) si factible,
          -log10(1 + max(Phi, 1e-6)/1e-6) si infactible.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class EvalResult:
    """Resultado inmutable de la evaluación de un candidato (x, phi).

    Atributos:
        f: Valor de fitness escalar para el GA (maximizar).
           f > 0 para individuos factibles (mu*(epsilon)),
           f in [-6, 0) para individuos infactibles con gradiente logarítmico.
        mu: Conectividad mínima alcanzada mu = min_a d_a (None si infactible).
        Phi: Métrica elástica de infactibilidad en [0, 1] de la Fase I (0 = factible).
        P: Matriz de densidades de enlace óptima encontrada por el solver (None si no factible).
        status: Estado de convergencia del solver ("OPTIMAL", "INFEASIBLE", "NUMERICAL_ERROR").
        solver_info: Diccionario con métricas de diagnóstico (tiempo, iteraciones, residuos).
    """

    f: float
    mu: float | None
    Phi: float
    P: tuple[tuple[float, ...], ...] | None = None
    status: str = "PENDING"
    solver_info: dict[str, float | str | int] = field(default_factory=dict)
