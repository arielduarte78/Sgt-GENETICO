"""Interfaz base para evaluadores de plantillas de Kuramoto.

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección IV (convex SDP para P dados x y phi).
    - Plan de Diseño del Algoritmo Genético, Sección 4 (Evaluación de Aptitud SDP Interno):
        Sección 4.1: Fase I elástica (siempre factible, acotada en [0, 1]).
        Sección 4.2: Fase II (maximización de conectividad mu*(epsilon) eliminando
                     el equilibrio por álgebra A_F p = 0 => p = N z vía SVD).
        Sección 4.3: Función de fitness escalar f(g).

NOTA DE ESQUELETO:
    PROHIBIDO escribir lógica SDP en esta etapa. Este módulo define exclusivamente
    el contrato abstracto del evaluador.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from kuramoto_templates.evaluators.results import EvalResult


class Evaluator(ABC):
    """Clase base abstracta para evaluadores de plantillas (x, phi)."""

    @abstractmethod
    def evaluate(self, x: tuple[float, ...], phi: tuple[float, ...]) -> EvalResult:
        """Evalúa un par (x, phi) resolviendo el problema de aptitud.

        Args:
            x: Tupla con tamaños de cluster x_a > 0, sum x_a = 1.
            phi: Tupla con fases de cluster phi_a en (-pi, pi].

        Returns:
            Instancia de EvalResult con fitness escalar f, conectividad mu y diagnósticos.
        """
        pass
