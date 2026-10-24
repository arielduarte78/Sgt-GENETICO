"""Configuración tipada e inmutable de parámetros de evaluación SDP.

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección IV (convex SDP).
    - Regla de Diseño: Parámetros NUNCA hardcodeados dentro de funciones.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class EvalConfig:
    """Configuración inmutable para evaluadores de plantillas.

    Atributos:
        eps: Margen espectral mínimo para estabilidad de modos de cluster: M(p) >> eps * I.
        eps_c: Cota inferior para la fuerza de estabilidad intra-cluster: c_a(p) >= eps_c.
        theta_tol: Tolerancia máxima permitida en el residuo de equilibrio: |F_a(p)| <= theta_tol.
        phase1_tol: Umbral en Phi de la Fase I para considerar factible e iniciar Fase II.
        solver_tol: Tolerancia absoluta y relativa configurada para el solver (default 1e-10).
        solver_name: Nombre del solver de cvxpy a utilizar ("CLARABEL" o "SCS").
    """

    eps: float = 1e-6
    eps_c: float = 1e-4
    theta_tol: float = 1e-4
    phase1_tol: float = 1e-6
    solver_tol: float = 1e-10
    solver_name: str = "CLARABEL"

    def __post_init__(self) -> None:
        """Valida la positividad y coherencia de los parámetros."""
        if self.eps <= 0.0:
            raise ValueError(f"eps debe ser estrictamente positivo, recibido {self.eps}.")
        if self.eps_c <= 0.0:
            raise ValueError(f"eps_c debe ser estrictamente positivo, recibido {self.eps_c}.")
        if self.theta_tol <= 0.0:
            raise ValueError(
                f"theta_tol debe ser estrictamente positivo, recibido {self.theta_tol}."
            )
        if self.phase1_tol <= 0.0:
            raise ValueError(
                f"phase1_tol debe ser estrictamente positivo, recibido {self.phase1_tol}."
            )
        if self.solver_tol <= 0.0:
            raise ValueError(
                f"solver_tol debe ser estrictamente positivo, recibido {self.solver_tol}."
            )
        valid_solvers = {"CLARABEL", "SCS"}
        if self.solver_name.upper() not in valid_solvers:
            raise ValueError(
                f"Solver '{self.solver_name}' no soportado. Solvers válidos: {valid_solvers}."
            )
