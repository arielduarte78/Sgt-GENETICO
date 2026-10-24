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
        solver_tol: Tolerancia absoluta y relativa configurada para el solver (default 1e-8).
        solver_name: Nombre del solver de cvxpy a utilizar ("CLARABEL" o "SCS").
    """

    eps: float = 1e-6
    eps_c: float = 1e-4
    theta_tol: float = 1e-4
    phase1_tol: float = 1e-6
    solver_tol: float = 1e-8
    solver_name: str = "CLARABEL"
    mu_t: float = 0.6929077
    rho: float = 1.0
    s_tol: float = 1e-6
    lam_floor: float = 1e-7
    # lam_floor_hat: Piso calibrado para evaluador escalado (default 1e-5).
    # Justificación de auditoría: en 200 individuos, los desacuerdos > 1e-4 provienen exclusivamente
    # de fallos de convergencia de SCS (OPTIMAL_INACCURATE). El desacuerdo Clarabel(1e-8) vs
    # Clarabel(1e-9) es <= 7.73e-8 y con SCS convergido es <= 5.04e-6. El valor 1e-5 representa
    # >100x el desacuerdo de Clarabel y >2x el de SCS convergido, preservando S7 (~3.2e-4).
    lam_floor_hat: float = 1e-5

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
        if self.mu_t <= 0.0:
            raise ValueError(
                f"mu_t debe ser estrictamente positivo, recibido {self.mu_t}."
            )
        if self.rho <= 0.0:
            raise ValueError(
                f"rho debe ser estrictamente positivo, recibido {self.rho}."
            )
        if self.s_tol <= 0.0:
            raise ValueError(
                f"s_tol debe ser estrictamente positivo, recibido {self.s_tol}."
            )
        if self.lam_floor <= 0.0:
            raise ValueError(
                f"lam_floor debe ser estrictamente positivo, recibido {self.lam_floor}."
            )
        if self.lam_floor_hat <= 0.0:
            raise ValueError(
                f"lam_floor_hat debe ser estrictamente positivo, recibido {self.lam_floor_hat}."
            )
        valid_solvers = {"CLARABEL", "SCS"}
        if self.solver_name.upper() not in valid_solvers:
            raise ValueError(
                f"Solver '{self.solver_name}' no soportado. Solvers válidos: {valid_solvers}."
            )
