"""Evaluador de plantillas mediante optimización semidefinida convexa (SDP).

Referencias Matemáticas:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección IV (convex SDP para P dados x y phi).
    - Plan de Diseño del Algoritmo Genético, Sección 4 (Evaluación de Aptitud SDP Interno):
        Sección 4.1: Fase I elástica (siempre factible, acotada en [0, 1]).
        Sección 4.2: Fase II (eliminación de F vía espacio nulo A_F p = 0 => p = N z).
        Sección 4.3: Función de fitness escalar f.
"""

from __future__ import annotations

import math
import time
from typing import Any

import clarabel  # type: ignore[import-untyped]
import cvxpy as cp
import numpy as np
import scipy.linalg  # type: ignore[import-untyped]

from kuramoto_templates.domain.conditions import (
    build_ac_matrix,
    build_ad_matrix,
    build_af_matrix,
    build_m_basis,
    get_pair_indices,
    p_vector_to_matrix,
    verify_solution,
)
from kuramoto_templates.evaluators.base import Evaluator
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.results import EvalResult


class SDPEvaluator(Evaluator):
    """Evaluador de plantillas mediante programación semidefinida (SDP).

    Resuelve la evaluación de aptitud en dos fases:
      - Fase I: Problema elástico de factibilidad (minimiza residuo tau y déficit delta).
      - Fase II: Maximización de conectividad mu*(eps) sobre el espacio nulo de equilibrio.
    """

    def __init__(self, config: EvalConfig | None = None) -> None:
        """Inicializa el evaluador con su configuración de tolerancias y parámetros.

        Args:
            config: Instancia de EvalConfig (si es None, se usan los valores por defecto).
        """
        self.config = config if config is not None else EvalConfig()

    def evaluate_phase1(
        self, x: tuple[float, ...], phi: tuple[float, ...]
    ) -> EvalResult:
        """Ejecuta exclusivamente la Fase I elástica sobre un candidato (x, phi).

        Args:
            x: Tamaños de cluster x_a > 0, sum x_a = 1.
            phi: Fases de cluster phi_a en (-pi, pi].

        Returns:
            EvalResult con la métrica Phi en [0, 1] y diagnósticos del solver.
        """
        return self._solve_phase1(x, phi)

    def evaluate_phase2(
        self,
        x: tuple[float, ...],
        phi: tuple[float, ...],
        eps: float | None = None,
    ) -> EvalResult:
        """Resuelve directamente la Fase II (maximización de mu*) sobre el espacio nulo de A_F.

        Args:
            x: Tamaños de cluster x_a > 0, sum x_a = 1.
            phi: Fases de cluster phi_a en (-pi, pi].
            eps: Margen espectral a exigir (si es None se usa self.config.eps; permite eps=0.0).

        Returns:
            EvalResult con mu*, matriz P y diagnósticos de verify_solution.
        """
        eps_target = self.config.eps if eps is None else eps
        return self._solve_phase2_core(
            x, phi, eps=eps_target, phi_phase1=0.0, base_solver_info={}
        )

    def evaluate(self, x: tuple[float, ...], phi: tuple[float, ...]) -> EvalResult:
        """Evalúa un candidato (x, phi) resolviendo Fase I y, si es factible, Fase II.

        Args:
            x: Tamaños de cluster x_a > 0, sum x_a = 1.
            phi: Fases de cluster phi_a en (-pi, pi].

        Returns:
            EvalResult con fitness escalar f, conectividad mu* y diagnósticos completos.
        """
        # 1. Ejecutar Fase I elástica
        res_phase1 = self._solve_phase1(x, phi)
        Phi = res_phase1.Phi

        # Si Fase I no es factible bajo el umbral configurado, retornar fitness penalizado
        if Phi > self.config.phase1_tol or res_phase1.status not in (
            "OPTIMAL",
            "OPTIMAL_INACCURATE",
        ):
            f_penalized = -math.log10(1.0 + max(Phi, 1e-6) / 1e-6)
            return EvalResult(
                f=f_penalized,
                mu=None,
                Phi=Phi,
                P=None,
                status=res_phase1.status,
                solver_info=res_phase1.solver_info,
            )

        # 2. Ejecutar Fase II: Maximización de conectividad mu*
        return self._solve_phase2(x, phi, res_phase1)

    def _get_solver_kwargs(self) -> dict[str, Any]:
        """Configura los argumentos y tolerancias según el solver seleccionado."""
        tol = self.config.solver_tol
        solver_name = self.config.solver_name.upper()

        if solver_name == "CLARABEL":
            return {
                "solver": cp.CLARABEL,
                "tol_gap_abs": tol,
                "tol_gap_rel": tol,
                "tol_feas": tol,
                "tol_infeas_abs": tol,
                "tol_infeas_rel": tol,
                "verbose": False,
            }
        elif solver_name == "SCS":
            return {
                "solver": cp.SCS,
                "eps_abs": tol,
                "eps_rel": tol,
                "verbose": False,
            }
        else:
            return {"solver": solver_name, "verbose": False}

    def _solve_phase1(
        self, x: tuple[float, ...], phi: tuple[float, ...]
    ) -> EvalResult:
        """Formula y resuelve la Fase I elástica en cvxpy."""
        x_arr = np.asarray(x, dtype=float)
        phi_arr = np.asarray(phi, dtype=float)
        m = len(x_arr)
        pairs = get_pair_indices(m)
        K = len(pairs)

        A_F = build_af_matrix(x_arr, phi_arr)
        A_c = build_ac_matrix(x_arr, phi_arr)
        Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
        M_basis = build_m_basis(x_arr, phi_arr, Q)

        eps = self.config.eps
        eps_c = self.config.eps_c
        theta_tol = self.config.theta_tol

        # Variables de optimización
        p = cp.Variable(K)
        tau = cp.Variable(nonneg=True)
        delta = cp.Variable(nonneg=True)

        M_expr = sum(p[k] * M_basis[k] for k in range(K))
        obj = cp.Minimize(tau / (eps * theta_tol) + delta / eps)

        constraints = [
            A_F @ p <= tau,
            -A_F @ p <= tau,
            M_expr >> (eps - delta) * np.eye(m - 1),
            x_arr + A_c @ p >= eps_c,
            p >= 0.0,
            p <= 1.0,
        ]

        prob = cp.Problem(obj, constraints)
        solver_kwargs = self._get_solver_kwargs()

        t0 = time.perf_counter()
        try:
            prob.solve(**solver_kwargs)  # type: ignore[no-untyped-call]
            solve_time = time.perf_counter() - t0
            raw_status = prob.status
        except Exception as e:
            solve_time = time.perf_counter() - t0
            return EvalResult(
                f=-math.log10(1.0 + 1.0 / 1e-6),
                mu=None,
                Phi=1.0,
                P=None,
                status="SOLVER_ERROR",
                solver_info={
                    "error": str(e),
                    "solve_time_phase1": solve_time,
                    "solver": self.config.solver_name,
                    "cvxpy_version": getattr(cp, "__version__", "unknown"),
                    "solver_version": getattr(clarabel, "__version__", "unknown"),
                },
            )

        status_norm = raw_status.upper()
        is_optimal = status_norm in ("OPTIMAL", "OPTIMAL_INACCURATE")

        Phi_val = float(prob.value) if (is_optimal and prob.value is not None) else 1.0
        # Asegurar cota matemática numérica
        Phi_val = max(0.0, min(1.0, Phi_val))

        f_val = -math.log10(1.0 + max(Phi_val, 1e-6) / 1e-6)

        P_opt: tuple[tuple[float, ...], ...] | None = None
        if is_optimal and p.value is not None:
            P_mat = p_vector_to_matrix(np.clip(p.value, 0.0, 1.0), m)
            P_opt = tuple(tuple(float(v) for v in row) for row in P_mat)

        solver_info: dict[str, float | str | int] = {
            "cvxpy_version": getattr(cp, "__version__", "unknown"),
            "solver": self.config.solver_name,
            "solver_version": getattr(clarabel, "__version__", "unknown"),
            "solve_time_phase1": solve_time,
            "status_phase1": status_norm,
            "optimal_inaccurate": int(status_norm == "OPTIMAL_INACCURATE"),
            "tau": float(tau.value) if tau.value is not None else -1.0,
            "delta": float(delta.value) if delta.value is not None else -1.0,
        }

        return EvalResult(
            f=f_val,
            mu=None,
            Phi=Phi_val,
            P=P_opt,
            status=status_norm,
            solver_info=solver_info,
        )

    def _solve_phase2(
        self,
        x: tuple[float, ...],
        phi: tuple[float, ...],
        res_phase1: EvalResult,
    ) -> EvalResult:
        """Invoca la Fase II usando la configuración estándar y metadatos de Fase I."""
        return self._solve_phase2_core(
            x,
            phi,
            eps=self.config.eps,
            phi_phase1=res_phase1.Phi,
            base_solver_info=res_phase1.solver_info,
        )

    def _solve_phase2_core(
        self,
        x: tuple[float, ...],
        phi: tuple[float, ...],
        eps: float,
        phi_phase1: float,
        base_solver_info: dict[str, float | str | int],
    ) -> EvalResult:
        """Formula y resuelve la Fase II en cvxpy eliminando F mediante N = null_space(A_F)."""
        x_arr = np.asarray(x, dtype=float)
        phi_arr = np.asarray(phi, dtype=float)
        m = len(x_arr)
        pairs = get_pair_indices(m)
        K = len(pairs)

        A_F = build_af_matrix(x_arr, phi_arr)
        A_c = build_ac_matrix(x_arr, phi_arr)
        A_d = build_ad_matrix(x_arr)
        Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
        M_basis = build_m_basis(x_arr, phi_arr, Q)

        eps_c = self.config.eps_c

        N = scipy.linalg.null_space(A_F)
        d = N.shape[1]
        if d == 0:
            f_penalized = -math.log10(1.0 + max(phi_phase1, 1e-6) / 1e-6)
            return EvalResult(
                f=f_penalized,
                mu=None,
                Phi=phi_phase1,
                P=None,
                status="NO_NULLSPACE_DEGREES_OF_FREEDOM",
                solver_info=base_solver_info,
            )

        z = cp.Variable(d)
        t_var = cp.Variable()
        p_z = N @ z

        M_z_expr = sum(p_z[k] * M_basis[k] for k in range(K))

        constraints = [
            x_arr + A_d @ p_z >= t_var,
            p_z >= 0.0,
            p_z <= 1.0,
            M_z_expr >> eps * np.eye(m - 1),
            x_arr + A_c @ p_z >= eps_c,
        ]

        prob = cp.Problem(cp.Maximize(t_var), constraints)
        solver_kwargs = self._get_solver_kwargs()

        t0 = time.perf_counter()
        try:
            prob.solve(**solver_kwargs)  # type: ignore[no-untyped-call]
            solve_time = time.perf_counter() - t0
            raw_status = prob.status
        except Exception as e:
            solve_time = time.perf_counter() - t0
            f_penalized = -math.log10(1.0 + max(phi_phase1, 1e-6) / 1e-6)
            info = dict(base_solver_info)
            info.update({"error_phase2": str(e), "solve_time_phase2": solve_time})
            return EvalResult(
                f=f_penalized,
                mu=None,
                Phi=phi_phase1,
                P=None,
                status="SOLVER_ERROR_PHASE2",
                solver_info=info,
            )

        status_norm = raw_status.upper()
        is_optimal = status_norm in ("OPTIMAL", "OPTIMAL_INACCURATE")

        solver_info = dict(base_solver_info)
        solver_info.update(
            {
                "solve_time_phase2": solve_time,
                "status_phase2": status_norm,
                "optimal_inaccurate_phase2": int(status_norm == "OPTIMAL_INACCURATE"),
                "nullspace_dim": d,
            }
        )

        if not is_optimal or t_var.value is None or z.value is None:
            f_penalized = -math.log10(1.0 + max(phi_phase1, 1e-6) / 1e-6)
            return EvalResult(
                f=f_penalized,
                mu=None,
                Phi=phi_phase1,
                P=None,
                status=status_norm,
                solver_info=solver_info,
            )

        # Extraer solución y reconstruir matriz P
        p_opt_vec = np.clip(N @ z.value, 0.0, 1.0)
        P_mat = p_vector_to_matrix(p_opt_vec, m)
        P_opt = tuple(tuple(float(v) for v in row) for row in P_mat)

        # Verificación independiente en numpy puro antes de reportar factible
        verif = verify_solution(
            x_arr,
            phi_arr,
            P_mat,
            eps=eps,
            eps_c=self.config.eps_c,
            theta_tol=self.config.theta_tol,
        )

        solver_info.update(
            {
                "verified": int(verif["is_valid"]),
                "verif_max_F": float(verif["max_F"]),
                "verif_min_c": float(verif["min_c"]),
                "verif_min_eig_M": float(verif["min_eig_M"]),
                "verif_mu": float(verif["mu"]),
            }
        )

        if not verif["is_valid"]:
            f_penalized = -math.log10(1.0 + max(phi_phase1, 1e-6) / 1e-6)
            solver_info["verification_reasons"] = "; ".join(verif["reasons"])
            return EvalResult(
                f=f_penalized,
                mu=None,
                Phi=phi_phase1,
                P=P_opt,
                status="VERIFICATION_FAILED",
                solver_info=solver_info,
            )

        mu_star = float(verif["mu"])
        return EvalResult(
            f=mu_star,
            mu=mu_star,
            Phi=phi_phase1,
            P=P_opt,
            status=status_norm,
            solver_info=solver_info,
        )

