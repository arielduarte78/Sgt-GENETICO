"""Evaluador SDP Escalado para plantillas de clusters (ScaledSDPEvaluator).

Formulación Matemática y de Rescaling:
    1. Agrupamiento en hojas (sheets):
       - Si existe genoma: por anclas discretas k[a] en {0, ..., K_ref - 1}.
       - Si solo hay (x, phi): respaldo por cercanía a múltiplos de 2*pi/K_ref
         con tolerancia sheet_tol.
       - Hoja de referencia j_ref: hoja de mayor masa sum_{a in C_j} x_a.
    2. Escala efectiva h_eff:
       - Raíz del promedio cuadrático (ponderado por x) de los desvíos respecto
         a la media circular de cada hoja.
       - Rango admisible: [1e-4, 0.1]. Fuera de rango -> status NOT_SCALABLE
         y fallback a SDPEvaluator plano.
    3. Base reducida T (m x (m-1)):
       - Columnas eta: indicadores de las n_s - 1 hojas no referencia.
       - Columnas zeta: para hojas con n_j >= 2, base ortonormal (QR) de suma
         cero en la hoja, multiplicada por h_eff.
       - Cumple rank([T, 1_m]) == m y span(T) complementa 1_m.
    4. Operador escalado:
       - G_hat(p) = T^T L_w(p) T / h_eff^2, lineal en p con matrices base G_k.
       - Por la Ley de Inercia de Sylvester, G_hat(p) y M(p) = Q^T L_w(p) Q
         tienen exactamente la misma inercia espectral.
"""

from __future__ import annotations

import math
import time
from typing import Any

import cvxpy as cp
import mpmath  # type: ignore[import-untyped]
import numpy as np
import scipy.linalg  # type: ignore[import-untyped]

from kuramoto_templates.domain.conditions import (
    build_ac_matrix,
    build_ad_matrix,
    build_af_matrix,
    build_lw_matrix,
    get_pair_indices,
    p_matrix_to_vector,
    p_vector_to_matrix,
)
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.fitness import compute_continuous_fitness
from kuramoto_templates.evaluators.results import EvalResult
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def compute_sheets(
    x: tuple[float, ...] | np.ndarray,
    phi: tuple[float, ...] | np.ndarray,
    genome: Any | None = None,
    K_ref: int = 4,
    sheet_tol: float = 0.2,
) -> tuple[dict[int, list[int]], int]:
    """Agrupa los clusters en hojas según el genoma (anclas k) o por cercanía angular (respaldo).

    Args:
        x: Vector de masas x_a > 0.
        phi: Vector de fases phi_a en (-pi, pi].
        genome: Genoma opcional del individuo (si existe, sus anclas determinan las hojas).
        K_ref: Número de anclas del retículo de referencia.
        sheet_tol: Tolerancia angular en radianes para el respaldo por cercanía.

    Returns:
        Tupla (sheets, ref_sheet) donde sheets mapea id_hoja -> lista de índices de cluster,
        y ref_sheet es el id de la hoja de mayor masa total.
    """
    x_arr = np.asarray(x, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)
    m = len(x_arr)

    sheets: dict[int, list[int]] = {}

    if genome is not None and hasattr(genome, "k") and len(genome.k) == m:
        for idx, anchor in enumerate(genome.k):
            sheets.setdefault(int(anchor), []).append(idx)
    else:
        two_pi = 2.0 * math.pi
        for idx, ph in enumerate(phi_arr):
            best_k = 0
            best_dist = float("inf")
            for k in range(K_ref):
                target_ang = (two_pi * k / float(K_ref) + math.pi) % two_pi - math.pi
                diff = abs((ph - target_ang + math.pi) % two_pi - math.pi)
                if diff < best_dist:
                    best_dist = diff
                    best_k = k

            if best_dist <= sheet_tol:
                sheets.setdefault(best_k, []).append(idx)
            else:
                # Si excede la tolerancia, asigna a un cluster singleton separado
                sheets.setdefault(1000 + idx, []).append(idx)

    # Identificar la hoja de mayor masa total (hoja de referencia)
    masses = {s: sum(x_arr[i] for i in idxs) for s, idxs in sheets.items()}
    ref_sheet = max(masses, key=lambda s: (masses[s], -s))
    return sheets, ref_sheet


def compute_h_eff(
    x: tuple[float, ...] | np.ndarray,
    phi: tuple[float, ...] | np.ndarray,
    sheets: dict[int, list[int]],
) -> float:
    """Calcula la escala efectiva h_eff como RMS ponderado de desvíos respecto a la media circular.

    Args:
        x: Vector de masas x_a.
        phi: Vector de fases phi_a.
        sheets: Diccionario de hojas con sus índices de cluster.

    Returns:
        Escala h_eff >= 0.0.
    """
    x_arr = np.asarray(x, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)
    two_pi = 2.0 * math.pi

    sum_sq_devs = 0.0
    for idxs in sheets.values():
        if len(idxs) == 0:
            continue
        sin_sum = sum(x_arr[i] * math.sin(phi_arr[i]) for i in idxs)
        cos_sum = sum(x_arr[i] * math.cos(phi_arr[i]) for i in idxs)
        mean_phi = math.atan2(sin_sum, cos_sum)

        for i in idxs:
            delta = (phi_arr[i] - mean_phi + math.pi) % two_pi - math.pi
            sum_sq_devs += float(x_arr[i]) * (delta**2)

    return float(math.sqrt(max(0.0, sum_sq_devs)))


def build_scaled_basis_t(
    x: tuple[float, ...] | np.ndarray,
    phi: tuple[float, ...] | np.ndarray,
    sheets: dict[int, list[int]],
    ref_sheet: int,
    h_eff: float,
) -> np.ndarray | None:
    """Construye la base reducida T de dimensiones m x (m-1).

    Args:
        x: Vector de masas.
        phi: Vector de fases.
        sheets: Agrupamiento de clusters por hoja.
        ref_sheet: Hoja de referencia no incluida en columnas eta.
        h_eff: Escala efectiva h_eff.

    Returns:
        Matriz numpy m x (m-1) con rank([T, 1_m]) == m, o None si no es reducible a base válida.
    """
    m = len(x)
    eta_cols: list[np.ndarray] = []

    # 1. Columnas eta: indicadores de hojas no referencia
    for s, idxs in sorted(sheets.items(), key=lambda it: it[0]):
        if s == ref_sheet:
            continue
        col = np.zeros(m, dtype=float)
        col[idxs] = 1.0
        eta_cols.append(col)

    # 2. Columnas zeta: base ortonormal de suma cero en hojas con n_j >= 2, multiplicadas por h_eff
    zeta_cols: list[np.ndarray] = []
    for _s, idxs in sorted(sheets.items(), key=lambda it: it[0]):
        n_j = len(idxs)
        if n_j >= 2:
            # Base del espacio nulo de ones(1, n_j)
            Z = scipy.linalg.null_space(np.ones((1, n_j), dtype=float))
            for k in range(Z.shape[1]):
                col_z = Z[:, k]
                col = np.zeros(m, dtype=float)
                col[idxs] = col_z * h_eff
                zeta_cols.append(col)

    all_cols = eta_cols + zeta_cols
    if len(all_cols) != m - 1:
        return None

    T = np.column_stack(all_cols)
    # Validar rango completo con 1_m
    aug = np.column_stack([T, np.ones(m, dtype=float)])
    if np.linalg.matrix_rank(aug) < m:
        return None

    return T


def build_g_hat_basis(
    x: np.ndarray,
    phi: np.ndarray,
    T: np.ndarray,
    h_eff: float,
) -> list[np.ndarray]:
    """Precalcula las matrices base simétricas G_k tal que G_hat(p) = sum_k p[k] * G_k.

    Args:
        x: Vector de masas (longitud m).
        phi: Vector de fases (longitud m).
        T: Matriz base m x (m-1).
        h_eff: Escala efectiva.

    Returns:
        Lista de K matrices simétricas (m-1) x (m-1).
    """
    m = len(x)
    pairs = get_pair_indices(m)
    h_eff_sq = h_eff**2

    G_basis: list[np.ndarray] = []
    for a, b in pairs:
        d_ab = np.zeros(m, dtype=float)
        d_ab[a] = 1.0
        d_ab[b] = -1.0
        Td = T.T @ d_ab
        weight = (x[a] * x[b] * math.cos(phi[a] - phi[b])) / h_eff_sq
        Gk = weight * np.outer(Td, Td)
        Gk = 0.5 * (Gk + Gk.T)
        G_basis.append(Gk)

    return G_basis


class ScaledSDPEvaluator:
    """Evaluador SDP de plantillas con reescalado analítico en base T y operador G_hat."""

    def __init__(
        self,
        config: EvalConfig | None = None,
        lam_floor_hat: float | None = None,
        h_eff_min: float | None = None,
    ) -> None:
        self.config = config or EvalConfig()
        self.plain_evaluator = SDPEvaluator(self.config)
        # lam_floor_hat calibrado: por defecto self.config.lam_floor_hat (1e-5) o especificado
        self.lam_floor_hat = (
            float(lam_floor_hat)
            if lam_floor_hat is not None
            else getattr(self.config, "lam_floor_hat", 1e-5)
        )
        self.h_eff_min = (
            float(h_eff_min)
            if h_eff_min is not None
            else getattr(self.config, "h_eff_min", 1e-4)
        )

    def evaluate_phase1(
        self,
        x: tuple[float, ...],
        phi: tuple[float, ...],
        genome: Any | None = None,
    ) -> EvalResult:
        """Formula y resuelve la Fase I escalada sobre G_hat(p).

        Si h_eff cae fuera de [1e-4, 0.1], retorna status NOT_SCALABLE
        y calcula f con evaluador plano.
        """
        x_arr = np.asarray(x, dtype=float)
        phi_arr = np.asarray(phi, dtype=float)
        m = len(x_arr)
        pairs = get_pair_indices(m)
        K = len(pairs)

        sheets, ref_sheet = compute_sheets(
            x_arr, phi_arr, genome=genome, K_ref=4, sheet_tol=0.2
        )
        h_eff = compute_h_eff(x_arr, phi_arr, sheets)

        # Validación de rango admisible de escala
        if h_eff < self.h_eff_min or h_eff > 0.1:
            res_plain = self.plain_evaluator.evaluate_phase1(x, phi)
            s_plain = float(res_plain.s) if res_plain.s is not None else 1.0
            lam_plain = float(res_plain.lam) if res_plain.lam is not None else float("-inf")
            is_stable_macro = (
                res_plain.status in ("OPTIMAL", "OPTIMAL_INACCURATE")
                and lam_plain >= self.config.lam_floor
                and s_plain <= self.config.s_tol
            )
            if is_stable_macro:
                f_val = lam_plain
                lam_scaled = lam_plain
                lam_norm_scaled = (
                    float(res_plain.lam_normalized)
                    if res_plain.lam_normalized is not None
                    else lam_plain
                )
            else:
                lam_scaled = min(lam_plain / 1e-8, -10.0)
                norm_base = (
                    float(res_plain.lam_normalized)
                    if res_plain.lam_normalized is not None
                    else lam_plain
                )
                lam_norm_scaled = min(norm_base / 1e-8, -10.0) - s_plain
                f_val = min(lam_plain / 1e-8 - self.lam_floor_hat, -10.0) - s_plain

            return EvalResult(
                f=f_val,
                mu=None,
                Phi=res_plain.Phi,
                lam=lam_scaled,
                s=s_plain,
                lam_normalized=lam_norm_scaled,
                P=res_plain.P,
                status="NOT_SCALABLE",
                solver_info={
                    "h_eff": h_eff,
                    "reason": "h_eff outside [1e-4, 0.1]",
                    "plain_status": res_plain.status,
                    "lam_plain": lam_plain,
                },
            )

        T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, h_eff)
        if T is None:
            res_plain = self.plain_evaluator.evaluate_phase1(x, phi)
            s_plain = float(res_plain.s) if res_plain.s is not None else 1.0
            lam_plain = float(res_plain.lam) if res_plain.lam is not None else float("-inf")
            is_stable_macro = (
                res_plain.status in ("OPTIMAL", "OPTIMAL_INACCURATE")
                and lam_plain >= self.config.lam_floor
                and s_plain <= self.config.s_tol
            )
            if is_stable_macro:
                f_val = lam_plain
                lam_scaled = lam_plain
                lam_norm_scaled = (
                    float(res_plain.lam_normalized)
                    if res_plain.lam_normalized is not None
                    else lam_plain
                )
            else:
                lam_scaled = min(lam_plain / 1e-8, -10.0)
                norm_base = (
                    float(res_plain.lam_normalized)
                    if res_plain.lam_normalized is not None
                    else lam_plain
                )
                lam_norm_scaled = min(norm_base / 1e-8, -10.0) - s_plain
                f_val = min(lam_plain / 1e-8 - self.lam_floor_hat, -10.0) - s_plain

            return EvalResult(
                f=f_val,
                mu=None,
                Phi=res_plain.Phi,
                lam=lam_scaled,
                s=s_plain,
                lam_normalized=lam_norm_scaled,
                P=res_plain.P,
                status="NOT_SCALABLE",
                solver_info={
                    "h_eff": h_eff,
                    "reason": "Could not form full rank basis T",
                    "plain_status": res_plain.status,
                    "lam_plain": lam_plain,
                },
            )

        G_basis = build_g_hat_basis(x_arr, phi_arr, T, h_eff)

        A_F = build_af_matrix(x_arr, phi_arr)
        A_c = build_ac_matrix(x_arr, phi_arr)
        A_d = build_ad_matrix(x_arr)

        N = scipy.linalg.null_space(A_F)
        d = N.shape[1]
        if d == 0:
            s_default = max(0.0, float(self.config.mu_t - np.min(x_arr)))
            f_val = compute_continuous_fitness(
                lam=0.0,
                s=s_default,
                tol=self.config.s_tol,
                lam_floor=self.lam_floor_hat,
            )
            return EvalResult(
                f=f_val,
                mu=None,
                Phi=s_default,
                lam=0.0,
                s=s_default,
                lam_normalized=0.0,
                P=None,
                status="NO_NULLSPACE_DEGREES_OF_FREEDOM",
                solver_info={"h_eff": h_eff, "res_inf": 0.0},
            )

        # Refinamiento iterativo para garantizar ||A_F N||_inf <= 1e-14
        res_inf = float(np.max(np.abs(A_F @ N)))
        if res_inf > 1e-14:
            corr, _, _, _ = np.linalg.lstsq(A_F, A_F @ N, rcond=None)
            N = N - corr
            N, _ = np.linalg.qr(N)
            res_inf = float(np.max(np.abs(A_F @ N)))

        z = cp.Variable(d)
        lam_var = cp.Variable()
        s_var = cp.Variable(nonneg=True)
        p_z = N @ z

        G_z_expr = sum(p_z[k] * G_basis[k] for k in range(K))
        obj = cp.Maximize(lam_var - self.config.rho * s_var)

        constraints = [
            p_z >= 0.0,
            p_z <= 1.0,
            G_z_expr >> lam_var * np.eye(m - 1),
            x_arr + A_c @ p_z >= self.config.eps_c,
            x_arr + A_d @ p_z >= self.config.mu_t - s_var,
            s_var >= 0.0,
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
            s_fail = 1.0
            return EvalResult(
                f=-s_fail,
                mu=None,
                Phi=s_fail,
                lam=None,
                s=s_fail,
                lam_normalized=None,
                P=None,
                status="SOLVER_ERROR",
                solver_info={
                    "error": str(e),
                    "solve_time_phase1": solve_time,
                    "h_eff": h_eff,
                },
            )

        status_norm = raw_status.upper()
        if (
            status_norm in ("OPTIMAL", "OPTIMAL_INACCURATE")
            and lam_var.value is not None
            and s_var.value is not None
            and z.value is not None
        ):
            lam_val = float(lam_var.value)
            s_val = float(s_var.value)

            f_val = compute_continuous_fitness(
                lam=lam_val,
                s=s_val,
                tol=self.config.s_tol,
                lam_floor=self.lam_floor_hat,
            )

            p_opt_vec = np.clip(N @ z.value, 0.0, 1.0)
            P_mat = p_vector_to_matrix(p_opt_vec, m)
            P_opt = tuple(tuple(float(v) for v in row) for row in P_mat)

            return EvalResult(
                f=f_val,
                mu=None,
                Phi=s_val,
                lam=lam_val,
                s=s_val,
                lam_normalized=lam_val,
                P=P_opt,
                status=status_norm,
                solver_info={
                    "solver": self.config.solver_name,
                    "solve_time_phase1": solve_time,
                    "status_phase1": status_norm,
                    "h_eff": h_eff,
                    "lam_hat": lam_val,
                    "s": s_val,
                    "af_res_inf": res_inf,
                    "lam_floor_hat": self.lam_floor_hat,
                },
            )
        else:
            s_fallback = 1.0
            return EvalResult(
                f=-s_fallback,
                mu=None,
                Phi=s_fallback,
                lam=None,
                s=s_fallback,
                lam_normalized=None,
                P=None,
                status=status_norm,
                solver_info={"h_eff": h_eff, "status_phase1": status_norm},
            )

    def evaluate_phase2(
        self,
        x: tuple[float, ...],
        phi: tuple[float, ...],
        eps: float | None = None,
        genome: Any | None = None,
    ) -> EvalResult:
        """Fase II escalada: max mu s.a. G_hat(p) >> lam_floor_hat*I, exacto por p = N z."""
        x_arr = np.asarray(x, dtype=float)
        phi_arr = np.asarray(phi, dtype=float)
        m = len(x_arr)
        pairs = get_pair_indices(m)
        K = len(pairs)

        sheets, ref_sheet = compute_sheets(
            x_arr, phi_arr, genome=genome, K_ref=4, sheet_tol=0.2
        )
        h_eff = compute_h_eff(x_arr, phi_arr, sheets)

        if h_eff < self.h_eff_min or h_eff > 0.1:
            return self.plain_evaluator.evaluate_phase2(x, phi, eps=eps)

        T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, h_eff)
        if T is None:
            return self.plain_evaluator.evaluate_phase2(x, phi, eps=eps)

        G_basis = build_g_hat_basis(x_arr, phi_arr, T, h_eff)
        A_F = build_af_matrix(x_arr, phi_arr)
        A_c = build_ac_matrix(x_arr, phi_arr)
        A_d = build_ad_matrix(x_arr)

        N = scipy.linalg.null_space(A_F)
        d = N.shape[1]
        if d == 0:
            return EvalResult(
                f=-1.0,
                mu=None,
                Phi=1.0,
                P=None,
                status="NO_NULLSPACE_DEGREES_OF_FREEDOM",
                solver_info={"h_eff": h_eff},
            )

        res_inf = float(np.max(np.abs(A_F @ N)))
        if res_inf > 1e-14:
            corr, _, _, _ = np.linalg.lstsq(A_F, A_F @ N, rcond=None)
            N = N - corr
            N, _ = np.linalg.qr(N)

        z = cp.Variable(d)
        t_var = cp.Variable()
        p_z = N @ z

        floor_margin = self.lam_floor_hat if eps is None else float(eps)
        G_z_expr = sum(p_z[k] * G_basis[k] for k in range(K))

        constraints = [
            x_arr + A_d @ p_z >= t_var,
            p_z >= 0.0,
            p_z <= 1.0,
            G_z_expr >> floor_margin * np.eye(m - 1),
            x_arr + A_c @ p_z >= self.config.eps_c,
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
            return EvalResult(
                f=-1.0,
                mu=None,
                Phi=1.0,
                P=None,
                status="SOLVER_ERROR_PHASE2",
                solver_info={"error": str(e), "solve_time_phase2": solve_time},
            )

        status_norm = raw_status.upper()
        if (
            status_norm in ("OPTIMAL", "OPTIMAL_INACCURATE")
            and t_var.value is not None
            and z.value is not None
        ):
            mu_val = float(t_var.value)
            p_opt_vec = np.clip(N @ z.value, 0.0, 1.0)
            P_mat = p_vector_to_matrix(p_opt_vec, m)
            P_opt = tuple(tuple(float(v) for v in row) for row in P_mat)

            return EvalResult(
                f=mu_val,
                mu=mu_val,
                Phi=0.0,
                P=P_opt,
                status=status_norm,
                solver_info={
                    "solver": self.config.solver_name,
                    "status_phase2": status_norm,
                    "solve_time_phase2": solve_time,
                    "mu": mu_val,
                    "h_eff": h_eff,
                },
            )
        else:
            return EvalResult(
                f=-1.0,
                mu=None,
                Phi=1.0,
                P=None,
                status=status_norm,
                solver_info={"h_eff": h_eff, "status_phase2": status_norm},
            )

    def evaluate(
        self,
        x: tuple[float, ...],
        phi: tuple[float, ...],
        genome: Any | None = None,
    ) -> EvalResult:
        """Evalúa un candidato ejecutando Fase I escalada y, si s* <= tol, Fase II escalada."""
        res_p1 = self.evaluate_phase1(x, phi, genome=genome)
        if (
            res_p1.s is None
            or res_p1.s > self.config.s_tol
            or res_p1.status not in ("OPTIMAL", "OPTIMAL_INACCURATE")
        ):
            return res_p1

        res_p2 = self.evaluate_phase2(x, phi, genome=genome)
        if res_p2.mu is not None:
            return EvalResult(
                f=res_p1.f,
                mu=res_p2.mu,
                Phi=res_p1.Phi,
                lam=res_p1.lam,
                s=res_p1.s,
                lam_normalized=res_p1.lam_normalized,
                P=res_p2.P,
                status=res_p2.status,
                solver_info=res_p2.solver_info,
            )
        return res_p2

    def verify_candidate_triplet(
        self,
        x: tuple[float, ...],
        phi: tuple[float, ...],
        P: tuple[tuple[float, ...], ...],
        genome: Any | None = None,
    ) -> dict[str, Any]:
        """Cuádruple verificación independiente (Antidoping Robusto):
        a) Re-evaluación con Clarabel (tol=1e-9) dando lam_hat >= lam_floor_hat.
        b) Verificación en NumPy float64 sobre P: max |F| <= 1e-12, min c_a >= eps_c,
           min d_a >= mu_t - s_tol, y eigvalsh(G_hat) >= lam_floor_hat.
        c) Verificación en mpmath a 50 dígitos (mp.eigsy) sobre G_hat(P) >= lam_floor_hat.
        d) Segundo solver SCS con max_iters=20000 (si da OPTIMAL_INACCURATE, el veredicto
           depende de Clarabel 1e-9 + NumPy + mpmath 50 dígitos).
        """
        x_arr = np.asarray(x, dtype=float)
        phi_arr = np.asarray(phi, dtype=float)
        m = len(x_arr)
        P_arr = np.asarray(P, dtype=float)
        p_vec = p_matrix_to_vector(P_arr)

        sheets, ref_sheet = compute_sheets(x_arr, phi_arr, genome=genome)
        h_eff = compute_h_eff(x_arr, phi_arr, sheets)
        T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, h_eff)

        L_w = build_lw_matrix(x_arr, phi_arr, p_vec)
        Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
        M = Q.T @ L_w @ Q

        # 1. NumPy float64 sobre P
        A_F = build_af_matrix(x_arr, phi_arr)
        max_abs_F = float(np.max(np.abs(A_F @ p_vec)))
        A_c = build_ac_matrix(x_arr, phi_arr)
        min_c_a = float(np.min(x_arr + A_c @ p_vec))
        A_d = build_ad_matrix(x_arr)
        min_d_a = float(np.min(x_arr + A_d @ p_vec))

        lam_min_M = float(np.min(np.linalg.eigvalsh(M)))
        if T is not None and h_eff >= self.h_eff_min:
            G_hat = (T.T @ L_w @ T) / (h_eff**2)
            lam_min_G = float(np.min(np.linalg.eigvalsh(G_hat)))
        else:
            G_hat = M
            lam_min_G = lam_min_M

        floor_target = (
            self.lam_floor_hat
            if (T is not None and h_eff >= self.h_eff_min)
            else self.config.lam_floor
        )

        numpy_ok = bool(
            max_abs_F <= 1e-12
            and min_c_a >= self.config.eps_c
            and min_d_a >= (self.config.mu_t - self.config.s_tol)
            and lam_min_G >= floor_target
        )

        # 2. mpmath 50 dígitos con mp.eigsy sobre G_hat
        mpmath.mp.dps = 50
        dim = G_hat.shape[0]
        M_mp = mpmath.mp.matrix(dim, dim)
        for r in range(dim):
            for c in range(dim):
                M_mp[r, c] = mpmath.mp.mpf(str(G_hat[r, c]))
        eigs_mp = mpmath.mp.eigsy(M_mp, eigvals_only=True)
        min_eig_mp = float(min(eigs_mp))
        mpmath_ok = bool(min_eig_mp >= floor_target)

        # 3. Re-evaluación con Clarabel (tol=1e-9)
        cl9_cfg = EvalConfig(
            mu_t=self.config.mu_t,
            solver_tol=1e-9,
            solver_name="CLARABEL",
            lam_floor=self.config.lam_floor,
            lam_floor_hat=self.lam_floor_hat,
            h_eff_min=self.h_eff_min,
        )
        cl9_ev = ScaledSDPEvaluator(
            config=cl9_cfg,
            lam_floor_hat=self.lam_floor_hat,
            h_eff_min=self.h_eff_min,
        )
        res_cl9 = cl9_ev.evaluate_phase1(x, phi, genome=genome)
        lam_cl9 = float(res_cl9.lam) if res_cl9.lam is not None else float("-inf")
        cl9_ok = bool(
            res_cl9.status in ("OPTIMAL", "OPTIMAL_INACCURATE")
            and lam_cl9 >= floor_target
            and res_cl9.s is not None
            and res_cl9.s <= self.config.s_tol
        )

        # 4. Solver SCS con max_iters=20000
        scs_cfg = EvalConfig(
            mu_t=self.config.mu_t,
            solver_tol=1e-8,
            solver_name="SCS",
            lam_floor=self.config.lam_floor,
            lam_floor_hat=self.lam_floor_hat,
            h_eff_min=self.h_eff_min,
        )
        scs_ev = ScaledSDPEvaluator(
            config=scs_cfg,
            lam_floor_hat=self.lam_floor_hat,
            h_eff_min=self.h_eff_min,
        )
        res_scs = scs_ev.evaluate_phase1(x, phi, genome=genome)
        lam_scs = float(res_scs.lam) if res_scs.lam is not None else float("-inf")

        if res_scs.status == "OPTIMAL":
            scs_ok = bool(lam_scs >= floor_target)
        elif res_scs.status == "OPTIMAL_INACCURATE":
            # Si da OPTIMAL_INACCURATE, veredicto depende de Cl(1e-9) + NumPy + mpmath(50d)
            scs_ok = True
        else:
            scs_ok = False

        is_valid = bool(cl9_ok and numpy_ok and mpmath_ok and scs_ok)

        return {
            "is_valid": is_valid,
            "lam_cl9": lam_cl9,
            "lam_np": lam_min_G,
            "lam_mp": min_eig_mp,
            "lam_scs": lam_scs,
            "lam_min_M": lam_min_M,
            "max_abs_F": max_abs_F,
            "min_c_a": min_c_a,
            "min_d_a": min_d_a,
            "status_scs": res_scs.status,
            "status_cl9": res_cl9.status,
            "cl9_ok": cl9_ok,
            "numpy_ok": numpy_ok,
            "mpmath_ok": mpmath_ok,
            "scs_ok": scs_ok,
        }

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
                "max_iters": 20000,
                "verbose": False,
            }
        else:
            return {"solver": solver_name, "verbose": False}
