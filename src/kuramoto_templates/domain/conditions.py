"""Condiciones de equilibrio y estabilidad para plantillas de clusters.

Referencias Matemáticas:
    O'Keeffe, K. P. (2026). "Beyond 11/16: dense Kuramoto networks that resist synchrony",
    arXiv:2610.21271v1, Secciones II y III:
      - Grado y Conectividad: d_a = sum_b P_ab x_b, mu(T) = min_a d_a (Sec. III).
      - Equilibrio de Fases: F_a = sum_b P_ab x_b sin(phi_a - phi_b) = 0 (Ec. 6).
      - Estabilidad intra-cluster: c_a = sum_b P_ab x_b cos(phi_a - phi_b) > 0 (Ec. 7).
      - Estabilidad inter-cluster: Q^T L_w(P) Q > 0, con w_ab = x_a x_b P_ab cos(phi_a - phi_b)
        y Q base ortonormal de 1_perp (Ec. 8 y Sec. II).

Todas las construcciones de matrices y operadores se implementan como funciones PURAS en numpy,
independientes de cualquier framework de optimización (cvxpy).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import scipy.linalg  # type: ignore[import-untyped]

from kuramoto_templates.domain.template import Template


def get_pair_indices(m: int) -> list[tuple[int, int]]:
    """Genera la lista canónica de pares ordenados (a, b) con 0 <= a < b < m.

    Args:
        m: Número de clusters (m >= 2).

    Returns:
        Lista de tuplas (a, b) con orden lexicográfico.
        Longitud K = m * (m - 1) // 2.
    """
    return [(a, b) for a in range(m) for b in range(a + 1, m)]


def build_af_matrix(x: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Construye el operador lineal A_F tal que A_F @ p = F (residuo de equilibrio).

    Referencia: O'Keeffe (2026), Ecuación (6):
        F_a(p) = sum_{b != a} x_b sin(phi_a - phi_b) P_ab

    Args:
        x: Vector 1D de tamaños de cluster x_a > 0 (longitud m).
        phi: Vector 1D de fases phi_a (longitud m).

    Returns:
        Matriz numpy de forma (m, K), donde K = m*(m-1)/2.
    """
    m = len(x)
    pairs = get_pair_indices(m)
    K = len(pairs)
    A_F = np.zeros((m, K), dtype=float)

    for k, (a, b) in enumerate(pairs):
        delta = phi[a] - phi[b]
        sin_val = math.sin(delta)
        A_F[a, k] = x[b] * sin_val
        A_F[b, k] = -x[a] * sin_val

    return A_F


def build_ac_matrix(x: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Construye el operador lineal A_c tal que c(p) = x + A_c @ p (estabilidad intra-cluster).

    Referencia: O'Keeffe (2026), Ecuación (7):
        c_a(p) = x_a * P_aa * cos(0) + sum_{b != a} x_b cos(phi_a - phi_b) P_ab
               = x_a + sum_{b != a} x_b cos(phi_a - phi_b) P_ab  (con P_aa = 1)

    Args:
        x: Vector 1D de tamaños de cluster (longitud m).
        phi: Vector 1D de fases (longitud m).

    Returns:
        Matriz numpy de forma (m, K).
    """
    m = len(x)
    pairs = get_pair_indices(m)
    K = len(pairs)
    A_c = np.zeros((m, K), dtype=float)

    for k, (a, b) in enumerate(pairs):
        cos_val = math.cos(phi[a] - phi[b])
        A_c[a, k] = x[b] * cos_val
        A_c[b, k] = x[a] * cos_val

    return A_c


def build_ad_matrix(x: np.ndarray) -> np.ndarray:
    """Construye el operador lineal A_d tal que d(p) = x + A_d @ p (grados de cluster).

    Referencia: O'Keeffe (2026), Sección III:
        d_a(p) = x_a * P_aa + sum_{b != a} x_b P_ab = x_a + sum_{b != a} x_b P_ab (con P_aa = 1)

    Args:
        x: Vector 1D de tamaños de cluster (longitud m).

    Returns:
        Matriz numpy de forma (m, K).
    """
    m = len(x)
    pairs = get_pair_indices(m)
    K = len(pairs)
    A_d = np.zeros((m, K), dtype=float)

    for k, (a, b) in enumerate(pairs):
        A_d[a, k] = x[b]
        A_d[b, k] = x[a]

    return A_d


def build_lw_matrix(x: np.ndarray, phi: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Construye el Laplaciano ponderado L_w(p) para los modos de cluster.

    Referencia: O'Keeffe (2026), Ecuación (8):
        L_w(p) = sum_{a < b} x_a x_b cos(phi_a - phi_b) P_ab (e_a - e_b)(e_a - e_b)^T

    Args:
        x: Vector de tamaños (longitud m).
        phi: Vector de fases (longitud m).
        p: Vector de densidades fuera de la diagonal (longitud K).

    Returns:
        Matriz simétrica m x m con L_w @ 1 = 0.
    """
    m = len(x)
    pairs = get_pair_indices(m)
    L_w = np.zeros((m, m), dtype=float)

    for k, (a, b) in enumerate(pairs):
        w_ab = x[a] * x[b] * math.cos(phi[a] - phi[b]) * float(p[k])
        L_w[a, a] += w_ab
        L_w[b, b] += w_ab
        L_w[a, b] -= w_ab
        L_w[b, a] -= w_ab

    return L_w


def build_m_basis(x: np.ndarray, phi: np.ndarray, Q: np.ndarray) -> list[np.ndarray]:
    """Construye la base matricial M_k tal que M(p) = sum_k p[k] * M_k.

    Referencia:
        M(p) = Q^T L_w(p) Q, con Q base ortonormal del complemento ortogonal de 1_m.
        Como L_w(p) = sum_k p[k] * B_k, cada M_k = Q^T B_k Q es una matriz constante (m-1) x (m-1).

    Args:
        x: Vector de tamaños (longitud m).
        phi: Vector de fases (longitud m).
        Q: Matriz ortonormal m x (m - 1) que satisface Q^T Q = I y Q^T 1 = 0.

    Returns:
        Lista de K matrices simétricas de forma (m - 1, m - 1).
    """
    m = len(x)
    pairs = get_pair_indices(m)
    M_basis: list[np.ndarray] = []

    for a, b in pairs:
        diff = np.zeros(m, dtype=float)
        diff[a] = 1.0
        diff[b] = -1.0
        q_diff = Q.T @ diff
        scale = x[a] * x[b] * math.cos(phi[a] - phi[b])
        Mk = scale * np.outer(q_diff, q_diff)
        Mk = np.asarray(0.5 * (Mk + Mk.T), dtype=float)
        M_basis.append(Mk)

    return M_basis


def build_m_matrix(x: np.ndarray, phi: np.ndarray, p: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Calcula la matriz proyectada M(p) = Q^T L_w(p) Q directamente.

    Args:
        x: Vector de tamaños (longitud m).
        phi: Vector de fases (longitud m).
        p: Vector de densidades (longitud K).
        Q: Matriz ortonormal m x (m - 1).

    Returns:
        Matriz simétrica (m - 1) x (m - 1).
    """
    L_w = build_lw_matrix(x, phi, p)
    M = Q.T @ L_w @ Q
    return np.asarray(0.5 * (M + M.T), dtype=float)


def p_matrix_to_vector(P: np.ndarray | tuple[tuple[float, ...], ...]) -> np.ndarray:
    """Extrae el vector p de densidades de pares (a < b) a partir de la matriz P.

    Args:
        P: Matriz m x m.

    Returns:
        Vector 1D de longitud K con entradas P_ab para a < b.
    """
    P_arr = np.asarray(P, dtype=float)
    m = P_arr.shape[0]
    pairs = get_pair_indices(m)
    return np.array([P_arr[a, b] for a, b in pairs], dtype=float)


def p_vector_to_matrix(p: np.ndarray, m: int) -> np.ndarray:
    """Reconstruye la matriz simétrica P m x m a partir de p, fijando diagonal P_aa = 1.

    Args:
        p: Vector de densidades (longitud K = m*(m-1)/2).
        m: Número de clusters.

    Returns:
        Matriz simétrica m x m con diagonal 1.0.
    """
    pairs = get_pair_indices(m)
    P = np.eye(m, dtype=float)
    for k, (a, b) in enumerate(pairs):
        val = float(p[k])
        P[a, b] = val
        P[b, a] = val
    return P


def degrees(template: Template) -> tuple[float, ...]:
    """Calcula los grados normalizados d_a de cada cluster.

    Referencia: O'Keeffe (2026), Sección III:
        d_a = sum_{b=1}^m P_ab x_b

    Args:
        template: Plantilla evaluada (debe tener matriz P definida).

    Returns:
        Tupla con los grados d_a para a = 1, ..., m.

    Raises:
        ValueError: Si la plantilla no contiene matriz P.
    """
    if template.P is None:
        raise ValueError("No se pueden calcular grados sin una matriz de densidades P.")

    x_arr = np.asarray(template.x, dtype=float)
    P_arr = np.asarray(template.P, dtype=float)
    d = P_arr @ x_arr
    return tuple(float(v) for v in d)


def connectivity(template: Template) -> float:
    """Calcula la conectividad mu(T) = min_a d_a de la plantilla.

    Referencia: O'Keeffe (2026), Sección III:
        mu(T) = min_a d_a

    Args:
        template: Plantilla evaluada (debe tener matriz P definida).

    Returns:
        Valor escalar de la conectividad mínima mu(T).

    Raises:
        ValueError: Si la plantilla no contiene matriz P.
    """
    if template.P is None:
        raise ValueError("No se puede calcular conectividad sin matriz P.")

    d = degrees(template)
    return float(min(d))


def equilibrium_residual(template: Template) -> tuple[float, ...]:
    """Calcula el residuo de equilibrio F_a para cada cluster.

    Referencia: O'Keeffe (2026), Sección III, Ecuación (6):
        F_a = sum_{b=1}^m P_ab x_b sin(phi_a - phi_b) = 0 para todo a

    Args:
        template: Plantilla con fases y densidades definidas.

    Returns:
        Tupla con los residuos (F_1, ..., F_m).

    Raises:
        ValueError: Si la plantilla no contiene matriz P.
    """
    if template.P is None:
        raise ValueError("No se puede calcular el equilibrio F_a sin matriz P.")

    x_arr = np.asarray(template.x, dtype=float)
    phi_arr = np.asarray(template.phi, dtype=float)
    p_vec = p_matrix_to_vector(template.P)
    A_F = build_af_matrix(x_arr, phi_arr)
    F = A_F @ p_vec
    return tuple(float(v) for v in F)


def c_a(template: Template) -> tuple[float, ...]:
    """Calcula la fuerza de estabilidad intra-cluster c_a para cada cluster.

    Referencia: O'Keeffe (2026), Sección III, Ecuación (7):
        c_a = sum_{b=1}^m P_ab x_b cos(phi_a - phi_b) > 0

    Args:
        template: Plantilla con fases y densidades definidas.

    Returns:
        Tupla con los coeficientes c_a para a = 1, ..., m.

    Raises:
        ValueError: Si la plantilla no contiene matriz P.
    """
    if template.P is None:
        raise ValueError("No se puede calcular c_a sin matriz P.")

    x_arr = np.asarray(template.x, dtype=float)
    phi_arr = np.asarray(template.phi, dtype=float)
    p_vec = p_matrix_to_vector(template.P)
    A_c = build_ac_matrix(x_arr, phi_arr)
    c_vec = x_arr + A_c @ p_vec
    return tuple(float(v) for v in c_vec)


def projected_laplacian(template: Template) -> tuple[tuple[float, ...], ...]:
    """Calcula la forma Laplaciana proyectada sobre 1_perp: Q^T L_w Q.

    Referencias:
        - O'Keeffe (2026), Sección III, Ecuación (8):
            L(y) = sum_{a < b} x_a x_b P_ab cos(phi_a - phi_b) (y_a - y_b)^2
        - Deflación de rotación global: Q base ortonormal m x (m-1) ortogonal al vector 1.

    Args:
        template: Plantilla con tamaños, fases y densidades definidas.

    Returns:
        Matriz (m-1) x (m-1) proyectada Q^T L_w Q como tupla de tuplas.

    Raises:
        ValueError: Si la plantilla no contiene matriz P.
    """
    if template.P is None:
        raise ValueError("No se puede construir L_w sin matriz P.")

    m = template.m
    x_arr = np.asarray(template.x, dtype=float)
    phi_arr = np.asarray(template.phi, dtype=float)
    p_vec = p_matrix_to_vector(template.P)

    Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
    M = build_m_matrix(x_arr, phi_arr, p_vec, Q)
    return tuple(tuple(float(v) for v in row) for row in M)


def verify_solution(
    x: tuple[float, ...] | np.ndarray,
    phi: tuple[float, ...] | np.ndarray,
    P: tuple[tuple[float, ...], ...] | np.ndarray,
    eps: float = 1e-6,
    eps_c: float = 1e-4,
    theta_tol: float = 1e-4,
    h_eff: float | None = None,
) -> dict[str, Any]:
    """Verificación independiente en numpy puro (sin cvxpy) de una solución candidata (x, phi, P).

    Todo resultado del solver SDP debe ser validado por esta función antes de
    reportarse como factible o de otorgar aptitud positiva en Fase II.

    Args:
        x: Tamaños de cluster (longitud m).
        phi: Fases de cluster (longitud m).
        P: Matriz de densidades simétrica m x m (con diagonal 1.0).
        eps: Margen espectral requerido para M (M >> eps * I).
        eps_c: Cota inferior para c_a (c_a >= eps_c).
        theta_tol: Tolerancia máxima en residuo de equilibrio (|F_a| <= theta_tol).
        h_eff: Escala efectiva opcional. Si se provee, el umbral de estabilidad
            para lambda_min(M) se escala como eps_efectivo = eps * (h_eff**2).

    Returns:
        Diccionario con las métricas recalculadas y bandera booleana `is_valid`:
          - `is_valid`: True si satisface todas las condiciones físicas de estabilidad.
          - `mu`: Conectividad recalculada min_a d_a.
          - `max_F`: Residuo máximo de equilibrio max_a |F_a|.
          - `min_c`: Coeficiente mínimo de estabilidad intra-cluster min_a c_a.
          - `min_eig_M`: Autovalor mínimo de M(p) = Q^T L_w Q.
          - `eps_effective`: Umbral efectivo aplicado para min_eig_M.
          - `p_bounds_ok`: True si todas las entradas de P caen en [-1e-9, 1.0 + 1e-9].
          - `reasons`: Lista de violaciones detectadas (vacía si is_valid es True).
    """
    x_arr = np.asarray(x, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)
    P_arr = np.asarray(P, dtype=float)
    m = len(x_arr)

    p_vec = p_matrix_to_vector(P_arr)
    A_F = build_af_matrix(x_arr, phi_arr)
    A_c = build_ac_matrix(x_arr, phi_arr)
    A_d = build_ad_matrix(x_arr)

    F = A_F @ p_vec
    c_vec = x_arr + A_c @ p_vec
    d_vec = x_arr + A_d @ p_vec
    mu_val = float(np.min(d_vec))

    Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
    M = build_m_matrix(x_arr, phi_arr, p_vec, Q)
    eigvals = scipy.linalg.eigvalsh(M)
    min_eig_M = float(np.min(eigvals))

    max_F = float(np.max(np.abs(F)))
    min_c = float(np.min(c_vec))

    min_p = float(np.min(p_vec))
    max_p = float(np.max(p_vec))
    p_bounds_ok = (min_p >= -1e-8) and (max_p <= 1.0 + 1e-8)

    eps_effective = eps * (h_eff**2) if h_eff is not None else eps
    slack = 1e-8 if h_eff is None else min(1e-8, 1e-3 * eps_effective)

    reasons: list[str] = []
    if max_F > theta_tol:
        reasons.append(f"Residuo de equilibrio excedido: max|F|={max_F:.2e} > {theta_tol:.2e}")
    if min_c < eps_c - 1e-8:
        reasons.append(f"Estabilidad intra-cluster insuficiente: min(c)={min_c:.2e} < {eps_c:.2e}")
    if min_eig_M < eps_effective - slack:
        reasons.append(
            f"Margen espectral insuficiente: lambda_min(M)={min_eig_M:.2e} < {eps_effective:.2e}"
        )
    if not p_bounds_ok:
        reasons.append(f"Densidades fuera de [0, 1]: min={min_p:.2e}, max={max_p:.2e}")

    is_valid = len(reasons) == 0

    return {
        "is_valid": is_valid,
        "mu": mu_val,
        "max_F": max_F,
        "min_c": min_c,
        "min_eig_M": min_eig_M,
        "eps_effective": eps_effective,
        "p_bounds_ok": p_bounds_ok,
        "reasons": reasons,
    }
