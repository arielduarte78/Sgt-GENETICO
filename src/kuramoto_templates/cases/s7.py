"""Familia S7: Deformación del estado de cuatro clusters (O'Keeffe 2026).

Referencias:
    - O'Keeffe, K. P. (2026), "Beyond 11/16: dense Kuramoto networks that resist synchrony",
      arXiv:2610.21271v1, Sección IV.A (Table II), Sección IV.C, Sección V y Apéndice A.
      [Estado: CITADO (no verificado) - preprint sin arbitraje externo].

Constantes Asintóticas Citadas (h -> 0):
    - q* ≈ 0.6505691669... (raíz única en (0, 1) del polinomio de grado 8, Ec. 3).
    - w* ≈ 1.3982781340... (razón de apertura central dada por Ec. 11).
    - mu_0 = (4 + 2*q*) / (7 + q*) ≈ 0.6929077063... (Ec. 2 y 10).
    - a* = 1 / (7 + q*) ≈ 0.130709... (tamaño de cada mitad dividida).
    - x_W = 1 - 6*a* ≈ 0.215744... (tamaño del 4to cluster entero).

Estructura Límite (Tabla II):
    - Clusters: (Am, Ap, Bm, Bp, Cm, Cp, W)
    - Fases: (-h, +h, pi/2 - w*h, pi/2 + w*h, pi - h, pi + h, 3*pi/2)
    - Enlaces en el límite:
        - Completos (P=1) internamente y entre vecinos sin cruce.
        - Cuatro enlaces cruzados atenuados con densidad q.
        - Enlaces entre opuestos nulos (P=0).
"""

from __future__ import annotations

import math

import mpmath  # type: ignore[import-untyped]

from kuramoto_templates.domain.template import Template

# Constantes asintóticas citadas de O'Keeffe (2026, Sec. IV.A) - [CITADO (no verificado)]
Q_STAR_CITADO: float = 0.6505691669447988
W_STAR_CITADO: float = 1.3982781340204612
MU_ZERO_CITADO: float = 0.6929077063695864


def compute_q_star_and_w_star(dps: int = 50) -> tuple[float, float]:
    """Calcula q* y w* resolviendo el sistema analítico (11) del preprint en alta precisión.

    Referencias:
        O'Keeffe (2026), arXiv:2610.21271v1, Ecuaciones (3) y (11):
          w^2 = (q^2 + 6q + 1) / (q + 1)^2
          4(1 - q^2) w = 2 q^2 w^2 + (3q - 1)(q + 1)
          Polinomio irreducible de grado 8:
            9q^8 + 72q^7 + 288q^6 + 448q^5 + 470q^4 + 88q^3 - 208q^2 - 128q - 15 = 0

    Args:
        dps: Número de dígitos decimales de precisión en aritmética mpmath (default 50).

    Returns:
        Tupla (q*, w*) como números flotantes de doble precisión estándar.
    """
    old_dps = mpmath.mp.dps
    try:
        mpmath.mp.dps = dps
        coeffs = [9, 72, 288, 448, 470, 88, -208, -128, -15]

        def poly(q_val: mpmath.mpf) -> mpmath.mpf:
            res = mpmath.mpf(0)
            for power, coeff in enumerate(reversed(coeffs)):
                res += coeff * (q_val**power)
            return res

        # Raíz en el intervalo admisible (0, 1)
        q_root = mpmath.findroot(poly, (mpmath.mpf("0.5"), mpmath.mpf("0.8")))
        w_sq = (q_root**2 + 6 * q_root + 1) / ((q_root + 1) ** 2)
        w_root = mpmath.sqrt(w_sq)

        return float(q_root), float(w_root)
    finally:
        mpmath.mp.dps = old_dps


def get_s7_template(h: float = 0.002, t: float = 0.0) -> Template:
    """Construye un miembro de la familia S7 a partir del patrón de enlaces de la Tabla II.

    Referencias:
        O'Keeffe (2026), arXiv:2610.21271v1, Sec. IV.A (Tabla II) y Sec. IV.C.
        - Tamaños: x_a = a = 1 / (7 + q) para las 6 mitades divididas, x_W = 1 - 6a.
        - Fases: phi_c = j_c * pi / 2 + h * s_c con j = (0, 0, 1, 1, 2, 2, 3) y
          s = (-1, 1, -w, w, -1, 1, 0).
        - Matriz P: patrón de enlaces de la Tabla II evaluado con densidad q = q* - t.

    Args:
        h: Ángulo de splitting positivo (pequeña deformación respecto a C4).
        t: Holgura respecto a la densidad crítica q = q* - t, t in [0, 0.05].

    Returns:
        Template inmutable con m=7 clusters y densidades P conformes a Tabla II.
    """
    q_star, w_star = compute_q_star_and_w_star(dps=50)
    q = q_star - t
    w = w_star

    a = 1.0 / (7.0 + q)
    x_w = 1.0 - 6.0 * a
    x = (a, a, a, a, a, a, x_w)

    # Orden de clusters de la Tabla II: (Am, Ap, Bm, Bp, Cm, Cp, W)
    j = (0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0)
    s = (-1.0, 1.0, -w, w, -1.0, 1.0, 0.0)
    phi = tuple(float(j[k] * (math.pi / 2.0) + h * s[k]) for k in range(7))

    # Matriz P de la Tabla II (Sec. IV.A):
    # Filas: Am, Ap, Bm, Bp, Cm, Cp, W
    P: tuple[tuple[float, ...], ...] = (
        (1.0, 1.0, 1.0, q, 0.0, 0.0, 1.0),  # Am
        (1.0, 1.0, q, 1.0, 0.0, 0.0, 1.0),  # Ap
        (1.0, q, 1.0, 1.0, 1.0, q, 0.0),  # Bm
        (q, 1.0, 1.0, 1.0, q, 1.0, 0.0),  # Bp
        (0.0, 0.0, 1.0, q, 1.0, 1.0, 1.0),  # Cm
        (0.0, 0.0, q, 1.0, 1.0, 1.0, 1.0),  # Cp
        (1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0),  # W
    )

    return Template(x=x, phi=phi, P=P)


def get_s7_appendix_a_template(
    h: float = 0.002, t: float = 0.0, kappa: float = 1.9
) -> Template:
    """Construye la plantilla S7 según la expansión analítica del Apéndice A.

    Referencias:
        O'Keeffe (2026), arXiv:2610.21271v1, Apéndice A (Sec. A.1, Eqs. A1-A3, Lemma 3).
        [Estado: CITADO (no verificado)]
        Densidades: P(h) = p + h * p' + Delta(h), donde:
          - p es el patrón de la Tabla II.
          - p' es la corrección de primer orden dependiente de kappa, u13, A1, A2.
          - Delta(h) = O(h^2) está soportada en el camino Am-Bm-Cm-W-Cp-Bp-Ap
            y resuelve F(P(h)) = 0 exactamente.

    Args:
        h: Ángulo de splitting h > 0.
        t: Holgura respecto a la densidad crítica q = q* - t, t in [0, 0.05].
        kappa: Parámetro de adelgazamiento uniforme (por defecto 19/10 = 1.9).

    Returns:
        Template inmutable con m=7 clusters y densidades P(h) de equilibrio exacto.
    """
    import numpy as np

    q_star, w_star = compute_q_star_and_w_star(dps=50)
    q = q_star - t
    w = w_star

    a = 1.0 / (7.0 + q)
    x_w = 1.0 - 6.0 * a
    x = (a, a, a, a, a, a, x_w)

    j = (0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0)
    s = (-1.0, 1.0, -w, w, -1.0, 1.0, 0.0)
    phi = tuple(float(j[k] * (math.pi / 2.0) + h * s[k]) for k in range(7))

    # Parámetros analíticos del Apéndice A.1
    u13 = 0.5 * (
        2.0 * (q**2) * (w**2)
        + 2.0 * (q**2) * w
        - (q**2)
        - 6.0 * q
        + 4.0 * (w**2)
        - 2.0 * w
        + 3.0
    )
    u = u13
    a1 = w * (w + 1.0) - u / 2.0
    a2 = u / 2.0 - w * (w - 1.0)

    # Patrón límite p (Tabla II)
    p_mat = np.zeros((7, 7), dtype=float)
    np.fill_diagonal(p_mat, 1.0)
    p_mat[0, 1] = p_mat[1, 0] = 1.0  # Am, Ap
    p_mat[2, 3] = p_mat[3, 2] = 1.0  # Bm, Bp
    p_mat[4, 5] = p_mat[5, 4] = 1.0  # Cm, Cp
    p_mat[0, 2] = p_mat[2, 0] = 1.0  # Am, Bm
    p_mat[1, 3] = p_mat[3, 1] = 1.0  # Ap, Bp
    p_mat[2, 4] = p_mat[4, 2] = 1.0  # Bm, Cm
    p_mat[3, 5] = p_mat[5, 3] = 1.0  # Bp, Cp
    p_mat[0, 6] = p_mat[6, 0] = 1.0  # Am, W
    p_mat[1, 6] = p_mat[6, 1] = 1.0  # Ap, W
    p_mat[4, 6] = p_mat[6, 4] = 1.0  # Cm, W
    p_mat[5, 6] = p_mat[6, 5] = 1.0  # Cp, W
    p_mat[0, 3] = p_mat[3, 0] = q  # Am, Bp
    p_mat[1, 2] = p_mat[2, 1] = q  # Ap, Bm
    p_mat[2, 5] = p_mat[5, 2] = q  # Bm, Cp
    p_mat[3, 4] = p_mat[4, 3] = q  # Bp, Cm

    # Corrección de primer orden p' (Eq. A1)
    p_prime = np.zeros((7, 7), dtype=float)
    p_prime[0, 2] = p_prime[2, 0] = a1 - kappa
    p_prime[3, 5] = p_prime[5, 3] = a1 - kappa

    p_prime[1, 2] = p_prime[2, 1] = a2 - kappa
    p_prime[3, 4] = p_prime[4, 3] = a2 - kappa

    p_prime[1, 3] = p_prime[3, 1] = -kappa
    p_prime[2, 4] = p_prime[4, 2] = -kappa
    p_prime[0, 3] = p_prime[3, 0] = -kappa
    p_prime[2, 5] = p_prime[5, 2] = -kappa

    p_prime[0, 6] = p_prime[6, 0] = (a / x_w) * (a1 + 2.0 - 2.0 * kappa)
    p_prime[5, 6] = p_prime[6, 5] = (a / x_w) * (a1 + 2.0 - 2.0 * kappa)

    p_prime[1, 6] = p_prime[6, 1] = (a / x_w) * (a2 - 2.0 - 2.0 * kappa)
    p_prime[4, 6] = p_prime[6, 4] = (a / x_w) * (a2 - 2.0 - 2.0 * kappa)

    p_lin = p_mat + h * p_prime

    # Camino de soporte de Delta(h): Am-Bm-Cm-W-Cp-Bp-Ap
    path = [(0, 2), (2, 4), (4, 6), (6, 5), (5, 3), (3, 1)]

    # Residuo de equilibrio F(p_lin)
    x_arr = np.asarray(x, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)
    f_lin = np.zeros(7, dtype=float)
    for i in range(7):
        f_lin[i] = float(np.sum(x_arr * np.sin(phi_arr[i] - phi_arr) * p_lin[i, :]))

    # Matriz jacobiana J(h) de 6x6 para las 6 aristas del camino
    j_mat = np.zeros((6, 6), dtype=float)
    for k, (c, d) in enumerate(path):
        diff = phi_arr[c] - phi_arr[d]
        s_val = float(np.sin(diff))
        if c < 6:
            j_mat[c, k] += x_arr[d] * s_val
        if d < 6:
            j_mat[d, k] -= x_arr[c] * s_val

    delta = -np.linalg.solve(j_mat, f_lin[:6])

    p_full = np.copy(p_lin)
    for k, (c, d) in enumerate(path):
        p_full[c, d] += delta[k]
        p_full[d, c] += delta[k]

    p_tuple = tuple(tuple(float(v) for v in row) for row in p_full)
    return Template(x=x, phi=phi, P=p_tuple)
