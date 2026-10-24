"""Pruebas unitarias de valores e identidades algebraicas exactas del dominio (Suite A).

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección III (A Cluster Model).
    - Kassabov, M., Strogatz, S. H., Townsend, A. (2021), Chaos 31, 073135.
"""

from __future__ import annotations

import math

import numpy as np

from kuramoto_templates.cases.c4 import get_c4_template
from kuramoto_templates.domain.conditions import (
    build_lw_matrix,
    c_a,
    connectivity,
    equilibrium_residual,
    get_pair_indices,
    projected_laplacian,
)
from kuramoto_templates.domain.template import Template


def test_equilibrium_antisymmetry_identity() -> None:
    """Valida la identidad algebraica fundamental de conservación de fase sum_a x_a F_a == 0.

    Qué valida:
        Para cualquier configuración (x, phi, P), la suma ponderada por los tamaños de cluster
        de los residuos de equilibrio F_a se anula idénticamente (rtol <= 1e-12):
            sum_{a=1}^m x_a F_a = sum_{a < b} x_a x_b P_ab [sin(phi_a-phi_b) + sin(phi_b-phi_a)] = 0

    De dónde viene el valor esperado:
        Identidad analítica exacta por antisimetría de la función seno: sin(-theta) = -sin(theta),
        válida para cualquier grafo y plantilla de Kuramoto (O'Keeffe 2026, Sección III, Ec. 6).
    """
    rng = np.random.default_rng(seed=42)
    m = 5
    raw_x = rng.exponential(scale=1.0, size=m)
    x = tuple(float(v) for v in raw_x / np.sum(raw_x))
    phi = tuple(float(v) for v in rng.uniform(-math.pi, math.pi, size=m))

    raw_p = rng.uniform(0.0, 1.0, size=(m, m))
    p_sym = 0.5 * (raw_p + raw_p.T)
    np.fill_diagonal(p_sym, 1.0)
    P = tuple(tuple(float(v) for v in row) for row in p_sym)

    template = Template(x=x, phi=phi, P=P)
    F = equilibrium_residual(template)

    weighted_sum = sum(x[a] * F[a] for a in range(m))
    assert math.isclose(weighted_sum, 0.0, abs_tol=1e-12)


def test_laplacian_ones_nullspace_identity() -> None:
    """Valida que el vector constante ones pertenezca al espacio nulo del Laplaciano L_w.

    Qué valida:
        L_w(p) @ ones == 0 exactamente para cualquier plantilla y vector de densidades p.
        Las sumas de cada fila y columna de L_w se anulan idénticamente.

    De dónde viene el valor esperado:
        Propiedad fundamental de los grafos y del operador Laplaciano ponderado:
            L_w(p) = sum_{a < b} w_ab (e_a - e_b)(e_a - e_b)^T.
        Como (e_a - e_b)^T @ 1 = 1 - 1 = 0 para cada arista (a, b), se tiene
        L_w @ 1 = sum_{a < b} w_ab (e_a - e_b) * 0 = 0
        (O'Keeffe 2026, Sección III, Ec. 8).
    """
    rng = np.random.default_rng(seed=123)
    m = 6
    raw_x = rng.exponential(scale=1.0, size=m)
    x_arr = raw_x / np.sum(raw_x)
    phi_arr = rng.uniform(-math.pi, math.pi, size=m)

    pairs = get_pair_indices(m)
    K = len(pairs)
    p_vec = rng.uniform(0.0, 1.0, size=K)

    L_w = build_lw_matrix(x_arr, phi_arr, p_vec)
    ones = np.ones(m)
    res = L_w @ ones

    assert np.allclose(res, 0.0, atol=1e-12)


def test_c4_exact_theoretical_properties() -> None:
    """Valida todas las propiedades exactas demostradas del estado C4.

    Qué valida:
        Para el estado de cuatro clusters C4 (x_a = 1/4, phi_a = a*pi/2, vecinos P=1, opuestos P=0):
          1. F_a == 0 para todo cluster a in {0, 1, 2, 3}.
          2. c_a == 1/4 para todo cluster a in {0, 1, 2, 3}.
          3. mu == 3/4 (conectividad mínima).
          4. M == 0 (matriz proyectada nula 3x3, autovalores idénticamente cero).

    De dónde viene el valor esperado:
        Demostrado analíticamente por cálculo directo en:
        - Kassabov, Strogatz & Townsend (2021), Chaos 31, 073135 (Sec. II, Fig. 1).
        - O'Keeffe (2026), arXiv:2610.21271v1, Sección III.a.
        Cálculo directo:
          - F_a = (1/4)[sin(pi/2) - sin(pi/2)] = 0.
          - c_a = x_a + sum_{vecinos} (1/4) cos(+-pi/2) = 1/4 + 0 = 1/4.
          - d_a = 1/4 * (1 + 1 + 1) = 3/4 => mu = 3/4.
          - cos(phi_a - phi_b) es 0 para vecinos (distancia pi/2) y P_ab = 0 para opuestos.
            Por tanto todas las aristas con P_ab > 0 tienen cos = 0, luego L_w = 0 y M = 0.
    """
    c4 = get_c4_template()

    # 1. Residuo de equilibrio F_a == 0
    F = equilibrium_residual(c4)
    for a in range(4):
        assert math.isclose(F[a], 0.0, abs_tol=1e-12)

    # 2. Estabilidad intra-cluster c_a == 1/4
    c = c_a(c4)
    for a in range(4):
        assert math.isclose(c[a], 0.25, abs_tol=1e-12)

    # 3. Conectividad mu == 3/4
    mu = connectivity(c4)
    assert math.isclose(mu, 0.75, abs_tol=1e-12)

    # 4. Hessiano proyectado M == 0 (3x3)
    M = projected_laplacian(c4)
    M_arr = np.array(M)
    assert M_arr.shape == (3, 3)
    assert np.allclose(M_arr, 0.0, atol=1e-12)
