"""Pruebas de regresión para la construcción analítica del Apéndice A (CITADO no verificado).

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Apéndice A (Proof of Theorem 1):
      P(h) = p + h p' + Delta(h), con p' (A1), A1, A2, u13 (A3), kappa=19/10,
      y Delta en el camino Am-Bm-Cm-W-Cp-Bp-Ap resuelto para F=0.
"""

from __future__ import annotations

import pytest

from kuramoto_templates.cases.s7 import (
    compute_q_star_and_w_star,
    get_s7_appendix_a_template,
)
from kuramoto_templates.domain.conditions import connectivity, equilibrium_residual


@pytest.mark.cited
def test_s7_appendix_a_perturbation_construction() -> None:
    """Valida la construcción analítica P(h) = p + h p' + Delta(h) del Apéndice A.

    Qué valida:
        Para t = 0 y h = 0.002:
          1. F(P(h)) ~ 0 (residuo de equilibrio nulo a precisión de máquina < 1e-12).
          2. P(h) en [0, 1] elemento a elemento.
          3. Conectividad mu >= mu0(q) - 0.9956 * h (cota explícita del Remark de la pág. 11).

    De dónde viene el valor esperado:
        O'Keeffe (2026), arXiv:2610.21271v1, Apéndice A, Eqs. (A1)-(A3), Lemma 3 y Remark:
        - Lemma 3(a): Delta(h) = -J(h)^(-1) F(p + h p') es O(h^2) y anula F exactamente.
        - Remark (pág. 11): para h <= 1/8 y t en [0, 0.05], |Delta| <= 1.135 h^2 y
          mu >= (4 + 2q)/(7 + q) - 0.9956 h.
    """
    h = 0.002
    t = 0.0
    kappa = 19.0 / 10.0

    q_star_mp, _ = compute_q_star_and_w_star(dps=50)
    q_val = float(q_star_mp)
    mu0 = (4.0 + 2.0 * q_val) / (7.0 + q_val)
    theoretical_bound = mu0 - 0.9956 * h

    template = get_s7_appendix_a_template(h=h, t=t, kappa=kappa)

    # 1. Verificar que F(P(h)) sea nulo a precisión numérica
    F = equilibrium_residual(template)
    max_F = max(abs(v) for v in F)
    assert max_F < 1e-12, f"Residuo de equilibrio no nulo en P(h): max|F| = {max_F:.2e}"

    # 2. Verificar que las densidades P(h) permanezcan en [0, 1]
    assert template.P is not None
    for i, row in enumerate(template.P):
        for j, val in enumerate(row):
            assert -1e-14 <= val <= 1.0 + 1e-14, (
                f"Densidad P[{i}, {j}] = {val} fuera de [0, 1]."
            )

    # 3. Verificar cota de conectividad mu >= mu0 - 0.9956 * h
    mu_obtained = connectivity(template)
    assert mu_obtained >= theoretical_bound - 1e-12, (
        f"Conectividad {mu_obtained:.6f} violó la cota del Remark {theoretical_bound:.6f}."
    )

    print(
        f"\n[APPENDIX A S7 h={h}] max|F|={max_F:.2e}, "
        f"mu={mu_obtained:.6f}, Cota Remark={theoretical_bound:.6f}, "
        f"Margen={mu_obtained - theoretical_bound:+.2e}"
    )
