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

ESTADO DEL MÓDULO: STUB CON PARÁMETROS (h, t)
    A splitting finito h > 0, las densidades reciben correcciones:
        P(h) = p + h*p' + Delta(h)
    donde p' proviene de la Ec. (A1) y Delta(h) = O(h^2) es una corrección sobre el
    camino de 6 aristas Am-Bm-Cm-W-Cp-Bp-Ap para forzar el equilibrio exacto F(P)=0.
    Los valores numéricos certificados para miembros particulares (como h=0.002, mu ≈ 0.69272,
    o Proposition 1 con mu* ≈ 0.6929050 y D ≈ 3.06e42) residen en el Supplementary Material
    (certificate.json), el cual no está disponible en este repositorio.
"""

from __future__ import annotations

from kuramoto_templates.domain.template import Template

# Constantes asintóticas citadas de O'Keeffe (2026, Sec. IV.A) - [CITADO (no verificado)]
Q_STAR_CITADO: float = 0.6505691669447988
W_STAR_CITADO: float = 1.3982781340204612
MU_ZERO_CITADO: float = 0.6929077063695864


def get_s7_template(h: float = 0.002, t: float = 0.0) -> Template:
    """Construye un miembro de la familia S7 dado el ángulo de splitting h y holgura t.

    Args:
        h: Ángulo de splitting positivo (pequeña deformación respecto a C4).
        t: Holgura respecto a la densidad crítica q = q* - t, t in [0, 0.05].

    Returns:
        Template con m=7 clusters y densidades P(h) equilibradas.

    TODO:
        1. Resolver el sistema lineal de caminos J(h) Delta = -F(p + h*p') descrito
           en el Lema 3 (Apéndice A) para computar Delta(h) = O(h^2) con equilibrio exacto,
           o bien importar `certificate.json` si se obtiene acceso al Supplementary Material.
        2. Mantener la verificación de que 0 <= P(h) <= 1 entrada a entrada.

    Raises:
        NotImplementedError: Stub a la espera de la resolución del sistema de corrección Delta(h).
    """
    raise NotImplementedError(
        f"TODO: Implementar corrección analítica/numérica Delta(h) para S7 (h={h}, t={t}) "
        "conforme a O'Keeffe (2026, Apéndice A.1) sin inventar datos no verificados."
    )
