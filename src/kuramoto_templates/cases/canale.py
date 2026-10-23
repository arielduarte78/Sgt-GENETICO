"""Estado de Canale de ocho clusters (Cota 11/16).

Referencias:
    - Canale, E. A. (2022). "From weighted to unweighted graphs in synchronizing graph theory",
      arXiv:2209.06362. [Estado: CITADO]
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección III.c (Canale's state)
      y Figura 1b. [Estado: CITADO (no verificado)]

Propiedades Citadas de la Literatura:
    - Configuración: m = 8 clusters iguales en dos retículos de cuartos entrelazados,
      uno rotado respecto al otro por arctan(1/2) ≈ 26.6°.
    - Cada cluster está acoplado a 5.5 de los 8 clusters (densidades en {0, 0.5, 1}).
    - Conectividad: mu -> 11/16 = 0.6875.
    - Estabilidad estricta demostrada analíticamente en la literatura justo por debajo de 11/16.

ESTADO DEL MÓDULO: STUB
    No se disponen en el repositorio de los coeficientes exactos de la matriz P 8x8 de Canale.
    Conforme a las Reglas de Honestidad Epistémica, queda prohibido fabricar densidades.
"""

from __future__ import annotations

from kuramoto_templates.domain.template import Template


def get_canale_template() -> Template:
    """Retorna la plantilla de 8 clusters de Canale.

    TODO:
        1. Consultar el paper original de Canale (2022, arXiv:2209.06362, Sec. 3)
           para transcribir la matriz exacta P (8x8) y el desfase angular exacto arctan(1/2).
        2. Determinar si los enlaces de densidad 0.5 requieren un ajuste de perturbación epsilon
           para alcanzar estabilidad estricta con conectividad próxima a 11/16.

    Raises:
        NotImplementedError: Stub a la espera de los datos exactos del paper de Canale.
    """
    raise NotImplementedError(
        "TODO: Implementar get_canale_template() extrayendo la matriz de enlaces 8x8 "
        "directamente de Canale (2022, arXiv:2209.06362) sin inventar densidades."
    )
