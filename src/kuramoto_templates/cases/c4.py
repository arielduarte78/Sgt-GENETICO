"""Estado C4 (Cuatro clusters a cuartos de vuelta - Degenerado).

Referencias:
    - Kassabov, M., Strogatz, S. H., Townsend, A. (2021).
      "Sufficiently dense Kuramoto networks are globally synchronizing",
      Chaos 31, 073135. [Estado: DEMOSTRADO]
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección III.a (The four-cluster state C4)
      y Figura 1a. [Estado: CITADO]

Propiedades Teóricas Demostradas:
    - Tamaños: x_a = 1/4 para a in {0, 1, 2, 3}.
    - Fases: phi_a = a * pi / 2.
    - Densidades de enlace: completas internamente y entre vecinos (|a-b|=1 mod 4),
      nulas entre opuestos (|a-b|=2 mod 4).
    - Conectividad: mu(T) = 3/4 = 0.75.
    - Residuo de equilibrio: F_a = 0 idénticamente.
    - Estabilidad intra-cluster: c_a = 1/4 = 0.25 > 0.
    - Estabilidad inter-cluster: L_w = 0 idénticamente porque cos(+-pi/2) = 0
      y los pares opuestos tienen P_ab = 0. Estado situado en filo de navaja (razor's edge).
"""

from __future__ import annotations

import math

from kuramoto_templates.domain.template import Template


def get_c4_template() -> Template:
    """Construye la plantilla analítica exacta para el estado C4.

    Returns:
        Template con m=4 clusters, conectividad mu=3/4 y Hessiano idénticamente nulo L=0.
    """
    x = (0.25, 0.25, 0.25, 0.25)
    phi = (0.0, math.pi / 2.0, math.pi, -math.pi / 2.0)

    # Matriz P: completa con vecinos y consigo mismo, 0 con opuestos
    # Indices 0: 0 rad, 1: pi/2, 2: pi, 3: -pi/2 (o 3pi/2)
    P: tuple[tuple[float, ...], ...] = (
        (1.0, 1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0, 0.0),
        (0.0, 1.0, 1.0, 1.0),
        (1.0, 0.0, 1.0, 1.0),
    )

    return Template(x=x, phi=phi, P=P)
