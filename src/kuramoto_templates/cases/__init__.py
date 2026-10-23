"""Casos de prueba con procedencia y trazabilidad científica.

Incluye estados de referencia teórica:
    - c4: Estado degenerado de cuatro clusters (Kassabov et al. 2021).
    - canale: Stub de la construcción de 8 clusters (Canale 2022).
    - s7: Stub de la familia split S7 con parámetros h y t (O'Keeffe 2026).
"""

from kuramoto_templates.cases.c4 import get_c4_template
from kuramoto_templates.cases.canale import get_canale_template
from kuramoto_templates.cases.s7 import (
    MU_ZERO_CITADO,
    Q_STAR_CITADO,
    W_STAR_CITADO,
    get_s7_template,
)

__all__ = [
    "get_c4_template",
    "get_canale_template",
    "get_s7_template",
    "Q_STAR_CITADO",
    "W_STAR_CITADO",
    "MU_ZERO_CITADO",
]
