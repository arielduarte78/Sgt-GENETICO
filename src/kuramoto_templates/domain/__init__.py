"""Módulo de Dominio Matemático Puro.

Este submódulo define las estructuras matemáticas y contratos algebraicos
para el modelo de Kuramoto sin dependencias de optimización numérica (cvxpy)
ni generadores estocásticos.
"""

from kuramoto_templates.domain.conditions import (
    c_a,
    connectivity,
    degrees,
    equilibrium_residual,
    projected_laplacian,
)
from kuramoto_templates.domain.gauge import canonicalize_gauge
from kuramoto_templates.domain.template import Template

__all__ = [
    "Template",
    "equilibrium_residual",
    "c_a",
    "projected_laplacian",
    "degrees",
    "connectivity",
    "canonicalize_gauge",
]
