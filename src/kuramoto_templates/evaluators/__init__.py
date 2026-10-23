"""Módulo de evaluadores de aptitud.

Contiene las interfaces abstractas y estructuras de datos para resolver el
problema SDP interno sobre plantillas de clusters.
"""

from kuramoto_templates.evaluators.base import Evaluator
from kuramoto_templates.evaluators.results import EvalResult

__all__ = ["Evaluator", "EvalResult"]
