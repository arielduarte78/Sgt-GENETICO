"""Módulo de evaluadores de aptitud.

Contiene las interfaces abstractas, configuraciones y solucionadores para el
problema SDP interno sobre plantillas de clusters.
"""

from kuramoto_templates.evaluators.base import Evaluator
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.fitness import compute_continuous_fitness
from kuramoto_templates.evaluators.results import EvalResult
from kuramoto_templates.evaluators.sdp import SDPEvaluator

__all__ = [
    "Evaluator",
    "EvalResult",
    "EvalConfig",
    "SDPEvaluator",
    "compute_continuous_fitness",
]
