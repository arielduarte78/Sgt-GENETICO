"""Registro estructurado y auditable de experimentos (JSONL + Manifiesto).

Referencias:
    - Plan de Diseño del Algoritmo Genético, Secciones 8 y 9 (Trazabilidad y Reproducibilidad).
    - Reglas de Honestidad Epistémica: Todo experimento debe asociar parámetros exactos,
      semilla aleatoria, versión de software, commit de git y tolerancias del solver.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


class RunLogger:
    """Gestiona el registro auditable de corridas de búsqueda y evaluación.

    Produce dos artefactos principales por corrida:
      1. `manifest.json`: Metadatos del entorno, versión de git, configuración completa y semilla.
      2. `metrics.jsonl`: Flujo de eventos y métricas por generación e isla en formato JSON Lines.
    """

    def __init__(self, run_dir: Path | str, config: dict[str, Any]) -> None:
        """Inicializa el logger y asegura la existencia del directorio de corrida.

        Args:
            run_dir: Directorio destino (usualmente results/runs/<run_id>).
            config: Diccionario con la configuración completa del experimento.
        """
        self.run_dir = Path(run_dir)
        self.config = config
        self.manifest_path = self.run_dir / "manifest.json"
        self.log_path = self.run_dir / "metrics.jsonl"

    def write_manifest(self, metadata: dict[str, Any] | None = None) -> Path:
        """Escribe el archivo manifest.json con metadatos del entorno y configuración.

        Args:
            metadata: Diccionario opcional con metadatos del sistema (git commit, fecha).

        Returns:
            Path al manifiesto generado.

        Raises:
            NotImplementedError: Pendiente de integración con metadatos de git/plataforma.
        """
        raise NotImplementedError(
            "TODO: Implementar escritura de manifest.json con hash de git y configuración."
        )

    def log_step(
        self,
        generation: int,
        island_id: int,
        metrics: dict[str, float | int | str],
    ) -> None:
        """Registra un evento generacional en formato JSON Lines.

        Args:
            generation: Número de generación.
            island_id: Identificador de la isla ejecutada.
            metrics: Diccionario de métricas (mejor mu, residuo, margen espectral, tiempo).

        Raises:
            NotImplementedError: Pendiente de implementación de streaming a metrics.jsonl.
        """
        raise NotImplementedError(
            "TODO: Implementar serialización a metrics.jsonl en modo append."
        )
