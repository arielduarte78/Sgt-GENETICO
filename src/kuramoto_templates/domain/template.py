"""Estructura de datos para plantillas de clusters en el modelo de Kuramoto.

Referencia Matemática:
    O'Keeffe, K. P. (2026). "Beyond 11/16: dense Kuramoto networks that resist synchrony",
    arXiv:2610.21271v1, Sección III (A Cluster Model), Definición 1.

    "Definición 1. Una plantilla con m clusters es una terna (x, phi, P) que consta de:
     (a) tamaños x_a > 0 con sum_a x_a = 1: fracción de osciladores en el cluster a;
     (b) fases phi_a: fase común del cluster a;
     (c) densidades P_ab = P_ba en [0, 1]: fracción de enlaces posibles presentes entre
         clusters a y b, o dentro del cluster a cuando b = a."
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Template:
    """Representación inmutable de una plantilla de clusters (x, phi, P).

    Atributos:
        x: Tupla de fracciones de osciladores (x_a > 0, sum x_a = 1).
        phi: Tupla de fases comunes de cada cluster (phi_a en (-pi, pi]).
        P: Matriz simétrica de densidades P_ab en [0, 1], o None si aún no se han
           resuelto las densidades de enlace.
    """

    x: tuple[float, ...]
    phi: tuple[float, ...]
    P: tuple[tuple[float, ...], ...] | None = None

    def __post_init__(self) -> None:
        """Valida las condiciones fundamentales de la Definición 1 (O'Keeffe 2026, Sec. III)."""
        m = len(self.x)
        if m < 2:
            raise ValueError(f"Una plantilla requiere al menos m=2 clusters, recibido m={m}.")

        if len(self.phi) != m:
            raise ValueError(
                f"Discrepancia dimensional: len(x)={m} no coincide con len(phi)={len(self.phi)}."
            )

        for a, val in enumerate(self.x):
            if val <= 0.0:
                raise ValueError(f"Tamaño de cluster inválido x[{a}]={val}: debe ser x_a > 0.")

        total_size = sum(self.x)
        if not math.isclose(total_size, 1.0, rel_tol=1e-7, abs_tol=1e-7):
            raise ValueError(
                f"La suma de tamaños de clusters debe ser 1.0, suma obtenida: {total_size}."
            )

        if self.P is not None:
            if len(self.P) != m:
                raise ValueError(
                    f"Matriz P con {len(self.P)} filas; se esperaban {m} clusters."
                )

            for a, row in enumerate(self.P):
                if len(row) != m:
                    raise ValueError(
                        f"Fila P[{a}] tiene longitud {len(row)}, se esperaba longitud {m}."
                    )
                for b, val in enumerate(row):
                    if not (0.0 <= val <= 1.0 + 1e-9):
                        raise ValueError(
                            f"Densidad P[{a}][{b}]={val} fuera del intervalo admisible [0, 1]."
                        )
                    # Simetría P_ab = P_ba
                    if not math.isclose(val, self.P[b][a], abs_tol=1e-8):
                        raise ValueError(
                            f"Asimetría en densidades de enlace: P[{a}][{b}]={val} != "
                            f"P[{b}][{a}]={self.P[b][a]}."
                        )

    @property
    def m(self) -> int:
        """Número de clusters de la plantilla."""
        return len(self.x)
