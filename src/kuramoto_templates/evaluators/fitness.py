"""Definición de la función de fitness continua para el evaluador SDP.

Especificación:
    si s* > tol:  f = -s*     (distancia a la conectividad objetivo, <= 0)
    si s* <= tol: f = +lam*   (margen espectral, >= 0)
Sin pesos 1e10 ni log10 de penalización.
"""

from __future__ import annotations


def compute_continuous_fitness(
    lam: float,
    s: float,
    tol: float = 1e-6,
    lam_floor: float = 1e-7,
) -> float:
    """Calcula la función de fitness continua sin discontinuidades ni escalas arbitrarias.

    Args:
        lam: Margen espectral alcanzado lambda* (autovalor mínimo de M(p)).
        s: Holgura de conectividad s* (déficit respecto a la conectividad objetivo mu_t).
        tol: Tolerancia numérica para considerar que el objetivo de conectividad se cumple.
        lam_floor: Piso espectral mínimo para considerar estrictamente factible y estable.

    Returns:
        - f = -s* si s* > tol (penalización lineal por déficit de conectividad, <= 0)
        - f = +lam* si s* <= tol y lam* >= lam_floor (premio por margen espectral, > 0)
        - f = lam* - lam_floor si s* <= tol y lam* < lam_floor (penalización por debajo del piso, < 0)
    """
    if tol <= 0.0:
        raise ValueError(f"tol debe ser estrictamente positiva, recibido {tol}.")
    if lam_floor <= 0.0:
        raise ValueError(f"lam_floor debe ser estrictamente positiva, recibido {lam_floor}.")

    if s > tol:
        return -float(s)
    if lam >= lam_floor:
        return float(lam)
    return float(lam) - float(lam_floor)
