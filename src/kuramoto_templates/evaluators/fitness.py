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

    Especificación ordenada (3 niveles):
    - a mu_t dado:
      * si lam >= lam_floor y s <= tol: f = +lam (factible estable, > 0)
      * si lam >= lam_floor y s > tol:  f = -s (estable, falta conectividad, <= 0)
      * si lam < lam_floor: f = -1.0 + lam_clipped (clip >= -1), peor que cualquier estable.

    Args:
        lam: Margen espectral alcanzado lambda* o lambda_hat (autovalor mínimo de M(p) o G_hat(z)).
        s: Holgura de conectividad s* (déficit respecto a la conectividad objetivo mu_t).
        tol: Tolerancia numérica para considerar que el objetivo de conectividad se cumple.
        lam_floor: Piso espectral mínimo para considerar estrictamente factible y estable.

    Returns:
        Valor de fitness escalar f.
    """
    if tol <= 0.0:
        raise ValueError(f"tol debe ser estrictamente positiva, recibido {tol}.")
    if lam_floor <= 0.0:
        raise ValueError(f"lam_floor debe ser estrictamente positiva, recibido {lam_floor}.")

    lam_val = float(lam)
    s_val = float(s)

    if lam_val >= lam_floor:
        if s_val <= tol:
            return lam_val
        return -s_val
    else:
        # lam < lam_floor: clip a >= -1, SIEMPRE peor que cualquier estable
        lam_clipped = max(lam_val, -1.0)
        if lam_clipped > 0.0:
            lam_clipped = lam_clipped - float(lam_floor)
        return -1.0 + lam_clipped

