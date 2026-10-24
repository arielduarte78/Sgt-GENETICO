"""Estructura del genoma, decodificación a fenotipo (x, phi) y canonicalización de gauge.

Referencias:
    - Plan de Diseño del Algoritmo Genético, Sección 5.1 (Genotipo de clusters):
      g = (m, Kref, k, sgn, ell, u), con m en [2, 10], Kref en {4, 8, 12, 16}.
      x = xmin + (1 - m * xmin) * softmax(u), xmin = 1e-3.
      phi = wrap_(-pi, pi](2 * pi * k / Kref + sgn * 10**ell), ell en [-6, -0.5].
    - Gauge: orden por fase, mayor cluster anclado a k=0, y mínimo lexicográfico con su reflejo.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

VALID_K_REF: tuple[int, ...] = (4, 8, 12, 16)
X_MIN: float = 1e-3
ELL_MIN: float = -6.0
ELL_MAX: float = -0.5


def wrap_minus_pi_to_pi(val: float) -> float:
    """Envuelve un ángulo estrictamente en el intervalo semiarbierto (-pi, pi]."""
    two_pi = 2.0 * math.pi
    res = (val + math.pi) % two_pi - math.pi
    if res <= -math.pi:
        res = math.pi
    return float(res)


def decode_x(u: tuple[float, ...], xmin: float = X_MIN) -> tuple[float, ...]:
    """Decodifica el vector continuo u a tamaños x en el símplex abierto con x_a >= xmin.

    Garantiza:
        sum_a x_a = 1.0 exactamente, y x_a >= xmin > 0 para todo a.
    """
    m = len(u)
    if m < 2 or m > 10:
        raise ValueError(f"m debe estar en [2, 10], recibido m={m}.")
    if xmin <= 0.0 or m * xmin >= 1.0:
        raise ValueError(f"xmin={xmin} inválido para m={m}.")

    max_u = max(u)
    exp_u = [math.exp(v - max_u) for v in u]
    sum_exp = sum(exp_u)
    softmax_u = [v / sum_exp for v in exp_u]

    slack = 1.0 - m * xmin
    raw_x = [xmin + slack * sm for sm in softmax_u]
    total = sum(raw_x)
    return tuple(float(v / total) for v in raw_x)


def decode_phi(
    K_ref: int,
    k: tuple[int, ...],
    sgn: tuple[int, ...],
    ell: tuple[float, ...],
) -> tuple[float, ...]:
    """Decodifica los parámetros angulares discretos y continuos a fases phi en (-pi, pi]."""
    m = len(k)
    two_pi = 2.0 * math.pi
    phases: list[float] = []
    for a in range(m):
        base_angle = two_pi * (k[a] % K_ref) / float(K_ref)
        offset = float(sgn[a]) * (10.0 ** ell[a])
        phases.append(wrap_minus_pi_to_pi(base_angle + offset))
    return tuple(phases)


@dataclass(frozen=True)
class Genome:
    """Genoma inmutable que define la parametrización de una plantilla de clusters.

    Atributos:
        m: Número de clusters en [2, 10].
        K_ref: Tamaño del retículo de referencia en {4, 8, 12, 16}.
        k: Tupla de longitud m con índices de ancla en {0, ..., K_ref - 1}.
        sgn: Tupla de longitud m con signos de desvío en {-1, 0, 1}.
        ell: Tupla de longitud m con exponente de desvío log10 en [-6.0, -0.5].
        u: Tupla de longitud m con logits para los tamaños de cluster x_a.
    """

    m: int
    K_ref: int
    k: tuple[int, ...]
    sgn: tuple[int, ...]
    ell: tuple[float, ...]
    u: tuple[float, ...]

    def __post_init__(self) -> None:
        """Valida rigurosamente las invariantes y cotas del genoma."""
        if not (2 <= self.m <= 10):
            raise ValueError(f"m debe estar en [2, 10], recibido m={self.m}.")
        if self.K_ref not in VALID_K_REF:
            raise ValueError(
                f"K_ref debe pertenecer a {VALID_K_REF}, recibido {self.K_ref}."
            )
        for name, tup in [
            ("k", self.k),
            ("sgn", self.sgn),
            ("ell", self.ell),
            ("u", self.u),
        ]:
            if len(tup) != self.m:
                raise ValueError(
                    f"Longitud de '{name}' ({len(tup)}) no coincide con m={self.m}."
                )

        for k_val in self.k:
            if not (0 <= k_val < self.K_ref):
                raise ValueError(
                    f"Índice de ancla k={k_val} fuera de rango [0, {self.K_ref - 1}]."
                )
        for s_val in self.sgn:
            if s_val not in (-1, 0, 1):
                raise ValueError(f"Signo sgn={s_val} debe pertenecer a {{-1, 0, 1}}.")
        for ell_val in self.ell:
            if not (ELL_MIN <= ell_val <= ELL_MAX):
                raise ValueError(
                    f"ell={ell_val} fuera de rango [{ELL_MIN}, {ELL_MAX}]."
                )


def to_phenotype(genome: Genome) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Decodifica el genoma a la tupla fenotípica (x, phi)."""
    x = decode_x(genome.u, xmin=X_MIN)
    phi = decode_phi(genome.K_ref, genome.k, genome.sgn, genome.ell)
    return x, phi


def _normalize_clusters(
    m: int,
    K_ref: int,
    k_list: list[int],
    sgn_list: list[int],
    ell_list: list[float],
    u_list: list[float],
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[float, ...], tuple[float, ...]]:
    """Ordena los clusters por fase y rota el cluster mayor a ancla k=0."""
    phases = decode_phi(K_ref, tuple(k_list), tuple(sgn_list), tuple(ell_list))
    sizes = decode_x(tuple(u_list))

    # 1. Ordenar por fase (con desempate por tamaño x)
    indices = list(range(m))
    indices.sort(key=lambda idx: (phases[idx], sizes[idx]))

    sorted_k = [k_list[i] for i in indices]
    sorted_sgn = [sgn_list[i] for i in indices]
    sorted_ell = [ell_list[i] for i in indices]
    sorted_u = [u_list[i] for i in indices]
    sorted_sizes = [sizes[i] for i in indices]

    # 2. Encontrar el cluster con mayor tamaño x (desempate por primer índice)
    max_idx = max(range(m), key=lambda idx: sorted_sizes[idx])
    anchor_shift = sorted_k[max_idx]

    # Rotar anclas para que el cluster mayor quede en k=0
    shifted_k = [(val - anchor_shift) % K_ref for val in sorted_k]

    # Reordenar nuevamente por fase resultante para preservar monotonicidad
    recomputed_phases = decode_phi(
        K_ref, tuple(shifted_k), tuple(sgn_list), tuple(ell_list)
    )
    final_indices = list(range(m))
    final_indices.sort(key=lambda idx: (recomputed_phases[idx], sorted_sizes[idx]))

    final_k = tuple(shifted_k[i] for i in final_indices)
    final_sgn = tuple(sorted_sgn[i] for i in final_indices)
    final_ell = tuple(sorted_ell[i] for i in final_indices)
    final_u = tuple(sorted_u[i] for i in final_indices)

    return final_k, final_sgn, final_ell, final_u


def canonicalize(genome: Genome) -> Genome:
    """Aplica la reducción de gauge canónico a un genoma (idempotente).

    1. Ordena los clusters por fase angular decodificada.
    2. Rota el genoma para anclar el cluster de mayor tamaño a k=0.
    3. Construye el genoma reflejado (phi -> -phi) y toma el mínimo lexicográfico.
    """
    m = genome.m
    K_ref = genome.K_ref

    # 1 y 2. Representación hacia adelante
    k_fwd, sgn_fwd, ell_fwd, u_fwd = _normalize_clusters(
        m, K_ref, list(genome.k), list(genome.sgn), list(genome.ell), list(genome.u)
    )

    # 3. Representación reflejada (phi -> -phi: k -> -k mod K_ref, sgn -> -sgn)
    refl_k = [(-k_val) % K_ref for k_val in genome.k]
    refl_sgn = [-s_val for s_val in genome.sgn]
    k_rev, sgn_rev, ell_rev, u_rev = _normalize_clusters(
        m, K_ref, refl_k, refl_sgn, list(genome.ell), list(genome.u)
    )

    # Comparación lexicográfica sobre tuplas de parámetros
    ell_fwd_round = tuple(round(v, 9) for v in ell_fwd)
    u_fwd_round = tuple(round(v, 9) for v in u_fwd)
    fwd_key = (k_fwd, sgn_fwd, ell_fwd_round, u_fwd_round)

    ell_rev_round = tuple(round(v, 9) for v in ell_rev)
    u_rev_round = tuple(round(v, 9) for v in u_rev)
    rev_key = (k_rev, sgn_rev, ell_rev_round, u_rev_round)

    if fwd_key <= rev_key:
        return Genome(m=m, K_ref=K_ref, k=k_fwd, sgn=sgn_fwd, ell=ell_fwd, u=u_fwd)
    else:
        return Genome(m=m, K_ref=K_ref, k=k_rev, sgn=sgn_rev, ell=ell_rev, u=u_rev)
