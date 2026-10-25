"""Benchmark T9: Tiempo medio por evaluación nueva (sin caché), por m, plano vs escalado."""

from __future__ import annotations

import sys
import time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.evaluators.scaled import ScaledSDPEvaluator
from kuramoto_templates.search.genome import Genome, canonicalize, decode_x, decode_phi


def run_benchmark_t9():
    m_values = [4, 5, 6, 7, 8]
    num_samples = 15
    rng = np.random.default_rng(999)

    cfg = EvalConfig(mu_t=0.690, solver_tol=1e-8)
    ev_plain = SDPEvaluator(config=cfg)
    ev_scaled = ScaledSDPEvaluator(config=cfg)

    results = {}

    print("=" * 70)
    print("BENCHMARK T9: TIEMPO MEDIO POR EVALUACIÓN NUEVA (SIN CACHÉ)")
    print(f"Muestras por m: {num_samples}")
    print("=" * 70)

    for m in m_values:
        # Generate random individuals near C4
        candidates = []
        for _ in range(num_samples):
            kref = 4
            k = tuple(int(i % kref) for i in range(m))
            sgn = tuple(int(rng.choice([-1, 1])) for _ in range(m))
            ell = tuple(float(rng.uniform(-3.0, -1.3)) for _ in range(m))
            u = tuple(float(rng.normal(0.0, 0.2)) for _ in range(m))
            genome = canonicalize(Genome(m=m, K_ref=kref, k=k, sgn=sgn, ell=ell, u=u))
            x = decode_x(genome.u)
            phi = decode_phi(genome.K_ref, genome.k, genome.sgn, genome.ell)
            candidates.append((x, phi, genome))

        # 1. Plain
        t0 = time.perf_counter()
        for x, phi, _ in candidates:
            _ = ev_plain.evaluate_phase1(x, phi)
        time_plain = (time.perf_counter() - t0) / num_samples

        # 2. Scaled
        t0 = time.perf_counter()
        for x, phi, g in candidates:
            _ = ev_scaled.evaluate_phase1(x, phi, genome=g)
        time_scaled = (time.perf_counter() - t0) / num_samples

        speedup = time_plain / time_scaled if time_scaled > 0 else 1.0
        results[m] = {
            "plain_sec": time_plain,
            "scaled_sec": time_scaled,
            "speedup": speedup,
        }
        print(
            f"m={m:2d} | Plano: {time_plain*1000:6.2f} ms | Escalado: {time_scaled*1000:6.2f} ms | Ratio: {speedup:5.2f}x"
        )

    print("=" * 70)
    return results


if __name__ == "__main__":
    run_benchmark_t9()
