"""Benchmark evaluation time without cache by m in [2, 10]."""

import math
import time
import numpy as np

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def main():
    config = EvalConfig(mu_t=0.6929077, rho=1.0, solver_tol=1e-8)
    evaluator = SDPEvaluator(config=config)
    rng = np.random.default_rng(seed=42)

    n_samples = 10
    results = {}

    print(f"{'m':>3} | {'Samples':>7} | {'Mean (ms)':>10} | {'Std (ms)':>10} | {'Min (ms)':>10} | {'Max (ms)':>10}")
    print("-" * 65)

    for m in range(2, 11):
        times_ms = []
        for _ in range(n_samples):
            # Generate random distinct (x, phi)
            raw_x = rng.exponential(scale=1.0, size=m)
            x_arr = np.maximum(raw_x / np.sum(raw_x), 2.0 * config.eps_c)
            x_arr = x_arr / np.sum(x_arr)
            x = tuple(float(v) for v in x_arr)
            phi = tuple(float(v) for v in rng.uniform(-math.pi, math.pi, size=m))

            t0 = time.perf_counter()
            res = evaluator.evaluate_phase1(x, phi)
            t1 = time.perf_counter()

            times_ms.append((t1 - t0) * 1000.0)

        mean_t = float(np.mean(times_ms))
        std_t = float(np.std(times_ms))
        min_t = float(np.min(times_ms))
        max_t = float(np.max(times_ms))
        results[m] = {
            "mean_ms": mean_t,
            "std_ms": std_t,
            "min_ms": min_t,
            "max_ms": max_t,
        }
        print(f"{m:>3} | {n_samples:>7} | {mean_t:>10.2f} | {std_t:>10.2f} | {min_t:>10.2f} | {max_t:>10.2f}")

    return results


if __name__ == "__main__":
    main()
