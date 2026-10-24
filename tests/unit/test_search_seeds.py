"""Unit tests for seeds generation, S7 roundtrip measurement, and direct injection."""

import math

from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.fitness import compute_continuous_fitness
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.search.seeds import (
    create_c4_seed,
    create_s7_direct_seed,
    measure_s7_roundtrip_error,
)


def test_c4_seed_generation() -> None:
    """C4 seed produces a valid 4-cluster genome."""
    seed_c4 = create_c4_seed()
    assert seed_c4.genome.m == 4
    assert len(seed_c4.phenotype_x) == 4
    assert len(seed_c4.phenotype_phi) == 4
    assert math.isclose(float(sum(seed_c4.phenotype_x)), 1.0, rel_tol=1e-7)


def test_s7_seed_and_roundtrip_measurement() -> None:
    """S7 seed roundtrip error must be measured and cluster W discrepancy identified."""
    metrics = measure_s7_roundtrip_error(h=0.002, t=0.0)

    assert "max_mass_error" in metrics
    assert "max_phase_error" in metrics
    assert "lam_direct" in metrics
    assert "lam_genome" in metrics
    assert "lam_abs_diff" in metrics

    # The phase error should be around 1e-6 due to cluster W having 0 deviation vs ell >= -6
    assert metrics["max_phase_error"] <= 1.5e-6
    assert metrics["max_phase_error"] >= 5e-7


def test_s7_acceptance_criterion_b() -> None:
    """CONTROL S7 BAJO PISO ESPECTRAL (Item 4):

    Semilla S7 evaluada sola con mu_t=0.690 y lam_floor=1e-7:
    Dado que lam* = 4.59e-9 < 1e-7, el fitness queda negativo f = lam* - 1e-7 < 0.
    """
    s7_seed = create_s7_direct_seed(h=0.002, t=0.0)
    evaluator = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-8, lam_floor=1e-7))

    res = evaluator.evaluate_phase1(s7_seed.phenotype_x, s7_seed.phenotype_phi)
    assert res.s is not None and res.lam is not None
    f = compute_continuous_fitness(lam=res.lam, s=res.s, tol=1e-8, lam_floor=1e-7)

    print("\nS7 MEASUREMENT WITH lam_floor=1e-7:")
    print(f"  status: {res.status}")
    print(f"  s_star: {res.s:.12e}")
    print(f"  lam_star: {res.lam:.12e}")
    print(f"  fitness f: {f:.12e}")

    # S7 da 4.59e-9 < 1e-7: queda estrictamente negativo con el piso
    assert f < 0.0, f"S7 debió dar f < 0 con lam_floor=1e-7, dio f = {f}"
