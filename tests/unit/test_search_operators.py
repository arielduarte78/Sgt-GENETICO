"""Property-based tests with Hypothesis for search genetic operators."""

import math

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from kuramoto_templates.search.genome import (
    Genome,
    to_phenotype,
)
from kuramoto_templates.search.operators.anchor_jump import anchor_jump
from kuramoto_templates.search.operators.crossover import crossover
from kuramoto_templates.search.operators.dispatcher import (
    OperatorRegistry,
    select_and_apply_operator,
)
from kuramoto_templates.search.operators.gauss_u import gauss_u
from kuramoto_templates.search.operators.merge import merge
from kuramoto_templates.search.operators.phase_dev import phase_dev
from kuramoto_templates.search.operators.split import split


@st.composite
def genome_strategy(draw: st.DrawFn, min_m: int = 2, max_m: int = 10) -> Genome:
    """Hypothesis strategy generating valid Genomes."""
    m = draw(st.integers(min_value=min_m, max_value=max_m))
    K_ref = draw(st.sampled_from([4, 8, 12, 16]))
    k = draw(
        st.lists(st.integers(min_value=0, max_value=K_ref - 1), min_size=m, max_size=m)
    )
    sgn = draw(st.lists(st.sampled_from([-1, 0, 1]), min_size=m, max_size=m))
    ell = draw(
        st.lists(
            st.floats(min_value=-6.0, max_value=-0.5, allow_nan=False),
            min_size=m,
            max_size=m,
        )
    )
    u = draw(
        st.lists(
            st.floats(min_value=-5.0, max_value=5.0, allow_nan=False),
            min_size=m,
            max_size=m,
        )
    )
    return Genome(
        m=m,
        K_ref=K_ref,
        k=tuple(k),
        sgn=tuple(sgn),
        ell=tuple(ell),
        u=tuple(u),
    )


def assert_valid_phenotype(genome: Genome) -> None:
    """Verify that decoded phenotype respects domain constraints."""
    assert 2 <= genome.m <= 10
    x, phi = to_phenotype(genome)
    assert len(x) == genome.m
    assert len(phi) == genome.m
    # Simplex: positive masses and sum == 1
    assert all(val >= 1e-3 - 1e-12 for val in x)
    assert math.isclose(float(sum(x)), 1.0, rel_tol=1e-6, abs_tol=1e-6)
    # Phases in (-pi, pi]
    assert all(val > -math.pi - 1e-12 for val in phi)
    assert all(val <= math.pi + 1e-12 for val in phi)


def assert_genome_unmodified(
    original: Genome, before_k: tuple[int, ...], before_u: tuple[float, ...]
) -> None:
    """Ensure original genome was not mutated in-place."""
    assert original.k == before_k
    assert original.u == before_u


# --- Tests for immutability and validity of gauss_u ---


@given(genome=genome_strategy(), seed=st.integers(min_value=0, max_value=2**31 - 1))
@settings(max_examples=30)
def test_gauss_u_properties(genome: Genome, seed: int) -> None:
    rng = np.random.default_rng(seed)
    before_k = genome.k
    before_u = genome.u

    mutant = gauss_u(genome, rng)

    assert_genome_unmodified(genome, before_k, before_u)
    assert mutant.m == genome.m
    assert_valid_phenotype(mutant)
    # Phases should remain unchanged
    _, orig_phi = to_phenotype(genome)
    _, mut_phi = to_phenotype(mutant)
    assert np.allclose(orig_phi, mut_phi)


# --- Tests for phase_dev ---


@given(genome=genome_strategy(), seed=st.integers(min_value=0, max_value=2**31 - 1))
@settings(max_examples=30)
def test_phase_dev_properties(genome: Genome, seed: int) -> None:
    rng = np.random.default_rng(seed)
    before_k = genome.k
    before_u = genome.u

    mutant = phase_dev(genome, rng)

    assert_genome_unmodified(genome, before_k, before_u)
    assert mutant.m == genome.m
    assert all(val >= -6.0 - 1e-12 for val in mutant.ell)
    assert all(val <= -0.5 + 1e-12 for val in mutant.ell)
    assert all(val in (-1, 0, 1) for val in mutant.sgn)
    assert_valid_phenotype(mutant)


# --- Tests for anchor_jump ---


@given(genome=genome_strategy(), seed=st.integers(min_value=0, max_value=2**31 - 1))
@settings(max_examples=30)
def test_anchor_jump_properties(genome: Genome, seed: int) -> None:
    rng = np.random.default_rng(seed)
    before_k = genome.k
    before_u = genome.u

    mutant = anchor_jump(genome, rng)

    assert_genome_unmodified(genome, before_k, before_u)
    assert mutant.m == genome.m
    assert all(0 <= val < mutant.K_ref for val in mutant.k)
    assert_valid_phenotype(mutant)


# --- Tests for split ---


@given(
    genome=genome_strategy(min_m=2, max_m=10),
    seed=st.integers(min_value=0, max_value=2**31 - 1),
)
@settings(max_examples=30)
def test_split_properties(genome: Genome, seed: int) -> None:
    rng = np.random.default_rng(seed)
    before_k = genome.k
    before_u = genome.u

    result = split(genome, rng)

    assert_genome_unmodified(genome, before_k, before_u)
    assert_valid_phenotype(result)

    if genome.m < 10:
        assert result.m == genome.m + 1
        # Mass must be conserved
        x_orig, _ = to_phenotype(genome)
        x_res, _ = to_phenotype(result)
        assert math.isclose(float(sum(x_orig)), float(sum(x_res)), rel_tol=1e-6)
    else:
        assert result.m == 10


# --- Tests for merge ---


@given(
    genome=genome_strategy(min_m=2, max_m=10),
    seed=st.integers(min_value=0, max_value=2**31 - 1),
)
@settings(max_examples=30)
def test_merge_properties(genome: Genome, seed: int) -> None:
    rng = np.random.default_rng(seed)
    before_k = genome.k
    before_u = genome.u

    result = merge(genome, rng)

    assert_genome_unmodified(genome, before_k, before_u)
    assert_valid_phenotype(result)

    if genome.m > 2:
        assert result.m == genome.m - 1
        # Mass must be conserved
        x_orig, _ = to_phenotype(genome)
        x_res, _ = to_phenotype(result)
        assert math.isclose(float(sum(x_orig)), float(sum(x_res)), rel_tol=1e-6)
    else:
        assert result.m == 2


# --- Tests for crossover ---


@given(
    g1=genome_strategy(),
    g2=genome_strategy(),
    seed=st.integers(min_value=0, max_value=2**31 - 1),
)
@settings(max_examples=30)
def test_crossover_properties(g1: Genome, g2: Genome, seed: int) -> None:
    rng = np.random.default_rng(seed)
    before_k1 = g1.k
    before_u1 = g1.u
    before_k2 = g2.k
    before_u2 = g2.u

    child = crossover(g1, g2, rng)

    assert_genome_unmodified(g1, before_k1, before_u1)
    assert_genome_unmodified(g2, before_k2, before_u2)
    assert 2 <= child.m <= 10
    assert_valid_phenotype(child)


# --- Tests for operator dispatcher (single rng.choice call) ---


def test_operator_dispatcher_uses_single_rng_choice() -> None:
    """Dispatcher must select operator using exactly one call to rng.choice."""
    genome = Genome(
        m=4,
        K_ref=8,
        k=(0, 2, 4, 6),
        sgn=(1, 1, 1, 1),
        ell=(-2.0, -2.0, -2.0, -2.0),
        u=(0.0, 0.0, 0.0, 0.0),
    )
    base_rng = np.random.default_rng(42)

    class TrackingRNG:
        def __init__(self, base: np.random.Generator) -> None:
            self._base = base
            self.choice_calls = 0

        def choice(self, *args: object, **kwargs: object) -> object:
            self.choice_calls += 1
            return getattr(self._base, "choice")(*args, **kwargs)

        def __getattr__(self, item: str) -> object:
            return getattr(self._base, item)

    rng = TrackingRNG(base_rng)

    registry = OperatorRegistry()
    child, op_name = select_and_apply_operator(genome, genome, registry, rng)  # type: ignore[arg-type]

    assert rng.choice_calls == 1
    assert op_name in registry.probabilities
    assert_valid_phenotype(child)
