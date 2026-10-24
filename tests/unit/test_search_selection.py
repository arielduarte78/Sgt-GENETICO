"""Unit tests for selection with fitness noise, tie-breaking, and deterministic crowding."""

import numpy as np

from kuramoto_templates.search.genome import Genome
from kuramoto_templates.search.selection import (
    Individual,
    compare_individuals,
    compute_f_sel,
    tournament_selection,
)


def _make_ind(
    ind_id: int,
    m: int,
    s_star: float,
    lam_star: float,
    lam_norm: float,
    tol: float = 1e-8,
) -> Individual:
    genome = Genome(
        m=m,
        K_ref=4,
        k=tuple(i % 4 for i in range(m)),
        sgn=tuple([1] * m),
        ell=tuple([-2.0] * m),
        u=tuple([0.0] * m),
    )
    f_sel = compute_f_sel(s_star, lam_norm, tol=tol)
    return Individual(
        id=ind_id,
        genome=genome,
        phenotype_x=tuple([1.0 / m] * m),
        phenotype_phi=tuple([0.0] * m),
        s_star=s_star,
        lam_star=lam_star,
        lam_normalized=lam_norm,
        f=f_sel,
        f_sel=f_sel,
        parent_ids=(),
        operator_name="init",
    )


def test_compute_f_sel() -> None:
    """f_sel follows the 3-tier ordered fitness hierarchy."""
    # Slack > tol, stable
    assert compute_f_sel(s_star=1e-4, lam_norm=0.5, tol=1e-8, lam_floor=1e-7) == -1e-4
    # Slack <= tol, stable
    assert compute_f_sel(s_star=1e-9, lam_norm=0.5, tol=1e-8, lam_floor=1e-7) == 0.5
    # Unstable state (lam_norm < lam_floor): f_sel <= -1.0
    f_unstable = compute_f_sel(s_star=1e-9, lam_norm=-0.33, tol=1e-8, lam_floor=1e-5)
    assert f_unstable == -1.33
    # S7 (stable, s* > tol) beats unstable state with tiny slack
    f_s7 = compute_f_sel(s_star=1.8e-3, lam_norm=1.36e-4, tol=1e-8, lam_floor=1e-5)
    assert f_s7 == -1.8e-3
    assert f_s7 > f_unstable



def test_compare_individuals_noise_and_tie_breaking() -> None:
    """When |f_sel_i - f_sel_j| < fitness_noise, break ties by smaller m, then smaller id."""
    # Distinct fitness: larger f_sel wins
    ind1 = _make_ind(ind_id=1, m=5, s_star=0.0, lam_star=0.1, lam_norm=0.1)
    ind2 = _make_ind(ind_id=2, m=4, s_star=0.0, lam_star=0.2, lam_norm=0.2)
    assert compare_individuals(ind2, ind1, fitness_noise=1e-9) > 0  # ind2 wins

    # Indistinguishable fitness (|diff| < 1e-9): smaller m wins
    ind3 = _make_ind(ind_id=3, m=4, s_star=0.0, lam_star=0.1, lam_norm=0.1 + 1e-11)
    ind4 = _make_ind(ind_id=4, m=5, s_star=0.0, lam_star=0.1, lam_norm=0.1)
    # ind3 has m=4 < ind4 m=5
    assert compare_individuals(ind3, ind4, fitness_noise=1e-9) > 0  # ind3 wins

    # Indistinguishable fitness and same m: smaller id wins
    ind5 = _make_ind(ind_id=5, m=4, s_star=0.0, lam_star=0.1, lam_norm=0.1)
    ind6 = _make_ind(ind_id=6, m=4, s_star=0.0, lam_star=0.1, lam_norm=0.1 + 1e-11)
    assert compare_individuals(ind5, ind6, fitness_noise=1e-9) > 0  # ind5 wins (id 5 < 6)


def test_tournament_selection() -> None:
    """Tournament of 3 returns the best individual according to comparison rule."""
    ind1 = _make_ind(ind_id=1, m=6, s_star=0.0, lam_star=0.05, lam_norm=0.05)
    ind2 = _make_ind(ind_id=2, m=5, s_star=0.0, lam_star=0.10, lam_norm=0.10)
    ind3 = _make_ind(ind_id=3, m=4, s_star=0.0, lam_star=0.02, lam_norm=0.02)

    pop = [ind1, ind2, ind3]
    rng = np.random.default_rng(123)

    winner = tournament_selection(pop, k=3, fitness_noise=1e-9, rng=rng)
    assert winner.id == 2  # ind2 has highest f_sel
