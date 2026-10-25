"""Seeds module for S7 (Appendix A) and C4 templates with roundtrip analysis."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from kuramoto_templates.cases.c4 import get_c4_template
from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.search.genome import (
    Genome,
    canonicalize,
    to_phenotype,
    wrap_minus_pi_to_pi,
)
from kuramoto_templates.search.operators.split import _encode_phase_to_genes


@dataclass(frozen=True)
class SeedIndividual:
    """Genetic algorithm seed holding genome and decoded or exact phenotype."""

    genome: Genome
    phenotype_x: tuple[float, ...]
    phenotype_phi: tuple[float, ...]
    name: str
    is_direct_injection: bool = False


def _encode_phenotype_to_genome(
    x: tuple[float, ...],
    phi: tuple[float, ...],
    K_ref: int = 4,
    canonical: bool = True,
) -> Genome:
    """Encode continuous phenotype (x, phi) into Genome representation."""
    m = len(x)
    xmin = 1e-3
    x_arr = np.array(x, dtype=np.float64)
    p = np.maximum(x_arr - xmin, 1e-12) / (1.0 - m * xmin)
    u_arr = np.log(p)
    u_arr = u_arr - np.mean(u_arr)

    k_list: list[int] = []
    sgn_list: list[int] = []
    ell_list: list[float] = []

    for phase in phi:
        k_val, sgn_val, ell_val = _encode_phase_to_genes(phase, K_ref)
        k_list.append(k_val)
        sgn_list.append(sgn_val)
        ell_list.append(ell_val)

    raw_genome = Genome(
        m=m,
        K_ref=K_ref,
        k=tuple(k_list),
        sgn=tuple(sgn_list),
        ell=tuple(ell_list),
        u=tuple(float(v) for v in u_arr),
    )
    return canonicalize(raw_genome) if canonical else raw_genome


def create_c4_seed() -> SeedIndividual:
    """Generate C4 seed template (4 orthogonal/antipodal clusters)."""
    c4_template = get_c4_template()
    genome = _encode_phenotype_to_genome(
        c4_template.x,
        c4_template.phi,
        K_ref=4,
    )
    x_dec, phi_dec = to_phenotype(genome)
    return SeedIndividual(
        genome=genome,
        phenotype_x=x_dec,
        phenotype_phi=phi_dec,
        name="C4",
        is_direct_injection=False,
    )


def create_s7_direct_seed(h: float = 0.002, t: float = 0.0) -> SeedIndividual:
    """Generate exact S7 seed from Appendix A with direct (x, phi) injection.

    Direct injection is required because cluster W has exactly zero deviation (offset = 0)
    which is not representable under ell in [-6.0, -0.5].
    """
    s7_template = get_s7_appendix_a_template(h=h, t=t)
    genome = _encode_phenotype_to_genome(
        s7_template.x,
        s7_template.phi,
        K_ref=4,
        canonical=False,
    )
    return SeedIndividual(
        genome=genome,
        phenotype_x=s7_template.x,
        phenotype_phi=s7_template.phi,
        name=f"S7_direct_h{h}_t{t}",
        is_direct_injection=True,
    )


def create_s7_genome_seed(h: float = 0.002, t: float = 0.0) -> SeedIndividual:
    """Generate S7 seed encoded and decoded strictly through Genome representation."""
    s7_template = get_s7_appendix_a_template(h=h, t=t)
    genome = _encode_phenotype_to_genome(
        s7_template.x,
        s7_template.phi,
        K_ref=4,
    )
    x_dec, phi_dec = to_phenotype(genome)
    return SeedIndividual(
        genome=genome,
        phenotype_x=x_dec,
        phenotype_phi=phi_dec,
        name="S7_genome",
        is_direct_injection=False,
    )


def create_champion_3178_seed() -> SeedIndividual:
    """Champion from with_seeds_mu_0.6920_seed_42 gen 50 (ID 3178, mu*=0.692569)."""
    genome = Genome(
        m=7,
        K_ref=4,
        k=(2, 3, 3, 0, 1, 1, 2),
        sgn=(-1, 1, -1, 1, 1, -1, 1),
        ell=(
            -2.854406433700391,
            -3.000000000000048,
            -3.000000000000048,
            -6.0,
            -3.000000000000048,
            -3.000000000000048,
            -2.854406433700391,
        ),
        u=(
            -0.06968680856854073,
            -0.06810554678401977,
            -0.07273613388036303,
            0.43449156693986546,
            -0.07034086205272089,
            -0.06991127518092961,
            -0.06850848195733376,
        ),
    )
    phenotype_x = (
        0.1307282343595373,
        0.13093353092993723,
        0.13033325328979067,
        0.2157812638662674,
        0.13064341289776796,
        0.13069911797019027,
        0.1308811866865093,
    )
    phenotype_phi = (
        3.1401943754557724,
        -1.5697963267948971,
        -1.571796326794896,
        1.000000000139778e-06,
        1.571796326794896,
        1.5697963267948971,
        -3.1401943754557724,
    )
    return SeedIndividual(
        genome=genome,
        phenotype_x=phenotype_x,
        phenotype_phi=phenotype_phi,
        name="champion_id3178",
        is_direct_injection=True,
    )



def measure_s7_roundtrip_error(
    h: float = 0.002,
    t: float = 0.0,
    results_dir: Path | str = "results",
) -> dict[str, float]:
    """Measure roundtrip error S7 -> Genome -> S7 and evaluate impact on lam*."""
    s7_exact = get_s7_appendix_a_template(h=h, t=t)
    x_exact = np.array(s7_exact.x)
    phi_exact = np.array(s7_exact.phi)

    # Encode to genome without gauge rotation to measure gene-by-gene reconstruction error
    genome = _encode_phenotype_to_genome(
        s7_exact.x, s7_exact.phi, K_ref=4, canonical=False
    )
    x_dec_tup, phi_dec_tup = to_phenotype(genome)
    x_dec = np.array(x_dec_tup)
    phi_dec = np.array(phi_dec_tup)

    max_mass_error = float(np.max(np.abs(x_exact - x_dec)))

    # Circular difference for phases
    phase_diffs = [
        abs(wrap_minus_pi_to_pi(float(p1 - p2)))
        for p1, p2 in zip(phi_exact, phi_dec, strict=False)
    ]
    max_phase_error = float(max(phase_diffs))

    # Evaluate both using SDPEvaluator at mu_t = 0.690
    evaluator = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-8))
    res_direct = evaluator.evaluate_phase1(s7_exact.x, s7_exact.phi)
    res_genome = evaluator.evaluate_phase1(x_dec_tup, phi_dec_tup)

    lam_direct = float(res_direct.lam) if res_direct.lam is not None else 0.0
    lam_genome = float(res_genome.lam) if res_genome.lam is not None else 0.0
    lam_abs_diff = float(abs(lam_direct - lam_genome))

    metrics: dict[str, float] = {
        "h": h,
        "t": t,
        "max_mass_error": max_mass_error,
        "max_phase_error": max_phase_error,
        "s_star_direct": float(res_direct.s) if res_direct.s is not None else 0.0,
        "lam_direct": lam_direct,
        "s_star_genome": float(res_genome.s) if res_genome.s is not None else 0.0,
        "lam_genome": lam_genome,
        "lam_abs_diff": lam_abs_diff,
    }

    # Save to results directory if exists
    out_dir = Path(results_dir)
    if out_dir.exists():
        with (out_dir / "s7_roundtrip_error.json").open("w", encoding="utf-8") as f:
            json.dump(metrics, f, indent=2)

    return metrics
