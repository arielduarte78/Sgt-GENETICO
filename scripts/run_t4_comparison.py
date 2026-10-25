"""Run T4 Solver Agreement Improvement Benchmark."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import ScaledSDPEvaluator
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def run_t4():
    tmpl = get_s7_appendix_a_template(h=0.002, t=0.0)

    # 1. Plain Evaluator
    lams_plain = {}
    for key, (s_name, tol) in [
        ("CL8", ("CLARABEL", 1e-8)),
        ("CL9", ("CLARABEL", 1e-9)),
        ("SCS", ("SCS", 1e-8)),
    ]:
        ev = SDPEvaluator(EvalConfig(mu_t=0.690, solver_name=s_name, solver_tol=tol))
        res = ev.evaluate_phase1(tmpl.x, tmpl.phi)
        lams_plain[key] = res.lam

    diff_plain_cl = abs(lams_plain["CL8"] - lams_plain["CL9"])
    diff_plain_scs = abs(lams_plain["CL8"] - lams_plain["SCS"])
    max_diff_plain = max(diff_plain_cl, diff_plain_scs)
    rel_diff_plain = max_diff_plain / abs(lams_plain["CL8"])

    # 2. Scaled Evaluator
    lams_scaled = {}
    for key, (s_name, tol) in [
        ("CL8", ("CLARABEL", 1e-8)),
        ("CL9", ("CLARABEL", 1e-9)),
        ("SCS", ("SCS", 1e-8)),
    ]:
        ev = ScaledSDPEvaluator(
            EvalConfig(mu_t=0.690, solver_name=s_name, solver_tol=tol)
        )
        res = ev.evaluate_phase1(tmpl.x, tmpl.phi)
        lams_scaled[key] = res.lam

    diff_scaled_cl = abs(lams_scaled["CL8"] - lams_scaled["CL9"])
    diff_scaled_scs = abs(lams_scaled["CL8"] - lams_scaled["SCS"])
    max_diff_scaled = max(diff_scaled_cl, diff_scaled_scs)
    rel_diff_scaled = max_diff_scaled / abs(lams_scaled["CL8"])

    ratio = rel_diff_plain / rel_diff_scaled
    print("=" * 70)
    print("RESULTADOS TEST T4: MEJORA DE ACUERDO ENTRE SOLVERS")
    print("=" * 70)
    print("EVALUADOR PLANO:")
    print(f"  Clarabel (1e-8): {lams_plain['CL8']:.12e}")
    print(f"  Clarabel (1e-9): {lams_plain['CL9']:.12e}")
    print(f"  SCS (1e-8):      {lams_plain['SCS']:.12e}")
    print(f"  Max diff abs:    {max_diff_plain:.12e}")
    print(f"  Rel diff:        {rel_diff_plain:.12e}")

    print("\nEVALUADOR ESCALADO:")
    print(f"  Clarabel (1e-8): {lams_scaled['CL8']:.12e}")
    print(f"  Clarabel (1e-9): {lams_scaled['CL9']:.12e}")
    print(f"  SCS (1e-8):      {lams_scaled['SCS']:.12e}")
    print(f"  Max diff abs:    {max_diff_scaled:.12e}")
    print(f"  Rel diff:        {rel_diff_scaled:.12e}")

    print(f"\nMEJORA DE ACUERDO RELATIVO (Razón): {ratio:.2f}x")
    print("=" * 70)


if __name__ == "__main__":
    run_t4()
