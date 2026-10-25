"""Diagnose step 0a and 0b as requested."""

import math
import sys
from pathlib import Path
import numpy as np
import mpmath

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.cases.s7 import (
    get_s7_template,
    get_s7_appendix_a_template,
    compute_q_star_and_w_star,
)
from kuramoto_templates.search.seeds import create_s7_direct_seed
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.domain.conditions import verify_solution


def run_0a():
    print("=" * 70)
    print("PASO 0a: INCONSISTENCIA S7 (h=0.02, t=0)")
    print("=" * 70)
    h = 0.02
    t = 0.0

    seed_direct = create_s7_direct_seed(h=h, t=t)
    tmpl_s7 = get_s7_template(h=h, t=t)
    tmpl_app_a = get_s7_appendix_a_template(h=h, t=t)

    # 1. Compare (x, phi)
    diff_x = np.max(np.abs(np.array(seed_direct.phenotype_x) - np.array(tmpl_s7.x)))
    diff_phi = np.max(np.abs(np.array(seed_direct.phenotype_phi) - np.array(tmpl_s7.phi)))

    print(f"(i) Comparación create_s7_direct_seed vs get_s7_template:")
    print(f"    diff_max(x):   {diff_x:.16e}")
    print(f"    diff_max(phi): {diff_phi:.16e}")
    print(f"    ¿Mismos (x, phi)?: {diff_x == 0.0 and diff_phi == 0.0}")

    # 2. Evaluate Phase 1 and Phase 2
    ev = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-8, lam_floor=1e-7))
    res_p1 = ev.evaluate_phase1(seed_direct.phenotype_x, seed_direct.phenotype_phi)
    res_p2 = ev.evaluate_phase2(tmpl_s7.x, tmpl_s7.phi, eps=0.0)

    print(f"\n(ii) Comparación Fase 1 vs Fase 2:")
    print(f"    Fase 1 (mu_t=0.690):")
    print(f"      s*:             {res_p1.s:.16e}" if res_p1.s is not None else "      s*: None")
    print(f"      lam*:           {res_p1.lam:.16e}" if res_p1.lam is not None else "      lam*: None")
    print(f"      status:         {res_p1.status}")
    if res_p1.P is not None:
        deg_p1 = np.array(res_p1.P) @ np.array(seed_direct.phenotype_x)
        mu_p1 = float(np.min(deg_p1))
        print(f"      mu real de P1:  {mu_p1:.16e}")

    print(f"    Fase 2 (eps=0.0):")
    print(f"      mu*:            {res_p2.mu:.16e}" if res_p2.mu is not None else "      mu*: None")
    print(f"      Phi (Phase 1):  {res_p2.Phi:.16e}" if res_p2.Phi is not None else "      Phi: None")
    print(f"      status:         {res_p2.status}")
    print(f"      solver_info:    {res_p2.solver_info}")
    if res_p2.P is not None:
        deg_p2 = np.array(res_p2.P) @ np.array(tmpl_s7.x)
        mu_p2 = float(np.min(deg_p2))
        print(f"      mu real de P2:  {mu_p2:.16e}")

    print(f"\n    Appendix A P analytic degrees:")
    if tmpl_app_a.P is not None:
        deg_app_a = np.array(tmpl_app_a.P) @ np.array(tmpl_app_a.x)
        print(f"      mu analitico:   {float(np.min(deg_app_a)):.16e}")

    # Check for discrepancy: mu* < 0.690 pero s* ≈ 0
    mu_p2_val = res_p2.mu if res_p2.mu is not None else 0.0
    s_val = res_p1.s if res_p1.s is not None else 0.0
    if mu_p2_val < 0.690 and s_val < 1e-4:
        print(f"\n    [ALERTA BUG DETECTADO]: mu* (Fase 2) = {mu_p2_val} < 0.690, pero s* (Fase 1) = {s_val} ≈ 0!")


if __name__ == "__main__":
    run_0a()
