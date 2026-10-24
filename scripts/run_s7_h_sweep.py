"""Barrido de h in {0.02, 0.01, 0.005, 0.002} con t=0 para la familia S7 (Item 5).

Reporta:
  - lam* (Fase I sin tau)
  - s* y lam*_norm
  - mu* (Fase II con eps=0)
  - Comparación con los valores de la Tabla de O'Keeffe (Sec. IV.C):
      h=0.02: 0.6906
      h=0.01: 0.6918
      h=0.005: 0.6924
      h=0.002: 0.69272
  - Brecha exacta sin explicaciones sin evidencia.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from kuramoto_templates.cases.s7 import get_s7_template
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def run_h_sweep() -> dict[str, Any]:
    h_values = [0.02, 0.01, 0.005, 0.002]
    table_cited = {
        0.02: 0.6906,
        0.01: 0.6918,
        0.005: 0.6924,
        0.002: 0.69272,
    }

    config = EvalConfig(mu_t=0.6929077, rho=1.0, solver_tol=1e-8)
    evaluator = SDPEvaluator(config=config)

    sweep_results: list[dict[str, Any]] = []

    print("=" * 80)
    print("BARRIDO S7 EN h in {0.02, 0.01, 0.005, 0.002} (t=0, eps=0.0)")
    print("=" * 80)

    for h in h_values:
        template = get_s7_template(h=h, t=0.0)

        # 1. Fase I sin tau
        res_p1 = evaluator.evaluate_phase1(template.x, template.phi)
        lam_star = res_p1.lam
        s_star = res_p1.s
        lam_norm = res_p1.lam_normalized

        # 2. Fase II con eps=0.0
        res_p2 = evaluator.evaluate_phase2(template.x, template.phi, eps=0.0)
        mu_sdp = res_p2.mu

        cited_mu = table_cited[h]
        gap = (mu_sdp - cited_mu) if mu_sdp is not None else None

        item: dict[str, Any] = {
            "h": h,
            "t": 0.0,
            "cited_mu_preprint": cited_mu,
            "phase1": {
                "lam_star": lam_star,
                "s_star": s_star,
                "lam_normalized": lam_norm,
                "status": res_p1.status,
                "max_m_basis_norm": res_p1.solver_info.get("max_m_basis_norm"),
            },
            "phase2_eps_0": {
                "mu_star": mu_sdp,
                "status": res_p2.status,
                "gap_vs_cited": gap,
            },
        }
        sweep_results.append(item)

        gap_str = f"{gap:+.6e}" if gap is not None else "N/A"
        mu_str = f"{mu_sdp:.6f}" if mu_sdp is not None else "N/A"
        lam_str = f"{lam_star:+.6e}" if lam_star is not None else "N/A"
        s_str = f"{s_star:.2e}" if s_star is not None else "N/A"
        print(
            f"h={h:6.4f} | lam*={lam_str} | s*={s_str} | "
            f"mu*(eps=0)={mu_str} | citado={cited_mu:.5f} | brecha={gap_str}"
        )

    out_data = {
        "description": (
            "Barrido de h para S7 a t=0 comparando Fase I (lam*, s*) y "
            "Fase II (mu* eps=0) vs preprint Table"
        ),
        "results": sweep_results,
    }

    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)
    json_path = results_dir / "s7_h_sweep.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out_data, f, indent=2)

    print(f"\nResultados guardados exitosamente en: {json_path}")
    return out_data


if __name__ == "__main__":
    run_h_sweep()
