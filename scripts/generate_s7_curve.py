"""Script para calcular y persistir la curva de Fase I y Fase II para S7."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from kuramoto_templates.cases.s7 import compute_q_star_and_w_star, get_s7_template
from kuramoto_templates.domain.conditions import verify_solution
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def main() -> None:
    t = get_s7_template(h=0.002, t=0.0)
    q_star, _ = compute_q_star_and_w_star(50)
    q_float = float(q_star)
    mu_analytical = (4.0 + 2.0 * q_float) / (7.0 + q_float)
    print(f"S7 Analytical mu: {mu_analytical:.10f}")

    eps_values = [1e-4, 1e-5, 1e-6, 1e-7]
    curve_data: dict[str, Any] = {}

    for eps in eps_values:
        cfg = EvalConfig(eps=eps, solver_tol=1e-8)
        ev = SDPEvaluator(cfg)
        r1 = ev.evaluate_phase1(t.x, t.phi)
        r2 = ev.evaluate_phase2(t.x, t.phi, eps=eps)

        entry: dict[str, Any] = {
            "eps": eps,
            "phase1": {
                "Phi": r1.Phi,
                "status": r1.status,
                "tau": r1.solver_info.get("tau"),
                "delta": r1.solver_info.get("delta"),
                "solve_time": r1.solver_info.get("solve_time_phase1"),
            },
            "phase2": {
                "mu": r2.mu,
                "status": r2.status,
                "solve_time": r2.solver_info.get("solve_time_phase2"),
            },
        }

        if r2.P is not None:
            P_mat = np.array(r2.P)
            v = verify_solution(
                t.x, t.phi, P_mat, eps=eps, eps_c=cfg.eps_c, theta_tol=cfg.theta_tol
            )
            entry["phase2"]["verification"] = {
                "max_F": float(v["max_F"]),
                "min_c": float(v["min_c"]),
                "min_eig_M": float(v["min_eig_M"]),
                "mu": float(v["mu"]),
                "is_valid": bool(v["is_valid"]),
                "reasons": v["reasons"],
            }
        else:
            entry["phase2"]["verification"] = None

        tau_val = r1.solver_info.get("tau")
        delta_val = r1.solver_info.get("delta")
        tau_str = f"{tau_val:.2e}" if isinstance(tau_val, (int | float)) else str(tau_val)
        delta_str = (
            f"{delta_val:.2e}" if isinstance(delta_val, (int | float)) else str(delta_val)
        )
        print(f"\n--- eps = {eps} ---")
        print(
            f"Phase I: Phi = {r1.Phi:.8e}, status = {r1.status}, "
            f"tau = {tau_str}, delta = {delta_str}"
        )
        print(f"Phase II: mu = {r2.mu}, status = {r2.status}")
        if entry["phase2"]["verification"]:
            ver = entry["phase2"]["verification"]
            print(
                f"Verify: max_F = {ver['max_F']:.2e}, min_c = {ver['min_c']:.6f}, "
                f"min_eig_M = {ver['min_eig_M']:.2e}, mu = {ver['mu']:.6f}, "
                f"is_valid = {ver['is_valid']}"
            )

        curve_data[str(eps)] = entry

    # Also evaluate Phase II with eps=0.0 (semidefinite relaxation M >= 0)
    ev0 = SDPEvaluator(EvalConfig(eps=1e-6, solver_tol=1e-8))
    r2_0 = ev0.evaluate_phase2(t.x, t.phi, eps=0.0)
    print("\n--- Phase II with eps = 0.0 (M >= 0) ---")
    print(f"Phase II (eps=0.0): mu = {r2_0.mu}, status = {r2_0.status}")
    if r2_0.P is not None:
        P_mat0 = np.array(r2_0.P)
        v0 = verify_solution(t.x, t.phi, P_mat0, eps=0.0, eps_c=1e-4, theta_tol=1e-4)
        print(
            f"Verify (eps=0.0): max_F = {v0['max_F']:.2e}, min_c = {v0['min_c']:.6f}, "
            f"min_eig_M = {v0['min_eig_M']:.2e}, mu = {v0['mu']:.6f}, "
            f"is_valid = {v0['is_valid']}"
        )
        curve_data["0.0_relaxation"] = {
            "eps": 0.0,
            "phase2": {
                "mu": r2_0.mu,
                "status": r2_0.status,
                "verification": {
                    "max_F": float(v0["max_F"]),
                    "min_c": float(v0["min_c"]),
                    "min_eig_M": float(v0["min_eig_M"]),
                    "mu": float(v0["mu"]),
                    "is_valid": bool(v0["is_valid"]),
                },
            },
        }

    # Save to results/s7_phase1_curve.json
    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)
    out_file = results_dir / "s7_phase1_curve.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(
            {
                "s7_analytical_mu": mu_analytical,
                "curve": curve_data,
            },
            f,
            indent=2,
        )
    print(f"\nSaved curve to {out_file}")


if __name__ == "__main__":
    main()
