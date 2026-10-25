"""Diagnose why verify_solution reports is_valid=False for overnight champions."""

from __future__ import annotations

import json
from kuramoto_templates.domain.conditions import verify_solution
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import ScaledSDPEvaluator

with open("results/overnight_12h_verified_report.json", "r", encoding="utf-8") as f:
    report = json.load(f)

for run in report["runs"]:
    name = run["run_name"]
    cid = run.get("id")
    print("=" * 80)
    print(f"Corrida: {name} (ID: {cid})")
    x = tuple(run["phenotype_x"])
    phi = tuple(run["phenotype_phi"])
    heff = run.get("h_eff", 1e-3)

    ev = ScaledSDPEvaluator(
        EvalConfig(
            mu_t=0.690,
            solver_tol=1e-9,
            lam_floor_hat=1e-5,
            h_eff_min=heff,
        )
    )
    res = ev.evaluate_phase1(x, phi)
    if res.P is not None:
        v_default = verify_solution(x, phi, res.P, eps=1e-6, eps_c=1e-4, theta_tol=1e-4)
        print("Parametros y umbrales por defecto en verify_solution:")
        print("  Umbral theta_tol: 1e-4")
        print("  Umbral eps_c    : 1e-4")
        print("  Umbral eps (M)  : 1e-6")
        print("Valores calculados:")
        print(f"  max_F           : {v_default['max_F']:.6e}  (pasa si <= 1e-4)")
        print(f"  min_c           : {v_default['min_c']:.6e}  (pasa si >= 1e-4)")
        print(f"  min_eig_M       : {v_default['min_eig_M']:.6e}  (pasa si >= 1e-6)")
        print(f"  p_bounds_ok     : {v_default['p_bounds_ok']}  (pasa si en [-1e-8, 1+1e-8])")
        print(f"  is_valid        : {v_default['is_valid']}")
        print(f"  reasons         : {v_default.get('reasons', [])}")

        eps_phys = 1e-5 * (heff**2)
        v_scaled = verify_solution(x, phi, res.P, eps=eps_phys)
        print(f"Con umbral fisico de escala eps = lam_floor_hat * h_eff^2 = {eps_phys:.6e}:")
        print(f"  is_valid        : {v_scaled['is_valid']}")
        print(f"  reasons         : {v_scaled.get('reasons', [])}")
