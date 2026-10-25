"""Final Fair Comparison based on absolute physical margin (unscaled).

Computes Phase II mu* at fixed absolute physical margins:
  lam_min(M) = 1e-10
  lam_min(M) = 1e-9

Using dynamically scaled floor margin:
  lam_floor_hat = lam_min_M / (h_eff**2)
Evaluator: ScaledSDPEvaluator with CLARABEL and solver_tol=1e-9.
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    compute_h_eff,
    compute_sheets,
)

REPORT_PATH = Path("results/overnight_12h_verified_report.json")
with REPORT_PATH.open("r", encoding="utf-8") as f:
    report_data = json.load(f)

# 1. Load the 3 seeded overnight champions
runs = report_data["runs"]
champ_r1 = runs[0]  # ID 125678
champ_r2 = runs[1]  # ID 2
champ_r3 = runs[2]  # ID 4

# 2. Generate analytical S7 template (Appendix A)
s7_template = get_s7_appendix_a_template(h=0.002, t=0.0)

candidates = [
    {
        "name": f"Campeón Corrida 1 (ID {champ_r1['id']}, mu_t=0.6923)",
        "x": tuple(champ_r1["phenotype_x"]),
        "phi": tuple(champ_r1["phenotype_phi"]),
    },
    {
        "name": f"Campeón Corrida 2 (ID {champ_r2['id']}, mu_t=0.6925)",
        "x": tuple(champ_r2["phenotype_x"]),
        "phi": tuple(champ_r2["phenotype_phi"]),
    },
    {
        "name": f"Campeón Corrida 3 (ID {champ_r3['id']}, mu_t=0.6926)",
        "x": tuple(champ_r3["phenotype_x"]),
        "phi": tuple(champ_r3["phenotype_phi"]),
    },
    {
        "name": "Semilla S7 Analítica (Apéndice A, h=0.002, t=0)",
        "x": s7_template.x,
        "phi": s7_template.phi,
    },
]

# 3. Setup ScaledSDPEvaluator (Clarabel, solver_tol=1e-9)
eval_cfg = EvalConfig(
    solver_name="Clarabel",
    solver_tol=1e-9,
    h_eff_min=1e-4,
)
evaluator = ScaledSDPEvaluator(eval_cfg, h_eff_min=1e-4)

# Targets for absolute physical margin on M:
lam_targets = [1e-10, 1e-9]

results = []
for cand in candidates:
    x_arr = np.asarray(cand["x"], dtype=float)
    phi_arr = np.asarray(cand["phi"], dtype=float)
    sheets, _ = compute_sheets(x_arr, phi_arr)
    heff = float(compute_h_eff(x_arr, phi_arr, sheets))

    row = {
        "Candidato": cand["name"],
        "h_eff": heff,
    }

    for lam_target in lam_targets:
        floor_hat = lam_target / (heff**2)
        res_p2 = evaluator.evaluate_phase2(cand["x"], cand["phi"], eps=floor_hat)
        mu_star = res_p2.mu if res_p2.mu is not None else float("nan")
        row[f"mu_{lam_target}"] = mu_star

    results.append(row)

# 4. Format and print clean Markdown Table
print("\n### Tabla de Comparación Justa por Margen Físico Absoluto en $M$\n")
print("| Candidato | $h_{\\text{eff}}$ | $\\mu^* (\\lambda_M = 10^{-10})$ | $\\mu^* (\\lambda_M = 10^{-9})$ |")
print("| :--- | :---: | :---: | :---: |")
for r in results:
    c_name = r["Candidato"]
    heff_str = f"{r['h_eff']:.6e}"
    mu_10 = f"{r['mu_1e-10']:.8f}"
    mu_9 = f"{r['mu_1e-09']:.8f}"
    print(f"| {c_name} | `{heff_str}` | **{mu_10}** | **{mu_9}** |")

# Also save JSON artifact for records
out_path = Path("results/final_fair_comparison.json")
with out_path.open("w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)
print(f"\nResultados estructurados guardados en: {out_path}")
