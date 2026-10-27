"""Cheap empirical experiments to validate/challenge the thesis:

1. S7 parametric sweep: h in {5e-4, 1e-3, 2e-3}
2. mu*(lambda_M) curve for lambda_M in {1e-12, 1e-11, 1e-10, 1e-9, 1e-8}
3. Spectral properties of T: sigma(T), kappa(T), and realized lambda_min(M)
4. Analysis of gap(h) and implications of h_min lock
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import scipy.linalg

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.domain.conditions import (
    build_lw_matrix,
    p_matrix_to_vector,
)
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)

# 1. Load overnight champions
with open("results/overnight_12h_verified_report.json", "r", encoding="utf-8") as f:
    report = json.load(f)

champs = [
    ("Campeon R1 (ID 125678)", tuple(report["runs"][0]["phenotype_x"]), tuple(report["runs"][0]["phenotype_phi"])),
    ("Campeon R2 (ID 2)", tuple(report["runs"][1]["phenotype_x"]), tuple(report["runs"][1]["phenotype_phi"])),
    ("Campeon R3 (ID 4)", tuple(report["runs"][2]["phenotype_x"]), tuple(report["runs"][2]["phenotype_phi"])),
]

# 2. Add S7 for h in {5e-4, 1e-3, 2e-3}
h_values = [5e-4, 1e-3, 2e-3]
s7_cases = []
for h in h_values:
    tmpl = get_s7_appendix_a_template(h=h, t=0.0)
    s7_cases.append((f"S7 (h={h:.1e})", tmpl.x, tmpl.phi))

all_candidates = champs + s7_cases

evaluator = ScaledSDPEvaluator(
    EvalConfig(solver_name="Clarabel", solver_tol=1e-9, h_eff_min=1e-5),
    h_eff_min=1e-5,
)

# Physical lambda_M targets
lam_M_targets = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8]

experiment_results = []

for name, x, phi in all_candidates:
    x_arr = np.asarray(x, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)
    m = len(x_arr)
    sheets, ref_sheet = compute_sheets(x_arr, phi_arr)
    heff = float(compute_h_eff(x_arr, phi_arr, sheets))
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, heff)

    if T is not None:
        sv = scipy.linalg.svdvals(T)
        sigma_max = float(sv[0])
        sigma_min = float(sv[-1])
        kappa_T = float(sigma_max / sigma_min) if sigma_min > 0 else float("inf")
    else:
        sv = []
        sigma_max = sigma_min = kappa_T = float("nan")

    # Sweep across lambda_M
    mu_by_lam_M = {}
    realized_lam_M = {}

    for lam_target in lam_M_targets:
        floor_hat = lam_target / (heff**2)
        res_p2 = evaluator.evaluate_phase2(x, phi, eps=floor_hat)
        mu_val = res_p2.mu if res_p2.mu is not None else float("nan")
        mu_by_lam_M[f"{lam_target:.0e}"] = mu_val

        # Compute realized lambda_min(M) in Phase II solution
        if res_p2.P is not None:
            p_vec = p_matrix_to_vector(np.asarray(res_p2.P, dtype=float))
            L_w = build_lw_matrix(x_arr, phi_arr, p_vec)
            Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
            M = Q.T @ L_w @ Q
            lam_min_real = float(np.min(np.linalg.eigvalsh(M)))
            realized_lam_M[f"{lam_target:.0e}"] = lam_min_real
        else:
            realized_lam_M[f"{lam_target:.0e}"] = float("nan")

    experiment_results.append({
        "name": name,
        "h_eff": heff,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "kappa_T": kappa_T,
        "mu_by_lam_M": mu_by_lam_M,
        "realized_lam_M": realized_lam_M,
    })

# Save JSON
out_file = Path("results/cheap_experiments_results.json")
out_file.parent.mkdir(parents=True, exist_ok=True)
with out_file.open("w", encoding="utf-8") as f:
    json.dump(experiment_results, f, indent=2)

print("EXPERIMENTOS COMPLETADOS:")
print(f"Resultados guardados en {out_file}\n")

# Print summary tables
print("=== TABLA 1: PROPIEDADES ESPECTRALES DE LA BASE T ===")
print(f"{'Candidato':<26} | {'h_eff':<10} | {'sigma_max':<10} | {'sigma_min':<10} | {'kappa(T)':<10}")
print("-" * 75)
for r in experiment_results:
    print(f"{r['name']:<26} | {r['h_eff']:<10.4e} | {r['sigma_max']:<10.4f} | {r['sigma_min']:<10.4e} | {r['kappa_T']:<10.2f}")

print("\n=== TABLA 2: CURVA mu*(lambda_M) DE FASE II ===")
headers = [f"{lm:.0e}" for lm in lam_M_targets]
header_str = " | ".join(f"{h:<10}" for h in headers)
print(f"{'Candidato':<26} | {'h_eff':<10} | {header_str}")
print("-" * (40 + 13 * len(headers)))
for r in experiment_results:
    vals = [f"{r['mu_by_lam_M'][h]:<10.8f}" for h in headers]
    val_str = " | ".join(vals)
    print(f"{r['name']:<26} | {r['h_eff']:<10.4e} | {val_str}")

print("\n=== TABLA 3: lambda_min(M) REALIZADO ===")
print(f"{'Candidato':<26} | {'lambda_target':<14} | {'lambda_realizado':<18} | {'Ratio Real/Target':<18}")
print("-" * 80)
for r in experiment_results:
    for h in ["1e-11", "1e-10", "1e-09"]:
        real = r['realized_lam_M'][h]
        targ = float(h)
        ratio = real / targ if targ > 0 else float("nan")
        print(f"{r['name']:<26} | {h:<14} | {real:<18.6e} | {ratio:<18.2f}")
