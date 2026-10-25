"""Verify champions with SCS in Phase II and compute single summary table at floor 1e-4."""

from __future__ import annotations

import json
from kuramoto_templates.cases.s7 import get_s7_template, get_s7_appendix_a_template
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import ScaledSDPEvaluator
from kuramoto_templates.search.genome import Genome

with open("results/overnight_12h_verified_report.json", "r", encoding="utf-8") as f:
    report = json.load(f)

print("=" * 80)
print("PARTE 2: VERIFICACION CON SCS DE FASE II EN CAMPEONES DE CORRIDAS 1 Y 3")
print("=" * 80)

scs_runs = [report["runs"][0], report["runs"][2]]  # Corrida 1 y Corrida 3

for run_data in scs_runs:
    rname = run_data["run_name"]
    cid = run_data["id"]
    x = tuple(run_data["phenotype_x"])
    phi = tuple(run_data["phenotype_phi"])
    heff = run_data["h_eff"]
    print(f"\n--- {rname} (ID: {cid}, h_eff: {heff:.6e}) ---")

    # Clarabel vs SCS en varios pisos de tolerancia
    floors = [1e-5, 5e-5, 1e-4, 1e-3]
    print(f"{'Piso (eps)':<12} | {'mu* (Clarabel)':<18} | {'mu* (SCS 1e-6)':<18} | {'mu* (SCS 1e-5)':<18} | {'Diferencia Clar-SCS':<20}")
    print("-" * 92)

    ev_clarabel = ScaledSDPEvaluator(
        EvalConfig(solver_name="Clarabel", solver_tol=1e-9, h_eff_min=heff),
        h_eff_min=heff,
    )
    ev_scs_1e6 = ScaledSDPEvaluator(
        EvalConfig(solver_name="SCS", solver_tol=1e-6, h_eff_min=heff),
        h_eff_min=heff,
    )
    ev_scs_1e5 = ScaledSDPEvaluator(
        EvalConfig(solver_name="SCS", solver_tol=1e-5, h_eff_min=heff),
        h_eff_min=heff,
    )

    for fl in floors:
        res_c = ev_clarabel.evaluate_phase2(x, phi, eps=fl)
        res_s6 = ev_scs_1e6.evaluate_phase2(x, phi, eps=fl)
        res_s5 = ev_scs_1e5.evaluate_phase2(x, phi, eps=fl)

        mu_c = res_c.mu if res_c.mu is not None else float("nan")
        mu_s6 = res_s6.mu if res_s6.mu is not None else float("nan")
        mu_s5 = res_s5.mu if res_s5.mu is not None else float("nan")
        diff = abs(mu_c - mu_s6)

        print(f"{fl:<12.1e} | {mu_c:<18.8f} | {mu_s6:<18.8f} | {mu_s5:<18.8f} | {diff:<20.2e}")


print("\n" + "=" * 80)
print("PARTE 3: TABLA UNICA A PISO 1e-4 CON EL MISMO EVALUADOR (ScaledSDPEvaluator)")
print("=" * 80)

# Configuramos un unico evaluador para todos con eps=1e-4
evaluator = ScaledSDPEvaluator(
    EvalConfig(solver_name="Clarabel", solver_tol=1e-9, lam_floor_hat=1e-4, h_eff_min=1e-4),
    lam_floor_hat=1e-4,
    h_eff_min=1e-4,
)

# Candidatos a evaluar:
# 1. Campeon Corrida 1
# 2. Campeon Corrida 2
# 3. Campeon Corrida 3
# 4. Campeon Corrida 4
# 5. S7 directa (h=0.002, t=0.0) de Tabla II
# 6. S7 Apendice A (h=0.002, t=0.0)

s7_direct = get_s7_template(h=0.002, t=0.0)
s7_app_a = get_s7_appendix_a_template(h=0.002, t=0.0)

items = [
    ("Corrida 1 (ID 125678, seeded)", tuple(report["runs"][0]["phenotype_x"]), tuple(report["runs"][0]["phenotype_phi"]), report["runs"][0]["h_eff"]),
    ("Corrida 2 (ID 2, seeded)", tuple(report["runs"][1]["phenotype_x"]), tuple(report["runs"][1]["phenotype_phi"]), report["runs"][1]["h_eff"]),
    ("Corrida 3 (ID 4, seeded)", tuple(report["runs"][2]["phenotype_x"]), tuple(report["runs"][2]["phenotype_phi"]), report["runs"][2]["h_eff"]),
    ("Corrida 4 (ID 346681, unseeded)", tuple(report["runs"][3]["phenotype_x"]), tuple(report["runs"][3]["phenotype_phi"]), report["runs"][3]["h_eff"]),
    ("S7 directa (h=0.002, t=0)", s7_direct.x, s7_direct.phi, None),
    ("S7 Apendice A (h=0.002, t=0)", s7_app_a.x, s7_app_a.phi, None),
]

table_rows = []
for label, x_cand, phi_cand, heff_hint in items:
    # Phase I at floor 1e-4
    res_p1 = evaluator.evaluate_phase1(x_cand, phi_cand)
    # Phase II at floor 1e-4
    res_p2 = evaluator.evaluate_phase2(x_cand, phi_cand, eps=1e-4)

    h_eff_val = res_p1.solver_info.get("h_eff") if res_p1.solver_info else None
    lam_hat_val = res_p1.lam
    s_star_val = res_p1.s
    status_p1 = res_p1.status
    mu_star_p2 = res_p2.mu

    table_rows.append({
        "label": label,
        "m": len(x_cand),
        "h_eff": h_eff_val,
        "status_p1": status_p1,
        "s_star": s_star_val,
        "lam_hat": lam_hat_val,
        "mu_star_1e4": mu_star_p2,
    })

print(f"{'Candidato / Estado':<35} | {'m':<2} | {'h_eff':<10} | {'Status P1':<10} | {'s* (slack)':<12} | {'lam_hat':<12} | {'mu* (piso 1e-4)':<15}")
print("-" * 115)
for r in table_rows:
    heff_s = f"{r['h_eff']:.4e}" if r['h_eff'] is not None else "N/A"
    s_s = f"{r['s_star']:.2e}" if r['s_star'] is not None else "N/A"
    lam_s = f"{r['lam_hat']:+.4e}" if r['lam_hat'] is not None else "N/A"
    mu_s = f"{r['mu_star_1e4']:.8f}" if r['mu_star_1e4'] is not None else "Infeasible"
    print(f"{r['label']:<35} | {r['m']:<2} | {heff_s:<10} | {r['status_p1']:<10} | {s_s:<12} | {lam_s:<12} | {mu_s:<15}")

# Save results for artifact/output
with open("results/comparison_piso_1e4.json", "w", encoding="utf-8") as f:
    json.dump(table_rows, f, indent=2)
print("\nResultados guardados en results/comparison_piso_1e4.json")
