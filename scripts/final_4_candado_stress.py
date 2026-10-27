"""Script 4: Estres Empirico del Candado de Escala.

Auditoria probabilistica de mutaciones phase_dev:
- Carga S7 con h = 1e-3.
- Aplica 1000 mutaciones independientes con phase_dev.
- Evalua h_eff resultante para cada mutacion.
- Reporta porcentaje de reduccion de h_eff, tasa de NOT_SCALABLE (< 1e-3) y media de h_eff.
"""

from __future__ import annotations

import numpy as np

from kuramoto_templates.evaluators.scaled import compute_h_eff, compute_sheets
from kuramoto_templates.search.genome import to_phenotype
from kuramoto_templates.search.operators.phase_dev import phase_dev
from kuramoto_templates.search.seeds import create_s7_genome_seed


def main() -> None:
    print("=" * 80)
    print("SCRIPT 4: ESTRES EMPIRICO DEL CANDADO DE ESCALA (1000 MUTACIONES)")
    print("Operador: phase_dev (sigma_ell=0.5, flip_sgn_prob=0.2) sobre S7(h=1e-3)")
    print("=" * 80)

    # 1. Cargar S7 nominal
    seed_ind = create_s7_genome_seed(h=0.001)
    base_genome = seed_ind.genome
    x_base, phi_base = to_phenotype(base_genome)
    sheets_base, _ = compute_sheets(x_base, phi_base, genome=base_genome)
    h_eff_base = compute_h_eff(x_base, phi_base, sheets_base)

    print(f"\n1. ESTADO BASE DE S7(h=1e-3):")
    print(f"   h_eff Nominal : {h_eff_base:.8f}")
    print(f"   Umbral Candado: h_eff_min = 0.001 (1e-3)")

    # 2. 1000 Mutaciones independientes
    rng = np.random.default_rng(42)
    n_trials = 1000
    below_1e3 = 0
    reduced = 0
    h_eff_list = []

    for _ in range(n_trials):
        mut_genome = phase_dev(base_genome, rng)
        x_m, phi_m = to_phenotype(mut_genome)
        sheets_m, _ = compute_sheets(x_m, phi_m, genome=mut_genome)
        h_m = compute_h_eff(x_m, phi_m, sheets_m)
        h_eff_list.append(h_m)
        if h_m < h_eff_base:
            reduced += 1
        if h_m < 0.001:
            below_1e3 += 1

    pct_reduced = (reduced / n_trials) * 100
    pct_locked = (below_1e3 / n_trials) * 100
    mean_h_eff = float(np.mean(h_eff_list))
    std_h_eff = float(np.std(h_eff_list))

    print(f"\n2. RESULTADOS DEL EXPERIMENTO (N = {n_trials:,} mutaciones):")
    print(f"   Mutaciones que redujeron h_eff : {reduced:,} / {n_trials:,} ({pct_reduced:.1f}%)")
    print(f"   Mutaciones que activan BLOQUEO : {below_1e3:,} / {n_trials:,} ({pct_locked:.1f}%) [h_eff < 1e-3]")
    print(f"   h_eff promedio resultante     : {mean_h_eff:.6f} +/- {std_h_eff:.6f}")
    print(f"   h_eff minimo alcanzado        : {min(h_eff_list):.6f}")
    print(f"   h_eff maximo alcanzado        : {max(h_eff_list):.6f}")

    print("\n3. COMPARACION CON EL REGISTRO DE 12 HORAS (run.jsonl):")
    print(f"   Tasa NOT_SCALABLE observada en run.jsonl (phase_dev) : 35.16% (41,594 / 118,306)")
    print(f"   Tasa NOT_SCALABLE en experimento controlado (S7)     : {pct_locked:.2f}% ({below_1e3} / {n_trials})")
    print(f"   Discrepancia entre prediccion y medicion            : {abs(pct_locked - 35.16):.2f}%")
    print("=" * 80)


if __name__ == "__main__":
    main()
