"""Script 1: Autopsia del GA y el Candado de Escala.

Auditoria del run nocturno de 12 horas:
- Conteo exhaustivo de eventos totales vs NOT_SCALABLE para cada operador.
- Snippet exacto del evaluador donde se activa el candado.
"""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    print("=" * 80)
    print("SCRIPT 1: AUTOPSIA DEL GA Y EL CANDADO DE ESCALA (run.jsonl)")
    print("=" * 80)

    run_path = Path("runs/night12h_seeded_mu06923_s42/run.jsonl")
    if not run_path.exists():
        print(f"ERROR: No se encontro el archivo {run_path}")
        return

    operators: dict[str, dict[str, int]] = {}
    total_evals = 0
    total_not_scalable = 0

    with open(run_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            ev = json.loads(line)
            op = ev.get("operator")
            if not op:
                continue
            total_evals += 1
            if op not in operators:
                operators[op] = {"total": 0, "not_scalable": 0}
            operators[op]["total"] += 1
            if ev.get("solver_status") == "NOT_SCALABLE":
                operators[op]["not_scalable"] += 1
                total_not_scalable += 1

    print(f"\n1. RESUMEN GLOBAL DE EVALUACIONES:")
    print(f"   Evaluaciones Totales : {total_evals:,}")
    print(f"   Total NOT_SCALABLE   : {total_not_scalable:,} ({total_not_scalable / total_evals * 100:.2f}%)")

    print(f"\n2. DISTRIBUCION POR OPERADOR GENETICO:")
    header = f"{'Operador':<30} | {'Intentos':<10} | {'NOT_SCALABLE':<12} | {'Tasa Rechazo':<12}"
    print(header)
    print("-" * len(header))

    for op, counts in sorted(operators.items(), key=lambda x: -x[1]["total"]):
        tot = counts["total"]
        ns = counts["not_scalable"]
        pct = (ns / tot * 100) if tot > 0 else 0.0
        print(f"{op:<30} | {tot:<10,d} | {ns:<12,d} | {pct:10.2f}%")

    print("-" * len(header))

    print(f"\n3. LOCALIZACION DEL CANDADO EN EL CODIGO FUENTE (scaled.py: L274-320):")
    snippet = """
    # src/kuramoto_templates/evaluators/scaled.py:
    sheets, ref_sheet = compute_sheets(
        x_arr, phi_arr, genome=genome, K_ref=4, sheet_tol=0.2
    )
    h_eff = compute_h_eff(x_arr, phi_arr, sheets)

    # Validacion de rango admisible de escala
    if h_eff < self.h_eff_min or h_eff > 0.1:
        res_plain = self.plain_evaluator.evaluate_phase1(x, phi)
        ...
        return EvalResult(..., status="NOT_SCALABLE")

    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, h_eff)
    if T is None:
        res_plain = self.plain_evaluator.evaluate_phase1(x, phi)
        ...
        return EvalResult(..., status="NOT_SCALABLE")
    """
    print(snippet.strip())
    print("\n   [EVALUACION]: Se evalua ANTES de armar la LMI escalada, pero")
    print("   al fallar ejecuta self.plain_evaluator.evaluate_phase1(x, phi),")
    print("   consumiendo un solve completo de Clarabel para asignar penalizacion.")
    print("=" * 80)


if __name__ == "__main__":
    main()
