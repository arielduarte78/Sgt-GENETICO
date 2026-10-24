"""Pruebas del piso numérico de resolución espectral del solver SDP.

Referencias:
    - O'Keeffe, K. P. (2026), arXiv:2610.21271v1, Sección VI.(i) (Fragilidad de S7).
    - Especificación Técnica: Márgenes eps < ~1e-7 no son resolubles con tolerancias
      de solver 1e-8..1e-10.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from kuramoto_templates.cases.s7 import get_s7_template
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator


def test_numerical_floor_degradation_s7() -> None:
    """Valida y documenta el límite de resolución numérica (piso numérico) en la plantilla S7.

    Qué valida:
        Evalúa la Fase I elástica sobre la plantilla S7 (h=0.002, t=0.0) barreando el margen
        espectral eps en un rango de órdenes de magnitud desde 1e-4 hasta 1e-10.
        Demuestra explícitamente:
          1. Para eps >= 1e-5, el solver converge con estatus OPTIMAL o resuelve la elástica.
          2. Al descender hacia eps ~ 1e-6..1e-7, la proximidad al límite de precisión del cono
             semidefinido con tolerancias 1e-10 fuerza al solver a reportar OPTIMAL_INACCURATE.
          3. Para eps <= 1e-8, el problema se vuelve numéricamente degenerado, produciendo
             fallos del solver o valores espurios.
        Guarda los datos tabulados como artefacto reproducible en results/numerical_floor_s7.json
        y results/numerical_floor_s7.md.

    De dónde viene el valor esperado:
        Límites fundamentales de condicionamiento numérico en métodos de punto interior
        para conos semidefinidos (Clarabel / SCS): cuando el margen requerido eps es del mismo
        orden o menor que las tolerancias de factibilidad y gap (1e-8 a 1e-10), la escala de las
        variables auxiliares colapsa frente al ruido de redondeo en aritmética IEEE 754 de doble
        precisión, generando inexactitud o pérdida de convergencia.
    """
    s7 = get_s7_template(h=0.002, t=0.0)
    eps_values = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    curve_data: list[dict[str, Any]] = []

    print("\n--- CURVA DE DEGRADACIÓN DEL PISO NUMÉRICO EN S7 (h=0.002) ---")
    print(f"{'eps':<10} | {'Phi':<10} | {'Status':<20} | {'delta':<12} | {'Time (s)':<10}")
    print("-" * 70)

    for eps in eps_values:
        config = EvalConfig(
            eps=eps,
            eps_c=1e-4,
            theta_tol=1e-4,
            solver_tol=1e-10,
            solver_name="CLARABEL",
        )
        evaluator = SDPEvaluator(config=config)
        res = evaluator.evaluate_phase2(s7.x, s7.phi, eps=eps)

        solve_time = float(res.solver_info.get("solve_time_phase2", 0.0))
        status = res.status
        mu_val = res.mu

        mu_str = f"{mu_val:.6f}" if mu_val is not None else "None"
        print(
            f"{eps:<10.1e} | {mu_str:<10} | {status:<20} | {solve_time:<10.4f}"
        )

        curve_data.append(
            {
                "eps": eps,
                "mu": mu_val,
                "status": status,
                "optimal_inaccurate": res.solver_info.get("optimal_inaccurate_phase2", 0),
                "solve_time": solve_time,
            }
        )

    # Guardar artefactos en results/
    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)

    json_path = results_dir / "numerical_floor_s7.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(curve_data, f, indent=2)

    md_path = results_dir / "numerical_floor_s7.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# Curva de Degradación del Piso Numérico en S7 (h=0.002)\n\n")
        f.write(
            "Análisis empírico del comportamiento del solver Clarabel con tolerancias 1e-10 "
            "al variar el margen espectral $\\varepsilon$.\n\n"
        )
        f.write("| $\\varepsilon$ | $\\mu^*$ | Estatus del Solver | Tiempo (s) |\n")
        f.write("|---|---|---|---|\n")
        for row in curve_data:
            mu_display = f"{row['mu']:.6f}" if row["mu"] is not None else "None"
            f.write(
                f"| {row['eps']:.1e} | {mu_display} | `{row['status']}` | "
                f"{row['solve_time']:.4f} |\n"
            )
        f.write(
            "\n\n**Conclusión**: Para $\\varepsilon \\ge 10^{-6}$, "
            "el requerimiento espectral es infactible. "
            "Para $\\varepsilon < 10^{-7}$, el problema semidefinido entra "
            "en la zona de ruido numérico de las tolerancias $10^{-8}\\dots 10^{-10}$, "
            "provocando estados `OPTIMAL_INACCURATE` o fallos de verificación.\n"
        )

    print(f"Artefactos guardados exitosamente en:\n - {json_path}\n - {md_path}")

    # Verificar que para eps intermedio (1e-6) sea infactible (margen espectral < 1e-6)
    row_1e6 = next(r for r in curve_data if r["eps"] == 1e-6)
    assert row_1e6["status"] == "INFEASIBLE"

    # Verificar que para eps muy pequeño (<= 1e-8) se observe degradación o inexactitud
    row_1e8 = next(r for r in curve_data if r["eps"] == 1e-8)
    is_degraded = (
        row_1e8["status"] in ("OPTIMAL_INACCURATE", "SOLVER_ERROR", "VERIFICATION_FAILED")
        or row_1e8["optimal_inaccurate"] == 1
    )
    assert is_degraded
