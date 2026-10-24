"""Run mandatory experiments for mu_t in {0.690, 0.6929077, 0.6935, 0.695} and output results."""

import json
from pathlib import Path
import sys

# Ensure src is on path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.engine import GAEngine

MU_T_VALUES = [0.690, 0.6929077, 0.6935, 0.695]


def main() -> None:
    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)
    runs_dir = Path("runs")
    runs_dir.mkdir(parents=True, exist_ok=True)

    summary_records: list[dict] = []

    for mu_t in MU_T_VALUES:
        print(f"\n=======================================================")
        print(f"Executing mandatory experiment: mu_t = {mu_t}")
        print(f"=======================================================")

        cfg = GAConfig(
            mu_t=mu_t,
            num_islands=2,
            island_pop_size=10,
            num_generations=5,
            use_seeds=True,
            seed=42,
            checkpoint_interval=5,
        )

        out_run_dir = runs_dir / f"exp_mu_{mu_t}"
        engine = GAEngine(config=cfg, output_dir=out_run_dir, allow_dirty=True)
        res = engine.run()

        # Best per m
        m_records = {}
        for m in range(2, 11):
            if m in engine.hall_of_fame_per_m:
                ind = engine.hall_of_fame_per_m[m]
                val_reported = ind.lam_star if ind.s_star <= cfg.tol else -ind.s_star
                metric_label = "lam*" if ind.s_star <= cfg.tol else "-s*"
                m_records[m] = {
                    "reported_value": val_reported,
                    "metric_label": metric_label,
                    "s_star": ind.s_star,
                    "lam_star": ind.lam_star,
                    "lam_normalized": ind.lam_normalized,
                    "f": ind.f,
                    "id": ind.id,
                }
            else:
                m_records[m] = None

        record = {
            "mu_t": mu_t,
            "best_global_f": res["best_f"],
            "best_global_m": res["best_m"],
            "best_global_lam": res["best_lam"],
            "best_global_s": res["best_s"],
            "total_evaluations": res["total_evaluations"],
            "cache_hits": res["cache_hits"],
            "per_m": m_records,
        }
        summary_records.append(record)

    # Save to JSON
    with (results_dir / "mandatory_experiments.json").open("w", encoding="utf-8") as f:
        json.dump(summary_records, f, indent=2)

    # Format Markdown Table labeled "evidencia numérica"
    md_lines = [
        "# Evidencia Numérica: Experimentos Obligatorios por $\\mu_t$ y por $m$",
        "",
        "> [!IMPORTANT]",
        "> **Evidencia numérica**: Ningún resultado aquí tabulado se declara cota inferior.",
        "> Todas las métricas representan valores observados por optimización numérica con tolerancia $10^{-8}$.",
        "",
        "## Resumen Global por $\\mu_t$",
        "",
        "| $\\mu_t$ | Mejor $f$ | Mejor $m$ | $\\lambda^*$ crudo | $s^*$ | Evals | Cache Hits |",
        "|---|---|---|---|---|---|---|",
    ]

    for rec in summary_records:
        md_lines.append(
            f"| {rec['mu_t']} | {rec['best_global_f']:.6e} | {rec['best_global_m']} | {rec['best_global_lam']:.6e} | {rec['best_global_s']:.6e} | {rec['total_evaluations']} | {rec['cache_hits']} |"
        )

    md_lines.append("")
    md_lines.append("## Mejor Métrica por cada $m \\in [2, 10]$ (rotulado evidencia numérica)")
    md_lines.append("")
    md_lines.append("| $\\mu_t$ | m=2 | m=3 | m=4 | m=5 | m=6 | m=7 | m=8 | m=9 | m=10 |")
    md_lines.append("|---|---|---|---|---|---|---|---|---|---|---|")

    for rec in summary_records:
        row = [str(rec["mu_t"])]
        for m in range(2, 11):
            m_info = rec["per_m"].get(m)
            if m_info is not None:
                val = m_info["reported_value"]
                lbl = m_info["metric_label"]
                row.append(f"{lbl}={val:.2e}")
            else:
                row.append("-")
        md_lines.append("| " + " | ".join(row) + " |")

    with (results_dir / "mandatory_experiments.md").open("w", encoding="utf-8") as f:
        f.write("\n".join(md_lines) + "\n")

    print("\nAll mandatory experiments completed successfully!")
    print(f"Results saved to {results_dir / 'mandatory_experiments.json'} and .md")


if __name__ == "__main__":
    main()
