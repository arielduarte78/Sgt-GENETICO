"""Orchestrator for the 12-hour overnight parallel GA run.

Runs 4 independent GA instances in parallel:
- Run 1 (runs/night12h_seeded_mu06923_s42): use_seeds: true, mu_t: 0.6923, lam_floor_hat: 5e-5, h_eff_min: 1e-3, seed: 42
- Run 2 (runs/night12h_seeded_mu06925_s123): use_seeds: true, mu_t: 0.6925, lam_floor_hat: 5e-5, h_eff_min: 1e-3, seed: 123
- Run 3 (runs/night12h_seeded_mu06926_s2026): use_seeds: true, mu_t: 0.6926, lam_floor_hat: 1e-5, h_eff_min: 8e-4, seed: 2026
- Run 4 (runs/night12h_unseeded_m5_8_s42): use_seeds: false, mu_t: 0.6900, lam_floor_hat: 5e-5, h_eff_min: 1e-3, m in [5, 8], seed: 42

Monitors execution every 30s, maintains runs/night12h_live_dashboard.txt,
terminates cleanly at 12 hours, performs multi-solver independent verification,
and saves the final report in results/overnight_12h_verified_report.json.
"""

from __future__ import annotations

import argparse
import ctypes
from datetime import datetime
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.domain.conditions import (
    build_lw_matrix,
    p_matrix_to_vector,
    verify_solution,
)
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)
from kuramoto_templates.search.config import GAConfig
from kuramoto_templates.search.genome import Genome, wrap_minus_pi_to_pi


RUN_SPECS: list[dict[str, Any]] = [
    {
        "name": "night12h_seeded_mu06923_s42",
        "mu_t": 0.6923,
        "use_seeds": True,
        "lam_floor_hat": 5.0e-5,
        "h_eff_min": 1.0e-3,
        "seed": 42,
        "min_m": 2,
        "max_m": 10,
    },
    {
        "name": "night12h_seeded_mu06925_s123",
        "mu_t": 0.6925,
        "use_seeds": True,
        "lam_floor_hat": 5.0e-5,
        "h_eff_min": 1.0e-3,
        "seed": 123,
        "min_m": 2,
        "max_m": 10,
    },
    {
        "name": "night12h_seeded_mu06926_s2026",
        "mu_t": 0.6926,
        "use_seeds": True,
        "lam_floor_hat": 1.0e-5,
        "h_eff_min": 8.0e-4,
        "seed": 2026,
        "min_m": 2,
        "max_m": 10,
    },
    {
        "name": "night12h_unseeded_m5_8_s42",
        "mu_t": 0.6900,
        "use_seeds": False,
        "lam_floor_hat": 5.0e-5,
        "h_eff_min": 1.0e-3,
        "seed": 42,
        "min_m": 5,
        "max_m": 8,
    },
]

COMMON_CONFIG = {
    "num_generations": 3500,
    "island_pop_size": 48,
    "num_islands": 6,
    "checkpoint_interval": 5,
    "checkpoint_keep_count": 3,
    "evaluator": "scaled",
}


def prevent_windows_sleep() -> None:
    """Keep system awake during the overnight run."""
    if sys.platform == "win32":
        try:
            # ES_CONTINUOUS (0x80000000) | ES_SYSTEM_REQUIRED (0x00000001)
            ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
            print("[INFO] Windows sleep prevention activated (ES_CONTINUOUS | ES_SYSTEM_REQUIRED).")
        except Exception as e:
            print(f"[WARN] Failed to set execution state: {e}")


def restore_windows_sleep() -> None:
    """Restore normal Windows power management."""
    if sys.platform == "win32":
        try:
            ctypes.windll.kernel32.SetThreadExecutionState(0x80000000)
            print("[INFO] Windows sleep prevention restored to default.")
        except Exception as e:
            print(f"[WARN] Failed to restore execution state: {e}")


def compute_gauge_aligned_s7_dist(
    x: tuple[float, ...],
    phi: tuple[float, ...],
    x_s7: tuple[float, ...],
    phi_s7: tuple[float, ...],
) -> float:
    """Compute gauge-aligned min euclidean distance in (x, phi) relative to S7."""
    from scipy.optimize import linear_sum_assignment  # type: ignore[import-untyped]

    if len(x) != len(x_s7):
        return float("inf")

    m = len(x)
    x_arr = np.asarray(x, dtype=np.float64)
    phi_arr = np.asarray(phi, dtype=np.float64)
    best_dist = float("inf")

    for i in range(m):
        alpha = wrap_minus_pi_to_pi(phi_s7[-1] - phi_arr[i])
        phi_rot = np.array([wrap_minus_pi_to_pi(p + alpha) for p in phi_arr])
        cost = np.zeros((m, m))
        for r in range(m):
            for c in range(m):
                cost[r, c] = (x_arr[r] - x_s7[c]) ** 2 + (wrap_minus_pi_to_pi(phi_rot[r] - phi_s7[c])) ** 2
        row_ind, col_ind = linear_sum_assignment(cost)
        total_sq = sum(cost[r, c] for r, c in zip(row_ind, col_ind))
        dist = float(np.sqrt(total_sq))
        if dist < best_dist:
            best_dist = dist
    return best_dist


def find_latest_checkpoint(out_dir: Path) -> Path | None:
    """Find the latest checkpoint file with highest generation number."""
    ckpt_dir = out_dir / "checkpoints"
    if not ckpt_dir.exists():
        return None
    ckpts = list(ckpt_dir.glob("checkpoint_gen_*.json"))
    if not ckpts:
        return None

    def _extract_gen(p: Path) -> int:
        try:
            return int(p.stem.split("_")[-1])
        except ValueError:
            return -1

    return max(ckpts, key=_extract_gen)


def count_false_positives(log_file: Path) -> int:
    """Count false positive rejections recorded in evaluations.jsonl."""
    if not log_file.exists():
        return 0
    count = 0
    try:
        with log_file.open("r", encoding="utf-8") as f:
            for line in f:
                if "falso_positivo" in line or "false_positive" in line:
                    count += 1
    except Exception:
        pass
    return count


def read_run_status(spec: dict[str, Any], runs_dir: Path) -> dict[str, Any]:
    """Read current generation, metrics and HoF best from status and checkpoints."""
    out_dir = runs_dir / spec["name"]
    status_file = out_dir / "status.json"
    evals_file = out_dir / "evaluations.jsonl"

    gen = 0
    evals = 0
    best_f = -1.0
    best_mu: float | None = None
    h_eff_val: float | None = None
    status_str = "INITIALIZING"

    if status_file.exists():
        try:
            with status_file.open("r", encoding="utf-8") as f:
                data = json.load(f)
            gen = int(data.get("generation", 0))
            evals = int(data.get("total_evaluations", 0))
            best_f = float(data.get("best_f", -1.0))
            status_str = str(data.get("status", "running")).upper()
            if data.get("best_mu_star") is not None:
                best_mu = float(data["best_mu_star"])
        except Exception:
            pass

    # Read latest checkpoint for best individual details
    latest_ckpt = find_latest_checkpoint(out_dir)
    if latest_ckpt:
        try:
            with latest_ckpt.open("r", encoding="utf-8") as f:
                ckpt_data = json.load(f)
            hof = ckpt_data.get("hall_of_fame_global", [])
            if hof:
                best_ind = hof[0]
                if best_ind.get("mu_star") is not None:
                    best_mu = float(best_ind["mu_star"])
                px = best_ind.get("phenotype_x")
                pphi = best_ind.get("phenotype_phi")
                g_dict = best_ind.get("genome", {})
                if px and pphi:
                    g_obj = Genome(
                        m=g_dict["m"],
                        K_ref=g_dict["K_ref"],
                        k=tuple(g_dict["k"]),
                        sgn=tuple(g_dict["sgn"]),
                        ell=tuple(g_dict["ell"]),
                        u=tuple(g_dict["u"]),
                    )
                    sheets, _ = compute_sheets(px, pphi, genome=g_obj)
                    h_eff_val = compute_h_eff(px, pphi, sheets)
        except Exception:
            pass

    fp_count = count_false_positives(evals_file)

    return {
        "name": spec["name"],
        "status": status_str,
        "generation": gen,
        "total_evaluations": evals,
        "best_f": best_f,
        "best_mu_star": best_mu,
        "h_eff": h_eff_val,
        "false_positives": fp_count,
    }


def write_live_dashboard(
    dashboard_file: Path,
    statuses: list[dict[str, Any]],
    start_time: float,
    duration_hours: float,
) -> None:
    """Format and write live dashboard to runs/night12h_live_dashboard.txt."""
    elapsed_sec = time.perf_counter() - start_time
    total_target_sec = duration_hours * 3600.0
    pct = min(100.0, (elapsed_sec / total_target_sec) * 100.0)

    elapsed_h = int(elapsed_sec // 3600)
    elapsed_m = int((elapsed_sec % 3600) // 60)
    elapsed_s = int(elapsed_sec % 60)

    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "=" * 115,
        "                     LOTE NOCTURNO 12 HORAS - KURAMOTO TEMPLATES GENETIC ALGORITHM",
        "=" * 115,
        f"Inicio: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time() - elapsed_sec))} | "
        f"Transcurrido: {elapsed_h:02d}h {elapsed_m:02d}m {elapsed_s:02d}s / {int(duration_hours):02d}h 00m 00s ({pct:5.1f}%)",
        f"Estatus Global: RUNNING (PID: {ctypes.windll.kernel32.GetCurrentProcessId() if sys.platform == 'win32' else ''}) | "
        "Anti-Suspension: ACTIVATED (ES_CONTINUOUS | ES_SYSTEM_REQUIRED)",
        "-" * 115,
        f"{'Corrida':<36} | {'Gen/Max':<10} | {'Evals/s':<8} | {'Mejor f':<13} | {'Mejor mu*':<12} | {'h_eff':<10} | {'FP'}",
        "-" * 115,
    ]

    for st in statuses:
        c_name = st["name"]
        gen_str = f"{st['generation']}/3500"
        evals_sec = st["total_evaluations"] / max(1.0, elapsed_sec)
        evals_str = f"{evals_sec:.1f}"

        bf = st["best_f"]
        bf_str = f"{bf:+.4e}" if bf > -0.99 else f"{bf:.4f}"

        bmu = st["best_mu_star"]
        bmu_str = f"{bmu:.8f}" if bmu is not None else "None"

        heff = st["h_eff"]
        heff_str = f"{heff:.3e}" if heff is not None else "N/A"

        fp = st["false_positives"]

        lines.append(
            f"{c_name:<36} | {gen_str:<10} | {evals_str:<8} | {bf_str:<13} | {bmu_str:<12} | {heff_str:<10} | {fp}"
        )

    lines.append("=" * 115)
    lines.append(f"Ultima actualizacion: {now_str} | Monitoreo en vivo: python scripts/watch_night.py")
    lines.append("")

    dashboard_file.parent.mkdir(parents=True, exist_ok=True)
    temp_file = dashboard_file.with_suffix(".tmp")
    with temp_file.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    temp_file.replace(dashboard_file)


def verify_run_champion(spec: dict[str, Any], runs_dir: Path) -> dict[str, Any]:
    """Execute complete post-run verification on the champion of a run."""
    import mpmath  # type: ignore[import-untyped]
    import scipy.linalg  # type: ignore[import-untyped]

    out_dir = runs_dir / spec["name"]
    latest_ckpt = find_latest_checkpoint(out_dir)
    if not latest_ckpt:
        return {"error": f"No checkpoint found in {out_dir}"}

    with latest_ckpt.open("r", encoding="utf-8") as f:
        ckpt_data = json.load(f)

    hof = ckpt_data.get("hall_of_fame_global", [])
    if not hof:
        return {"error": "Hall of Fame is empty"}

    champ = hof[0]
    x_tup = tuple(champ["phenotype_x"])
    phi_tup = tuple(champ["phenotype_phi"])
    g_dict = champ["genome"]
    g_obj = Genome(
        m=g_dict["m"],
        K_ref=g_dict["K_ref"],
        k=tuple(g_dict["k"]),
        sgn=tuple(g_dict["sgn"]),
        ell=tuple(g_dict["ell"]),
        u=tuple(g_dict["u"]),
    )

    x_arr = np.asarray(x_tup, dtype=float)
    phi_arr = np.asarray(phi_tup, dtype=float)
    m = len(x_arr)

    # 1. Scale metrics
    sheets, ref_sheet = compute_sheets(x_arr, phi_arr, genome=g_obj)
    h_eff = compute_h_eff(x_arr, phi_arr, sheets)
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, h_eff)

    # 2. Phase I evaluation with Clarabel 1e-9
    ev_scaled = ScaledSDPEvaluator(
        EvalConfig(
            mu_t=spec["mu_t"],
            solver_tol=1e-9,
            lam_floor_hat=spec["lam_floor_hat"],
            h_eff_min=spec["h_eff_min"],
        ),
        lam_floor_hat=spec["lam_floor_hat"],
        h_eff_min=spec["h_eff_min"],
    )
    res_p1 = ev_scaled.evaluate_phase1(x_tup, phi_tup, genome=g_obj)

    # 3. Pure NumPy & mpmath 50 digits on G_hat and M
    lam_numpy_M = None
    lam_mpmath_M = None
    lam_numpy_G = None
    lam_mpmath_G = None
    verif_res = None

    if res_p1.P is not None:
        p_vec = p_matrix_to_vector(np.asarray(res_p1.P, dtype=float))
        L_w = build_lw_matrix(x_arr, phi_arr, p_vec)
        Q = scipy.linalg.null_space(np.ones((1, m), dtype=float))
        M = Q.T @ L_w @ Q
        lam_numpy_M = float(np.min(np.linalg.eigvalsh(M)))

        mpmath.mp.dps = 50
        dim_M = M.shape[0]
        M_mp = mpmath.mp.matrix(dim_M, dim_M)
        for r in range(dim_M):
            for c in range(dim_M):
                M_mp[r, c] = mpmath.mp.mpf(str(M[r, c]))
        eigs_M_mp = mpmath.mp.eigsy(M_mp, eigvals_only=True)
        lam_mpmath_M = float(min(eigs_M_mp))

        if T is not None and h_eff >= spec["h_eff_min"]:
            G_hat = (T.T @ L_w @ T) / (h_eff**2)
            lam_numpy_G = float(np.min(np.linalg.eigvalsh(G_hat)))
            dim_G = G_hat.shape[0]
            G_mp = mpmath.mp.matrix(dim_G, dim_G)
            for r in range(dim_G):
                for c in range(dim_G):
                    G_mp[r, c] = mpmath.mp.mpf(str(G_hat[r, c]))
            eigs_G_mp = mpmath.mp.eigsy(G_mp, eigvals_only=True)
            lam_mpmath_G = float(min(eigs_G_mp))

        verif_res = verify_solution(x_tup, phi_tup, res_p1.P, eps=spec["lam_floor_hat"], eps_c=1e-4)

    # 4. Phase II curve across floors: [1e-5, 5e-5, 1e-4, 1e-3]
    floors = [1e-5, 5e-5, 1e-4, 1e-3]
    mu_vs_floor = {}
    for fl in floors:
        ev_fl = ScaledSDPEvaluator(
            EvalConfig(
                mu_t=spec["mu_t"],
                solver_tol=1e-8,
                lam_floor_hat=fl,
                h_eff_min=spec["h_eff_min"],
            ),
            lam_floor_hat=fl,
            h_eff_min=spec["h_eff_min"],
        )
        r_p2 = ev_fl.evaluate_phase2(x_tup, phi_tup, eps=fl, genome=g_obj)
        mu_vs_floor[str(fl)] = r_p2.mu

    # 5. Distance to S7
    s7_ref = get_s7_appendix_a_template(h=0.002, t=0.0)
    dist_s7 = compute_gauge_aligned_s7_dist(x_tup, phi_tup, s7_ref.x, s7_ref.phi)

    return {
        "run_name": spec["name"],
        "checkpoint": str(latest_ckpt),
        "id": champ.get("id"),
        "m": m,
        "h_eff": h_eff,
        "f": champ.get("f"),
        "s_star": champ.get("s_star"),
        "lam_hat": champ.get("lam_normalized"),
        "clarabel_1e9_lam": res_p1.lam_normalized,
        "clarabel_1e9_s": res_p1.s,
        "clarabel_1e9_status": res_p1.status,
        "lam_numpy_M": lam_numpy_M,
        "lam_mpmath_M": lam_mpmath_M,
        "lam_numpy_G_hat": lam_numpy_G,
        "lam_mpmath_G_hat": lam_mpmath_G,
        "verify_solution": {
            "is_valid": verif_res["is_valid"] if verif_res else False,
            "max_F": float(verif_res["max_F"]) if verif_res else None,
            "min_c": float(verif_res["min_c"]) if verif_res else None,
            "mu": float(verif_res["mu"]) if verif_res else None,
        } if verif_res else None,
        "mu_vs_floor": mu_vs_floor,
        "dist_s7": dist_s7,
        "phenotype_x": list(x_tup),
        "phenotype_phi": list(phi_tup),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Overnight 12-hour GA Orchestrator")
    parser.add_argument("--duration-hours", type=float, default=12.0, help="Duration in hours (default 12.0)")
    parser.add_argument("--poll-interval", type=int, default=30, help="Polling interval in seconds (default 30)")
    args = parser.parse_args()

    duration_sec = args.duration_hours * 3600.0
    runs_dir = Path("runs")
    runs_dir.mkdir(parents=True, exist_ok=True)
    dashboard_file = runs_dir / "night12h_live_dashboard.txt"
    report_file = Path("results/overnight_12h_verified_report.json")
    report_file.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print(f"INICIANDO ORQUESTADOR DE CORRIDA NOCTURNA ({args.duration_hours:.1f} HORAS)")
    print("=" * 80)

    # 1. Prevent Windows sleep
    prevent_windows_sleep()

    # 2. Launch 4 parallel processes
    processes: list[subprocess.Popen[Any]] = []
    configs: list[GAConfig] = []

    for spec in RUN_SPECS:
        out_dir = runs_dir / spec["name"]
        out_dir.mkdir(parents=True, exist_ok=True)

        cfg = GAConfig(
            mu_t=spec["mu_t"],
            lam_floor_hat=spec["lam_floor_hat"],
            h_eff_min=spec["h_eff_min"],
            seed=spec["seed"],
            use_seeds=spec["use_seeds"],
            min_m=spec["min_m"],
            max_m=spec["max_m"],
            num_generations=COMMON_CONFIG["num_generations"],
            island_pop_size=COMMON_CONFIG["island_pop_size"],
            num_islands=COMMON_CONFIG["num_islands"],
            checkpoint_interval=COMMON_CONFIG["checkpoint_interval"],
            checkpoint_keep_count=COMMON_CONFIG["checkpoint_keep_count"],
            evaluator=COMMON_CONFIG["evaluator"],
        )
        configs.append(cfg)
        cfg_path = out_dir / "config.yaml"
        cfg.to_yaml(cfg_path)

        cmd = [
            sys.executable,
            "-m",
            "kuramoto_templates.search.cli",
            "run",
            str(cfg_path),
            "--output-dir",
            str(out_dir),
            "--allow-dirty",
        ]
        if not spec["use_seeds"]:
            cmd.append("--no-seeds")

        log_out = open(out_dir / "stdout.log", "w", encoding="utf-8")
        p = subprocess.Popen(cmd, stdout=log_out, stderr=subprocess.STDOUT)
        processes.append(p)
        print(f"Lanzada {spec['name']} (PID {p.pid}) -> {out_dir}")

    start_time = time.perf_counter()

    try:
        while True:
            elapsed = time.perf_counter() - start_time
            if elapsed >= duration_sec:
                print(f"\n[TIEMPO CUMPLIDO] Alcanzadas {args.duration_hours:.2f} horas. Iniciando frenado limpio...")
                break

            # Poll statuses
            statuses = [read_run_status(spec, runs_dir) for spec in RUN_SPECS]
            write_live_dashboard(dashboard_file, statuses, start_time, args.duration_hours)

            # Check if all processes have completed on their own
            if all(p.poll() is not None for p in processes):
                print("\n[TODAS COMPLETADAS] Todas las corridas alcanzaron 3500 generaciones.")
                break

            time.sleep(args.poll_interval)

    except KeyboardInterrupt:
        print("\n[INTERRUPCION] Deteniendo orquestador y subprocesos...")

    finally:
        # Cleanly stop processes
        for i, p in enumerate(processes):
            if p.poll() is None:
                print(f"Terminando subproceso {RUN_SPECS[i]['name']} (PID {p.pid})...")
                p.terminate()
                try:
                    p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    p.kill()

        # Update dashboard final
        statuses = [read_run_status(spec, runs_dir) for spec in RUN_SPECS]
        write_live_dashboard(dashboard_file, statuses, start_time, args.duration_hours)

        # Execute full post-run verification on all champions
        print("\n" + "=" * 80)
        print("EJECUTANDO AUTO-VERIFICACION COMPLETA DE GANADORES NOCTURNOS...")
        print("=" * 80)

        verified_results = []
        for spec in RUN_SPECS:
            print(f"Verificando campeon de {spec['name']}...")
            verif = verify_run_champion(spec, runs_dir)
            verified_results.append(verif)

        final_report = {
            "timestamp": datetime.now().isoformat(),
            "duration_hours_target": args.duration_hours,
            "duration_hours_elapsed": (time.perf_counter() - start_time) / 3600.0,
            "runs": verified_results,
        }

        with report_file.open("w", encoding="utf-8") as f:
            json.dump(final_report, f, indent=2)

        print(f"[REPORTE GUARDADO] Reporte verificado guardado en: {report_file}")

        # Restore Windows sleep
        restore_windows_sleep()
        print("Orquestador nocturno finalizado.")


if __name__ == "__main__":
    main()
