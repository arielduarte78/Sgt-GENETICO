"""Generate exact numerical data for final closure report:
- T4 & T6 raw metrics for S7 (h=0.002, t=0, mu_t=0.690) across Clarabel(1e-8), Clarabel(1e-9), SCS(1e-8), NumPy float64 and mpmath 50d.
- T9 benchmark table for m=2..10 without cache.
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path
import numpy as np
import scipy.linalg
import mpmath

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.domain.conditions import (
    build_af_matrix,
    build_ac_matrix,
    build_ad_matrix,
    build_lw_matrix,
    p_matrix_to_vector,
)
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    compute_sheets,
    compute_h_eff,
    build_scaled_basis_t,
)
from kuramoto_templates.search.genome import Genome, canonicalize, decode_x, decode_phi


def get_s7_metrics():
    tmpl = get_s7_appendix_a_template(h=0.002, t=0.0)
    mu_t = 0.690
    x = tmpl.x
    phi = tmpl.phi

    # 1. Scaled Clarabel(1e-8)
    ev_cl8 = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, solver_name="CLARABEL"))
    res_cl8 = ev_cl8.evaluate_phase1(x, phi)

    # 2. Scaled Clarabel(1e-9)
    ev_cl9 = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-9, solver_name="CLARABEL"))
    res_cl9 = ev_cl9.evaluate_phase1(x, phi)

    # 3. Scaled SCS(1e-8)
    ev_scs = ScaledSDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, solver_name="SCS"))
    res_scs = ev_scs.evaluate_phase1(x, phi)

    # 4. Plain Clarabel(1e-8) and SCS(1e-8) for T4 comparison
    ev_p_cl8 = SDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, solver_name="CLARABEL"))
    res_p_cl8 = ev_p_cl8.evaluate_phase1(x, phi)

    ev_p_cl9 = SDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-9, solver_name="CLARABEL"))
    res_p_cl9 = ev_p_cl9.evaluate_phase1(x, phi)

    ev_p_scs = SDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, solver_name="SCS"))
    res_p_scs = ev_p_scs.evaluate_phase1(x, phi)

    # 5. NumPy & mpmath verification on P from Clarabel 1e-9
    P_opt = res_cl9.P
    x_arr = np.array(x, dtype=float)
    phi_arr = np.array(phi, dtype=float)
    p_vec = p_matrix_to_vector(np.array(P_opt))

    sheets, ref_sheet = compute_sheets(x_arr, phi_arr)
    h_eff = compute_h_eff(x_arr, phi_arr, sheets)
    T = build_scaled_basis_t(x_arr, phi_arr, sheets, ref_sheet, h_eff)

    L_w = build_lw_matrix(x_arr, phi_arr, p_vec)
    G_hat = (T.T @ L_w @ T) / (h_eff**2)

    # NumPy
    eigs_G_np = np.linalg.eigvalsh(G_hat)
    lam_np = float(np.min(eigs_G_np))

    # mpmath 50 digits
    mpmath.mp.dps = 50
    dim = G_hat.shape[0]
    G_mp = mpmath.mp.matrix(dim, dim)
    for r in range(dim):
        for c in range(dim):
            G_mp[r, c] = mpmath.mp.mpf(str(G_hat[r, c]))
    eigs_mp = mpmath.mp.eigsy(G_mp, eigvals_only=True)
    lam_mp = float(min(eigs_mp))

    # NumPy conditions
    A_F = build_af_matrix(x_arr, phi_arr)
    max_F = float(np.max(np.abs(A_F @ p_vec)))
    A_c = build_ac_matrix(x_arr, phi_arr)
    min_c = float(np.min(x_arr + A_c @ p_vec))
    A_d = build_ad_matrix(x_arr)
    min_d = float(np.min(x_arr + A_d @ p_vec))

    # T4 ratio
    diff_plain = max(abs(res_p_cl8.lam - res_p_cl9.lam), abs(res_p_cl8.lam - res_p_scs.lam))
    rel_plain = diff_plain / abs(res_p_cl8.lam)

    diff_scaled = max(abs(res_cl8.lam - res_cl9.lam), abs(res_cl8.lam - res_scs.lam))
    rel_scaled = diff_scaled / abs(res_cl8.lam)
    t4_ratio = rel_plain / rel_scaled

    return {
        "h_eff": h_eff,
        "cl8": {"lam": res_cl8.lam, "s": res_cl8.s, "f": res_cl8.f, "status": res_cl8.status},
        "cl9": {"lam": res_cl9.lam, "s": res_cl9.s, "f": res_cl9.f, "status": res_cl9.status},
        "scs": {"lam": res_scs.lam, "s": res_scs.s, "f": res_scs.f, "status": res_scs.status},
        "np": {"lam": lam_np, "max_F": max_F, "min_c": min_c, "min_d": min_d},
        "mp50": {"lam": lam_mp},
        "plain": {
            "cl8": res_p_cl8.lam,
            "cl9": res_p_cl9.lam,
            "scs": res_p_scs.lam,
            "rel_diff": rel_plain,
        },
        "scaled": {
            "cl8": res_cl8.lam,
            "cl9": res_cl9.lam,
            "scs": res_cl9.lam,
            "rel_diff": rel_scaled,
        },
        "t4_ratio": t4_ratio,
    }


def run_t9_benchmark():
    print("Running T9 benchmark for m=2..10 (15 individuals per m, no cache)...")
    rng = np.random.default_rng(2026)
    cfg = EvalConfig(mu_t=0.690, solver_tol=1e-8)
    ev_plain = SDPEvaluator(cfg)
    ev_scaled = ScaledSDPEvaluator(cfg)

    results_by_m = {}

    for m in range(2, 11):
        times_plain = []
        times_scaled = []
        samples_count = 15

        for _ in range(samples_count):
            kref = 4 if m <= 4 else 8
            k = tuple(int(i % kref) for i in range(m))
            sgn = tuple(int(rng.choice([-1, 1])) for _ in range(m))
            ell = tuple(float(rng.uniform(-3.0, -1.3)) for _ in range(m))
            u = tuple(float(rng.normal(0.0, 0.2)) for _ in range(m))

            g = canonicalize(Genome(m=m, K_ref=kref, k=k, sgn=sgn, ell=ell, u=u))
            x = decode_x(g.u)
            phi = decode_phi(g.K_ref, g.k, g.sgn, g.ell)

            # Plain
            t0 = time.perf_counter()
            ev_plain.evaluate_phase1(x, phi)
            times_plain.append((time.perf_counter() - t0) * 1000.0)

            # Scaled
            t0 = time.perf_counter()
            ev_scaled.evaluate_phase1(x, phi, genome=g)
            times_scaled.append((time.perf_counter() - t0) * 1000.0)

        mean_p = float(np.mean(times_plain))
        mean_s = float(np.mean(times_scaled))
        ratio = mean_s / mean_p if mean_p > 0 else 1.0

        results_by_m[m] = {
            "plain_ms": mean_p,
            "scaled_ms": mean_s,
            "ratio": ratio,
        }

    return results_by_m


def main():
    s7_res = get_s7_metrics()
    print("=" * 70)
    print("VALORES CRUDOS T4 Y T6 PARA S7 (h=0.002, t=0, mu_t=0.690):")
    print("=" * 70)
    print(f"h_eff de S7: {s7_res['h_eff']:.8e}")
    print(f"Clarabel(1e-8): lam_hat={s7_res['cl8']['lam']:.12e}, s*={s7_res['cl8']['s']:.12e}, f={s7_res['cl8']['f']:.12e}, status={s7_res['cl8']['status']}")
    print(f"Clarabel(1e-9): lam_hat={s7_res['cl9']['lam']:.12e}, s*={s7_res['cl9']['s']:.12e}, f={s7_res['cl9']['f']:.12e}, status={s7_res['cl9']['status']}")
    print(f"SCS(1e-8):      lam_hat={s7_res['scs']['lam']:.12e}, s*={s7_res['scs']['s']:.12e}, f={s7_res['scs']['f']:.12e}, status={s7_res['scs']['status']}")
    print(f"NumPy (float64): lam_min(G_hat)={s7_res['np']['lam']:.12e}, max|F|={s7_res['np']['max_F']:.12e}, min c_a={s7_res['np']['min_c']:.6f}, min d_a={s7_res['np']['min_d']:.6f}")
    print(f"mpmath (50d):    lam_min(G_hat)={s7_res['mp50']['lam']:.12e}")
    print(f"T4 Mejora de Acuerdo: {s7_res['t4_ratio']:.2f}x (rel_plain={s7_res['plain']['rel_diff']:.6e}, rel_scaled={s7_res['scaled']['rel_diff']:.6e})")
    print("=" * 70)

    t9_res = run_t9_benchmark()
    print("=" * 70)
    print("TABLA DE TIEMPOS T9 (ms/evaluación nueva, sin caché):")
    print("=" * 70)
    print(f"{'m':2s} | {'Plano (ms)':12s} | {'Escalado (ms)':14s} | {'Ratio (Esc/Plano)':18s}")
    print("-" * 55)
    for m in range(2, 11):
        r = t9_res[m]
        print(f"{m:2d} | {r['plain_ms']:10.2f} ms | {r['scaled_ms']:12.2f} ms | {r['ratio']:16.2f}x")
    print("=" * 70)


if __name__ == "__main__":
    main()
