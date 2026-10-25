"""Analyze piso and outliers between Clarabel(1e-8), Clarabel(1e-9) and SCS(1e-8)."""

from __future__ import annotations

import math
import sys
from pathlib import Path
import numpy as np
import scipy.linalg
import cvxpy as cp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.domain.conditions import (
    build_af_matrix,
    build_ac_matrix,
    build_ad_matrix,
    get_pair_indices,
)
from kuramoto_templates.evaluators.scaled import (
    compute_sheets,
    compute_h_eff,
    build_scaled_basis_t,
    build_g_hat_basis,
)
from kuramoto_templates.search.genome import Genome, canonicalize, decode_x, decode_phi


def solve_full_info(x, phi, genome, solver_name, tol, mu_t=0.690, rho=1.0, eps_c=1e-4):
    sheets, ref_sheet = compute_sheets(x, phi, genome=genome)
    h_eff = compute_h_eff(x, phi, sheets)
    if h_eff < 1e-4 or h_eff > 0.1:
        return None

    T = build_scaled_basis_t(x, phi, sheets, ref_sheet, h_eff)
    if T is None:
        return None

    m = len(x)
    K = m * (m - 1) // 2
    G_basis = build_g_hat_basis(x, phi, T, h_eff)

    A_F = build_af_matrix(x, phi)
    A_c = build_ac_matrix(x, phi)
    A_d = build_ad_matrix(x)

    N = scipy.linalg.null_space(A_F)
    res_inf = np.max(np.abs(A_F @ N))
    if res_inf > 1e-14:
        corr, _, _, _ = np.linalg.lstsq(A_F, A_F @ N, rcond=None)
        N = N - corr
        N, _ = np.linalg.qr(N)

    d = N.shape[1]
    if d == 0:
        return None

    z = cp.Variable(d)
    lam_var = cp.Variable()
    s_var = cp.Variable(nonneg=True)
    p_z = N @ z

    G_z = sum(p_z[k] * G_basis[k] for k in range(K))
    obj = cp.Maximize(lam_var - rho * s_var)
    constraints = [
        p_z >= 0.0,
        p_z <= 1.0,
        G_z >> lam_var * np.eye(m - 1),
        x + A_c @ p_z >= eps_c,
        x + A_d @ p_z >= mu_t - s_var,
        s_var >= 0.0,
    ]
    prob = cp.Problem(obj, constraints)

    if solver_name == "CLARABEL":
        kwargs = {
            "solver": cp.CLARABEL,
            "tol_gap_abs": tol,
            "tol_gap_rel": tol,
            "tol_feas": tol,
            "tol_infeas_abs": tol,
            "tol_infeas_rel": tol,
            "verbose": False,
        }
    else:  # SCS
        kwargs = {
            "solver": cp.SCS,
            "eps_abs": tol,
            "eps_rel": tol,
            "verbose": False,
        }

    try:
        prob.solve(**kwargs)
        return {
            "status": prob.status,
            "lam": float(lam_var.value) if lam_var.value is not None else None,
            "s": float(s_var.value) if s_var.value is not None else None,
            "obj": float(prob.value) if prob.value is not None else None,
        }
    except Exception as e:
        return {"status": "ERROR", "error": str(e), "lam": None, "s": None, "obj": None}


def main():
    rng = np.random.default_rng(2026)
    valid_count = 0
    records = []

    print("Recolectando datos detallados de 200 individuos...")

    while valid_count < 200:
        m = int(rng.integers(4, 8))
        kref = 4
        k = tuple(int(i % kref) for i in range(m))
        sgn = tuple(int(rng.choice([-1, 1])) for _ in range(m))
        ell = tuple(float(rng.uniform(-3.0, -1.3)) for _ in range(m))
        u = tuple(float(rng.normal(0.0, 0.2)) for _ in range(m))

        genome = canonicalize(Genome(m=m, K_ref=kref, k=k, sgn=sgn, ell=ell, u=u))
        x = np.array(decode_x(genome.u))
        phi = np.array(decode_phi(genome.K_ref, genome.k, genome.sgn, genome.ell))

        sheets, _ = compute_sheets(x, phi, genome=genome)
        h_eff = compute_h_eff(x, phi, sheets)
        if not (1e-3 <= h_eff <= 5e-2):
            continue

        r_cl8 = solve_full_info(x, phi, genome, "CLARABEL", tol=1e-8)
        if r_cl8 is None or r_cl8["lam"] is None:
            continue
        r_cl9 = solve_full_info(x, phi, genome, "CLARABEL", tol=1e-9)
        if r_cl9 is None or r_cl9["lam"] is None:
            continue
        r_scs = solve_full_info(x, phi, genome, "SCS", tol=1e-8)
        if r_scs is None or r_scs["lam"] is None:
            continue

        diff_cl = abs(r_cl8["lam"] - r_cl9["lam"])
        diff_scs8 = abs(r_cl8["lam"] - r_scs["lam"])
        diff_scs9 = abs(r_cl9["lam"] - r_scs["lam"])
        max_d = max(diff_cl, diff_scs8, diff_scs9)

        rec = {
            "idx": valid_count,
            "m": m,
            "h_eff": h_eff,
            "x": x,
            "phi": phi,
            "genome": genome,
            "cl8": r_cl8,
            "cl9": r_cl9,
            "scs": r_scs,
            "diff_cl": diff_cl,
            "diff_scs8": diff_scs8,
            "max_d": max_d,
        }
        records.append(rec)
        valid_count += 1

    # 1. Piso
    overall_max_diff = max(r["max_d"] for r in records)
    lam_floor_hat_raw = 100.0 * overall_max_diff

    print("=" * 70)
    print("ITEM 1: PISO RECALCULADO SIN DESCARTAR OUTLIERS")
    print("=" * 70)
    print(f"Max desacuerdo absoluto overall: {overall_max_diff:.12e}")
    print(f"lam_floor_hat = 100 * max_desacuerdo: {lam_floor_hat_raw:.12e}")

    # Individuos con desacuerdo > 1e-5
    above_1e5 = [r for r in records if r["max_d"] > 1e-5]
    print(f"\nCantidad de individuos con desacuerdo > 1e-5: {len(above_1e5)} de {len(records)} ({len(above_1e5)/len(records)*100:.1f}%)")
    print("\nCaracterización de individuos con desacuerdo > 1e-5:")
    for r in above_1e5:
        print(f"  idx={r['idx']:3d} | m={r['m']} | h_eff={r['h_eff']:.6e} | diff_scs={r['diff_scs8']:.6e} | diff_cl={r['diff_cl']:.6e} | "
              f"status_cl8={r['cl8']['status']} | status_scs={r['scs']['status']}")

    # 2. Los 5 peores outliers (por diff_scs8 o max_d)
    print("\n" + "=" * 70)
    print("ITEM 2: EXPLICACIÓN DE LOS PEORES OUTLIERS (Clarabel vs SCS)")
    print("=" * 70)
    worst_5 = sorted(records, key=lambda r: r["diff_scs8"], reverse=True)[:5]
    for rank, r in enumerate(worst_5, 1):
        print(f"\n--- Outlier #{rank} (idx={r['idx']}, diff={r['diff_scs8']:.8e}) ---")
        print(f"  m: {r['m']}, h_eff: {r['h_eff']:.8e}")
        print(f"  x:   tuple({list(r['x'])})")
        print(f"  phi: tuple({list(r['phi'])})")
        print(f"  Clarabel(1e-8): lam={r['cl8']['lam']:.12e}, s={r['cl8']['s']:.12e}, obj={r['cl8']['obj']:.12e}, status={r['cl8']['status']}")
        print(f"  Clarabel(1e-9): lam={r['cl9']['lam']:.12e}, s={r['cl9']['s']:.12e}, obj={r['cl9']['obj']:.12e}, status={r['cl9']['status']}")
        print(f"  SCS(1e-8):      lam={r['scs']['lam']:.12e}, s={r['scs']['s']:.12e}, obj={r['scs']['obj']:.12e}, status={r['scs']['status']}")
        print(f"  Discrepancia Cl(1e-8) vs Cl(1e-9): {r['diff_cl']:.8e}")


if __name__ == "__main__":
    main()
