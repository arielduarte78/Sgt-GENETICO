"""Analyze the 20 genome vs fallback discrepancies."""

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


def solve_scaled_phase1_custom_sheets(x, phi, sheets, ref_sheet, tol=1e-8, mu_t=0.690, rho=1.0, eps_c=1e-4):
    h_eff = compute_h_eff(x, phi, sheets)
    if h_eff < 1e-4 or h_eff > 0.1:
        return {"h_eff": h_eff, "status": "NOT_SCALABLE", "lam": None, "s": None}

    T = build_scaled_basis_t(x, phi, sheets, ref_sheet, h_eff)
    if T is None:
        return {"h_eff": h_eff, "status": "NO_BASIS_T", "lam": None, "s": None}

    m = len(x)
    K = m * (m - 1) // 2
    G_basis = build_g_hat_basis(x, phi, T, h_eff)

    A_F = build_af_matrix(x, phi)
    A_c = build_ac_matrix(x, phi)
    A_d = build_ad_matrix(x)

    N = scipy.linalg.null_space(A_F)
    res_inf = float(np.max(np.abs(A_F @ N)))
    if res_inf > 1e-14:
        corr, _, _, _ = np.linalg.lstsq(A_F, A_F @ N, rcond=None)
        N = N - corr
        N, _ = np.linalg.qr(N)

    d = N.shape[1]
    if d == 0:
        return {"h_eff": h_eff, "status": "NO_NULLSPACE", "lam": None, "s": None}

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

    kwargs = {
        "solver": cp.CLARABEL,
        "tol_gap_abs": tol,
        "tol_gap_rel": tol,
        "tol_feas": tol,
        "tol_infeas_abs": tol,
        "tol_infeas_rel": tol,
        "verbose": False,
    }
    try:
        prob.solve(**kwargs)
        return {
            "h_eff": h_eff,
            "status": prob.status,
            "lam": float(lam_var.value) if lam_var.value is not None else None,
            "s": float(s_var.value) if s_var.value is not None else None,
        }
    except Exception as e:
        return {"h_eff": h_eff, "status": "ERROR", "error": str(e), "lam": None, "s": None}


def main():
    rng = np.random.default_rng(42)
    num_ind = 200

    discrepancies = []

    for idx in range(num_ind):
        m = int(rng.integers(2, 11))
        kref = int(rng.choice([4, 8, 12, 16]))
        k = tuple(int(rng.integers(0, kref)) for _ in range(m))
        sgn = tuple(int(rng.choice([-1, 1])) for _ in range(m))
        ell = tuple(float(rng.uniform(-6.0, -0.5)) for _ in range(m))
        u = tuple(float(rng.normal(0.0, 1.0)) for _ in range(m))

        raw_genome = Genome(m=m, K_ref=kref, k=k, sgn=sgn, ell=ell, u=u)
        genome = canonicalize(raw_genome)

        x = np.array(decode_x(genome.u))
        phi = np.array(decode_phi(genome.K_ref, genome.k, genome.sgn, genome.ell))

        sheets_g, ref_g = compute_sheets(x, phi, genome=genome, K_ref=genome.K_ref)
        sheets_fb, ref_fb = compute_sheets(x, phi, genome=None, K_ref=genome.K_ref)

        if sheets_g != sheets_fb:
            discrepancies.append({
                "idx": idx,
                "genome": genome,
                "x": x,
                "phi": phi,
                "sheets_g": sheets_g,
                "ref_g": ref_g,
                "sheets_fb": sheets_fb,
                "ref_fb": ref_fb,
            })

    print(f"Total discrepancias encontradas: {len(discrepancies)} de {num_ind}")
    print("=" * 80)
    print(f"{'#':2s} | {'idx':3s} | {'m':2s} | {'Kref':4s} | {'h_eff(gen)':11s} | {'h_eff(fb)':11s} | {'lam(gen)':13s} | {'lam(fb)':13s} | {'sign_change?':12s}")
    print("-" * 80)

    sign_changes = 0
    for num, item in enumerate(discrepancies, 1):
        res_g = solve_scaled_phase1_custom_sheets(item["x"], item["phi"], item["sheets_g"], item["ref_g"])
        res_fb = solve_scaled_phase1_custom_sheets(item["x"], item["phi"], item["sheets_fb"], item["ref_fb"])

        lam_g = res_g["lam"]
        lam_fb = res_fb["lam"]
        h_g = res_g["h_eff"]
        h_fb = res_fb["h_eff"]

        str_lam_g = f"{lam_g:+.5e}" if lam_g is not None else res_g["status"]
        str_lam_fb = f"{lam_fb:+.5e}" if lam_fb is not None else res_fb["status"]

        sign_change = False
        if lam_g is not None and lam_fb is not None:
            if (lam_g > 0 and lam_fb < 0) or (lam_g < 0 and lam_fb > 0):
                sign_change = True
                sign_changes += 1

        sign_str = "SI (CAMBIA!)" if sign_change else "no"
        m = item["genome"].m
        kref = item["genome"].K_ref
        idx = item["idx"]

        print(f"{num:2d} | {idx:3d} | {m:2d} | {kref:4d} | {h_g:11.4e} | {h_fb:11.4e} | {str_lam_g:13s} | {str_lam_fb:13s} | {sign_str:12s}")

    print("=" * 80)
    print(f"Total casos con cambio de signo en lam_hat: {sign_changes}")


if __name__ == "__main__":
    main()
