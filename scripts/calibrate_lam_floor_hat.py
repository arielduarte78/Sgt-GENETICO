"""Calibrate lam_floor_hat as requested:
Measure on 200 random individuals near C4 (h_eff in [1e-3, 5e-2]) the maximum
disagreement of lam_hat between Clarabel(1e-8), Clarabel(1e-9) and SCS(1e-8);
lam_floor_hat = 100 * max_disagreement.
"""

from __future__ import annotations

import math
import sys
import time
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
from kuramoto_templates.search.genome import Genome, canonicalize, decode_x, decode_phi


def compute_sheets(x, phi, genome=None, K_ref=4, sheet_tol=0.2):
    m = len(x)
    if genome is not None and hasattr(genome, "k") and len(genome.k) == m:
        sheets = {}
        for i, anchor in enumerate(genome.k):
            sheets.setdefault(int(anchor), []).append(i)
    else:
        sheets = {}
        for i, ph in enumerate(phi):
            dists = [abs((ph - (2 * np.pi * a / K_ref) + np.pi) % (2 * np.pi) - np.pi) for a in range(K_ref)]
            min_a = int(np.argmin(dists))
            if dists[min_a] <= sheet_tol:
                sheets.setdefault(min_a, []).append(i)
            else:
                sheets.setdefault(f"sing_{i}", []).append(i)
    masses = {s: sum(x[i] for i in idxs) for s, idxs in sheets.items()}
    ref_sheet = max(masses, key=masses.get)
    return sheets, ref_sheet


def compute_h_eff(x, phi, sheets):
    devs = []
    for s, idxs in sheets.items():
        if len(idxs) == 0:
            continue
        sin_m = sum(x[i] * np.sin(phi[i]) for i in idxs)
        cos_m = sum(x[i] * np.cos(phi[i]) for i in idxs)
        mean_ph = np.arctan2(sin_m, cos_m)
        for i in idxs:
            d = (phi[i] - mean_ph + np.pi) % (2 * np.pi) - np.pi
            devs.append(x[i] * (d**2))
    return float(np.sqrt(sum(devs)))


def build_scaled_basis_t(x, phi, sheets, ref_sheet, h_eff):
    m = len(x)
    eta_cols = []
    for s, idxs in sorted(sheets.items(), key=lambda item: str(item[0])):
        if s == ref_sheet:
            continue
        col = np.zeros(m)
        col[idxs] = 1.0
        eta_cols.append(col)
    zeta_cols = []
    for s, idxs in sorted(sheets.items(), key=lambda item: str(item[0])):
        if len(idxs) >= 2:
            Z = scipy.linalg.null_space(np.ones((1, len(idxs))))
            for col_z in Z.T:
                full_col = np.zeros(m)
                full_col[idxs] = col_z * h_eff
                zeta_cols.append(full_col)
    cols = eta_cols + zeta_cols
    if len(cols) != m - 1:
        return None
    T = np.column_stack(cols)
    if np.linalg.matrix_rank(np.column_stack([T, np.ones(m)])) < m:
        return None
    return T


def build_g_hat_basis(x, phi, T, h_eff):
    m = len(x)
    pairs = get_pair_indices(m)
    G_basis = []
    for a, b in pairs:
        d_ab = np.zeros(m)
        d_ab[a] = 1.0
        d_ab[b] = -1.0
        Td = T.T @ d_ab
        scale = (x[a] * x[b] * math.cos(phi[a] - phi[b])) / (h_eff**2)
        Gk = scale * np.outer(Td, Td)
        Gk = 0.5 * (Gk + Gk.T)
        G_basis.append(Gk)
    return G_basis


def solve_scaled_phase1(x, phi, genome, solver_name, tol, mu_t=0.690, rho=1.0, eps_c=1e-4):
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
        # Refinement step
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
        if prob.status in ("optimal", "optimal_inaccurate"):
            return float(lam_var.value)
    except Exception:
        pass
    return None


def run_calibration():
    print("Iniciando calibración de lam_floor_hat sobre 200 individuos válidos...")
    rng = np.random.default_rng(2026)
    valid_count = 0
    max_diff_overall = 0.0
    diffs = []

    t0 = time.time()
    while valid_count < 200:
        # Sample near C4: m in [4, 8], K_ref in [4, 8]
        m = int(rng.integers(4, 8))
        kref = 4
        # distribute clusters around anchors 0, 1, 2, 3 with small deviation
        k = tuple(int(i % kref) for i in range(m))
        sgn = tuple(int(rng.choice([-1, 1])) for _ in range(m))
        # ell to produce h_eff in [1e-3, 5e-2]: 10^ell in [1e-3, 5e-2] -> ell in [-3.0, -1.3]
        ell = tuple(float(rng.uniform(-3.0, -1.3)) for _ in range(m))
        u = tuple(float(rng.normal(0.0, 0.2)) for _ in range(m))

        genome = canonicalize(Genome(m=m, K_ref=kref, k=k, sgn=sgn, ell=ell, u=u))
        x = np.array(decode_x(genome.u))
        phi = np.array(decode_phi(genome.K_ref, genome.k, genome.sgn, genome.ell))

        sheets, _ = compute_sheets(x, phi, genome=genome)
        h_eff = compute_h_eff(x, phi, sheets)
        if not (1e-3 <= h_eff <= 5e-2):
            continue

        lam_cl8 = solve_scaled_phase1(x, phi, genome, "CLARABEL", tol=1e-8)
        if lam_cl8 is None:
            continue
        lam_cl9 = solve_scaled_phase1(x, phi, genome, "CLARABEL", tol=1e-9)
        if lam_cl9 is None:
            continue
        lam_scs = solve_scaled_phase1(x, phi, genome, "SCS", tol=1e-8)
        if lam_scs is None:
            continue

        diff1 = abs(lam_cl8 - lam_cl9)
        diff2 = abs(lam_cl8 - lam_scs)
        diff3 = abs(lam_cl9 - lam_scs)
        max_d = max(diff1, diff2, diff3)
        diffs.append(max_d)
        if max_d > max_diff_overall:
            max_diff_overall = max_d

        valid_count += 1
        if valid_count % 25 == 0:
            print(f"  [{valid_count}/200] max_diff actual: {max_diff_overall:.6e} (tiempo: {time.time()-t0:.1f}s)")

    lam_floor_hat = 100.0 * max_diff_overall
    print("=" * 60)
    print("RESULTADO DE CALIBRACIÓN DE lam_floor_hat:")
    print(f"  Muestras analizadas: {valid_count}")
    print(f"  Máximo desacuerdo Clarabel(1e-8) vs Clarabel(1e-9) vs SCS(1e-8): {max_diff_overall:.6e}")
    print(f"  Percentil 95 de desacuerdos: {np.percentile(diffs, 95):.6e}")
    print(f"  Percentil 99 de desacuerdos: {np.percentile(diffs, 99):.6e}")
    print(f"  lam_floor_hat calibrado (100 x max): {lam_floor_hat:.6e}")
    print("=" * 60)


if __name__ == "__main__":
    run_calibration()
