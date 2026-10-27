"""Script de Reproducibilidad y Medición del Piso de Ruido Numérico.

Captura del entorno experimental y evaluación plana del Óptimo Par bajo barrido
de tolerancias de convergencia del solver SDP Clarabel:
- Captura de hash Git y versiones de bibliotecas del entorno.
- Configuración evaluada: Óptimo Par a h = 1e-3.
- Barrido de tolerancias Clarabel: [1e-8, 1e-9, 1e-10, 1e-11, 1e-12].
- LMI reescalada por h_eff^2.
- Certificación primal mpmath a 50 dígitos sobre P recortada.
- Estimación empírica del piso de ruido numérico (desviación estándar).
"""

from __future__ import annotations

import math
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")

from typing import TypedDict

import cvxpy as cp
import mpmath
import numpy as np
import scipy.linalg

from kuramoto_templates.cases.s7 import compute_q_star_and_w_star
from kuramoto_templates.domain.conditions import (
    build_ac_matrix,
    build_ad_matrix,
    build_af_matrix,
    build_lw_matrix,
    get_pair_indices,
    p_matrix_to_vector,
    p_vector_to_matrix,
)
from kuramoto_templates.evaluators.scaled import (
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)

q_star, w_star = compute_q_star_and_w_star(dps=50)
a_star = float(1.0 / (7.0 + q_star))
x_W_base = float(1.0 - 6.0 * a_star)

h = 0.001
half_pi = math.pi / 2.0
pi = math.pi

phi = np.array([
    pi - w_star * h,      # 0
    -half_pi + h,         # 1
    -half_pi - h,         # 2
    0.0,                  # 3 (W fijo en 0)
    half_pi + h,          # 4
    half_pi - h,          # 5
    -pi + w_star * h,     # 6
], dtype=float)

# Masas del Optimo Par a h = 1e-3
dx_0_6 = 1.112791e-05
dx_1_5 = 9.353764e-05
dx_2_4 = -1.348106e-04
dx_W = -2.0 * (dx_0_6 + dx_1_5 + dx_2_4)

xc_opt = np.array([
    a_star + dx_0_6,
    a_star + dx_1_5,
    a_star + dx_2_4,
    x_W_base + dx_W,
    a_star + dx_2_4,
    a_star + dx_1_5,
    a_star + dx_0_6,
], dtype=float)


def compute_mpmath_M(xc: np.ndarray, phic: np.ndarray, P_mat: np.ndarray) -> float:
    """Calcula el autovalor minimo fisico lam_min(M) en mpmath a 50 digitos."""
    m_len = len(xc)
    p_v = p_matrix_to_vector(np.array(P_mat, dtype=float))
    L_w = build_lw_matrix(xc, phic, p_v)
    Q_np = scipy.linalg.null_space(np.ones((1, m_len), dtype=float))
    dim = m_len - 1
    mpmath.mp.dps = 50
    Q_mp = mpmath.mp.matrix(m_len, dim)
    for r in range(m_len):
        for c in range(dim):
            Q_mp[r, c] = mpmath.mp.mpf(str(Q_np[r, c]))
    for j in range(dim):
        v = Q_mp[:, j]
        for k in range(j):
            proj = mpmath.mp.fdot(Q_mp[:, k], v)
            for r in range(m_len):
                v[r] -= proj * Q_mp[r, k]
        norm_v = mpmath.sqrt(sum(v[r] ** 2 for r in range(m_len)))
        for r in range(m_len):
            Q_mp[r, j] = v[r] / norm_v
    L_mp = mpmath.mp.matrix(m_len, m_len)
    for r in range(m_len):
        for c in range(m_len):
            L_mp[r, c] = mpmath.mp.mpf(str(L_w[r, c]))
    M_mp = Q_mp.T * L_mp * Q_mp
    eigsy_fn = getattr(mpmath.mp, "eigsy")
    eigs = eigsy_fn(M_mp, eigvals_only=True)
    return float(min(eigs))


def solve_exact(
    xc: np.ndarray,
    phic: np.ndarray,
    lam_M: float = 1e-10,
    tol: float = 1e-10,
) -> tuple[float | None, float | None, str, dict[str, float]]:
    """Resuelve la LMI reescalada y certifica primalmente en mpmath."""
    if np.min(xc) <= 1e-4:
        return None, None, "INFEASIBLE_MASS", {}

    m = len(xc)
    pairs = get_pair_indices(m)
    K = len(pairs)
    sheets, ref_sheet = compute_sheets(xc, phic)
    h_eff = float(compute_h_eff(xc, phic, sheets))
    if h_eff < 1e-6:
        return None, None, "SCALE_TOO_SMALL", {}

    T = build_scaled_basis_t(xc, phic, sheets, ref_sheet, h_eff)
    if T is None:
        return None, None, "NOT_SCALABLE", {}

    P_perp = np.eye(m) - np.ones((m, m)) / float(m)
    T_prime = P_perp @ T
    metric_M_rescaled = (T_prime.T @ T_prime) / (h_eff**2)

    A_F = build_af_matrix(xc, phic)
    A_c = build_ac_matrix(xc, phic)
    A_d = build_ad_matrix(xc)
    N = scipy.linalg.null_space(A_F)
    if N.shape[1] == 0:
        return None, None, "NO_NULLSPACE", {}

    res_inf = float(np.max(np.abs(A_F @ N)))
    if res_inf > 1e-14:
        corr, _, _, _ = np.linalg.lstsq(A_F, A_F @ N, rcond=None)
        N = N - corr
        N, _ = np.linalg.qr(N)

    G_rescaled_basis = []
    for a, b in pairs:
        d_ab = np.zeros(m, dtype=float)
        d_ab[a] = 1.0
        d_ab[b] = -1.0
        Td = T.T @ d_ab
        weight = xc[a] * xc[b] * math.cos(phic[a] - phic[b])
        Gk = weight * np.outer(Td, Td) / (h_eff**2)
        G_rescaled_basis.append(0.5 * (Gk + Gk.T))

    z = cp.Variable(N.shape[1])
    t = cp.Variable()
    p = N @ z
    L_expr = p[0] * G_rescaled_basis[0]
    for k in range(1, K):
        L_expr = L_expr + p[k] * G_rescaled_basis[k]

    prob = cp.Problem(
        cp.Maximize(t),
        [
            xc + A_d @ p >= t,
            p >= 0.0,
            p <= 1.0,
            L_expr >> lam_M * metric_M_rescaled,
            xc + A_c @ p >= 1e-4,
        ],
    )
    try:
        prob.solve(
            solver=cp.CLARABEL,
            tol_feas=tol,
            tol_gap_abs=tol,
            tol_gap_rel=tol,
            verbose=False,
        )
    except Exception:
        return None, None, "SOLVER_EXCEPTION", {}

    extra_info: dict[str, float] = {}
    if hasattr(prob, "solver_stats") and prob.solver_stats is not None:
        extra_info["solve_time"] = float(getattr(prob.solver_stats, "solve_time", 0.0) or 0.0)
        extra_info["iterations"] = float(getattr(prob.solver_stats, "num_iters", 0.0) or 0.0)

    if prob.status != "optimal":
        return None, None, prob.status, extra_info

    p_clip = np.clip(N @ z.value, 0.0, 1.0)
    mu_cert = float(np.min(xc + A_d @ p_clip))
    P_mat = p_vector_to_matrix(p_clip, m)
    lam_mp = compute_mpmath_M(xc, phic, P_mat)
    return mu_cert, lam_mp, prob.status, extra_info


class TolResultItem(TypedDict):
    tol: float
    status: str
    mu_cert: float | None
    lam_mp: float | None
    solve_time: float
    iterations: float


def main() -> None:
    print("=" * 105)
    print("SCRIPT DE REPRODUCIBILIDAD Y MEDICION DEL PISO DE RUIDO NUMERICO")
    print("Evaluacion plana del Optimo Par a h = 1e-3 bajo barrido de tolerancias de Clarabel")
    print("=" * 105)

    # 1. Captura del Entorno Experimental
    try:
        git_hash = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        git_hash = "NO_GIT_REPO"

    try:
        git_branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        git_branch = "UNKNOWN"

    cp_version = getattr(cp, "__version__", "unknown")
    print("\n1. ENTORNO EXPERIMENTAL Y REPRODUCIBILIDAD:")
    print(f"  Git Commit Hash     : {git_hash}")
    print(f"  Git Branch          : {git_branch}")
    print(f"  Python Version      : {sys.version.split()[0]} ({sys.platform})")
    print(f"  CVXPY Version       : {cp_version}")
    print(f"  NumPy Version       : {np.__version__}")
    print(f"  SciPy Version       : {scipy.__version__}")
    print(f"  mpmath Version      : {mpmath.__version__}")
    print(f"  Solver LMI          : CLARABEL (Rust SDP Engine)")

    # 2. Parametros Evaluados
    print("\n2. CONFIGURACION EVALUADA (Optimo Par a h = 1e-3, delta_W = 0.0):")
    print(f"  dx_0_6              : {dx_0_6:+.6e}")
    print(f"  dx_1_5              : {dx_1_5:+.6e}")
    print(f"  dx_2_4              : {dx_2_4:+.6e}")
    print(f"  dx_W                : {dx_W:+.6e}")
    print(f"  Vector x*           : [{', '.join(f'{v:.8f}' for v in xc_opt)}]")
    print(f"  Suma sum(x*)        : {np.sum(xc_opt):.16f}")
    print(f"  Requisito LMI lam_M : 1.00e-10")

    # 3. Barrido de Tolerancias
    tolerances = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
    results: list[TolResultItem] = []

    for tol_val in tolerances:
        mu_c, lam_c, st_c, info = solve_exact(xc_opt, phi, lam_M=1e-10, tol=tol_val)
        results.append({
            "tol": tol_val,
            "status": st_c,
            "mu_cert": mu_c,
            "lam_mp": lam_c,
            "solve_time": info.get("solve_time", 0.0),
            "iterations": info.get("iterations", 0.0),
        })

    # 4. Tabla de Resultados
    print(f"\n3. BARRIDO DE TOLERANCIAS DE CONVERGENCIA CLARABEL:")
    print(f"{'Tolerancia':<12} | {'Estado Solver':<16} | {'mu_cert (mpmath 50d)':<22} | {'lam_min(M) Realizado':<22} | {'Iter':<6} | {'Tiempo (s)'}")
    print("-" * 95)
    valid_mus: list[float] = []
    for r in results:
        mu_str = f"{r['mu_cert']:.10f}" if r["mu_cert"] is not None else "N/A"
        lam_str = f"{r['lam_mp']:.4e}" if r["lam_mp"] is not None else "N/A"
        print(f"{r['tol']:<12.1e} | {r['status']:<16} | {mu_str:<22} | {lam_str:<22} | {r['iterations']:<6.0f} | {r['solve_time']:.4f}")
        if r["status"] == "optimal" and r["mu_cert"] is not None:
            valid_mus.append(r["mu_cert"])
    print("-" * 95)

    # 5. Estadisticas de Ruido Numerico (Barra de Error)
    if len(valid_mus) >= 2:
        mu_arr = np.array(valid_mus, dtype=np.float64)
        mean_mu = float(np.mean(mu_arr))
        std_mu = float(np.std(mu_arr, ddof=1))
        range_mu = float(np.max(mu_arr) - np.min(mu_arr))
        
        print("\n4. ANALISIS DEL PISO DE RUIDO NUMERICO (BARRA DE ERROR):")
        print(f"  Puntos optimos validos   : {len(valid_mus)} / {len(tolerances)}")
        print(f"  mu_cert Promedio         : {mean_mu:.12f}")
        print(f"  Desviacion Estandar (s)  : {std_mu:.4e} ({std_mu * 1e4:.4e}e-4)")
        print(f"  Rango Maximo - Minimo    : {range_mu:.4e} ({range_mu * 1e4:.4e}e-4)")
        print(f"  Incertidumbre Absoluta   : +/- {std_mu:.2e}")
        print(f"  Ganancia Reportada       : ~ +7.11e-4")
        print(f"  Relacion Senal / Ruido   : {((7.11e-4) / (std_mu if std_mu > 0 else 1e-15)):.1f}x")
    else:
        print("\n4. ANALISIS DEL PISO DE RUIDO: Puntos optimos insuficientes para desviacion estandar.")

    print("=" * 105)


if __name__ == "__main__":
    main()
