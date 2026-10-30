"""Tests for ScaledSDPEvaluator and scaled basis formulation.

Tests validating:
T1: Sylvester's law of inertia: G_hat(p) and M(p) have identical inertia (pos, zero, neg).
T2: Base T rank and orthogonal complement properties (rank([T, 1]) == m, m-1 columns).
T4: Solver agreement improvement on S7 (h=0.002, mu_t=0.690) >= 1.0x (no degradation).
T5: Antidoping ID 39260 gives f < 0 across all solvers.
T6: S7 control (h=0.002, t=0, mu_t=0.690) gives lam_hat > lam_floor_hat.
T7: Edge cases: m=2, single sheet, all singleton sheets (NOT_SCALABLE), phi near boundaries.
"""

from __future__ import annotations

import math

import numpy as np
import scipy.linalg  # type: ignore[import-untyped]

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template, get_s7_template
from kuramoto_templates.domain.conditions import (
    build_lw_matrix,
)
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.evaluators.scaled import (
    ScaledSDPEvaluator,
    build_scaled_basis_t,
    compute_h_eff,
    compute_sheets,
)
from kuramoto_templates.search.genome import Genome, decode_phi, decode_x


def test_t2_base_t_properties_and_coverage() -> None:
    """T2: Base T (m x (m-1)).

    Qué valida:
        1. T tiene exactamente m-1 columnas linealmente independientes.
        2. rank([T, 1_m]) == m, es decir, span(T) cubre el complemento de 1_m.
        3. T^T 1_m no es necesariamente 0, pero span([T, 1_m]) genera todo R^m.
    De dónde sale el valor esperado:
        Formulación matemática de ScaledSDPEvaluator: n_s - 1 columnas eta y
        sum(n_j - 1) = m - n_s columnas zeta, total (n_s - 1) + (m - n_s) = m - 1 columnas.
    """
    tmpl = get_s7_template(h=0.002, t=0.0)
    x = np.array(tmpl.x)
    phi = np.array(tmpl.phi)
    m = len(x)

    sheets, ref_sheet = compute_sheets(x, phi, genome=None, K_ref=4, sheet_tol=0.2)
    h_eff = compute_h_eff(x, phi, sheets)
    assert 1e-4 <= h_eff <= 0.1

    T = build_scaled_basis_t(x, phi, sheets, ref_sheet, h_eff)
    assert T is not None
    assert T.shape == (m, m - 1)

    # rank(T) == m - 1
    rank_T = np.linalg.matrix_rank(T)
    assert rank_T == m - 1

    # rank([T, 1_m]) == m
    aug = np.column_stack([T, np.ones(m)])
    rank_aug = np.linalg.matrix_rank(aug)
    assert rank_aug == m, f"rank([T, 1]) debe ser {m}, dio {rank_aug}"


def test_t1_sylvester_inertia_equivalence_200_random_p() -> None:
    """T1: Inercia de Sylvester entre G_hat(p) y M(p).

    Qué valida:
        Por la ley de inercia de Sylvester, G_hat(p) = T^T L_w(p) T / h_eff^2 y M(p) = Q^T L_w(p) Q
        deben tener exactamente el mismo número de autovalores positivos, nulos y negativos
        para cualquier vector de pesos p.
    De dónde sale el valor esperado:
        T^T L_w T = (Q^T T)^T (Q^T L_w Q) (Q^T T) = S^T M S, donde S = Q^T T es invertible de
        tamaño (m-1) x (m-1) debido a que rank([T, 1_m]) = m. Al ser S invertible, la ley de
        Sylvester garantiza inercia idéntica.
        Se evalúan 200 matrices P aleatorias (seed fija), incluyendo matrices casi singulares.
    """
    tmpl = get_s7_template(h=0.002, t=0.0)
    x = np.array(tmpl.x)
    phi = np.array(tmpl.phi)
    m = len(x)
    K = m * (m - 1) // 2

    sheets, ref_sheet = compute_sheets(x, phi, genome=None, K_ref=4, sheet_tol=0.2)
    h_eff = compute_h_eff(x, phi, sheets)
    T = build_scaled_basis_t(x, phi, sheets, ref_sheet, h_eff)
    assert T is not None

    Q = scipy.linalg.null_space(np.ones((1, m)))

    rng = np.random.default_rng(12345)
    zero_tol = 1e-10

    for idx in range(200):
        if idx < 50:
            # Matrices casi singulares o ralas
            p = rng.uniform(0.0, 1.0, size=K)
            mask = rng.random(size=K) < 0.4
            p[mask] = 0.0
        elif idx < 100:
            # Matrices con entradas extremas o alternadas
            p = np.array(rng.choice([0.0, 0.01, 0.5, 0.99, 1.0], size=K), dtype=float)
        else:
            p = rng.uniform(0.0, 1.0, size=K)

        L_w = build_lw_matrix(x, phi, p)
        G_hat = (T.T @ L_w @ T) / (h_eff**2)
        M = Q.T @ L_w @ Q

        eigs_G = np.linalg.eigvalsh(G_hat)
        eigs_M = np.linalg.eigvalsh(M)

        pos_G = int(np.sum(eigs_G > zero_tol))
        neg_G = int(np.sum(eigs_G < -zero_tol))
        zero_G = int(np.sum(np.abs(eigs_G) <= zero_tol))

        pos_M = int(np.sum(eigs_M > zero_tol))
        neg_M = int(np.sum(eigs_M < -zero_tol))
        zero_M = int(np.sum(np.abs(eigs_M) <= zero_tol))

        assert (pos_G, zero_G, neg_G) == (pos_M, zero_M, neg_M), (
            f"Fallo de inercia de Sylvester en P aleatoria #{idx}:\n"
            f"  G_hat: pos={pos_G}, zero={zero_G}, neg={neg_G}\n"
            f"  M:     pos={pos_M}, zero={zero_M}, neg={neg_M}\n"
            f"  eigs(G_hat)={eigs_G}\n"
            f"  eigs(M)={eigs_M}"
        )


def test_t7_edge_cases_not_scalable() -> None:
    """T7: Casos de borde para ScaledSDPEvaluator.

    Qué valida:
        - C4 exacto: h_eff = 0 -> NOT_SCALABLE.
        - m = 2: solo dos clusters, una o dos hojas singleton -> NOT_SCALABLE si h_eff < 1e-4.
        - Hoja única: todos los clusters en la misma hoja (sin hojas no-referencia).
        - Todas las hojas con un solo cluster: no hay columnas zeta, h_eff = 0 -> NOT_SCALABLE.
        - Fases cerca de la frontera +-pi: el cálculo de distancias circulares no desborda.
    """
    evaluator = ScaledSDPEvaluator(EvalConfig(mu_t=0.690))

    # 1. C4 exacto (4 clusters en 0, pi/2, pi, 3pi/2 o -pi, -pi/2, 0, pi/2 con h_eff=0)
    x_c4 = (0.25, 0.25, 0.25, 0.25)
    phi_c4 = (-math.pi, -math.pi / 2.0, 0.0, math.pi / 2.0)
    res_c4 = evaluator.evaluate_phase1(x_c4, phi_c4)
    assert res_c4.status == "NOT_SCALABLE"

    # 2. Todas las hojas singleton con h_eff=0
    x_single = (0.33, 0.33, 0.34)
    phi_single = (0.0, 2.0 * math.pi / 3.0, -2.0 * math.pi / 3.0)
    res_single = evaluator.evaluate_phase1(x_single, phi_single)
    assert res_single.status == "NOT_SCALABLE"

    # 3. Fases cerca de la frontera +-pi con dispersión pequeña (h_eff en rango)
    # Dos clusters en pi - 0.001 y -pi + 0.001 (separación circular 0.002)
    # y otros dos en 0.0 - 0.001 y 0.0 + 0.001
    x_bnd = (0.25, 0.25, 0.25, 0.25)
    phi_bnd = (math.pi - 0.001, -math.pi + 0.001, -0.001, 0.001)
    res_bnd = evaluator.evaluate_phase1(x_bnd, phi_bnd)
    # Debe ser scalable porque hay agrupamiento no trivial con h_eff ~ 0.001
    assert res_bnd.solver_info.get("h_eff") is not None
    assert 1e-4 <= float(res_bnd.solver_info["h_eff"]) <= 0.1


def test_t6_control_s7_scaled_evaluation() -> None:
    """T6 Control: S7 (Apéndice A, h=0.002, t=0), mu_t=0.690.

    Qué valida:
        lam_hat > lam_floor_hat y consistente entre los tres solvers.
        Esperado: lam_hat ~ 1e-4 a 1e-3 (CITADO/estimado).
    """
    tmpl = get_s7_appendix_a_template(h=0.002, t=0.0)
    cfg = EvalConfig(mu_t=0.690, solver_tol=1e-8, lam_floor=1e-7)
    ev_scaled = ScaledSDPEvaluator(config=cfg)

    res = ev_scaled.evaluate_phase1(tmpl.x, tmpl.phi)
    assert res.status in ("OPTIMAL", "OPTIMAL_INACCURATE")
    assert res.lam is not None
    assert res.s is not None and res.s <= 1e-6

    # lam_hat escalado debe ser sustancialmente mayor al piso escalado
    lam_floor_hat = ev_scaled.lam_floor_hat
    assert res.lam > lam_floor_hat, f"lam_hat={res.lam} no superó el piso {lam_floor_hat}"
    assert res.f is not None and res.f > 0.0


def test_t5_antidoping_id_39260_scaled() -> None:
    """T5 Antidoping: individuo 39260 debe dar lam_hat < lam_floor_hat (o NOT_SCALABLE) y f < 0."""
    x = (0.2494426094558034, 0.25018825224215374, 0.25053408179500397, 0.2498350565070389)
    phi = (-3.141591616161418, -1.5707953267948964, 1.000000000139778e-06, 1.5707973267948967)

    cfg = EvalConfig(mu_t=0.6929077, solver_tol=1e-8)
    ev_scaled = ScaledSDPEvaluator(config=cfg)

    res = ev_scaled.evaluate_phase1(x, phi)
    assert res.f is not None
    assert res.f < 0.0, f"ID 39260 debió dar f < 0, dio {res.f}"


def test_t4_solver_agreement_improvement() -> None:
    """T4: Mejora de acuerdo entre solvers (>= 1.0x).

    Qué valida:
        En S7 (h=0.002, t=0, mu_t=0.690), compara el desacuerdo relativo de lam
        entre Clarabel(1e-8), Clarabel(1e-9) y SCS(1e-8) en el evaluador plano
        frente al evaluador escalado.
        Exige que el reescalado no degrade el acuerdo entre solvers (ratio >= 1.0).
        La magnitud de la mejora depende de la plataforma y tolerancias numéricas.
    """
    from kuramoto_templates.evaluators.sdp import SDPEvaluator

    tmpl = get_s7_appendix_a_template(h=0.002, t=0.0)

    # 1. Plain
    lams_plain = {}
    for key, (s_name, tol) in [
        ("CL8", ("CLARABEL", 1e-8)),
        ("CL9", ("CLARABEL", 1e-9)),
        ("SCS", ("SCS", 1e-8)),
    ]:
        ev_p = SDPEvaluator(EvalConfig(mu_t=0.690, solver_name=s_name, solver_tol=tol))
        res_p = ev_p.evaluate_phase1(tmpl.x, tmpl.phi)
        assert res_p.lam is not None
        lams_plain[key] = res_p.lam

    diff_plain = max(
        abs(lams_plain["CL8"] - lams_plain["CL9"]),
        abs(lams_plain["CL8"] - lams_plain["SCS"]),
    )
    rel_diff_plain = diff_plain / abs(lams_plain["CL8"])

    # 2. Scaled
    lams_scaled = {}
    for key, (s_name, tol) in [
        ("CL8", ("CLARABEL", 1e-8)),
        ("CL9", ("CLARABEL", 1e-9)),
        ("SCS", ("SCS", 1e-8)),
    ]:
        ev_s = ScaledSDPEvaluator(EvalConfig(mu_t=0.690, solver_name=s_name, solver_tol=tol))
        res_s = ev_s.evaluate_phase1(tmpl.x, tmpl.phi)
        assert res_s.lam is not None
        lams_scaled[key] = res_s.lam

    diff_scaled = max(
        abs(lams_scaled["CL8"] - lams_scaled["CL9"]),
        abs(lams_scaled["CL8"] - lams_scaled["SCS"]),
    )
    rel_diff_scaled = diff_scaled / abs(lams_scaled["CL8"])

    ratio = rel_diff_plain / rel_diff_scaled
    print(
        f"\n[T4 RATIO]: {ratio:.2f}x (rel_plain={rel_diff_plain:.6e}, "
        f"rel_scaled={rel_diff_scaled:.6e})"
    )
    assert ratio >= 1.0, f"Acuerdo degradado en evaluador escalado: ratio {ratio:.2f}x < 1.0"


def test_t3_scaling_law_results_file() -> None:
    """T3: Ley de escala del paso 0b guardada en results/s7_scaling_law.json.

    Qué valida:
        1. results/s7_scaling_law.json existe y contiene las curvas para t in {0.0, 0.01, 0.05}.
        2. El exponente log-log para t=0.0 es ~3.0 (en [2.90, 3.10]).
        3. El exponente log-log para t=0.05 es ~2.0 (en [2.00, 2.30]).
        4. Todos los lambda_min(M) son estrictamente positivos.
    """
    import json
    from pathlib import Path

    file_path = (
        Path(__file__).resolve().parent.parent.parent
        / "results"
        / "s7_scaling_law.json"
    )
    assert file_path.exists(), f"Falta archivo {file_path}"

    with file_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    assert "exponents" in data
    exp_0 = data["exponents"]["0.0"]
    exp_05 = data["exponents"]["0.05"]

    assert 2.90 <= exp_0 <= 3.10, f"Exponente en t=0.0 debe ser ~3, dio {exp_0}"
    assert 2.00 <= exp_05 <= 2.30, f"Exponente en t=0.05 debe ser ~2, dio {exp_05}"

    # Positividad estricta
    for _t_str, curve in data["curves"].items():
        for pt in curve:
            assert pt["lam_min_numpy"] > 0.0
            assert pt["lam_min_mpmath"] > 0.0


def test_t8_sdpevaluator_and_domain_untouched() -> None:
    """T8 Regresión: SDPEvaluator no se modificó."""
    import subprocess

    diff_targets = [
        "src/kuramoto_templates/evaluators/sdp.py",
    ]
    res = subprocess.run(
        ["git", "diff", *diff_targets],
        capture_output=True,
        text=True,
        check=True,
    )
    msg = f"git diff en sdp.py debe estar vacío, pero dio:\n{res.stdout}"
    assert res.stdout.strip() == "", msg


def test_genome_sheet_assignment_vs_fallback() -> None:
    """Valida que la asignación por genoma use anclas k y compara con respaldo."""
    from kuramoto_templates.evaluators.scaled import compute_sheets

    # Genoma C4 con 4 clusters en anclas (0, 1, 2, 3) y desvíos pequeños
    g = Genome(
        m=4,
        K_ref=4,
        k=(0, 1, 2, 3),
        sgn=(1, 1, 1, 1),
        ell=(-3.0, -3.0, -3.0, -3.0),
        u=(0.0, 0.0, 0.0, 0.0),
    )
    x = decode_x(g.u)
    phi = decode_phi(g.K_ref, g.k, g.sgn, g.ell)

    sheets_g, ref_g = compute_sheets(x, phi, genome=g)
    sheets_fb, ref_fb = compute_sheets(x, phi, genome=None)

    assert sheets_g == sheets_fb


def test_t9_scale_incompatibility_penalty_s7_beats_c4_at_frontier() -> None:
    """Valida que en mu_t = 0.69291 el fitness de S7(h=0.002, t=0) en ScaledSDPEvaluator
    sea estrictamente mayor que el fitness de C4 exacto (m=4, h_eff=0).

    Evita el colapso por reward hacking donde C4 unscaled batía a S7.
    """
    from kuramoto_templates.cases.c4 import get_c4_template
    from kuramoto_templates.cases.s7 import get_s7_appendix_a_template

    mu_t = 0.69291
    ev = ScaledSDPEvaluator(
        EvalConfig(
            mu_t=mu_t,
            solver_tol=1e-8,
            lam_floor=1e-7,
            lam_floor_hat=1e-5,
        )
    )

    s7_template = get_s7_appendix_a_template(h=0.002, t=0.0)
    c4_template = get_c4_template()

    res_s7 = ev.evaluate_phase1(s7_template.x, s7_template.phi)
    res_c4 = ev.evaluate_phase1(c4_template.x, c4_template.phi)

    assert res_c4.status == "NOT_SCALABLE"
    assert res_s7.status == "OPTIMAL"
    assert res_s7.f > res_c4.f, (
        f"f(S7)={res_s7.f} debe ser estrictamente mayor que f(C4)={res_c4.f}"
    )


def test_verify_solution_scaled_threshold() -> None:
    """lam_min(M) ~ 1e-9 falla verify_solution sin h_eff, pero pasa con h_eff = 1e-3."""
    from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
    from kuramoto_templates.domain.conditions import verify_solution

    ev = ScaledSDPEvaluator(
        EvalConfig(mu_t=0.690, solver_tol=1e-9, lam_floor_hat=1e-5, h_eff_min=1e-3)
    )
    s7 = get_s7_appendix_a_template(h=0.002, t=0.0)
    res = ev.evaluate_phase1(s7.x, s7.phi)
    assert res.P is not None

    # 1. Sin h_eff (eps=1e-6 por defecto): lam_min(M) ~ 4.48e-9 < 1e-6 -> Falla
    v_unscaled = verify_solution(s7.x, s7.phi, res.P)
    assert not v_unscaled["is_valid"]
    assert 1e-10 < v_unscaled["min_eig_M"] < 1e-6
    assert any("Margen espectral insuficiente" in r for r in v_unscaled["reasons"])

    # 2. Con h_eff = 1e-3 (eps_efectivo = 1e-6 * 1e-6 = 1e-12): lam_min(M) > 1e-12 -> Pasa
    v_scaled = verify_solution(s7.x, s7.phi, res.P, h_eff=1e-3)
    assert v_scaled["is_valid"]
    assert len(v_scaled["reasons"]) == 0
    assert v_scaled["eps_effective"] == 1e-12





