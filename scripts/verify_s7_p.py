"""Verification of solution P for S7 at mu_t=0.690 with Clarabel and SCS."""

import sys
sys.path.insert(0, "src")

from kuramoto_templates.cases.s7 import get_s7_appendix_a_template
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.domain.conditions import verify_solution

def main():
    s7 = get_s7_appendix_a_template(h=0.002, t=0.0)
    x, phi = s7.x, s7.phi

    # 1. Clarabel tol=1e-8
    eval_cl_8 = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-8, solver_name="CLARABEL"))
    res_cl_8 = eval_cl_8.evaluate_phase1(x, phi)

    # 2. Clarabel tol=1e-9
    eval_cl_9 = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-9, solver_name="CLARABEL"))
    res_cl_9 = eval_cl_9.evaluate_phase1(x, phi)

    # 3. SCS tol=1e-8
    eval_scs = SDPEvaluator(EvalConfig(mu_t=0.690, solver_tol=1e-8, solver_name="SCS"))
    res_scs = eval_scs.evaluate_phase1(x, phi)

    print("=== SDP RESULTS FOR S7 (mu_t = 0.690) ===")
    print(f"1) Clarabel (tol=1e-8): lam* = {res_cl_8.lam:.12e}, s* = {res_cl_8.s:.12e}, status = {res_cl_8.status}")
    print(f"2) Clarabel (tol=1e-9): lam* = {res_cl_9.lam:.12e}, s* = {res_cl_9.s:.12e}, status = {res_cl_9.status}")
    print(f"3) SCS      (tol=1e-8): lam* = {res_scs.lam:.12e}, s* = {res_scs.s:.12e}, status = {res_scs.status}")

    # Verify solution on res_cl_8.P
    if res_cl_8.P is not None:
        v = verify_solution(x, phi, res_cl_8.P, eps=1e-8, eps_c=1e-4, theta_tol=1e-4)
        print("\n=== PURE NUMPY VERIFICATION (res_cl_8.P) ===")
        print(f"max |F|:             {v['max_F']:.12e}")
        print(f"min c_a:             {v['min_c']:.12e}")
        print(f"min eigval of M:     {v['min_eig_M']:.12e}")
        print(f"mu (min d_a):        {v['mu']:.12e}")
        print(f"p_bounds_ok:         {v['p_bounds_ok']}")
        print(f"is_valid:            {v['is_valid']}")
        print(f"mu_t:                0.690")
        print(f"deficit (0.690 - mu):{0.690 - v['mu']:.12e}")

    # Also verify solution on res_cl_9.P and res_scs.P
    if res_cl_9.P is not None:
        v9 = verify_solution(x, phi, res_cl_9.P, eps=1e-9, eps_c=1e-4, theta_tol=1e-4)
        print("\n=== PURE NUMPY VERIFICATION (res_cl_9.P) ===")
        print(f"min eigval of M:     {v9['min_eig_M']:.12e}")
        print(f"mu (min d_a):        {v9['mu']:.12e}")

    if res_scs.P is not None:
        v_scs = verify_solution(x, phi, res_scs.P, eps=1e-8, eps_c=1e-4, theta_tol=1e-4)
        print("\n=== PURE NUMPY VERIFICATION (res_scs.P) ===")
        print(f"min eigval of M:     {v_scs['min_eig_M']:.12e}")
        print(f"mu (min d_a):        {v_scs['mu']:.12e}")

if __name__ == "__main__":
    main()
