import sys
sys.path.insert(0, 'src')
import numpy as np
from kuramoto_templates.evaluators.sdp import SDPEvaluator
from kuramoto_templates.evaluators.config import EvalConfig
from kuramoto_templates.domain.conditions import verify_solution

# Los datos exactos del ID 39260 (del prompt del usuario)
x = (0.249443, 0.250188, 0.250534, 0.249835)
phi = (-3.141592, -1.570795, 1.0e-6, 1.570797)
mu_t = 0.6929077

print('--- Test Antidoping ID 39260 (m=4) [Valores del prompt] ---')
ev_cl8 = SDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, solver_name='CLARABEL'))
r_cl8 = ev_cl8.evaluate_phase1(x, phi)
print(f'Clarabel (1e-8): lam*={r_cl8.lam:.12e}, status={r_cl8.status}')

ev_cl9 = SDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-9, solver_name='CLARABEL'))
r_cl9 = ev_cl9.evaluate_phase1(x, phi)
print(f'Clarabel (1e-9): lam*={r_cl9.lam:.12e}, status={r_cl9.status}')

ev_scs = SDPEvaluator(EvalConfig(mu_t=mu_t, solver_tol=1e-8, solver_name='SCS'))
r_scs = ev_scs.evaluate_phase1(x, phi)
print(f'SCS (1e-8):      lam*={r_scs.lam:.12e}, status={r_scs.status}')

if r_cl8.P is not None:
    v = verify_solution(x, phi, r_cl8.P, eps=1e-8, eps_c=1e-4, theta_tol=1e-4)
    print('\n--- Verificacion Fisica (NumPy Puro) ---')
    print(f'max |F| (Equilibrio): {v["max_F"]:.12e}')
    print(f'min eigval M (Tension): {v["min_eig_M"]:.12e}')
    print(f'is_valid: {v["is_valid"]}')

# Adicional: Datos con precision completa float64 del checkpoint
import json
data = json.load(open('runs/prod_mu_0.6929077_8isl_mu64_gen200/checkpoints/checkpoint_gen_200.json'))
ind = data['hall_of_fame_global'][0]
x_full = tuple(ind['phenotype_x'])
phi_full = tuple(ind['phenotype_phi'])

print('\n--- Test Antidoping ID 39260 (m=4) [Precision completa float64 del checkpoint] ---')
r_cl8_full = ev_cl8.evaluate_phase1(x_full, phi_full)
print(f'Clarabel (1e-8): lam*={r_cl8_full.lam:.12e}, status={r_cl8_full.status}')
r_cl9_full = ev_cl9.evaluate_phase1(x_full, phi_full)
print(f'Clarabel (1e-9): lam*={r_cl9_full.lam:.12e}, status={r_cl9_full.status}')
r_scs_full = ev_scs.evaluate_phase1(x_full, phi_full)
print(f'SCS (1e-8):      lam*={r_scs_full.lam:.12e}, status={r_scs_full.status}')

if r_cl8_full.P is not None:
    v_full = verify_solution(x_full, phi_full, r_cl8_full.P, eps=1e-8, eps_c=1e-4, theta_tol=1e-4)
    print('\n--- Verificacion Fisica (NumPy Puro - float64) ---')
    print(f'max |F| (Equilibrio): {v_full["max_F"]:.12e}')
    print(f'min eigval M (Tension): {v_full["min_eig_M"]:.12e}')
    print(f'is_valid: {v_full["is_valid"]}')
