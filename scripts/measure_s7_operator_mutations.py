import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kuramoto_templates.search.seeds import create_s7_direct_seed
from kuramoto_templates.search.operators import (
    anchor_jump,
    phase_dev,
    gauss_u,
    split,
    merge,
)
from kuramoto_templates.search.genome import decode_x, decode_phi
from kuramoto_templates.evaluators.scaled import compute_sheets, compute_h_eff

s7_seed = create_s7_direct_seed(h=0.002, t=0.0).genome

sheets_s7, _ = compute_sheets(decode_x(s7_seed.u), decode_phi(s7_seed.K_ref, s7_seed.k, s7_seed.sgn, s7_seed.ell), genome=s7_seed)
h_eff_s7 = compute_h_eff(decode_x(s7_seed.u), decode_phi(s7_seed.K_ref, s7_seed.k, s7_seed.sgn, s7_seed.ell), sheets_s7)

print(f"S7 inicial: m={s7_seed.m}, K_ref={s7_seed.K_ref}, h_eff={h_eff_s7:.6e}")

N_MUT = 1000
rng = np.random.default_rng(2026)

operators = {
    "phase_dev": lambda g, r: phase_dev(g, r),
    "gauss_u": lambda g, r: gauss_u(g, r),
    "anchor_jump": lambda g, r: anchor_jump(g, r),
    "split": lambda g, r: split(g, r),
    "merge": lambda g, r: merge(g, r),
}



print("=" * 70)
print(f"FRACCIÓN DE NOT_SCALABLE (h_eff < 1e-4) TRAS {N_MUT} MUTACIONES DESDE S7")
print("=" * 70)

for name, op in operators.items():
    not_scalable_count = 0
    h_effs = []
    for _ in range(N_MUT):
        child = op(s7_seed, rng)
        x = decode_x(child.u)
        phi = decode_phi(child.K_ref, child.k, child.sgn, child.ell)
        sheets, _ = compute_sheets(x, phi, genome=child)
        h_eff = compute_h_eff(x, phi, sheets)
        h_effs.append(h_eff)
        if h_eff < 1e-4:
            not_scalable_count += 1
    pct = (not_scalable_count / N_MUT) * 100.0
    print(f"Operador: {name:24s} | NOT_SCALABLE: {not_scalable_count:4d}/{N_MUT} ({pct:5.2f}%) | "
          f"h_eff min={min(h_effs):.3e}, mean={np.mean(h_effs):.3e}, max={np.max(h_effs):.3e}")
