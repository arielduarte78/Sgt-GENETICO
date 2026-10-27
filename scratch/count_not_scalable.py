import json
from collections import Counter

target_ops = {"anchor_jump", "merge", "split", "gauss_u", "phase_dev", "crossover"}
counts_total = Counter()
counts_not_scalable = Counter()

with open("runs/night12h_seeded_mu06923_s42/run.jsonl", "r", encoding="utf-8") as f:
    for line in f:
        if '"type": "eval"' not in line:
            continue
        d = json.loads(line)
        op = d.get("operator")
        if op in target_ops:
            counts_total[op] += 1
            if d.get("solver_status") == "NOT_SCALABLE":
                counts_not_scalable[op] += 1

print("=" * 65)
print("RECUENTO DE ESTADO NOT_SCALABLE PARA gauss_u Y phase_dev")
print("=" * 65)
for op in ["gauss_u", "phase_dev"]:
    tot = counts_total[op]
    ns = counts_not_scalable[op]
    print(f"  {op:<15}: {ns:>6} / {tot:>6} ({ns/tot*100:6.2f}%)")

tot_comb = counts_total["gauss_u"] + counts_total["phase_dev"]
ns_comb = counts_not_scalable["gauss_u"] + counts_not_scalable["phase_dev"]
print(f"  {'TOTAL COMBINADO':<15}: {ns_comb:>6} / {tot_comb:>6} ({ns_comb/tot_comb*100:6.2f}%)")
