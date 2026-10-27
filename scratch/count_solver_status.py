import json
from collections import Counter

status_counts = Counter()
top_evals = []

with open("runs/night12h_seeded_mu06923_s42/run.jsonl", "r", encoding="utf-8") as f:
    for line in f:
        if not line.strip():
            continue
        try:
            d = json.loads(line)
        except Exception:
            continue
        if d.get("type") == "eval":
            st = d.get("solver_status", "UNKNOWN")
            status_counts[st] += 1
            lam = d.get("lam_norm") or d.get("lam_raw") or 0.0
            top_evals.append((lam, d.get("eval_id"), st, d.get("f"), d.get("generation")))

print(f"HISTOGRAMA DE SOLVER_STATUS EN run.jsonl:")
total = sum(status_counts.values())
for st, cnt in status_counts.most_common():
    print(f"  {st:<25}: {cnt:>8} ({cnt/total*100:6.2f}%)")
print(f"  TOTAL                    : {total:>8}")

# Top feasible (f > 0, which means s_star <= 1e-6 and f = lam_norm)
feasible = [e for e in top_evals if e[3] is not None and e[3] > 0]
feasible.sort(key=lambda x: x[3], reverse=True)

st_feas = Counter(item[2] for item in feasible)
print(f"\nDESGLOSE DE STATUS EN EVALUACIONES FACTIBLES (f > 0):")
for st, cnt in st_feas.most_common():
    print(f"  {st:<25}: {cnt:>8} ({cnt/len(feasible)*100:6.2f}%)")

print(f"\nTOP 20 POR FITNESS FACTIBLE f (lam_hat):")
print(f"  {'Rank':<5} | {'Eval ID':<10} | {'Gen':<6} | {'f (lam_hat)':<18} | {'Status'}")
print("  " + "-" * 60)
for r, item in enumerate(feasible[:20], 1):
    print(f"  {r:<5} | {item[1]:<10} | {item[4]:<6} | {item[3]:<18.10e} | {item[2]}")
