"""Trace R1 ancestry and compute progression curve vs evaluations."""

import json
from pathlib import Path

run_log = Path("runs/night12h_seeded_mu06923_s42/run.jsonl")

eval_map = {}
thresholds = [0.6915, 0.6920, 0.6921, 0.6923, 0.6925, 0.6928]
first_hit = {th: None for th in thresholds}

# Track running best
running_best_f = -1e9
running_best_eval = 0

print("Analizando run.jsonl (165 MB)...")

with open(run_log, "r", encoding="utf-8") as f:
    for line_idx, line in enumerate(f):
        if not line.strip():
            continue
        try:
            d = json.loads(line)
        except Exception:
            continue
        if d.get("type") != "eval":
            continue

        eid = d.get("eval_id")
        f_val = d.get("f", -1e9)
        mu_val = d.get("mu_eval") or d.get("mu_star")
        gen = d.get("generation")
        op = d.get("operator")
        parents = d.get("parent_ids", [])

        # Store minimal metadata for lineage tracing
        eval_map[eid] = {
            "eval_id": eid,
            "gen": gen,
            "parents": parents,
            "op": op,
            "f": f_val,
            "lam": d.get("lam_norm") or d.get("lam_raw"),
            "s": d.get("s_star"),
        }

        # Check if mu crossed thresholds
        if mu_val is not None:
            for th in thresholds:
                if first_hit[th] is None and mu_val >= th:
                    first_hit[th] = (eid, gen, mu_val)

        if f_val > running_best_f:
            running_best_f = f_val
            running_best_eval = eid

print(f"Total evaluaciones procesadas: {len(eval_map)}")
print(f"Mejor fitness alcanzado: {running_best_f:.6e} en eval_id {running_best_eval}")

# ==============================================================================
# 1. ANCESTRY TRACE OF R1 (eval_id 125678)
# ==============================================================================
print("\n" + "=" * 90)
print("1. CADENA COMPLETA DE ANCESTROS DE R1 (ID 125678)")
print("=" * 90)

target = 125678
chain = []
curr = target
visited = set()

while curr and curr in eval_map and curr not in visited:
    visited.add(curr)
    node = eval_map[curr]
    chain.append(node)
    parents = node["parents"]
    if parents:
        curr = parents[0]
    else:
        break

chain.reverse()

print(f"{'Gen':<6} | {'Eval ID':<10} | {'Parent IDs':<18} | {'Operador':<20} | {'Fitness f':<16} | {'lam_norm':<14}")
print("-" * 90)
for node in chain:
    p_str = str(node["parents"])
    print(f"{str(node['gen']):<6} | {node['eval_id']:<10} | {p_str:<18} | {str(node['op']):<20} | {node['f']:<16.6e} | {str(node['lam']):<14}")


# ==============================================================================
# 2. EVALUATION PROGRESSION FOR THRESHOLDS
# ==============================================================================
print("\n" + "=" * 90)
print("2. HITOS DE EVALUACION EN R1")
print("=" * 90)
print(f"{'Umbral mu*':<14} | {'Eval ID':<10} | {'Generacion':<12} | {'mu* Real':<14}")
print("-" * 55)
for th in thresholds:
    res = first_hit[th]
    if res is not None:
        print(f"{th:<14.4f} | {res[0]:<10} | {res[1]:<12} | {res[2]:<14.6f}")
    else:
        print(f"{th:<14.4f} | {'NUNCA':<10} | {'NUNCA':<12} | {'N/A':<14}")

# Save lineage and thresholds
out_data = {
    "chain": chain,
    "first_hit": first_hit,
    "total_evals": len(eval_map),
    "running_best_eval": running_best_eval,
    "running_best_f": running_best_f,
}
with open("results/r1_ancestry_and_hits.json", "w", encoding="utf-8") as f:
    json.dump(out_data, f, indent=2)

print("\nResultados guardados en results/r1_ancestry_and_hits.json")
