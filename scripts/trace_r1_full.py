"""Trace full ancestry chain of R1 (ID 125678) back to the initial seed."""

import json
from pathlib import Path

run_log = Path("runs/night12h_seeded_mu06923_s42/run.jsonl")

# We want to find ancestry chain starting from 125678
# First, scan run.jsonl to build an index of the needed ancestors.
# Since run.jsonl is large (165MB), we can do it efficiently.

needed_ids = {125678}
chain = {}

print("Buscando ancestros en run.jsonl...")

# Multiple passes backwards if needed, or index on the fly
# Since IDs are sequential (parent_id < child_id), reading backwards or reading all ancestors:
with open(run_log, "r", encoding="utf-8") as f:
    for line in f:
        # Quick check if any needed_id in line
        for nid in list(needed_ids):
            if f'"id": {nid}' in line or f'"id":{nid}' in line:
                d = json.loads(line)
                if d.get("id") == nid:
                    chain[nid] = {
                        "id": d.get("id"),
                        "eval_id": d.get("eval_id"),
                        "gen": d.get("generation"),
                        "parent_ids": d.get("parent_ids", []),
                        "op": d.get("operator_name") or d.get("operator"),
                        "f": d.get("f"),
                        "lam": d.get("lam_normalized") or d.get("lam_star"),
                        "mu_eval": d.get("mu_eval") or d.get("mu_star"),
                    }
                    for p in d.get("parent_ids", []):
                        if p not in chain and p > 0:
                            needed_ids.add(p)

print(f"Ancestros indexados: {len(chain)}")

# Now trace from 125678 down to seeds
curr = 125678
lineage = []
visited = set()
queue = [125678]

while queue:
    cid = queue.pop(0)
    if cid in visited or cid not in chain:
        continue
    visited.add(cid)
    node = chain[cid]
    lineage.append(node)
    for p in node["parent_ids"]:
        if p not in visited and p in chain:
            queue.append(p)

# Sort chronologically by gen/id
lineage.sort(key=lambda x: (x["gen"] if x["gen"] is not None else 0, x["id"]))

print("\n" + "=" * 90)
print(f"{'Gen':<6} | {'ID':<8} | {'Parent IDs':<18} | {'Operador':<14} | {'f (fitness)':<14} | {'lam_hat':<14}")
print("-" * 90)
for node in lineage:
    p_str = str(node["parent_ids"])
    print(f"{str(node['gen']):<6} | {node['id']:<8} | {p_str:<18} | {str(node['op']):<14} | {node['f']:<14.6e} | {str(node['lam']):<14}")

with open("results/r1_lineage.json", "w", encoding="utf-8") as f:
    json.dump(lineage, f, indent=2)

print("\nLinaje guardado en results/r1_lineage.json")
