import json

path = "runs/night12h_seeded_mu06923_s42/checkpoints/checkpoint_gen_1675.json"
with open(path, "r", encoding="utf-8") as f:
    data = json.load(f)

inds = []
seen = set()
for island in data.get("island_populations", []):
    for ind in island:
        iid = ind["id"]
        if iid not in seen:
            seen.add(iid)
            inds.append(ind)

# Also check hall of fame
for ind in data.get("hall_of_fame_global", []):
    iid = ind["id"]
    if iid not in seen:
        seen.add(iid)
        inds.append(ind)

inds.sort(key=lambda x: x.get("lam_normalized", x.get("f", -1e9)), reverse=True)

print(f"Total individuos unicos analizados: {len(inds)}")
print(f"{'Rank':<5} | {'ID':<8} | {'lam_hat (Fase I)':<20} | {'s_star':<14} | {'mu_star':<16} | {'Gen':<6} | {'Operator'}")
print("-" * 95)
for r, ind in enumerate(inds[:20], 1):
    lam = ind.get("lam_normalized", ind.get("f"))
    s = ind.get("s_star")
    mu = ind.get("mu_star")
    mu_str = f"{mu:.8f}" if mu is not None else "None"
    print(f"{r:<5} | {ind['id']:<8} | {lam:<20.10e} | {s:<14.4e} | {mu_str:<16} | {str(ind.get('generation', '?')):<6} | {ind.get('operator_name', '?')}")
