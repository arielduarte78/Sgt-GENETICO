import json

e_272 = None
e_274 = None

with open("runs/night12h_seeded_mu06923_s42/run.jsonl", "r", encoding="utf-8") as f:
    for line in f:
        if '"eval_id": 272008' in line:
            e_272 = json.loads(line)
        elif '"eval_id": 274309' in line:
            e_274 = json.loads(line)
        if e_272 and e_274:
            break

print("272008 (Parent):", e_272)
print("274309 (Child) :", e_274)

# Look in checkpoint for 274309
path = "runs/night12h_seeded_mu06923_s42/checkpoints/checkpoint_gen_1675.json"
with open(path, "r", encoding="utf-8") as f:
    data = json.load(f)

c274 = next((ind for ind in data.get("hall_of_fame_global", []) if ind["id"] == 274309), None)
c276 = next((ind for ind in data.get("hall_of_fame_global", []) if ind["id"] == 276921), None)

print("\n--- HoF 274309 ---")
print("Genome:", c274["genome"] if c274 else None)
print("x:", c274["phenotype_x"] if c274 else None)
print("phi:", c274["phenotype_phi"] if c274 else None)
print("f:", c274["f"] if c274 else None)

print("\n--- HoF 276921 (Sibling clone from 272008) ---")
print("Genome:", c276["genome"] if c276 else None)
print("x:", c276["phenotype_x"] if c276 else None)
print("phi:", c276["phenotype_phi"] if c276 else None)
print("f:", c276["f"] if c276 else None)

if c274 and c276:
    import numpy as np
    dx = np.array(c274["phenotype_x"]) - np.array(c276["phenotype_x"])
    dphi = np.array(c274["phenotype_phi"]) - np.array(c276["phenotype_phi"])
    print("Diff x (274309 vs 276921):", np.linalg.norm(dx))
    print("Diff phi (274309 vs 276921):", np.linalg.norm(dphi))
