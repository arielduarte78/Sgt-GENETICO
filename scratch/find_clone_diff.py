import json

# Checkpoint has hall_of_fame_global which contains 274309!
path = "runs/night12h_seeded_mu06923_s42/checkpoints/checkpoint_gen_1675.json"
with open(path, "r", encoding="utf-8") as f:
    data = json.load(f)

hof = data.get("hall_of_fame_global", [])
c274 = None
for ind in hof:
    if ind["id"] == 274309:
        c274 = ind
        break

print("274309 in HoF:", c274 is not None)
if c274:
    print("274309 parent_ids:", c274.get("parent_ids"))
    print("274309 operator:", c274.get("operator_name"))
    print("274309 f:", c274.get("f"))
    print("274309 genome:", c274.get("genome"))

# Now let's find the parent in run.jsonl or earlier checkpoints
# If parent is in run.jsonl, find it
pid = c274["parent_ids"][0] if c274 and c274.get("parent_ids") else None
print(f"Parent ID: {pid}")

p_ind = None
# Let's search checkpoints for pid first
import glob
for cp_file in sorted(glob.glob("runs/night12h_seeded_mu06923_s42/checkpoints/*.json")):
    with open(cp_file, "r", encoding="utf-8") as f:
        cp_d = json.load(f)
    for ind in cp_d.get("hall_of_fame_global", []):
        if ind["id"] == pid:
            p_ind = ind
            break
    if p_ind:
        break

if not p_ind and pid:
    # Search in run.jsonl
    with open("runs/night12h_seeded_mu06923_s42/run.jsonl", "r", encoding="utf-8") as f:
        for line in f:
            if f'"eval_id": {pid}' in line or f'"eval_id":{pid}' in line:
                d = json.loads(line)
                if d.get("eval_id") == pid:
                    p_ind = d
                    break

print("Parent found:", p_ind is not None)
if p_ind:
    print("Parent keys:", list(p_ind.keys()))
    print("Parent f:", p_ind.get("f"))
    print("Parent genome:", p_ind.get("genome"))
