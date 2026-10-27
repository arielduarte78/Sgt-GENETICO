import json

target_id = 125678
log_path = "runs/night12h_seeded_mu06923_s42/run.jsonl"

with open(log_path, "r", encoding="utf-8") as f:
    for line in f:
        if f'"id": {target_id}' in line:
            d = json.loads(line)
            print("Encontrado R1:")
            print(f"  id: {d.get('id')}")
            print(f"  eval_id: {d.get('eval_id')}")
            print(f"  generation: {d.get('generation')}")
            print(f"  parent_ids: {d.get('parent_ids')}")
            print(f"  operator: {d.get('operator')}")
            print(f"  mu_eval: {d.get('mu_eval')}")
            break
