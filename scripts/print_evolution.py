"""Print generation milestones and evolution from prod run."""

import json

with open("results/prod_run_0.6929077.json") as f:
    d = json.load(f)

print(f"Total generations: {len(d['generation_evolution'])}")
prev_f = None
for x in d['generation_evolution']:
    g = x['generation']
    f_val = x['best_f']
    m = x['best_m']
    lam = x['best_lam']
    if prev_f is None or f_val != prev_f or g in [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200]:
        print(f"Gen {g:3d} | best_f = {f_val:14.6e} | best_m = {m:2d} | best_lam = {lam:14.6e}")
        prev_f = f_val
