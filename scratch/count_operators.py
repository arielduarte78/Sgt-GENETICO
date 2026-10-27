import json
from collections import Counter

ops_early_all = Counter()
ops_late_all = Counter()
ops_early_imp = Counter()
ops_late_imp = Counter()

# Track running global best and island best
global_best = -1e9
ops_global_imp_early = Counter()
ops_global_imp_late = Counter()
island_best = {}

with open("runs/night12h_seeded_mu06923_s42/run.jsonl", "r", encoding="utf-8") as f:
    for line in f:
        if '"type": "eval"' not in line:
            continue
        d = json.loads(line)
        op = d.get("operator", "unknown")
        gen = d.get("generation", 0)
        isl = d.get("island", 0)
        f_val = d.get("f", -1e9)

        if gen <= 50:
            ops_early_all[op] += 1
        else:
            ops_late_all[op] += 1

        # Global improvement
        if f_val > global_best and f_val > 0:
            global_best = f_val
            if gen <= 50:
                ops_global_imp_early[op] += 1
            else:
                ops_global_imp_late[op] += 1

        # Island improvement
        if isl not in island_best:
            island_best[isl] = f_val
            is_imp = True
        elif f_val > island_best[isl]:
            island_best[isl] = f_val
            is_imp = True
        else:
            is_imp = False

        if is_imp and f_val > 0:
            if gen <= 50:
                ops_early_imp[op] += 1
            else:
                ops_late_imp[op] += 1

print("=" * 80)
print("HISTOGRAMA DE OPERADORES: EVALUACIONES TOTALES")
print("=" * 80)
all_ops = sorted(set(list(ops_early_all.keys()) + list(ops_late_all.keys())))
print(f"{'Operador':<20} | {'gen <= 50':<12} | {'gen > 50':<12} | {'Total':<12}")
print("-" * 62)
for op in all_ops:
    c1 = ops_early_all[op]
    c2 = ops_late_all[op]
    print(f"{op:<20} | {c1:<12} | {c2:<12} | {c1+c2:<12}")

print("\n" + "=" * 80)
print("HISTOGRAMA DE OPERADORES: MEJORAS FACTIBLES EN ISLA (f > 0)")
print("=" * 80)
print(f"{'Operador':<20} | {'gen <= 50':<12} | {'gen > 50':<12} | {'Total':<12}")
print("-" * 62)
for op in all_ops:
    c1 = ops_early_imp[op]
    c2 = ops_late_imp[op]
    print(f"{op:<20} | {c1:<12} | {c2:<12} | {c1+c2:<12}")

print("\n" + "=" * 80)
print("HISTOGRAMA DE OPERADORES: MEJORAS DEL RECORD GLOBAL (f > 0)")
print("=" * 80)
print(f"{'Operador':<20} | {'gen <= 50':<12} | {'gen > 50':<12} | {'Total':<12}")
print("-" * 62)
for op in all_ops:
    c1 = ops_global_imp_early[op]
    c2 = ops_global_imp_late[op]
    print(f"{op:<20} | {c1:<12} | {c2:<12} | {c1+c2:<12}")
