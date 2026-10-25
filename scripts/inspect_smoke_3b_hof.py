import json

ckpt = json.load(open("runs/smoke_frontier_scaled_rerun/checkpoints/checkpoint_gen_50.json"))

hof_global = ckpt.get("hall_of_fame_global", [])
hof_m = ckpt.get("hall_of_fame_per_m", {})
lam_floor_hat = 1e-5

print("==================================================")
print("AUDITORÍA HALL OF FAME: SMOKE TEST 3b (mu_t = 0.69291)")
print("==================================================")
print(f"Piso espectral de referencia: lam_floor_hat = {lam_floor_hat:.1e}\n")

print(f"--- 1. Hall of Fame Global (Total: {len(hof_global)}) ---")
ge_global = 0
for i, ind in enumerate(hof_global):
    ind_id = ind.get("id")
    m = ind.get("genome", {}).get("m")
    f_val = ind.get("f")
    s_star = ind.get("s_star")
    lam_star = ind.get("lam_star")
    lam_norm = ind.get("lam_normalized")
    is_ge = (lam_norm is not None and lam_norm >= lam_floor_hat) or (lam_star is not None and lam_star >= lam_floor_hat)
    if is_ge:
        ge_global += 1
    print(f"Rank {i+1:2d} | ID {ind_id:5d} | m={m} | f={f_val:+.8e} | s*={s_star:.6e} | lam_hat={lam_star:.6e} | >= piso: {is_ge}")

print(f"\nTotal en HoF Global con lam_hat >= lam_floor_hat: {ge_global}")

print(f"\n--- 2. Hall of Fame por m (Total: {len(hof_m)}) ---")
ge_m = 0
for m_val in sorted(hof_m.keys(), key=lambda x: int(x)):
    ind = hof_m[m_val]
    ind_id = ind.get("id")
    m = ind.get("genome", {}).get("m")
    f_val = ind.get("f")
    s_star = ind.get("s_star")
    lam_star = ind.get("lam_star")
    lam_norm = ind.get("lam_normalized")
    is_ge = (lam_norm is not None and lam_norm >= lam_floor_hat) or (lam_star is not None and lam_star >= lam_floor_hat)
    if is_ge:
        ge_m += 1
    print(f"m={m_val:2s} | ID {ind_id:5d} | m={m} | f={f_val:+.8e} | s*={s_star:.6e} | lam_hat={lam_star:.6e} | >= piso: {is_ge}")

print(f"\nTotal en HoF por m con lam_hat >= lam_floor_hat: {ge_m}")
