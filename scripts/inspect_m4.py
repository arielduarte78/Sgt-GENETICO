"""Inspect m=4 individuals and compare against C4."""

import json
from pathlib import Path
import numpy as np

from kuramoto_templates.cases.c4 import get_c4_template


def wrap_angle(theta):
    """Wrap angle to (-pi, pi]."""
    return (theta + np.pi) % (2 * np.pi) - np.pi


def main():
    c4 = get_c4_template()
    c4_x = np.array(c4.x)
    c4_phi = np.array(c4.phi)

    print(f"Theoretical C4:")
    print(f"  x:   {c4_x}")
    print(f"  phi: {c4_phi}")
    print()

    # Load results from mandatory_experiments.json
    with open("results/mandatory_experiments.json") as f:
        data = json.load(f)

    # Also search in the run jsonl logs in runs/exp_mu_*
    for exp_dir in sorted(Path("runs").glob("exp_mu_*")):
        log_file = exp_dir / "run.jsonl"
        manifest_file = exp_dir / "manifest.json"
        if not log_file.exists():
            continue

        print(f"=== Inspecting run directory: {exp_dir.name} ===")
        # Look for m=4 individuals
        with open(log_file, "r", encoding="utf-8") as f:
            for line in f:
                record = json.loads(line)
                if record.get("event") == "individual_evaluated":
                    ind_data = record["individual"]
                    if ind_data["genome"]["m"] == 4:
                        ind_id = ind_data["id"]
                        x = np.array(ind_data["phenotype_x"])
                        phi = np.array(ind_data["phenotype_phi"])
                        lam = ind_data["lam_star"]
                        s = ind_data["s_star"]
                        f_val = ind_data["f"]
                        op = ind_data["operator_name"]

                        # Check distance to C4
                        # Note: Gauge invariance allows global phase shift and permutation of clusters
                        # For canonicalized genomes, cluster sizes and phases are ordered
                        # Let's compute min distance under cyclic shifts and permutations
                        dist_x = np.max(np.abs(x - c4_x))
                        
                        # Phase distance modulo 2pi and gauge shift
                        # Try all gauge shifts where phi[0] aligns or best circular alignment
                        best_dphi = float("inf")
                        for shift in np.linspace(-np.pi, np.pi, 360):
                            dphi = np.max(np.abs(wrap_angle(phi - c4_phi - shift)))
                            if dphi < best_dphi:
                                best_dphi = dphi

                        # Also permutation check if needed
                        import itertools
                        perm_best_dx = float("inf")
                        perm_best_dphi = float("inf")
                        for p in itertools.permutations(range(4)):
                            xp = x[list(p)]
                            phip = phi[list(p)]
                            dx = np.max(np.abs(xp - c4_x))
                            # alignment
                            shift = phip[0] - c4_phi[0]
                            dph = np.max(np.abs(wrap_angle(phip - c4_phi - shift)))
                            if dx + dph < perm_best_dx + perm_best_dphi:
                                perm_best_dx = dx
                                perm_best_dphi = dph

                        print(f"  ID: {ind_id:3d} | op: {op:12s} | lam*: {lam:14.6e} | s*: {s:10.2e} | f: {f_val:14.6e}")
                        print(f"    x:   {np.round(x, 6)}")
                        print(f"    phi: {np.round(phi, 6)}")
                        print(f"    dist(x, C4_x) = {perm_best_dx:.6e}, dist(phi, C4_phi) = {perm_best_dphi:.6e}")


if __name__ == "__main__":
    main()
