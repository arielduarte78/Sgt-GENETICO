"""Demonstrate real background run with SIGTERM/SIGINT interception and clean ga resume."""

import sys
sys.path.insert(0, "src")
import json
from pathlib import Path
import subprocess
import sys
import time

def main():
    import time
    ts = int(time.time())
    out_dir = Path(f"runs/demo_sigterm_{ts}")
    out_dir.mkdir(parents=True, exist_ok=True)
    from kuramoto_templates.search.config import GAConfig
    test_cfg = GAConfig(
        mu_t=0.690,
        num_islands=2,
        island_pop_size=6,
        num_generations=100,
        checkpoint_interval=1,
    )
    cfg_file = out_dir / "demo_cfg.yaml"
    test_cfg.to_yaml(cfg_file)

    import os
    env = dict(os.environ)
    env["PYTHONPATH"] = "src"

    cmd_run = [
        sys.executable,
        "-m",
        "kuramoto_templates.search.cli",
        "run",
        str(cfg_file),
        "-o",
        str(out_dir),
        "--allow-dirty",
    ]

    creationflags = 0
    if sys.platform == "win32":
        creationflags = subprocess.CREATE_NEW_PROCESS_GROUP

    print(">>> 1. Launching ga run in background...")
    proc = subprocess.Popen(
        cmd_run,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        creationflags=creationflags,
    )

    # Wait 3.5s so population is initialized and runs a couple of generations
    time.sleep(3.5)
    print(f"Poll before signal: {proc.poll()}")
    print(">>> 2. Sending SIGINT/SIGBREAK/termination signal to running process...")
    import signal
    if sys.platform == "win32":
        proc.send_signal(signal.CTRL_BREAK_EVENT)
    else:
        proc.send_signal(signal.SIGTERM)

    stdout, stderr = proc.communicate(timeout=10)
    print(f"Process returncode: {proc.returncode}")
    print("--- Background Process STDOUT ---")
    print(stdout)
    print("--- Background Process STDERR ---")
    print(stderr)

    status_path = out_dir / "status.json"
    if status_path.exists():
        with status_path.open("r", encoding="utf-8") as f:
            print("--- status.json after termination ---")
            print(f.read())

    # Find latest checkpoint
    cp_dir = out_dir / "checkpoints"
    checkpoints = sorted(cp_dir.glob("*.json"))
    print(f"Checkpoints available: {[cp.name for cp in checkpoints]}")
    if not checkpoints:
        print("No checkpoints found to resume.")
        return

    latest_cp = checkpoints[-1]
    print(f"\n>>> 3. Executing ga resume on {latest_cp}...")
    cmd_resume = [
        sys.executable,
        "-m",
        "kuramoto_templates.search.cli",
        "resume",
        str(latest_cp),
        "-g",
        "5",
    ]
    res_proc = subprocess.run(cmd_resume, capture_output=True, text=True, env=env)
    print("--- Resume STDOUT ---")
    print(res_proc.stdout)
    print("--- Resume STDERR ---")
    print(res_proc.stderr)

if __name__ == "__main__":
    main()
