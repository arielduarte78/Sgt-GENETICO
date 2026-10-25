"""Live monitoring script for the 12-hour overnight GA run.

Reads and displays runs/night12h_live_dashboard.txt in real time.
Usage:
    python scripts/watch_night.py          # Show current status once
    python scripts/watch_night.py --watch  # Refresh continuously every 5s
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import time


def main() -> None:
    parser = argparse.ArgumentParser(description="Watch overnight 12h run live dashboard")
    parser.add_argument("--watch", "-w", action="store_true", help="Refresh continuously every 5 seconds")
    parser.add_argument("--interval", "-i", type=float, default=5.0, help="Refresh interval in seconds (default 5)")
    args = parser.parse_args()

    dashboard_path = Path("runs/night12h_live_dashboard.txt")

    while True:
        if args.watch:
            # Clear console for live refresh
            os.system("cls" if os.name == "nt" else "clear")

        if dashboard_path.exists():
            try:
                content = dashboard_path.read_text(encoding="utf-8")
                print(content)
            except Exception as e:
                print(f"[WARN] Error leyendo dashboard: {e}")
        else:
            print("[INFO] Dashboard aun no disponible. Las corridas se estan inicializando...")

        if not args.watch:
            break

        try:
            time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\nMonitoreo finalizado.")
            break


if __name__ == "__main__":
    main()
