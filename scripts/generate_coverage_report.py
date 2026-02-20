#!/usr/bin/env python3
"""Generate coverage calibration report (MC vs UT vs STM)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from scripts.validate_coverage_short_arc import run_validation


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate short-arc coverage report.")
    parser.add_argument("--mc_particles", type=int, default=192)
    parser.add_argument("--duration_s", type=float, default=600.0)
    parser.add_argument("--dt_s", type=float, default=30.0)
    parser.add_argument("--seed", type=int, default=77)
    parser.add_argument("--out", default="results/reports/coverage_report.json")
    args = parser.parse_args()

    report = run_validation(
        mc_particles=int(args.mc_particles),
        duration_s=float(args.duration_s),
        dt_s=float(args.dt_s),
        seed=int(args.seed),
    )
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[coverage-report] wrote {out_path}")
    raise SystemExit(0 if bool(report.get("pass", False)) else 1)


if __name__ == "__main__":
    main()
