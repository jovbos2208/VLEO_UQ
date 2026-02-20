#!/usr/bin/env python3
"""Run automated validation/report bundle for SV-3 and output artifacts."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _run(cmd: list[str]) -> tuple[int, str]:
    env = dict(os.environ)
    extra = f"{ROOT}:{ROOT / 'python'}"
    if env.get("PYTHONPATH"):
        env["PYTHONPATH"] = extra + ":" + env["PYTHONPATH"]
    else:
        env["PYTHONPATH"] = extra
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    return int(p.returncode), p.stdout


def main() -> None:
    parser = argparse.ArgumentParser(description="Run validation/report bundle.")
    parser.add_argument("--outdir", default="results/reports")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    steps = [
        {
            "name": "tolerance_sweep",
            "cmd": [
                args.python,
                "scripts/validate_tolerance_sweep.py",
                "--out",
                str(outdir / "tolerance_sweep_report.json"),
            ],
        },
        {
            "name": "coverage_report",
            "cmd": [
                args.python,
                "scripts/generate_coverage_report.py",
                "--out",
                str(outdir / "coverage_report.json"),
            ],
        },
        {
            "name": "mc_convergence_report",
            "cmd": [
                args.python,
                "scripts/generate_mc_convergence_report.py",
                "--out",
                str(outdir / "mc_convergence_report.json"),
            ],
        },
        {
            "name": "reproducibility_gate",
            "cmd": [
                args.python,
                "scripts/validate_reproducibility_gate.py",
                "--out",
                str(outdir / "reproducibility_gate_report.json"),
            ],
        },
    ]

    rows = []
    all_pass = True
    for step in steps:
        rc, out = _run(step["cmd"])
        rows.append({"name": step["name"], "returncode": rc, "pass": bool(rc == 0), "command": step["cmd"]})
        (outdir / f"{step['name']}.log").write_text(out, encoding="utf-8")
        all_pass = bool(all_pass and rc == 0)

    summary = {"pass": bool(all_pass), "steps": rows}
    (outdir / "validation_bundle_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[validation-bundle] wrote {outdir / 'validation_bundle_summary.json'}")
    raise SystemExit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
