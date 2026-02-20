#!/usr/bin/env python3
"""Generate consolidated maneuver realism report from gate artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def generate_report(results_dir: Path) -> dict[str, Any]:
    pre_path = results_dir / "maneuver_gate_pre.json"
    post_path = results_dir / "maneuver_gate_post.json"
    summary_path = results_dir / "maneuver_gate_summary.json"

    pre = _load_json(pre_path) if pre_path.exists() else {"results": [], "all_pass": None}
    post = _load_json(post_path) if post_path.exists() else {"results": [], "all_pass": None}
    summary = _load_json(summary_path) if summary_path.exists() else {}

    pre_rows = pre.get("results", []) if isinstance(pre.get("results"), list) else []
    post_rows = post.get("results", []) if isinstance(post.get("results"), list) else []
    post_fail = []
    for row in post_rows:
        if not isinstance(row, dict):
            continue
        if bool(row.get("pass", False)):
            continue
        post_fail.append(
            {
                "scenario": row.get("scenario"),
                "reasons": row.get("reasons", []),
            }
        )

    report = {
        "results_dir": str(results_dir),
        "files": {
            "pre": str(pre_path) if pre_path.exists() else None,
            "post": str(post_path) if post_path.exists() else None,
            "summary": str(summary_path) if summary_path.exists() else None,
        },
        "pre_count": int(len(pre_rows)),
        "post_count": int(len(post_rows)),
        "pre_all_pass": pre.get("all_pass"),
        "post_all_pass": post.get("all_pass"),
        "summary_pass": summary.get("pass"),
        "failed_post_scenarios": post_fail,
    }
    report["pass"] = bool(report.get("summary_pass", False))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate maneuver realism report.")
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    report = generate_report(results_dir)
    out_path = Path(args.out) if args.out else (results_dir / "maneuver_realism_report.json")
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"[maneuver-report] wrote {out_path}")
    raise SystemExit(0 if bool(report.get("pass", False)) else 1)


if __name__ == "__main__":
    main()
