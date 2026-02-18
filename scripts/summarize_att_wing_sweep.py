#!/usr/bin/env python3
"""Summarize ATT wing-angle sweep outputs into CSV and Markdown tables."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


def scalar_or_none(v: Any) -> float | None:
    try:
        out = float(v)
    except Exception:
        return None
    if out != out:  # NaN
        return None
    return out


def fmt(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:.3f}"
    if v is None:
        return "-"
    return str(v)


def parse_reached_fraction(note: Any) -> float | None:
    text = "" if note is None else str(note)
    m = re.search(r"\((\d+)\s*/\s*(\d+)\)", text)
    if not m:
        return None
    den = int(m.group(2))
    if den <= 0:
        return None
    return float(int(m.group(1)) / den)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def write_md(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("No rows.\n", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    lines = [
        "| " + " | ".join(keys) + " |",
        "| " + " | ".join(["---"] * len(keys)) + " |",
    ]
    for r in rows:
        lines.append("| " + " | ".join(fmt(r[k]) for k in keys) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize ATT sweep times to 90/180 deg.")
    parser.add_argument("--config", required=True, help="sweep config JSON")
    parser.add_argument("--results_dir", required=True, help="run output directory")
    parser.add_argument("--out_csv", required=True, help="output CSV file")
    parser.add_argument("--out_md", required=True, help="output markdown file")
    args = parser.parse_args()

    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    scenarios = cfg.get("scenarios", [])
    if not isinstance(scenarios, list):
        raise ValueError("invalid config: scenarios is not a list")

    results_dir = Path(args.results_dir)
    rows: list[dict[str, Any]] = []
    for s in scenarios:
        name = str(s["name"])
        axis = str(s.get("maneuver_axis", "")).lower()
        wing_deg = scalar_or_none(s.get("wing_constant_deg"))
        summary_path = results_dir / name / "summary.json"
        if not summary_path.exists():
            rows.append(
                {
                    "scenario": name,
                    "axis": axis,
                    "wing_deg": wing_deg,
                    "in_recommended_range": None,
                    "n_mc_effective": None,
                    "time_to_90_s": None,
                    "time_to_90_status": "missing_summary",
                    "time_to_180_s": None,
                    "time_to_180_status": "missing_summary",
                    "reached_fraction_180": None,
                }
            )
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        core = summary.get("catalog_core_outputs", {})
        m90 = core.get(f"time_to_{axis}_90_deg", {})
        m180 = core.get(f"time_to_{axis}_180_deg", {})
        frac180 = parse_reached_fraction(m180.get("note"))
        rows.append(
            {
                "scenario": name,
                "axis": axis,
                "wing_deg": wing_deg,
                "in_recommended_range": bool(wing_deg is not None and abs(wing_deg) <= 60.0),
                "n_mc_effective": summary.get("n_mc_effective"),
                "time_to_90_s": scalar_or_none(m90.get("value")),
                "time_to_90_status": m90.get("status"),
                "time_to_180_s": scalar_or_none(m180.get("value")),
                "time_to_180_status": m180.get("status"),
                "reached_fraction_180": frac180,
            }
        )

    rows.sort(key=lambda r: (str(r["axis"]), float(r["wing_deg"]) if r["wing_deg"] is not None else 1e9))
    out_csv = Path(args.out_csv)
    out_md = Path(args.out_md)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    write_csv(out_csv, rows)
    write_md(out_md, rows)
    print(f"[att-sweep] wrote {out_csv}")
    print(f"[att-sweep] wrote {out_md}")


if __name__ == "__main__":
    main()
