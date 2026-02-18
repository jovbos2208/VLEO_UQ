#!/usr/bin/env python3
"""Recompute convergence metrics from existing det_stm outputs with a new reference dt."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.convergence_study import compare_type, summarize_rows, write_csv


def parse_types(raw: str) -> list[str]:
    vals = []
    for tok in raw.split(","):
        t = tok.strip().lower()
        if not t:
            continue
        vals.append(t)
    if not vals:
        raise ValueError("types is empty")
    valid = {"attitude", "mission_od"}
    bad = [t for t in vals if t not in valid]
    if bad:
        raise ValueError(f"invalid types: {bad}")
    return vals


def dt_from_config_name(path: Path) -> float:
    m = re.match(r"config_dt_([0-9]+p?[0-9]*)\.json$", path.name)
    if not m:
        raise ValueError(f"unexpected config filename: {path.name}")
    return float(m.group(1).replace("p", "."))


def gather_configs(type_dir: Path) -> dict[float, dict]:
    cfgs: dict[float, dict] = {}
    for cfg_path in sorted(type_dir.glob("config_dt_*.json")):
        dt = dt_from_config_name(cfg_path)
        cfgs[dt] = json.loads(cfg_path.read_text(encoding="utf-8"))
    if not cfgs:
        raise FileNotFoundError(f"no config_dt_*.json found in {type_dir}")
    return cfgs


def gather_run_dirs(type_dir: Path, dt_values: list[float]) -> dict[float, Path]:
    out = {}
    for dt in dt_values:
        tag = f"dt_{str(dt).replace('.', 'p')}"
        d = type_dir / tag
        if not d.is_dir():
            raise FileNotFoundError(f"missing run directory for dt={dt}: {d}")
        out[dt] = d
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Recompute convergence metrics from existing outputs.")
    parser.add_argument("--outdir", required=True, help="Path to results/convergence_*")
    parser.add_argument("--types", default="attitude,mission_od")
    parser.add_argument("--reference_dt", type=float, default=1.0)
    parser.add_argument("--mu", type=float, default=3.986004418e14)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if not outdir.is_dir():
        raise FileNotFoundError(outdir)

    types = parse_types(args.types)
    all_rows = []

    for type_name in types:
        type_dir = outdir / type_name
        if not type_dir.is_dir():
            raise FileNotFoundError(type_dir)

        config_by_dt = gather_configs(type_dir)
        dt_values = sorted(config_by_dt.keys())
        ref_dt = float(args.reference_dt)
        if ref_dt not in config_by_dt:
            raise ValueError(f"reference_dt={ref_dt} not found for {type_name}; available={dt_values}")

        out_by_dt = gather_run_dirs(type_dir, dt_values)
        rows = compare_type(
            type_name=type_name,
            config_by_dt=config_by_dt,
            out_by_dt=out_by_dt,
            ref_dt=ref_dt,
            mu=float(args.mu),
        )
        summary = summarize_rows(rows)
        (type_dir / "metrics.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
        (type_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        write_csv(type_dir / "metrics.csv", rows)
        write_csv(type_dir / "summary.csv", summary)
        print(f"[recompute] updated {type_name}: reference_dt_s={ref_dt} rows={len(rows)}")
        all_rows.extend(rows)

    (outdir / "metrics_all.json").write_text(json.dumps(all_rows, indent=2), encoding="utf-8")
    write_csv(outdir / "metrics_all.csv", all_rows)
    print(f"[recompute] updated aggregate files in {outdir}")


if __name__ == "__main__":
    main()
