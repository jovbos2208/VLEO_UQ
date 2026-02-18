#!/usr/bin/env python3
"""Build a reproducible conference config for A1/A2/C1 runs."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


DEFAULT_IDS = ",".join(
    [
        "ATT_A1_DETUMBLE_AERO_ONLY",
        "ATT_A2_NADIR_POINTING",
        "OD_C1_GNSS_CONTINUOUS_BESTCASE",
    ]
)


def run_catalog_converter(
    root: Path,
    catalog: Path,
    out: Path,
    ids: str,
    duration_scale: float,
    start_utc: str | None,
    no_env_sources: bool,
    dt_floor: float,
    default_particles: int,
) -> None:
    cmd = [
        sys.executable,
        str(root / "scripts" / "catalog_to_case_config.py"),
        "--catalog",
        str(catalog),
        "--out",
        str(out),
        "--ids",
        ids,
        "--duration_scale",
        str(duration_scale),
        "--dt_floor",
        str(dt_floor),
        "--default_particles",
        str(default_particles),
    ]
    if start_utc:
        cmd.extend(["--start_utc", start_utc])
    if no_env_sources:
        cmd.append("--no_env_sources")
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build A1/A2/C1 conference config.")
    parser.add_argument("--catalog", default="scenario_catalog.md")
    parser.add_argument("--out", default="configs/generated/paper_a1_a2_c1.json")
    parser.add_argument("--ids", default=DEFAULT_IDS)
    parser.add_argument("--duration_scale", type=float, default=1.0)
    parser.add_argument("--dt_floor", type=float, default=1.0)
    parser.add_argument("--default_particles", type=int, default=256)
    parser.add_argument("--particles_attitude", type=int, default=256)
    parser.add_argument("--particles_mission", type=int, default=256)
    parser.add_argument("--start_utc", default="2025-01-24T00:00:00")
    parser.add_argument("--no_env_sources", action="store_true")
    parser.add_argument("--pod_estimator", choices=("batch", "enkf"), default="enkf")
    parser.add_argument("--pod_enkf_members", type=int, default=128)
    parser.add_argument("--pod_enkf_inflation", type=float, default=1.02)
    parser.add_argument("--pod_enkf_use_carrier", action="store_true")
    parser.add_argument("--pod_arc_s", type=float, default=7200.0)
    parser.add_argument("--pod_overlap_s", type=float, default=0.0)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    catalog = Path(args.catalog)
    if not catalog.is_absolute():
        catalog = root / catalog
    out = Path(args.out)
    if not out.is_absolute():
        out = root / out
    out.parent.mkdir(parents=True, exist_ok=True)

    run_catalog_converter(
        root=root,
        catalog=catalog,
        out=out,
        ids=args.ids,
        duration_scale=args.duration_scale,
        start_utc=args.start_utc,
        no_env_sources=args.no_env_sources,
        dt_floor=args.dt_floor,
        default_particles=args.default_particles,
    )

    raw = json.loads(out.read_text(encoding="utf-8"))
    scenarios = raw.get("scenarios", [])
    if not scenarios:
        raise ValueError("no scenarios generated")

    for scenario in scenarios:
        name = str(scenario.get("name", ""))
        typ = str(scenario.get("type", "")).lower()
        if typ == "attitude":
            scenario["particles"] = int(args.particles_attitude)
            scenario["pod_uq"] = False
        elif typ == "mission":
            scenario["particles"] = int(args.particles_mission)
            scenario["pod_uq"] = True
            scenario["pod_estimator"] = str(args.pod_estimator)
            scenario["pod_enkf_members"] = int(args.pod_enkf_members)
            scenario["pod_enkf_inflation"] = float(args.pod_enkf_inflation)
            scenario["pod_enkf_use_carrier"] = bool(args.pod_enkf_use_carrier)
            scenario["pod_arc_s"] = float(args.pod_arc_s)
            scenario["pod_overlap_s"] = float(args.pod_overlap_s)
            if "od_c1" in name:
                scenario["pod_skip_slr"] = True
                scenario["pod_disable_gnss"] = False
        if args.start_utc:
            scenario["start_utc"] = args.start_utc
        scenario["use_env_sources"] = not bool(args.no_env_sources)

    out.write_text(json.dumps({"scenarios": scenarios}, indent=2), encoding="utf-8")
    print(f"[paper-config] wrote {len(scenarios)} scenarios -> {out}")


if __name__ == "__main__":
    main()
