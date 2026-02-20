#!/usr/bin/env python3
"""Run the Plan_progression Phase A/B/C scenario campaign end-to-end."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]

PHASE_SCENARIO_IDS: dict[str, list[str]] = {
    "phase_a": [
        "ATT_A1_DETUMBLE_AERO_ONLY",
        "ATT_A2_NADIR_POINTING",
        "OD_C1_GNSS_CONTINUOUS_BESTCASE",
        "OD_C2_GNSS_REALISTIC_OUTAGES",
        "OD_C3_GNSS_PLUS_ACCEL",
    ],
    "phase_b": [
        "ATT_A7_ECLIPSE_THERMAL_TRANSITION",
        "OD_C9_MEAS_ERROR_STRESS_TESTS",
        "OD_C8_STATION_WEATHER_AVAILABILITY",
        "FORM_B6_ROBUST_FORMATION_STORM",
    ],
    "phase_c": [
        "FORM_B2_ALONGTRACK_2SAT_DD_KEEPING",
        "FORM_B3_LVLH_BOX_CONSTRAINT_3SAT",
        "FORM_B5_CLOSE_APPROACH_RISK",
        "FORM_B7_CLUSTER_SPARSE_GROUND_CONTACT",
    ],
}

AREA_KEYS = ("att", "mis", "orb", "aero")


def _timestamp() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d_%H%M%S")


def _run(cmd: list[str]) -> tuple[int, str]:
    env = dict(os.environ)
    extra = f"{ROOT}:{ROOT / 'python'}"
    if env.get("PYTHONPATH"):
        env["PYTHONPATH"] = extra + ":" + env["PYTHONPATH"]
    else:
        env["PYTHONPATH"] = extra
    p = subprocess.run(
        cmd,
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return int(p.returncode), p.stdout


def _apply_config_overrides(
    config_path: Path,
    *,
    max_particles: int | None,
    disable_pod: bool,
) -> None:
    if max_particles is None and not disable_pod:
        return
    payload = json.loads(config_path.read_text(encoding="utf-8"))
    scenarios = payload.get("scenarios", [])
    for scenario in scenarios:
        if max_particles is not None:
            try:
                particles = int(scenario.get("particles", max_particles))
            except Exception:
                particles = max_particles
            scenario["particles"] = int(min(particles, max_particles))
        if disable_pod:
            scenario["pod_uq"] = False
    config_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _selected_phases(phase: str) -> list[str]:
    if phase in AREA_KEYS:
        return [phase]
    if phase == "all":
        return []
    return [phase]


def _catalog_area_ids(catalog_path: Path) -> dict[str, list[str]]:
    try:
        import yaml  # type: ignore
    except Exception:
        return {k: [] for k in AREA_KEYS}
    text = catalog_path.read_text(encoding="utf-8")
    blocks = []
    # Keep regex local here to avoid introducing a hard dependency on converter internals.
    import re

    pattern = re.compile(r"```yaml\s*\n(.*?)```", re.DOTALL)
    for m in pattern.finditer(text):
        raw = m.group(1).strip()
        if not raw:
            continue
        try:
            obj = yaml.safe_load(raw)
        except Exception:
            continue
        if isinstance(obj, dict):
            blocks.append(obj)
    out: dict[str, list[str]] = {k: [] for k in AREA_KEYS}
    for obj in blocks:
        sid = str(obj.get("scenario_id", "")).strip()
        area = str(obj.get("area", "")).strip().lower()
        if not sid or sid == "<string>":
            continue
        if area in out:
            out[area].append(sid)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Phase A/B/C campaign and publish standardized summaries.")
    parser.add_argument("--catalog", default="scenario_catalog.md")
    parser.add_argument(
        "--phase",
        choices=("all", "att", "mis", "orb", "aero", "phase_a", "phase_b", "phase_c"),
        default="all",
    )
    parser.add_argument(
        "--ids",
        default=None,
        help="Optional comma-separated scenario IDs; overrides phase/area selection.",
    )
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--default_particles", type=int, default=256)
    parser.add_argument("--max_particles", type=int, default=None)
    parser.add_argument("--disable_pod", action="store_true", help="Disable POD UQ in generated configs.")
    parser.add_argument("--duration_scale", type=float, default=1.0)
    parser.add_argument("--dt_floor", type=float, default=1.0)
    parser.add_argument("--start_utc", default=None)
    parser.add_argument("--use_env_sources", action="store_true")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--strict_campaign",
        action="store_true",
        help="Fail unless campaign report passes all phases.",
    )
    args = parser.parse_args()

    outdir = (
        Path(args.outdir)
        if args.outdir
        else (ROOT / "results" / f"run_campaign_{_timestamp()}")
    )
    outdir.mkdir(parents=True, exist_ok=True)
    cfg_dir = outdir / "configs"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    catalog_path = Path(args.catalog)
    area_ids = _catalog_area_ids(catalog_path)

    rows: list[dict[str, Any]] = []
    overall_ok = True

    manual_ids = None
    if args.ids:
        manual_ids = [x.strip() for x in str(args.ids).split(",") if x.strip()]

    selected = _selected_phases(args.phase)
    if args.phase == "all" and manual_ids is None:
        any_new = any(len(area_ids[k]) > 0 for k in AREA_KEYS)
        if any_new:
            selected = [k for k in AREA_KEYS if area_ids[k]]
        else:
            selected = ["phase_a", "phase_b", "phase_c"]

    if manual_ids is not None:
        selected = ["custom_ids"]

    for phase in selected:
        if phase in AREA_KEYS:
            ids = area_ids.get(phase, [])
        elif phase in PHASE_SCENARIO_IDS:
            ids = PHASE_SCENARIO_IDS[phase]
        elif phase == "custom_ids":
            ids = manual_ids or []
        else:
            ids = []
        if not ids:
            rows.append(
                {
                    "name": f"{phase}_config",
                    "command": None,
                    "returncode": 0,
                    "pass": True,
                    "note": "no scenarios selected for this group",
                }
            )
            continue
        config_path = cfg_dir / f"{phase}.json"
        convert_cmd = [
            args.python,
            "scripts/catalog_to_case_config.py",
            "--catalog",
            str(args.catalog),
            "--out",
            str(config_path),
            "--ids",
            ",".join(ids),
            "--default_particles",
            str(args.default_particles),
            "--duration_scale",
            str(args.duration_scale),
            "--dt_floor",
            str(args.dt_floor),
        ]
        if not args.use_env_sources:
            convert_cmd.append("--no_env_sources")
        if args.start_utc:
            convert_cmd.extend(["--start_utc", str(args.start_utc)])

        rc_cfg, out_cfg = _run(convert_cmd)
        (outdir / f"log_{phase}_config.txt").write_text(out_cfg, encoding="utf-8")
        rows.append(
            {
                "name": f"{phase}_config",
                "command": convert_cmd,
                "returncode": rc_cfg,
                "pass": bool(rc_cfg == 0),
            }
        )
        overall_ok = bool(overall_ok and rc_cfg == 0)
        if rc_cfg != 0:
            continue

        _apply_config_overrides(
            config_path,
            max_particles=args.max_particles,
            disable_pod=bool(args.disable_pod),
        )

        run_cmd = [
            args.python,
            "scripts/run_case_studies.py",
            "--config",
            str(config_path),
            "--threads",
            str(args.threads),
            "--seed",
            str(args.seed),
            "--outdir",
            str(outdir),
        ]
        rc_run, out_run = _run(run_cmd)
        (outdir / f"log_{phase}_run.txt").write_text(out_run, encoding="utf-8")
        rows.append(
            {
                "name": f"{phase}_run",
                "command": run_cmd,
                "returncode": rc_run,
                "pass": bool(rc_run == 0),
            }
        )
        overall_ok = bool(overall_ok and rc_run == 0)

    run_summary_cmd = [
        args.python,
        "scripts/aggregate_run_summary.py",
        "--run_dir",
        str(outdir),
        "--out",
        str(outdir / "run_summary.json"),
    ]
    rc_summary, out_summary = _run(run_summary_cmd)
    (outdir / "log_run_summary.txt").write_text(out_summary, encoding="utf-8")
    rows.append(
        {
            "name": "aggregate_run_summary",
            "command": run_summary_cmd,
            "returncode": rc_summary,
            "pass": bool(rc_summary == 0),
        }
    )
    overall_ok = bool(overall_ok and rc_summary == 0)

    campaign_cmd = [
        args.python,
        "scripts/generate_campaign_report.py",
        "--run_dir",
        str(outdir),
        "--out",
        str(outdir / "campaign_report.json"),
    ]
    rc_campaign, out_campaign = _run(campaign_cmd)
    (outdir / "log_campaign_report.txt").write_text(out_campaign, encoding="utf-8")
    rows.append(
        {
            "name": "campaign_report",
            "command": campaign_cmd,
            "returncode": rc_campaign,
            "pass": bool(rc_campaign == 0),
        }
    )
    if args.strict_campaign:
        overall_ok = bool(overall_ok and rc_campaign == 0)

    summary = {
        "run_dir": str(outdir),
        "phase": args.phase,
        "catalog": str(catalog_path),
        "seed": int(args.seed),
        "duration_scale": float(args.duration_scale),
        "max_particles": None if args.max_particles is None else int(args.max_particles),
        "disable_pod": bool(args.disable_pod),
        "use_env_sources": bool(args.use_env_sources),
        "strict_campaign": bool(args.strict_campaign),
        "pass": bool(overall_ok),
        "steps": rows,
    }
    summary_path = outdir / "phase_campaign_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[phase-campaign] wrote {summary_path}")
    raise SystemExit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
