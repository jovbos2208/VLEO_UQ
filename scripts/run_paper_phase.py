#!/usr/bin/env python3
"""Run one paper phase (A/B/C) and generate reports/plots."""

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

PHASE_IDS = {
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

PHASE_CASES = {k: [x.lower() for x in v] for k, v in PHASE_IDS.items()}
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


def _apply_config_overrides(config_path: Path, *, max_particles: int | None, disable_pod: bool) -> None:
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


def _catalog_area_ids(catalog_path: Path) -> dict[str, list[str]]:
    try:
        import yaml  # type: ignore
    except Exception:
        return {k: [] for k in AREA_KEYS}
    import re

    text = catalog_path.read_text(encoding="utf-8")
    pattern = re.compile(r"```yaml\s*\n(.*?)```", re.DOTALL)
    out: dict[str, list[str]] = {k: [] for k in AREA_KEYS}
    for m in pattern.finditer(text):
        raw = m.group(1).strip()
        if not raw:
            continue
        try:
            obj = yaml.safe_load(raw)
        except Exception:
            continue
        if not isinstance(obj, dict):
            continue
        sid = str(obj.get("scenario_id", "")).strip()
        area = str(obj.get("area", "")).strip().lower()
        if sid and sid != "<string>" and area in out:
            out[area].append(sid)
    return out


def _build_phase_report(run_dir: Path, *, phase: str, cases: list[str]) -> dict[str, Any]:
    rows = []
    present = 0
    for case in cases:
        path = run_dir / case / "summary.json"
        ok = path.exists()
        rows.append({"case": case, "present": bool(ok), "summary": str(path)})
        if ok:
            present += 1
    report = {
        "phase": phase,
        "run_dir": str(run_dir),
        "required_cases": int(len(cases)),
        "present_cases": int(present),
        "pass": bool(present == len(cases)),
        "cases": rows,
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Run one paper phase (A/B/C).")
    parser.add_argument(
        "--phase",
        choices=("phase_a", "phase_b", "phase_c", "att", "mis", "orb", "aero"),
        required=True,
    )
    parser.add_argument("--catalog", default="scenario_catalog.md")
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--config_out", default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--default_particles", type=int, default=256)
    parser.add_argument("--max_particles", type=int, default=None)
    parser.add_argument("--duration_scale", type=float, default=1.0)
    parser.add_argument("--dt_floor", type=float, default=1.0)
    parser.add_argument("--disable_pod", action="store_true")
    parser.add_argument("--use_env_sources", action="store_true")
    parser.add_argument("--start_utc", default=None)
    parser.add_argument("--omni_path", default="data/space_weather/omni/omni2_all_years.dat")
    parser.add_argument("--hwm14_lib", default="data/hwm14/libhwm14.so")
    parser.add_argument("--hwm14_data", default="data/hwm14")
    parser.add_argument("--env_interpolation", default="nearest")
    parser.add_argument("--pod_sp3", default=None)
    parser.add_argument("--pod_eop", default="data/eop/EOP-All.txt")
    parser.add_argument("--pod_use_sat_clock", action="store_true")
    parser.add_argument("--pod_estimator", choices=("batch", "enkf"), default="batch")
    parser.add_argument("--pod_enkf_members", type=int, default=128)
    parser.add_argument("--pod_enkf_inflation", type=float, default=1.02)
    parser.add_argument("--pod_enkf_use_carrier", action="store_true")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--plot_all_states", action="store_true")
    parser.add_argument("--make_phase_plots", action="store_true")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    stamp = _timestamp()
    outdir = Path(args.outdir) if args.outdir else (ROOT / "results" / f"paper_{args.phase}_{stamp}")
    outdir.mkdir(parents=True, exist_ok=True)
    config_out = (
        Path(args.config_out)
        if args.config_out
        else (ROOT / "configs" / "generated" / f"paper_{args.phase}_{stamp}.json")
    )
    config_out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    overall_ok = True
    area_ids = _catalog_area_ids(Path(args.catalog))
    if args.phase in PHASE_IDS:
        ids = PHASE_IDS[args.phase]
        cases_for_report = PHASE_CASES[args.phase]
    else:
        ids = area_ids.get(args.phase, [])
        cases_for_report = [x.lower() for x in ids]
    if not ids:
        summary = {"phase": args.phase, "run_dir": str(outdir), "pass": False, "steps": rows, "reason": "no scenario ids selected"}
        (outdir / "paper_phase_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        raise SystemExit(1)

    cfg_cmd = [
        args.python,
        "scripts/catalog_to_case_config.py",
        "--catalog",
        str(args.catalog),
        "--out",
        str(config_out),
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
        cfg_cmd.append("--no_env_sources")
    if args.start_utc:
        cfg_cmd.extend(["--start_utc", str(args.start_utc)])
    rc, out = _run(cfg_cmd)
    (outdir / "log_build_config.txt").write_text(out, encoding="utf-8")
    rows.append({"step": "build_config", "returncode": rc, "pass": bool(rc == 0), "command": cfg_cmd})
    overall_ok = bool(overall_ok and rc == 0)
    if rc != 0:
        summary = {"phase": args.phase, "run_dir": str(outdir), "pass": False, "steps": rows}
        (outdir / "paper_phase_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        raise SystemExit(1)

    _apply_config_overrides(config_out, max_particles=args.max_particles, disable_pod=bool(args.disable_pod))

    run_cmd = [
        args.python,
        "scripts/run_case_studies.py",
        "--config",
        str(config_out),
        "--threads",
        str(args.threads),
        "--seed",
        str(args.seed),
        "--outdir",
        str(outdir),
        "--pod_estimator",
        str(args.pod_estimator),
        "--pod_enkf_members",
        str(args.pod_enkf_members),
        "--pod_enkf_inflation",
        str(args.pod_enkf_inflation),
        "--pod_eop",
        str(args.pod_eop),
    ]
    if args.use_env_sources:
        run_cmd.extend(
            [
                "--use_env_sources",
                "--omni_path",
                str(args.omni_path),
                "--hwm14_lib",
                str(args.hwm14_lib),
                "--hwm14_data",
                str(args.hwm14_data),
                "--env_interpolation",
                str(args.env_interpolation),
            ]
        )
        if args.start_utc:
            run_cmd.extend(["--start_utc", str(args.start_utc)])
    if args.pod_sp3:
        run_cmd.extend(["--pod_sp3", str(args.pod_sp3)])
    if args.pod_use_sat_clock:
        run_cmd.append("--pod_use_sat_clock")
    if args.pod_enkf_use_carrier:
        run_cmd.append("--pod_enkf_use_carrier")
    if args.plot:
        run_cmd.append("--plot")
    if args.plot_all_states:
        run_cmd.append("--plot_all_states")

    rc, out = _run(run_cmd)
    (outdir / "log_run_case_studies.txt").write_text(out, encoding="utf-8")
    rows.append({"step": "run_case_studies", "returncode": rc, "pass": bool(rc == 0), "command": run_cmd})
    overall_ok = bool(overall_ok and rc == 0)

    agg_cmd = [
        args.python,
        "scripts/aggregate_run_summary.py",
        "--run_dir",
        str(outdir),
        "--out",
        str(outdir / "run_summary.json"),
    ]
    rc, out = _run(agg_cmd)
    (outdir / "log_aggregate_run_summary.txt").write_text(out, encoding="utf-8")
    rows.append({"step": "aggregate_run_summary", "returncode": rc, "pass": bool(rc == 0), "command": agg_cmd})
    overall_ok = bool(overall_ok and rc == 0)

    phase_report = _build_phase_report(outdir, phase=args.phase, cases=cases_for_report)
    (outdir / "phase_report.json").write_text(json.dumps(phase_report, indent=2), encoding="utf-8")
    overall_ok = bool(overall_ok and bool(phase_report.get("pass", False)))

    if args.make_phase_plots and args.phase in {"phase_a", "phase_b", "phase_c"}:
        plot_cmd = [
            args.python,
            "scripts/plot_phase_results.py",
            "--run_dir",
            str(outdir),
            "--phase",
            str(args.phase),
            "--outdir",
            str(outdir / "paper_plots"),
        ]
        rc, out = _run(plot_cmd)
        (outdir / "log_phase_plots.txt").write_text(out, encoding="utf-8")
        rows.append({"step": "phase_plots", "returncode": rc, "pass": bool(rc == 0), "command": plot_cmd})
        overall_ok = bool(overall_ok and rc == 0)

    summary = {
        "phase": args.phase,
        "run_dir": str(outdir),
        "config": str(config_out),
        "seed": int(args.seed),
        "threads": int(args.threads),
        "pass": bool(overall_ok),
        "steps": rows,
        "phase_report": phase_report,
    }
    summary_path = outdir / "paper_phase_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[paper-phase] wrote {summary_path}")
    raise SystemExit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
