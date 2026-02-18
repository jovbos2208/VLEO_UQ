#!/usr/bin/env python3
"""Convergence study for DET+STM across timestep floors by scenario type.

This script builds case configs from `scenario_catalog.md`, runs `run_det_stm_phase.py`
for multiple timestep floors, and compares each run against a reference timestep.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np


def parse_dt_values(raw: str) -> list[float]:
    vals = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        vals.append(float(token))
    if not vals:
        raise ValueError("dt_values is empty")
    return sorted(set(vals))


def parse_ids_by_type(catalog_path: Path) -> tuple[list[str], list[str]]:
    text = catalog_path.read_text(encoding="utf-8")
    ids = re.findall(r"^scenario_id:\s*([A-Za-z0-9_]+)\s*$", text, flags=re.MULTILINE)
    ids = [x for x in ids if x != "<string>"]
    att = [x for x in ids if x.startswith("ATT_A")]
    od = [x for x in ids if x.startswith("OD_C")]
    return att, od


def run_command(cmd: list[str]) -> None:
    print(f"[conv] run: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def build_config(
    *,
    python: str,
    catalog: Path,
    out_config: Path,
    ids: list[str],
    start_utc: str | None,
    dt_floor: float,
) -> dict:
    cmd = [
        python,
        "scripts/catalog_to_case_config.py",
        "--catalog",
        str(catalog),
        "--out",
        str(out_config),
        "--ids",
        ",".join(ids),
        "--dt_floor",
        str(dt_floor),
    ]
    if start_utc:
        cmd.extend(["--start_utc", start_utc])
    run_command(cmd)
    return json.loads(out_config.read_text(encoding="utf-8"))


def run_det_stm(
    *,
    python: str,
    config: Path,
    outdir: Path,
    use_env_sources: bool,
    start_utc: str | None,
    omni_path: str,
    hwm14_lib: str,
    hwm14_data: str,
    env_interpolation: str,
) -> None:
    cmd = [
        python,
        "scripts/run_det_stm_phase.py",
        "--config",
        str(config),
        "--outdir",
        str(outdir),
    ]
    if use_env_sources:
        cmd.append("--use_env_sources")
    if start_utc:
        cmd.extend(["--start_utc", start_utc])
    cmd.extend(
        [
            "--omni_path",
            omni_path,
            "--hwm14_lib",
            hwm14_lib,
            "--hwm14_data",
            hwm14_data,
            "--env_interpolation",
            env_interpolation,
        ]
    )
    run_command(cmd)


def interp_states(t_ref: np.ndarray, t_cmp: np.ndarray, x_cmp: np.ndarray) -> np.ndarray:
    out = np.zeros((t_ref.size, x_cmp.shape[1]), dtype=float)
    for j in range(x_cmp.shape[1]):
        out[:, j] = np.interp(t_ref, t_cmp, x_cmp[:, j])
    return out


def normalize_quat(q: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(q, axis=1, keepdims=True)
    n[n == 0.0] = 1.0
    return q / n


def quat_angle_error_deg(q_ref: np.ndarray, q_cmp: np.ndarray) -> np.ndarray:
    q_ref_n = normalize_quat(q_ref)
    q_cmp_n = normalize_quat(q_cmp)
    dots = np.abs(np.sum(q_ref_n * q_cmp_n, axis=1))
    dots = np.clip(dots, -1.0, 1.0)
    ang = 2.0 * np.arccos(dots)
    return np.rad2deg(ang)


def quat_wxyz_to_dcm_batch(q_wxyz: np.ndarray) -> np.ndarray:
    q = np.array(q_wxyz, dtype=float)
    n = np.linalg.norm(q, axis=1, keepdims=True)
    n[n == 0.0] = 1.0
    q = q / n
    w = q[:, 0]
    x = q[:, 1]
    y = q[:, 2]
    z = q[:, 3]
    R = np.zeros((q.shape[0], 3, 3), dtype=float)
    R[:, 0, 0] = 1.0 - 2.0 * (y * y + z * z)
    R[:, 0, 1] = 2.0 * (x * y - w * z)
    R[:, 0, 2] = 2.0 * (x * z + w * y)
    R[:, 1, 0] = 2.0 * (x * y + w * z)
    R[:, 1, 1] = 1.0 - 2.0 * (x * x + z * z)
    R[:, 1, 2] = 2.0 * (y * z - w * x)
    R[:, 2, 0] = 2.0 * (x * z - w * y)
    R[:, 2, 1] = 2.0 * (y * z + w * x)
    R[:, 2, 2] = 1.0 - 2.0 * (x * x + y * y)
    return R


def boresight_eci_from_quat(q_wxyz: np.ndarray, boresight_body: np.ndarray | None = None) -> np.ndarray:
    if boresight_body is None:
        boresight_body = np.array([0.0, 0.0, 1.0], dtype=float)
    b = np.array(boresight_body, dtype=float).reshape(3)
    R_BI = quat_wxyz_to_dcm_batch(q_wxyz)
    # Convert body-frame boresight to inertial frame: b_I = R_BI^T * b_B.
    return np.einsum("nij,j->ni", np.transpose(R_BI, (0, 2, 1)), b)


def vector_angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    aa = np.array(a, dtype=float)
    bb = np.array(b, dtype=float)
    an = np.linalg.norm(aa, axis=1, keepdims=True)
    bn = np.linalg.norm(bb, axis=1, keepdims=True)
    an[an == 0.0] = 1.0
    bn[bn == 0.0] = 1.0
    aa = aa / an
    bb = bb / bn
    dots = np.clip(np.sum(aa * bb, axis=1), -1.0, 1.0)
    return np.rad2deg(np.arccos(dots))


def specific_energy(mu: float, states: np.ndarray) -> np.ndarray:
    r = states[:, 0:3]
    v = states[:, 3:6]
    rn = np.linalg.norm(r, axis=1)
    vn2 = np.sum(v * v, axis=1)
    return 0.5 * vn2 - mu / rn


def load_det_stm(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dat = np.load(path)
    return dat["t_grid"], dat["det_states"], dat["stm_cov"]


def case_path(base: Path, scenario: dict) -> Path:
    name = scenario["name"]
    typ = scenario["type"]
    if typ == "formation":
        # This script is for ATT and OD(mission) convergence sweeps.
        # Keep explicit guard if a formation case is passed accidentally.
        raise ValueError(f"unsupported type for this study: {typ} ({name})")
    return base / name / "det_stm.npz"


def summarize_rows(rows: list[dict]) -> list[dict]:
    by_dt: dict[float, list[dict]] = {}
    for row in rows:
        by_dt.setdefault(float(row["dt_floor_s"]), []).append(row)
    out = []
    for dt_floor, chunk in sorted(by_dt.items()):
        def mean(key: str) -> float:
            vals = [float(x[key]) for x in chunk if key in x]
            return float(np.mean(vals)) if vals else float("nan")

        def vmax(key: str) -> float:
            vals = [float(x[key]) for x in chunk if key in x]
            return float(np.max(vals)) if vals else float("nan")

        out.append(
            {
                "dt_floor_s": dt_floor,
                "n_cases": len(chunk),
                "mean_pos_rms_m": mean("pos_rms_m"),
                "max_pos_rms_m": vmax("pos_rms_m"),
                "mean_pos_final_m": mean("pos_final_m"),
                "max_pos_final_m": vmax("pos_final_m"),
                "mean_vel_rms_mps": mean("vel_rms_mps"),
                "max_vel_rms_mps": vmax("vel_rms_mps"),
                "mean_energy_drift_rms": mean("energy_rel_drift_rms"),
                "max_energy_drift_rms": vmax("energy_rel_drift_rms"),
                "mean_stm_cov_final_rel_frob": mean("stm_cov_final_rel_frob"),
                "max_stm_cov_final_rel_frob": vmax("stm_cov_final_rel_frob"),
                "mean_attitude_rms_deg": mean("attitude_rms_deg"),
                "max_attitude_rms_deg": vmax("attitude_rms_deg"),
                "mean_boresight_rms_deg": mean("boresight_rms_deg"),
                "max_boresight_rms_deg": vmax("boresight_rms_deg"),
                "mean_boresight_p95_deg": mean("boresight_p95_deg"),
                "max_boresight_p95_deg": vmax("boresight_p95_deg"),
                "mean_rate_rms_deg_s": mean("rate_rms_deg_s"),
                "max_rate_rms_deg_s": vmax("rate_rms_deg_s"),
            }
        )
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = sorted({k for r in rows for k in r.keys()})
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def compare_type(
    *,
    type_name: str,
    config_by_dt: dict[float, dict],
    out_by_dt: dict[float, Path],
    ref_dt: float,
    mu: float,
) -> list[dict]:
    ref_cfg = config_by_dt[ref_dt]
    ref_scenarios = ref_cfg.get("scenarios", [])
    rows: list[dict] = []

    for scenario in ref_scenarios:
        ref_path = case_path(out_by_dt[ref_dt], scenario)
        t_ref, x_ref, p_ref = load_det_stm(ref_path)
        e_ref = specific_energy(mu, x_ref)
        e_ref_drift = e_ref - e_ref[0]
        ref_scale = max(float(np.linalg.norm(p_ref[-1], ord="fro")), 1e-16)

        for dt_floor, cfg in sorted(config_by_dt.items()):
            if dt_floor == ref_dt:
                continue
            cmp_path = case_path(out_by_dt[dt_floor], scenario)
            t_cmp, x_cmp_raw, p_cmp = load_det_stm(cmp_path)
            x_cmp = interp_states(t_ref, t_cmp, x_cmp_raw)

            dr = x_cmp[:, 0:3] - x_ref[:, 0:3]
            dv = x_cmp[:, 3:6] - x_ref[:, 3:6]
            pos_norm = np.linalg.norm(dr, axis=1)
            vel_norm = np.linalg.norm(dv, axis=1)

            e_cmp = specific_energy(mu, x_cmp)
            e_cmp_drift = e_cmp - e_cmp[0]
            de_drift = e_cmp_drift - e_ref_drift

            row = {
                "scenario": scenario["name"],
                "scenario_type": type_name,
                "dt_floor_s": float(dt_floor),
                "reference_dt_s": float(ref_dt),
                "pos_rms_m": float(np.sqrt(np.mean(pos_norm * pos_norm))),
                "pos_final_m": float(pos_norm[-1]),
                "vel_rms_mps": float(np.sqrt(np.mean(vel_norm * vel_norm))),
                "energy_rel_drift_rms": float(np.sqrt(np.mean(de_drift * de_drift))),
                "stm_cov_final_rel_frob": float(np.linalg.norm(p_cmp[-1] - p_ref[-1], ord="fro") / ref_scale),
            }

            if scenario.get("type") == "attitude" and x_ref.shape[1] >= 13:
                q_ref = x_ref[:, 6:10]
                q_cmp = x_cmp[:, 6:10]
                ang_err = quat_angle_error_deg(q_ref, q_cmp)
                w_ref = x_ref[:, 10:13]
                w_cmp = x_cmp[:, 10:13]
                w_err = np.linalg.norm(w_cmp - w_ref, axis=1)
                b_ref = boresight_eci_from_quat(q_ref)
                b_cmp = boresight_eci_from_quat(q_cmp)
                b_err = vector_angle_deg(b_ref, b_cmp)
                row["attitude_rms_deg"] = float(np.sqrt(np.mean(ang_err * ang_err)))
                row["boresight_rms_deg"] = float(np.sqrt(np.mean(b_err * b_err)))
                row["boresight_p95_deg"] = float(np.percentile(b_err, 95.0))
                row["rate_rms_deg_s"] = float(np.rad2deg(np.sqrt(np.mean(w_err * w_err))))

            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DET+STM convergence study by scenario type.")
    parser.add_argument("--catalog", default="scenario_catalog.md")
    parser.add_argument("--outdir", default=None, help="Default: results/convergence_<UTC>")
    parser.add_argument("--types", default="attitude,mission_od", help="Comma-separated: attitude,mission_od")
    parser.add_argument("--dt_values", default="1,2,5,10,30,60")
    parser.add_argument("--reference_dt", type=float, default=None, help="Default: smallest dt in dt_values")
    parser.add_argument("--max_scenarios_per_type", type=int, default=0, help="0 means all")
    parser.add_argument("--start_utc", default=None)
    parser.add_argument("--use_env_sources", action="store_true")
    parser.add_argument("--omni_path", default="data/space_weather/omni/omni2_all_years.dat")
    parser.add_argument("--hwm14_lib", default="data/hwm14/libhwm14.so")
    parser.add_argument("--hwm14_data", default="data/hwm14")
    parser.add_argument("--env_interpolation", default="nearest")
    args = parser.parse_args()

    catalog = Path(args.catalog)
    if not catalog.is_file():
        raise FileNotFoundError(catalog)

    dt_values = parse_dt_values(args.dt_values)
    ref_dt = float(args.reference_dt) if args.reference_dt is not None else float(min(dt_values))
    if ref_dt not in dt_values:
        raise ValueError("reference_dt must be one of dt_values")

    if args.outdir is None:
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d_%H%M%S")
        outdir = Path("results") / f"convergence_{stamp}"
    else:
        outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_requested = {x.strip().lower() for x in args.types.split(",") if x.strip()}
    alias = {
        "a": "attitude",
        "att": "attitude",
        "attitude": "attitude",
        "c": "mission_od",
        "od": "mission_od",
        "mission_od": "mission_od",
        "orbit_determination": "mission_od",
    }
    requested = {alias.get(x, x) for x in raw_requested}
    valid = {"attitude", "mission_od"}
    bad = requested - valid
    if bad:
        raise ValueError(f"unknown types: {sorted(bad)}")

    att_ids, od_ids = parse_ids_by_type(catalog)
    if args.max_scenarios_per_type and args.max_scenarios_per_type > 0:
        att_ids = att_ids[: args.max_scenarios_per_type]
        od_ids = od_ids[: args.max_scenarios_per_type]

    type_to_ids = {
        "attitude": att_ids,
        "mission_od": od_ids,
    }

    python = sys.executable
    mu = 3.986004418e14
    all_rows = []

    for type_name in ("attitude", "mission_od"):
        if type_name not in requested:
            continue
        ids = type_to_ids[type_name]
        if not ids:
            print(f"[conv] skip {type_name}: no scenario ids found", flush=True)
            continue

        type_dir = outdir / type_name
        type_dir.mkdir(parents=True, exist_ok=True)
        config_by_dt: dict[float, dict] = {}
        run_by_dt: dict[float, Path] = {}

        for dt_floor in dt_values:
            dt_tag = f"dt_{str(dt_floor).replace('.', 'p')}"
            cfg_path = type_dir / f"config_{dt_tag}.json"
            run_dir = type_dir / dt_tag
            run_dir.mkdir(parents=True, exist_ok=True)

            cfg = build_config(
                python=python,
                catalog=catalog,
                out_config=cfg_path,
                ids=ids,
                start_utc=args.start_utc,
                dt_floor=dt_floor,
            )
            config_by_dt[dt_floor] = cfg
            run_by_dt[dt_floor] = run_dir

            run_det_stm(
                python=python,
                config=cfg_path,
                outdir=run_dir,
                use_env_sources=args.use_env_sources,
                start_utc=args.start_utc,
                omni_path=args.omni_path,
                hwm14_lib=args.hwm14_lib,
                hwm14_data=args.hwm14_data,
                env_interpolation=args.env_interpolation,
            )

        rows = compare_type(
            type_name=type_name,
            config_by_dt=config_by_dt,
            out_by_dt=run_by_dt,
            ref_dt=ref_dt,
            mu=mu,
        )
        all_rows.extend(rows)
        summary = summarize_rows(rows)

        (type_dir / "metrics.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
        (type_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        write_csv(type_dir / "metrics.csv", rows)
        write_csv(type_dir / "summary.csv", summary)
        print(f"[conv] wrote {type_dir}/metrics.json and summary files", flush=True)

    (outdir / "metrics_all.json").write_text(json.dumps(all_rows, indent=2), encoding="utf-8")
    write_csv(outdir / "metrics_all.csv", all_rows)
    print(f"[conv] done -> {outdir}", flush=True)


if __name__ == "__main__":
    main()
