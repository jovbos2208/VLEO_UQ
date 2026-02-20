#!/usr/bin/env python3
"""Validate deterministic integrator tolerance convergence (VG-1)."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))


def _normalize_quat(q: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(q, axis=1, keepdims=True)
    n[n == 0.0] = 1.0
    return q / n


def _quat_angle_deg(q_ref: np.ndarray, q_cmp: np.ndarray) -> np.ndarray:
    qa = _normalize_quat(q_ref)
    qb = _normalize_quat(q_cmp)
    dots = np.clip(np.abs(np.sum(qa * qb, axis=1)), -1.0, 1.0)
    return np.rad2deg(2.0 * np.arccos(dots))


def _make_env(EnvInputs, t_grid: np.ndarray, eta1: np.ndarray | None, eta2: np.ndarray | None):
    env = []
    for i, _t in enumerate(t_grid):
        e = EnvInputs()
        e.density = 0.0
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        e.wind_I = np.zeros(3)
        if eta1 is not None:
            e.eta1_rad = float(eta1[i])
        if eta2 is not None:
            e.eta2_rad = float(eta2[i])
        env.append(e)
    return env


def _run_case(*, freeze_attitude: bool, rtol: float, atol: float, duration_s: float, dt_s: float):
    from vleo_uq import (
        AeroAdapter,
        DeterministicPropagator,
        EnvInputs,
        PropagatorConfig,
        VehicleParams,
        default_geometry,
    )

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()

    cfg = PropagatorConfig()
    cfg.rtol = float(rtol)
    cfg.atol = float(atol)
    cfg.max_step_s = 120.0
    cfg.freeze_attitude = bool(freeze_attitude)
    cfg.rho_fast_sigma = 0.0
    cfg.rho_bias_sigma = 0.0
    cfg.wind_sigma = 0.0
    prop = DeterministicPropagator(aero, vehicle, cfg)

    mu = cfg.mu_earth_m3_s2
    r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
    v0 = np.array([0.0, math.sqrt(mu / np.linalg.norm(r0)), 0.0], dtype=float)
    x0 = np.zeros(prop.state_size, dtype=float)
    x0[0:3] = r0
    x0[3:6] = v0
    x0[6] = 1.0
    if not freeze_attitude:
        x0[10:13] = np.deg2rad(np.array([0.2, -0.15, 0.1], dtype=float))

    t_grid = np.arange(0.0, duration_s + dt_s, dt_s, dtype=float)
    eta1 = eta2 = None
    if not freeze_attitude:
        eta1 = np.deg2rad(10.0) * np.sin(2.0 * np.pi * t_grid / 120.0)
        eta2 = np.deg2rad(8.0) * np.sin(2.0 * np.pi * t_grid / 150.0 + 0.4)
    env = _make_env(EnvInputs, t_grid, eta1=eta1, eta2=eta2)
    states = prop.propagate(x0, t_grid, env)
    return t_grid, states


def _compute_metrics(x_ref: np.ndarray, x_cmp: np.ndarray, t_grid: np.ndarray) -> dict:
    dr = x_cmp[:, 0:3] - x_ref[:, 0:3]
    dv = x_cmp[:, 3:6] - x_ref[:, 3:6]
    pos = np.linalg.norm(dr, axis=1)
    vel = np.linalg.norm(dv, axis=1)

    mid_idx = int(np.argmin(np.abs(t_grid - 0.5 * t_grid[-1])))
    out = {
        "pos_rms_m": float(np.sqrt(np.mean(pos * pos))),
        "pos_mid_m": float(pos[mid_idx]),
        "pos_final_m": float(pos[-1]),
        "vel_rms_mps": float(np.sqrt(np.mean(vel * vel))),
        "vel_mid_mps": float(vel[mid_idx]),
        "vel_final_mps": float(vel[-1]),
    }
    if x_ref.shape[1] >= 10 and x_cmp.shape[1] >= 10:
        ang = _quat_angle_deg(x_ref[:, 6:10], x_cmp[:, 6:10])
        out["att_rms_deg"] = float(np.sqrt(np.mean(ang * ang)))
        out["att_mid_deg"] = float(ang[mid_idx])
        out["att_final_deg"] = float(ang[-1])
    return out


def run_sweep(
    *,
    duration_s: float,
    dt_s: float,
    rtol_values: list[float],
    ref_rtol: float,
    ref_atol: float,
) -> dict:
    cases = [("mission", True), ("attitude", False)]
    out = {"cases": {}, "pass": True}
    for name, freeze_attitude in cases:
        t_ref, x_ref = _run_case(
            freeze_attitude=freeze_attitude,
            rtol=ref_rtol,
            atol=ref_atol,
            duration_s=duration_s,
            dt_s=dt_s,
        )

        metrics = []
        for rtol in rtol_values:
            atol = max(ref_atol, rtol * 1e-2)
            t_cmp, x_cmp = _run_case(
                freeze_attitude=freeze_attitude,
                rtol=rtol,
                atol=atol,
                duration_s=duration_s,
                dt_s=dt_s,
            )
            if t_cmp.shape != t_ref.shape or not np.allclose(t_cmp, t_ref):
                raise RuntimeError("time grids are inconsistent across sweep")
            m = _compute_metrics(x_ref, x_cmp, t_ref)
            m["rtol"] = float(rtol)
            m["atol"] = float(atol)
            metrics.append(m)

        metrics = sorted(metrics, key=lambda x: x["rtol"], reverse=True)
        coarse = metrics[0]
        fine = metrics[-1]
        # Convergence may appear as strict improvement or as saturation at machine precision.
        case_pass = (
            fine["pos_final_m"] <= coarse["pos_final_m"] + 1e-12
            and fine["vel_final_mps"] <= coarse["vel_final_mps"] + 1e-15
            and fine["pos_final_m"] <= 5.0
            and fine["vel_final_mps"] <= 0.02
        )
        if "att_final_deg" in fine and "att_final_deg" in coarse:
            case_pass = case_pass and (
                fine["att_final_deg"] <= coarse["att_final_deg"] + 1e-12
                and fine["att_final_deg"] <= 1e-2
            )
        out["cases"][name] = {"metrics": metrics, "pass": bool(case_pass)}
        out["pass"] = bool(out["pass"] and case_pass)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrator tolerance convergence validation.")
    parser.add_argument("--out", default=None, help="Optional JSON output path.")
    parser.add_argument("--duration_s", type=float, default=1200.0)
    parser.add_argument("--dt_s", type=float, default=20.0)
    parser.add_argument("--rtol_values", default="1e-5,1e-6,1e-7")
    parser.add_argument("--ref_rtol", type=float, default=1e-9)
    parser.add_argument("--ref_atol", type=float, default=1e-11)
    args = parser.parse_args()

    rtol_values = [float(x.strip()) for x in args.rtol_values.split(",") if x.strip()]
    if len(rtol_values) < 2:
        raise ValueError("Need at least two rtol values for a sweep")

    result = run_sweep(
        duration_s=float(args.duration_s),
        dt_s=float(args.dt_s),
        rtol_values=rtol_values,
        ref_rtol=float(args.ref_rtol),
        ref_atol=float(args.ref_atol),
    )
    payload = json.dumps(result, indent=2)
    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(payload, encoding="utf-8")
        print(f"[tol] wrote {out_path}")
    print(payload)
    raise SystemExit(0 if result.get("pass", False) else 1)


if __name__ == "__main__":
    main()
