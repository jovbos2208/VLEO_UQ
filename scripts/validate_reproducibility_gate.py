#!/usr/bin/env python3
"""VG-6 reproducibility gate: same seed should be invariant across thread counts."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))


def _build_env(EnvInputs, t_grid: np.ndarray) -> list:
    env = []
    for _ in t_grid:
        e = EnvInputs()
        e.density = 0.0
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        e.wind_I = np.zeros(3)
        env.append(e)
    return env


def _make_x0(state_size: int, mu_earth: float) -> np.ndarray:
    r0 = np.array([7000e3, 0.0, 0.0], dtype=float)
    v0 = np.array([0.0, math.sqrt(mu_earth / np.linalg.norm(r0)), 0.0], dtype=float)
    x0 = np.zeros(state_size, dtype=float)
    x0[0:3] = r0
    x0[3:6] = v0
    x0[6] = 1.0
    return x0


def _sample_mc_initial(x0: np.ndarray, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    X0 = np.repeat(x0[None, :], n, axis=0)
    X0[:, 0:3] += rng.normal(0.0, 2.0, size=(n, 3))
    X0[:, 3:6] += rng.normal(0.0, 0.02, size=(n, 3))
    X0[:, 6:10] = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    X0[:, 10:13] = 0.0
    return X0


def _array_compare(a: np.ndarray, b: np.ndarray, atol: float) -> dict:
    if a.shape != b.shape:
        return {"equal": False, "shape_equal": False, "max_abs_diff": None}
    diff = np.abs(a - b)
    max_abs = float(np.max(diff)) if diff.size > 0 else 0.0
    equal = bool(np.array_equal(a, b) if atol <= 0.0 else np.all(diff <= atol))
    return {"equal": equal, "shape_equal": True, "max_abs_diff": max_abs}


def run_gate(
    *,
    threads: list[int],
    mc_particles: int,
    duration_s: float,
    dt_s: float,
    seed: int,
    atol: float,
) -> dict:
    from vleo_uq import (
        AeroAdapter,
        DeterministicPropagator,
        EnsemblePropagatorMC,
        EnvInputs,
        PropagatorConfig,
        VehicleParams,
        default_geometry,
    )

    if len(threads) < 2:
        raise ValueError("need at least two thread counts for reproducibility gate")

    t_grid = np.arange(0.0, duration_s + dt_s, dt_s, dtype=float)
    ref_det = None
    ref_mc = None
    ref_thread = None
    comparisons = []

    prev_env = {
        "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
        "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
        "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
    }
    try:
        for th in threads:
            os.environ["OMP_NUM_THREADS"] = str(int(th))
            os.environ["OPENBLAS_NUM_THREADS"] = "1"
            os.environ["MKL_NUM_THREADS"] = "1"

            geom = default_geometry()
            aero = AeroAdapter()
            aero.init(geom)
            vehicle = VehicleParams()
            cfg = PropagatorConfig()
            cfg.rho_fast_sigma = 0.0
            cfg.rho_bias_sigma = 0.0
            cfg.wind_sigma = 0.0
            cfg.rng_seed = int(seed)
            cfg.rtol = 1e-7
            cfg.atol = 1e-9

            prop_det = DeterministicPropagator(aero, vehicle, cfg)
            prop_mc = EnsemblePropagatorMC(aero, vehicle, cfg)
            x0 = _make_x0(prop_det.state_size, cfg.mu_earth_m3_s2)
            env = _build_env(EnvInputs, t_grid)
            X0 = _sample_mc_initial(x0, int(mc_particles), int(seed))

            det = prop_det.propagate(x0, t_grid, env)
            mc = prop_mc.propagate(X0, t_grid, env)

            if ref_det is None:
                ref_det = det
                ref_mc = mc
                ref_thread = int(th)
                continue

            det_cmp = _array_compare(ref_det, det, atol=atol)
            mc_cmp = _array_compare(ref_mc, mc, atol=atol)
            comparisons.append(
                {
                    "reference_thread": int(ref_thread),
                    "thread": int(th),
                    "deterministic": det_cmp,
                    "mc": mc_cmp,
                }
            )
    finally:
        for k, v in prev_env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    all_pass = all(
        bool(row["deterministic"]["equal"] and row["mc"]["equal"])
        for row in comparisons
    )
    return {
        "pass": bool(all_pass),
        "threads": [int(t) for t in threads],
        "mc_particles": int(mc_particles),
        "duration_s": float(duration_s),
        "dt_s": float(dt_s),
        "seed": int(seed),
        "atol": float(atol),
        "comparisons": comparisons,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate reproducibility across thread counts (VG-6).")
    parser.add_argument("--threads", default="1,2", help="Comma-separated thread counts, e.g. 1,2,4")
    parser.add_argument("--mc_particles", type=int, default=64)
    parser.add_argument("--duration_s", type=float, default=300.0)
    parser.add_argument("--dt_s", type=float, default=30.0)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    threads = [int(x.strip()) for x in str(args.threads).split(",") if x.strip()]
    report = run_gate(
        threads=threads,
        mc_particles=int(args.mc_particles),
        duration_s=float(args.duration_s),
        dt_s=float(args.dt_s),
        seed=int(args.seed),
        atol=float(args.atol),
    )
    payload = json.dumps(report, indent=2)
    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(payload, encoding="utf-8")
        print(f"[repro-gate] wrote {out_path}")
    print(payload)
    raise SystemExit(0 if bool(report.get("pass", False)) else 1)


if __name__ == "__main__":
    main()
