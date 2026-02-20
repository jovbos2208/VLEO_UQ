#!/usr/bin/env python3
"""Generate MC convergence report versus a high-particle reference."""

from __future__ import annotations

import argparse
import csv
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


def _sample_states(x0: np.ndarray, P0: np.ndarray, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    X = rng.multivariate_normal(mean=x0, cov=P0, size=int(n))
    X[:, 6:10] = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    X[:, 10:13] = 0.0
    return X


def _mean_cov(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(samples, axis=0)
    cov = np.cov(samples, rowvar=False)
    return mean, cov


def run_report(
    *,
    particles: list[int],
    ref_particles: int,
    seed: int,
    duration_s: float,
    dt_s: float,
) -> dict:
    from vleo_uq import (
        AeroAdapter,
        EnsemblePropagatorMC,
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
    cfg.rho_fast_sigma = 0.0
    cfg.rho_bias_sigma = 0.0
    cfg.wind_sigma = 0.0
    cfg.rng_seed = int(seed)
    cfg.freeze_attitude = True
    cfg.rtol = 1e-7
    cfg.atol = 1e-9
    prop = EnsemblePropagatorMC(aero, vehicle, cfg)

    t_grid = np.arange(0.0, duration_s + dt_s, dt_s, dtype=float)
    env = _build_env(EnvInputs, t_grid)
    x0 = _make_x0(prop.state_size, cfg.mu_earth_m3_s2)
    P0 = np.diag(
        [
            9.0,
            9.0,
            9.0,
            1e-2,
            1e-2,
            1e-2,
            1e-8,
            1e-8,
            1e-8,
            1e-8,
            1e-10,
            1e-10,
            1e-10,
            1e-4,
            1e-4,
            1e-4,
            1e-4,
            1e-4,
        ]
    )

    X_ref = _sample_states(x0, P0, int(ref_particles), int(seed))
    Y_ref = prop.propagate(X_ref, t_grid, env)[-1]
    mean_ref, cov_ref = _mean_cov(Y_ref)
    norm_cov_ref = max(float(np.linalg.norm(cov_ref, ord="fro")), 1e-12)

    rows = []
    for n in sorted(set(int(x) for x in particles if int(x) > 0)):
        X = _sample_states(x0, P0, n, int(seed))
        Y = prop.propagate(X, t_grid, env)[-1]
        mean, cov = _mean_cov(Y)
        row = {
            "particles": int(n),
            "mean_err_norm": float(np.linalg.norm(mean - mean_ref)),
            "mean_err_pos_norm_m": float(np.linalg.norm(mean[0:3] - mean_ref[0:3])),
            "cov_rel_frob_err": float(np.linalg.norm(cov - cov_ref, ord="fro") / norm_cov_ref),
        }
        rows.append(row)

    rows.sort(key=lambda r: int(r["particles"]))
    monotone = True
    for i in range(1, len(rows)):
        if rows[i]["cov_rel_frob_err"] > rows[i - 1]["cov_rel_frob_err"] + 1e-12:
            monotone = False
            break
    pass_flag = bool(len(rows) >= 2 and rows[-1]["cov_rel_frob_err"] <= 0.35 and monotone)
    return {
        "pass": pass_flag,
        "seed": int(seed),
        "duration_s": float(duration_s),
        "dt_s": float(dt_s),
        "reference_particles": int(ref_particles),
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate MC convergence report.")
    parser.add_argument("--particles", default="32,64,128")
    parser.add_argument("--ref_particles", type=int, default=512)
    parser.add_argument("--seed", type=int, default=2030)
    parser.add_argument("--duration_s", type=float, default=600.0)
    parser.add_argument("--dt_s", type=float, default=30.0)
    parser.add_argument("--out", default="results/reports/mc_convergence_report.json")
    parser.add_argument("--out_csv", default=None)
    args = parser.parse_args()

    particles = [int(x.strip()) for x in str(args.particles).split(",") if x.strip()]
    report = run_report(
        particles=particles,
        ref_particles=int(args.ref_particles),
        seed=int(args.seed),
        duration_s=float(args.duration_s),
        dt_s=float(args.dt_s),
    )
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    csv_path = Path(args.out_csv) if args.out_csv else out_path.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["particles", "mean_err_norm", "mean_err_pos_norm_m", "cov_rel_frob_err"])
        w.writeheader()
        for row in report.get("rows", []):
            w.writerow(row)
    print(f"[mc-convergence] wrote {csv_path}")
    print(f"[mc-convergence] wrote {out_path}")
    raise SystemExit(0 if bool(report.get("pass", False)) else 1)


if __name__ == "__main__":
    main()
