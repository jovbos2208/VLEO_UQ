#!/usr/bin/env python3
"""Validate UT/STM coverage realism against MC (VG-3)."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


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


def _sample_mc_initial(x0: np.ndarray, P0: np.ndarray, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    X0 = rng.multivariate_normal(mean=x0, cov=P0, size=n)
    # Keep attitude deterministic for stable short-arc orbital coverage comparison.
    X0[:, 6:10] = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    X0[:, 10:13] = 0.0
    return X0


def _coverage_from_samples(
    samples: np.ndarray,
    mean_pred: np.ndarray,
    cov_pred: np.ndarray,
    levels: tuple[float, ...] = (1.0, 2.0, 3.0),
    indices: tuple[int, ...] = (0, 1, 2, 3, 4, 5),
) -> dict:
    out = {}
    for k in levels:
        vals = []
        for idx in indices:
            sigma = float(np.sqrt(max(cov_pred[idx, idx], 0.0)))
            if sigma <= 0.0:
                continue
            z = np.abs(samples[:, idx] - mean_pred[idx]) / sigma
            vals.append(float(np.mean(z <= k)))
        out[f"{int(k)}sigma"] = float(np.mean(vals)) if vals else 0.0
    return out


def _coverage_pass(coverage: dict) -> bool:
    c1 = float(coverage.get("1sigma", 0.0))
    c2 = float(coverage.get("2sigma", 0.0))
    c3 = float(coverage.get("3sigma", 0.0))
    # Broad acceptance against MC finite-sample variability.
    return (0.50 <= c1 <= 0.90) and (0.80 <= c2 <= 1.0) and (0.90 <= c3 <= 1.0) and (c1 <= c2 <= c3)


def run_validation(
    *,
    mc_particles: int = 192,
    duration_s: float = 600.0,
    dt_s: float = 30.0,
    seed: int = 77,
) -> dict:
    from vleo_uq import (
        AeroAdapter,
        EnsemblePropagatorMC,
        EnvInputs,
        PropagatorConfig,
        SigmaPointPropagatorUT,
        StmPropagator,
        VehicleParams,
        default_geometry,
    )

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()

    cfg = PropagatorConfig()
    cfg.freeze_attitude = True
    cfg.rho_fast_sigma = 0.0
    cfg.rho_bias_sigma = 0.0
    cfg.wind_sigma = 0.0
    cfg.rtol = 1e-7
    cfg.atol = 1e-9
    cfg.rng_seed = int(seed)

    prop_mc = EnsemblePropagatorMC(aero, vehicle, cfg)
    prop_ut = SigmaPointPropagatorUT(aero, vehicle, cfg)
    prop_stm = StmPropagator(aero, vehicle, cfg)

    x0 = _make_x0(prop_mc.state_size, cfg.mu_earth_m3_s2)
    P0 = np.diag(
        [
            10.0,
            10.0,
            10.0,
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
    X0 = _sample_mc_initial(x0, P0, int(mc_particles), int(seed))

    t_grid = np.arange(0.0, duration_s + dt_s, dt_s, dtype=float)
    env = _build_env(EnvInputs, t_grid)

    mc = prop_mc.propagate(X0, t_grid, env)
    ut_mean, ut_cov = prop_ut.propagate(x0, P0, t_grid, env)
    stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)

    samples_f = mc[-1]
    cov_ut = _coverage_from_samples(samples_f, ut_mean[-1], ut_cov[-1])
    cov_stm = _coverage_from_samples(samples_f, stm_mean[-1], stm_cov[-1])

    result = {
        "mc_particles": int(mc_particles),
        "duration_s": float(duration_s),
        "dt_s": float(dt_s),
        "ut_coverage": cov_ut,
        "stm_coverage": cov_stm,
        "ut_pass": _coverage_pass(cov_ut),
        "stm_pass": _coverage_pass(cov_stm),
    }
    result["pass"] = bool(result["ut_pass"] and result["stm_pass"])
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate UT/STM coverage realism against MC.")
    parser.add_argument("--mc_particles", type=int, default=192)
    parser.add_argument("--duration_s", type=float, default=600.0)
    parser.add_argument("--dt_s", type=float, default=30.0)
    parser.add_argument("--seed", type=int, default=77)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    result = run_validation(
        mc_particles=int(args.mc_particles),
        duration_s=float(args.duration_s),
        dt_s=float(args.dt_s),
        seed=int(args.seed),
    )
    payload = json.dumps(result, indent=2)
    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(payload, encoding="utf-8")
        print(f"[coverage] wrote {out_path}")
    print(payload)
    raise SystemExit(0 if result.get("pass", False) else 1)


if __name__ == "__main__":
    main()
