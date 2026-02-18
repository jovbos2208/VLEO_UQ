#!/usr/bin/env python3
"""Benchmark deterministic, MC, UT, and STM propagators with timing and stats."""

from __future__ import annotations

import argparse
import os
import time

import numpy as np


def build_env(t_grid: np.ndarray, EnvInputs, density: float, temperature_K: float, particle_mass_kg: float) -> list:
    env = []
    for _ in t_grid:
        e = EnvInputs()
        e.density = density
        e.temperature_K = temperature_K
        e.particles_mass_kg = particle_mass_kg
        e.wind_I = np.zeros(3)
        env.append(e)
    return env


def summarize(values: list[float]) -> tuple[float, float]:
    arr = np.array(values, dtype=float)
    return float(arr.mean()), float(arr.std(ddof=1)) if arr.size > 1 else 0.0


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark propagators.")
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--particles", type=int, default=32)
    parser.add_argument("--duration_s", type=float, default=600.0)
    parser.add_argument("--dt_s", type=float, default=10.0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--augment_process_noise", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))

    from vleo_uq import (
        AeroAdapter,
        DeterministicPropagator,
        EnsemblePropagatorMC,
        EnvInputs,
        PropagatorConfig,
        SigmaPointPropagatorUT,
        StmPropagator,
        VehicleParams,
        default_geometry,
    )

    rng = np.random.default_rng(args.seed)

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()
    config = PropagatorConfig()
    config.rtol = 1e-6
    config.atol = 1e-6
    config.max_step_s = 120.0
    if args.verbose:
        config.verbose = True
        config.progress_stride = 10
        config.particle_stride = 5
    if args.debug:
        config.debug_state = True
        config.debug_stride = 1
    config.rho_fast_sigma = 0.1
    config.rho_fast_tau_s = 600.0
    config.rho_bias_sigma = 0.02
    config.rho_bias_tau_s = 3600.0
    config.wind_sigma = 5.0
    config.wind_tau_s = 600.0

    prop_det = DeterministicPropagator(aero, vehicle, config)
    prop_mc = EnsemblePropagatorMC(aero, vehicle, config)
    prop_ut = SigmaPointPropagatorUT(aero, vehicle, config)
    prop_stm = StmPropagator(aero, vehicle, config)

    t_grid = np.arange(0.0, args.duration_s + args.dt_s, args.dt_s)
    env = build_env(t_grid, EnvInputs, density=1e-12, temperature_K=1000.0,
                    particle_mass_kg=28.0 * 1.6605390689252e-27)

    mu = config.mu_earth_m3_s2
    r0 = np.array([7000e3, 0.0, 0.0])
    v0 = np.array([0.0, np.sqrt(mu / np.linalg.norm(r0)), 0.0])
    x0 = np.zeros(prop_det.state_size)
    x0[0:3] = r0
    x0[3:6] = v0
    x0[6] = 1.0
    x0[10:13] = np.array([0.0, 0.0, 0.01])
    P0 = np.diag([10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [0.0] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3)

    det_times = []
    mc_times = []
    ut_times = []
    stm_times = []
    ut_mean_errs = []
    stm_mean_errs = []
    ut_cov_errs = []
    stm_cov_errs = []

    print("Propagator benchmark")
    print(f"Runs={args.runs}  Particles={args.particles}  dt={args.dt_s}s  duration={args.duration_s}s  threads={args.threads}")

    for run in range(args.runs):
        config.rng_seed = args.seed + run

        t0 = time.perf_counter()
        det_states = prop_det.propagate(x0, t_grid, env)
        det_times.append(time.perf_counter() - t0)

        X0 = np.tile(x0, (args.particles, 1))
        t0 = time.perf_counter()
        mc_states = prop_mc.propagate(X0, t_grid, env)
        mc_times.append(time.perf_counter() - t0)

        mc_mean = mc_states.mean(axis=1)
        mc_cov_final = np.cov(mc_states[-1], rowvar=False)

        t0 = time.perf_counter()
        ut_mean, ut_cov = prop_ut.propagate(x0, P0, t_grid, env, augment_process_noise=args.augment_process_noise)
        ut_times.append(time.perf_counter() - t0)

        t0 = time.perf_counter()
        stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)
        stm_times.append(time.perf_counter() - t0)

        ut_mean_errs.append(float(np.linalg.norm(mc_mean[-1] - ut_mean[-1])))
        stm_mean_errs.append(float(np.linalg.norm(mc_mean[-1] - stm_mean[-1])))
        ut_cov_errs.append(float(np.linalg.norm(mc_cov_final - ut_cov[-1]) / np.linalg.norm(mc_cov_final)))
        stm_cov_errs.append(float(np.linalg.norm(mc_cov_final - stm_cov[-1]) / np.linalg.norm(mc_cov_final)))

        print(
            f"Run {run+1}: det={det_times[-1]:.3f}s  mc={mc_times[-1]:.3f}s  "
            f"ut={ut_times[-1]:.3f}s  stm={stm_times[-1]:.3f}s  "
            f"ut_mean_err={ut_mean_errs[-1]:.3e}  stm_mean_err={stm_mean_errs[-1]:.3e}  "
            f"ut_cov_rel_err={ut_cov_errs[-1]:.3e}  stm_cov_rel_err={stm_cov_errs[-1]:.3e}"
        )

    det_mean, det_std = summarize(det_times)
    mc_mean_t, mc_std = summarize(mc_times)
    ut_mean_t, ut_std = summarize(ut_times)
    stm_mean_t, stm_std = summarize(stm_times)

    print("\nTiming summary (mean +/- std)")
    print(f"det={det_mean:.3f} +/- {det_std:.3f}s")
    print(f"mc ={mc_mean_t:.3f} +/- {mc_std:.3f}s")
    print(f"ut ={ut_mean_t:.3f} +/- {ut_std:.3f}s")
    print(f"stm={stm_mean_t:.3f} +/- {stm_std:.3f}s")

    print("\nMC vs UT/STM summary (mean +/- std)")
    ut_mean_m, ut_mean_s = summarize(ut_mean_errs)
    stm_mean_m, stm_mean_s = summarize(stm_mean_errs)
    ut_cov_m, ut_cov_s = summarize(ut_cov_errs)
    stm_cov_m, stm_cov_s = summarize(stm_cov_errs)
    print(f"ut_mean_err={ut_mean_m:.3e} +/- {ut_mean_s:.3e}")
    print(f"stm_mean_err={stm_mean_m:.3e} +/- {stm_mean_s:.3e}")
    print(f"ut_cov_rel_err={ut_cov_m:.3e} +/- {ut_cov_s:.3e}")
    print(f"stm_cov_rel_err={stm_cov_m:.3e} +/- {stm_cov_s:.3e}")


if __name__ == "__main__":
    main()
