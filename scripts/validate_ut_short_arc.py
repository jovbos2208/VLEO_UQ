#!/usr/bin/env python3
import argparse
import datetime as dt
from pathlib import Path

import numpy as np

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


def build_env(t_grid, EnvInputs, eta1_rad=None, eta2_rad=None):
    env = []
    for i, _ in enumerate(t_grid):
        e = EnvInputs()
        e.density = 1e-12
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        if eta1_rad is not None:
            e.eta1_rad = float(eta1_rad[i])
        if eta2_rad is not None:
            e.eta2_rad = float(eta2_rad[i])
        env.append(e)
    return env


def wing_profile(t_grid, amp_deg, period_s):
    amp_rad = np.deg2rad(amp_deg)
    omega = 2.0 * np.pi / period_s
    eta1 = amp_rad * np.sin(omega * t_grid)
    eta2 = amp_rad * np.cos(omega * t_grid)
    return eta1, eta2


def build_initial_state(mu):
    r0 = np.array([7000e3, 0.0, 0.0])
    v0 = np.array([0.0, np.sqrt(mu / np.linalg.norm(r0)), 0.0])
    x0 = np.zeros(DeterministicPropagator.state_size)
    x0[0:3] = r0
    x0[3:6] = v0
    x0[6] = 1.0
    x0[10:13] = np.array([0.0, 0.0, 0.01])
    return x0


def sample_mc_states(x0, P0, n, seed, freeze_attitude):
    rng = np.random.default_rng(seed)
    X0 = np.tile(x0, (n, 1))
    if n <= 1:
        return X0
    jitter = 0.0
    P = 0.5 * (P0 + P0.T)
    for _ in range(6):
        try:
            L = np.linalg.cholesky(P + jitter * np.eye(P.shape[0]))
            break
        except np.linalg.LinAlgError:
            diag_max = float(np.max(np.diag(P))) if P.size else 1.0
            jitter = (1e-12 if jitter == 0.0 else jitter * 10.0) * max(1.0, diag_max)
    else:
        raise ValueError("P0 is not positive definite")
    X0 = x0 + rng.standard_normal((n, x0.size)) @ L.T
    q_ref = x0[6:10].copy()
    if np.linalg.norm(q_ref) == 0.0:
        q_ref = np.array([1.0, 0.0, 0.0, 0.0])
    for i in range(n):
        q = X0[i, 6:10]
        norm = np.linalg.norm(q)
        if norm == 0.0:
            q = q_ref.copy()
        else:
            q = q / norm
        if float(np.dot(q, q_ref)) < 0.0:
            q = -q
        X0[i, 6:10] = q
    if freeze_attitude:
        X0[:, 6:10] = q_ref
        X0[:, 10:13] = x0[10:13]
    return X0


def summarize_case(name, mc_states, ut_mean, ut_cov, stm_mean, stm_cov):
    mc_mean = mc_states.mean(axis=1)
    mc_cov = np.cov(mc_states[-1], rowvar=False)
    ut_mean_err = float(np.linalg.norm(mc_mean[-1] - ut_mean[-1]))
    stm_mean_err = float(np.linalg.norm(mc_mean[-1] - stm_mean[-1]))
    ut_cov_err = float(np.linalg.norm(mc_cov - ut_cov[-1]) / np.linalg.norm(mc_cov))
    stm_cov_err = float(np.linalg.norm(mc_cov - stm_cov[-1]) / np.linalg.norm(mc_cov))
    print(f"[validate] {name} ut_mean_err={ut_mean_err:.3e} stm_mean_err={stm_mean_err:.3e} "
          f"ut_cov_rel_err={ut_cov_err:.3e} stm_cov_rel_err={stm_cov_err:.3e}")
    return {
        "case": name,
        "ut_mean_err": ut_mean_err,
        "stm_mean_err": stm_mean_err,
        "ut_cov_rel_err": ut_cov_err,
        "stm_cov_rel_err": stm_cov_err,
    }


def main():
    parser = argparse.ArgumentParser(description="Short-arc UT/MC/STM validation.")
    parser.add_argument("--mc_particles", type=int, default=256)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--ut_alpha", type=float, default=1.0)
    parser.add_argument("--ut_beta", type=float, default=2.0)
    parser.add_argument("--ut_kappa", type=float, default=0.0)
    args = parser.parse_args()

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()

    os_env = __import__("os").environ
    os_env.setdefault("OMP_NUM_THREADS", str(args.threads))

    config_mission = PropagatorConfig()
    config_mission.freeze_attitude = True
    config_att = PropagatorConfig()

    prop_ut_m = SigmaPointPropagatorUT(aero, vehicle, config_mission, args.ut_alpha, args.ut_beta, args.ut_kappa)
    prop_stm_m = StmPropagator(aero, vehicle, config_mission)
    prop_mc_m = EnsemblePropagatorMC(aero, vehicle, config_mission)

    prop_ut_a = SigmaPointPropagatorUT(aero, vehicle, config_att, args.ut_alpha, args.ut_beta, args.ut_kappa)
    prop_stm_a = StmPropagator(aero, vehicle, config_att)
    prop_mc_a = EnsemblePropagatorMC(aero, vehicle, config_att)

    mu = config_mission.mu_earth_m3_s2
    x0 = build_initial_state(mu)

    P0_mission = np.diag(
        [10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [0.0] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3
    )
    P0_att = np.diag(
        [10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [1e-4] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3
    )

    period = 2.0 * np.pi * np.sqrt(np.linalg.norm(x0[0:3]) ** 3 / mu)
    t_grid_m = np.linspace(0.0, 0.25 * period, 81)
    env_m = build_env(t_grid_m, EnvInputs)
    X0_m = sample_mc_states(x0, P0_mission, args.mc_particles, 42, True)
    mc_m = prop_mc_m.propagate(X0_m, t_grid_m, env_m)
    ut_mean_m, ut_cov_m = prop_ut_m.propagate(x0, P0_mission, t_grid_m, env_m)
    stm_mean_m, stm_cov_m = prop_stm_m.propagate(x0, P0_mission, t_grid_m, env_m)
    summary_m = summarize_case("mission_short_orbit", mc_m, ut_mean_m, ut_cov_m, stm_mean_m, stm_cov_m)

    t_grid_a = np.linspace(0.0, 600.0, 121)
    eta1, eta2 = wing_profile(t_grid_a, 10.0, 120.0)
    env_a = build_env(t_grid_a, EnvInputs, eta1, eta2)
    X0_a = sample_mc_states(x0, P0_att, args.mc_particles, 43, False)
    mc_a = prop_mc_a.propagate(X0_a, t_grid_a, env_a)
    ut_mean_a, ut_cov_a = prop_ut_a.propagate(x0, P0_att, t_grid_a, env_a)
    stm_mean_a, stm_cov_a = prop_stm_a.propagate(x0, P0_att, t_grid_a, env_a)
    summary_a = summarize_case("attitude_short", mc_a, ut_mean_a, ut_cov_a, stm_mean_a, stm_cov_a)

    outdir = Path("results") / f"ut_short_arc_{dt.datetime.now(dt.timezone.utc).strftime('%Y%m%d_%H%M%S')}"
    outdir.mkdir(parents=True, exist_ok=True)
    np.savez(outdir / "mission_short_orbit.npz", t_grid=t_grid_m, mc=mc_m, ut_mean=ut_mean_m, ut_cov=ut_cov_m,
             stm_mean=stm_mean_m, stm_cov=stm_cov_m)
    np.savez(outdir / "attitude_short.npz", t_grid=t_grid_a, mc=mc_a, ut_mean=ut_mean_a, ut_cov=ut_cov_a,
             stm_mean=stm_mean_a, stm_cov=stm_cov_a)
    (outdir / "summary.json").write_text(
        __import__("json").dumps({"mission": summary_m, "attitude": summary_a}, indent=2)
    )
    print(f"[validate] wrote {outdir}")


if __name__ == "__main__":
    main()
