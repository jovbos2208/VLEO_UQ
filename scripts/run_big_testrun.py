#!/usr/bin/env python3
"""Multi-core testrun: MC/UT/STM + POD UQ + attitude UQ/filters."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path

import numpy as np


def build_env(
    t_grid: np.ndarray,
    EnvInputs,
    density: float,
    temperature_K: float,
    particle_mass_kg: float,
) -> list:
    env = []
    for _ in t_grid:
        e = EnvInputs()
        e.density = density
        e.temperature_K = temperature_K
        e.particles_mass_kg = particle_mass_kg
        e.wind_I = np.zeros(3)
        env.append(e)
    return env


def synthetic_gnss_positions(
    t_grid: np.ndarray,
    mu: float,
    n_sats: int = 24,
    radius_m: float = 26560e3,
    inc_deg: float = 55.0,
) -> np.ndarray:
    inc = np.deg2rad(inc_deg)
    n = np.sqrt(mu / radius_m ** 3)
    raan_list = np.linspace(0.0, 2.0 * np.pi, n_sats, endpoint=False)
    m0_list = np.linspace(0.0, 2.0 * np.pi, n_sats, endpoint=False)

    def r1(angle: float) -> np.ndarray:
        c = np.cos(angle)
        s = np.sin(angle)
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])

    def r3(angle: float) -> np.ndarray:
        c = np.cos(angle)
        s = np.sin(angle)
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

    rot_i = r1(inc)
    positions = np.zeros((t_grid.size, n_sats, 3))
    for s_idx in range(n_sats):
        rot = r3(raan_list[s_idx]) @ rot_i
        for k, t in enumerate(t_grid):
            m = m0_list[s_idx] + n * t
            r_orb = np.array([radius_m * np.cos(m), radius_m * np.sin(m), 0.0])
            positions[k, s_idx] = rot @ r_orb
    return positions


def build_arcs(duration_s: float, dt_s: float, arc_s: float, overlap_s: float) -> list[tuple[int, int]]:
    arcs = []
    step = max(1, int(round((arc_s - overlap_s) / dt_s)))
    arc_len = max(2, int(round(arc_s / dt_s)))
    n = int(round(duration_s / dt_s)) + 1
    start = 0
    while start < n - 1:
        end = min(n, start + arc_len)
        arcs.append((start, end))
        start += step
        if end == n:
            break
    return arcs


def save_summary(path: Path, data: dict) -> None:
    def to_jsonable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        return obj

    payload = {k: to_jsonable(v) for k, v in data.items()}
    path.write_text(json.dumps(payload, indent=2))


def _timestamp() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")


def banner(title: str) -> None:
    line = "=" * max(8, len(title))
    print(f"\n{line}\n{title}\n{line}")


def status(msg: str) -> None:
    print(f"[{_timestamp()}] [vleo-uq] {msg}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a multi-core VLEO UQ testrun.")
    parser.add_argument("--threads", type=int, default=4, help="OpenMP threads.")
    parser.add_argument("--duration_s", type=float, default=3600.0, help="Orbit run duration.")
    parser.add_argument("--dt_s", type=float, default=30.0, help="Orbit timestep.")
    parser.add_argument("--mc_particles", type=int, default=32, help="MC ensemble size.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--outdir", default=None, help="Output directory.")
    parser.add_argument("--sp3", default=None, help="Optional SP3 file for GNSS.")
    parser.add_argument("--eop", default="data/eop/EOP-All.txt", help="Optional EOP file for SP3.")
    parser.add_argument("--skip_slr", action="store_true", help="Skip SLR simulation.")
    parser.add_argument("--attitude_duration_s", type=float, default=600.0, help="Attitude run duration.")
    parser.add_argument("--attitude_dt_s", type=float, default=1.0, help="Attitude timestep.")
    parser.add_argument("--fullstate_duration_s", type=float, default=600.0, help="Full-state filter duration.")
    parser.add_argument("--fullstate_dt_s", type=float, default=10.0, help="Full-state filter timestep.")
    parser.add_argument("--skip_fullstate_filter", action="store_true", help="Skip full-state filters.")
    parser.add_argument("--fullstate_enkf_members", type=int, default=64, help="EnKF members for full-state filter.")
    parser.add_argument("--fullstate_enkf_inflation", type=float, default=1.0, help="EnKF inflation for full-state filter.")
    parser.add_argument("--verbose", action="store_true", help="Enable propagator progress logs.")
    parser.add_argument("--debug", action="store_true", help="Enable debug state logging.")
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
        boresight_eci_from_states,
        boresight_eci_with_errors,
        dcm_to_quat_wxyz,
        default_geometry,
        default_ilrs_stations,
        geolocation_error_m,
        ground_points_from_los,
        parse_eop_all,
        parse_sp3,
        run_mekf_aero,
        run_ukf_aero,
        run_mekf_fullstate,
        run_ukf_fullstate,
        run_enkf_fullstate,
        simulate_attitude_ou,
        simulate_attitude_truth_aero,
        simulate_gnss_measurements,
        simulate_gyro_measurements,
        simulate_slr_measurements,
        simulate_star_tracker_measurements,
        slice_gnss_measurements,
        slice_slr_measurements,
        slr_availability_mask,
        smear_over_exposure_m,
        run_pod_uq_multi_arc,
        summarize_multi_arc,
    )

    rng = np.random.default_rng(args.seed)
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir or f"results/big_testrun_{timestamp}")
    outdir.mkdir(parents=True, exist_ok=True)

    banner("VLEO UQ Big Testrun")
    status(f"Threads: {args.threads}")
    status(f"Output dir: {outdir}")

    banner("Initialize Models")
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
        config.debug_stride = 10
    config.rho_fast_sigma = 0.1
    config.rho_fast_tau_s = 600.0
    config.rho_bias_sigma = 0.02
    config.rho_bias_tau_s = 3600.0
    config.wind_sigma = 5.0
    config.wind_tau_s = 600.0
    config.rng_seed = args.seed

    prop_det = DeterministicPropagator(aero, vehicle, config)
    prop_mc = EnsemblePropagatorMC(aero, vehicle, config)
    prop_ut = SigmaPointPropagatorUT(aero, vehicle, config)
    prop_stm = StmPropagator(aero, vehicle, config)
    status("Propagators ready (det/MC/UT/STM).")

    banner("Orbit Propagation: MC/UT/STM")
    t_grid = np.arange(0.0, args.duration_s + args.dt_s, args.dt_s)
    env = build_env(
        t_grid,
        EnvInputs,
        density=1e-12,
        temperature_K=1000.0,
        particle_mass_kg=28.0 * 1.6605390689252e-27,
    )
    status(f"Orbit grid: {t_grid.size} points, dt={args.dt_s}s, duration={args.duration_s}s")

    mu = config.mu_earth_m3_s2
    r0 = np.array([7000e3, 0.0, 0.0])
    v0 = np.array([0.0, np.sqrt(mu / np.linalg.norm(r0)), 0.0])
    x0 = np.zeros(prop_det.state_size)
    x0[0:3] = r0
    x0[3:6] = v0
    x0[6] = 1.0
    x0[10:13] = np.array([0.0, 0.0, 0.01])

    P0 = np.diag(
        [10.0 ** 2] * 3
        + [0.01 ** 2] * 3
        + [1e-6] * 4
        + [0.0] * 3
        + [0.05 ** 2]
        + [0.02 ** 2]
        + [1.0 ** 2] * 3
    )

    X0 = np.tile(x0, (args.mc_particles, 1))
    status(f"MC propagation: N={args.mc_particles}")
    mc_states = prop_mc.propagate(X0, t_grid, env)
    mc_mean = mc_states.mean(axis=1)
    mc_cov_final = np.cov(mc_states[-1], rowvar=False)

    status("UT propagation with process-noise augmentation")
    ut_mean, ut_cov = prop_ut.propagate(x0, P0, t_grid, env, augment_process_noise=True)
    status("STM propagation")
    stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)

    ut_mean_err = float(np.linalg.norm(mc_mean[-1] - ut_mean[-1]))
    stm_mean_err = float(np.linalg.norm(mc_mean[-1] - stm_mean[-1]))
    ut_cov_err = float(np.linalg.norm(mc_cov_final - ut_cov[-1]) / np.linalg.norm(mc_cov_final))
    stm_cov_err = float(np.linalg.norm(mc_cov_final - stm_cov[-1]) / np.linalg.norm(mc_cov_final))

    status("Saving MC/UT/STM outputs")
    np.savez(
        outdir / "mc_ut_stm.npz",
        t_grid=t_grid,
        mc_states=mc_states,
        mc_mean=mc_mean,
        mc_cov_final=mc_cov_final,
        ut_mean=ut_mean,
        ut_cov=ut_cov,
        stm_mean=stm_mean,
        stm_cov=stm_cov,
    )

    banner("POD UQ: GNSS/SLR")
    gnss_positions = None
    sat_clock_bias_m = None
    if args.sp3:
        sp3_path = Path(args.sp3)
        if sp3_path.exists():
            status(f"Loading SP3: {sp3_path}")
            eph = parse_sp3(str(sp3_path), constellation_prefix=("G",))
            t0 = eph.epochs[0]
            t_grid_dt = [t0 + dt.timedelta(seconds=float(sec)) for sec in t_grid]
            eop = None
            if args.eop and Path(args.eop).exists():
                status(f"Loading EOP: {args.eop}")
                eop = parse_eop_all(args.eop)
            gnss_positions = eph.sample(t_grid_dt, frame="eci", eop=eop)
            clocks = eph.sample_clock(t_grid_dt)
            if clocks is not None:
                sat_clock_bias_m = clocks * 299792458.0
    if gnss_positions is None:
        status("Using synthetic GNSS constellation")
        gnss_positions = synthetic_gnss_positions(t_grid, mu)

    truth_states = prop_det.propagate(x0, t_grid, env)
    status("Simulating GNSS measurements")
    gnss = simulate_gnss_measurements(
        truth_states,
        t_grid,
        gnss_positions,
        rng=rng,
        frequencies_hz=(1575.42e6, 1227.60e6),
        code_sigma_m=0.5,
        carrier_sigma_m=0.003,
        require_los=True,
        earth_margin_m=50e3,
        fov_half_angle_deg=80.0,
        dropout_prob=0.01,
        cycle_slip_gap_s=30.0,
        cycle_slip_prob_per_min=1e-4,
        sat_clock_bias_m=sat_clock_bias_m,
    )

    slr_sets = None
    if not args.skip_slr:
        status("Simulating SLR measurements")
        stations = default_ilrs_stations()
        slr_mask = slr_availability_mask(
            truth_states,
            t_grid,
            stations,
            min_elevation_deg=20.0,
            require_night=False,
            t0=None,
        )
        slr = simulate_slr_measurements(
            truth_states,
            t_grid,
            stations,
            sigma_m=0.01,
            cadence_s=20.0,
            availability_mask=slr_mask,
            rng=rng,
        )
    else:
        slr = None

    arcs = build_arcs(args.duration_s, args.dt_s, arc_s=1800.0, overlap_s=600.0)
    status(f"POD arcs: {len(arcs)} (arc=1800s, overlap=600s)")
    gnss_sets = []
    slr_sets = []
    for start, end in arcs:
        gnss_sets.append(slice_gnss_measurements(gnss, start, end))
        if slr is not None:
            try:
                slr_sets.append(slice_slr_measurements(slr, start, end))
            except ValueError:
                slr_sets.append(None)
        else:
            slr_sets.append(None)

    x0_guess = x0.copy()
    x0_guess[0:3] += np.array([50.0, -20.0, 10.0])
    P0_guess = np.diag(np.ones(prop_stm.state_size))

    status("Running batch OD per arc")
    pod_results = run_pod_uq_multi_arc(
        prop_det,
        prop_stm,
        x0,
        x0_guess,
        P0_guess,
        t_grid,
        env,
        arcs,
        gnss_sets=gnss_sets,
        slr_sets=None if all(s is None for s in slr_sets) else slr_sets,
        max_iter=6,
    )
    summaries = summarize_multi_arc(pod_results, horizons_s=(1800.0, 3600.0))
    status("Saving POD arc outputs")

    for i, result in enumerate(pod_results):
        np.savez(
            outdir / f"pod_arc_{i+1}.npz",
            t_grid=result.t_grid,
            rtn_error=result.rtn_error,
            rtn_sigma=result.rtn_sigma,
            radial_rms_m=result.radial_rms_m,
        )

    banner("Attitude UQ")
    att_t = np.arange(0.0, args.attitude_duration_s + args.attitude_dt_s, args.attitude_dt_s)
    att_states = np.zeros((att_t.size, 18))
    for i, t in enumerate(att_t):
        theta = np.sqrt(mu / np.linalg.norm(r0) ** 3) * t
        r = np.array([np.linalg.norm(r0) * np.cos(theta), np.linalg.norm(r0) * np.sin(theta), 0.0])
        v = np.array([-np.linalg.norm(v0) * np.sin(theta), np.linalg.norm(v0) * np.cos(theta), 0.0])
        r_hat = r / np.linalg.norm(r)
        h_hat = np.cross(r, v)
        h_hat /= np.linalg.norm(h_hat)
        z_b = -r_hat
        y_b = -h_hat
        x_b = np.cross(y_b, z_b)
        R_BI = np.vstack((x_b, y_b, z_b))
        q = dcm_to_quat_wxyz(R_BI)
        att_states[i, 0:3] = r
        att_states[i, 3:6] = v
        att_states[i, 6:10] = q
    att_err = simulate_attitude_ou(att_t, tau_s=200.0, sigma_rad=1e-4)
    boresight_nom = boresight_eci_from_states(att_states)
    boresight_err = boresight_eci_with_errors(att_states, att_err.errors_rad)
    ground_nom = ground_points_from_los(att_states[:, 0:3], boresight_nom)
    ground_err = ground_points_from_los(att_states[:, 0:3], boresight_err)
    geo_err = geolocation_error_m(ground_nom, ground_err)
    smear = smear_over_exposure_m(att_t, ground_err, exposure_s=1.0)

    status("Saving attitude UQ outputs")
    np.savez(
        outdir / "attitude_uq.npz",
        t_grid=att_t,
        geolocation_error_m=geo_err,
        smear_m=smear,
    )

    banner("Attitude Filters (Aero Torque)")
    att_filt_t = np.arange(0.0, args.attitude_duration_s + args.attitude_dt_s, args.attitude_dt_s)
    env_att = build_env(
        att_filt_t,
        EnvInputs,
        density=1e-12,
        temperature_K=1000.0,
        particle_mass_kg=28.0 * 1.6605390689252e-27,
    )
    v_I = np.zeros((att_filt_t.size, 3))
    v_I[:, 1] = 7500.0
    q0 = np.array([1.0, 0.0, 0.0, 0.0])
    w0 = np.array([0.0, 0.0, 0.01])
    truth_att = simulate_attitude_truth_aero(
        att_filt_t,
        q0,
        w0,
        v_I,
        env_att,
        vehicle,
        aero,
        bias_tau_s=200.0,
        bias_sigma_rad_s=1e-5,
    )
    gyro = simulate_gyro_measurements(truth_att, gyro_sigma_rad_s=1e-4)
    star_idx, star_q = simulate_star_tracker_measurements(truth_att, sigma_rad=5e-4, cadence_s=5.0)

    P0_att = np.diag([1e-6] * 9)
    status("Running MEKF (aero)")
    mekf = run_mekf_aero(
        att_filt_t,
        v_I,
        env_att,
        aero,
        vehicle,
        gyro,
        star_idx,
        star_q,
        q0,
        w0,
        np.zeros(3),
        P0_att,
        gyro_noise_std=1e-4,
        bias_rw_std=1e-6,
        meas_noise_std=5e-4,
    )
    status("Running UKF (aero)")
    ukf = run_ukf_aero(
        att_filt_t,
        v_I,
        env_att,
        aero,
        vehicle,
        gyro,
        star_idx,
        star_q,
        q0,
        w0,
        np.zeros(3),
        P0_att,
        gyro_noise_std=1e-4,
        bias_rw_std=1e-6,
        meas_noise_std=5e-4,
    )

    fullstate_metrics = {}
    if not args.skip_fullstate_filter:
        banner("Full-State Filters (Coupled 18-state)")
        fs_t = np.arange(0.0, args.fullstate_duration_s + args.fullstate_dt_s, args.fullstate_dt_s)
        env_fs = build_env(
            fs_t,
            EnvInputs,
            density=1e-12,
            temperature_K=1000.0,
            particle_mass_kg=28.0 * 1.6605390689252e-27,
        )
        truth_fs = prop_det.propagate(x0, fs_t, env_fs)
        star_idx_fs = np.arange(0, fs_t.size, 5, dtype=int)
        star_meas_fs = truth_fs[star_idx_fs, 6:10]
        gyro_meas_fs = truth_fs[:, 10:13]
        status("Running full-state MEKF")
        mekf_fs = run_mekf_fullstate(
            prop_stm,
            x0,
            P0,
            fs_t,
            env_fs,
            star_idx_fs,
            star_meas_fs,
            star_sigma=1e-6,
            gyro_meas=gyro_meas_fs,
            gyro_sigma=1e-6,
        )
        status("Running full-state UKF")
        ukf_fs = run_ukf_fullstate(
            prop_ut,
            x0,
            P0,
            fs_t,
            env_fs,
            star_idx_fs,
            star_meas_fs,
            star_sigma=1e-6,
            gyro_meas=gyro_meas_fs,
            gyro_sigma=1e-6,
        )
        status("Running full-state EnKF")
        enkf_fs = run_enkf_fullstate(
            prop_mc,
            x0,
            P0,
            fs_t,
            env_fs,
            star_idx_fs,
            star_meas_fs,
            star_sigma=1e-6,
            gyro_meas=gyro_meas_fs,
            gyro_sigma=1e-6,
            members=args.fullstate_enkf_members,
            seed=args.seed + 1000,
            inflation=args.fullstate_enkf_inflation,
        )
        fullstate_metrics = {
            "mekf_final_att_error_norm": float(np.linalg.norm(mekf_fs.states[-1, 6:10] - truth_fs[-1, 6:10])),
            "ukf_final_att_error_norm": float(np.linalg.norm(ukf_fs.states[-1, 6:10] - truth_fs[-1, 6:10])),
            "enkf_final_att_error_norm": float(np.linalg.norm(enkf_fs.states[-1, 6:10] - truth_fs[-1, 6:10])),
        }

    summary = {
        "threads": args.threads,
        "mc_particles": args.mc_particles,
        "fullstate_enkf_members": args.fullstate_enkf_members,
        "fullstate_enkf_inflation": args.fullstate_enkf_inflation,
        "duration_s": args.duration_s,
        "dt_s": args.dt_s,
        "gnss_source": "sp3" if args.sp3 and Path(args.sp3).exists() else "synthetic",
        "ut_mean_err": ut_mean_err,
        "stm_mean_err": stm_mean_err,
        "ut_cov_rel_err": ut_cov_err,
        "stm_cov_rel_err": stm_cov_err,
        "pod_summaries": summaries,
        "attitude_geo_rms_m": float(np.sqrt(np.mean(geo_err ** 2))),
        "attitude_smear_median_m": float(np.median(smear)),
        "attitude_mekf_final_err_norm": float(np.linalg.norm(mekf.q_wxyz[-1] - truth_att.q_wxyz[-1])),
        "attitude_ukf_final_err_norm": float(np.linalg.norm(ukf.q_wxyz[-1] - truth_att.q_wxyz[-1])),
    }
    summary.update(fullstate_metrics)
    status("Writing summary.json")
    save_summary(outdir / "summary.json", summary)
    banner("Done")
    print(f"Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
