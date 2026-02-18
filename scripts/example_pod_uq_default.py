#!/usr/bin/env python3
"""Offline POD UQ example using SP3 ephemeris and default station list."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np

from vleo_uq import (
    AeroAdapter,
    DeterministicPropagator,
    EnvInputs,
    PropagatorConfig,
    StmPropagator,
    VehicleParams,
    default_geometry,
    default_ilrs_stations,
    load_stations_csv,
    load_stations_itrf_csv,
    parse_antex,
    parse_eop_all,
    parse_sp3,
    simulate_gnss_measurements,
    simulate_slr_measurements,
    slr_availability_mask,
    slice_gnss_measurements,
    slice_slr_measurements,
    run_pod_uq_multi_arc,
    summarize_multi_arc,
)


def build_env(t_grid: np.ndarray) -> list[EnvInputs]:
    env = []
    for _ in t_grid:
        e = EnvInputs()
        e.density = 0.0
        e.temperature_K = 1000.0
        e.particles_mass_kg = 28.0 * 1.6605390689252e-27
        e.wind_I = np.zeros(3)
        env.append(e)
    return env


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


def main() -> None:
    parser = argparse.ArgumentParser(description="Offline POD UQ example.")
    parser.add_argument("--sp3", required=True, help="Path to SP3 file.")
    parser.add_argument("--eop", default="data/eop/EOP-All.txt", help="Path to EOP-All.txt.")
    parser.add_argument("--stations", default=None, help="CSV station list (name, lat_deg, lon_deg, alt_m).")
    parser.add_argument("--stations-itrf", default=None,
                        help="ITRF CSV (name,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,epoch_year).")
    parser.add_argument("--duration_s", type=float, default=2 * 3600.0, help="Total duration in seconds.")
    parser.add_argument("--dt_s", type=float, default=30.0, help="Propagation timestep in seconds.")
    parser.add_argument("--arc_s", type=float, default=3600.0, help="Arc length in seconds.")
    parser.add_argument("--overlap_s", type=float, default=600.0, help="Arc overlap in seconds.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--no-los", action="store_true", help="Disable line-of-sight mask.")
    parser.add_argument("--no-fov", action="store_true", help="Disable antenna field-of-view mask.")
    parser.add_argument("--use-sat-clock", action="store_true", help="Add SP3 satellite clock bias.")
    parser.add_argument("--antex", default=None, help="Optional ANTEX file for antenna metadata.")
    parser.add_argument("--skip-slr", action="store_true", help="Skip SLR simulation.")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    eph = parse_sp3(args.sp3, constellation_prefix=("G",))
    eop = None
    if args.eop and Path(args.eop).exists():
        eop = parse_eop_all(args.eop)

    t0 = eph.epochs[0]
    t_grid = np.arange(0.0, args.duration_s + args.dt_s, args.dt_s)
    t_grid_dt = [t0 + dt.timedelta(seconds=float(sec)) for sec in t_grid]

    sat_positions_eci = eph.sample(t_grid_dt, frame="eci", eop=eop)
    sat_clock_bias_m = None
    if args.use_sat_clock:
        clocks = eph.sample_clock(t_grid_dt)
        if clocks is not None:
            sat_clock_bias_m = clocks * 299792458.0

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    vehicle = VehicleParams()
    config = PropagatorConfig()
    config.rho_fast_sigma = 0.0
    config.rho_bias_sigma = 0.0
    config.wind_sigma = 0.0

    prop_det = DeterministicPropagator(aero, vehicle, config)
    prop_stm = StmPropagator(aero, vehicle, config)

    mu = config.mu_earth_m3_s2
    r0 = np.array([7000e3, 0.0, 0.0])
    v0 = np.array([0.0, np.sqrt(mu / np.linalg.norm(r0)), 0.0])

    x0_truth = np.zeros(prop_stm.state_size)
    x0_truth[0:3] = r0
    x0_truth[3:6] = v0
    x0_truth[6] = 1.0

    x0_guess = x0_truth.copy()
    x0_guess[0:3] += np.array([50.0, -20.0, 10.0])

    P0_guess = np.diag(np.ones(prop_stm.state_size))
    env = build_env(t_grid)

    truth_states = prop_det.propagate(x0_truth, t_grid, env)

    gnss = simulate_gnss_measurements(
        truth_states,
        t_grid,
        sat_positions_eci,
        rng=rng,
        frequencies_hz=(1575.42e6, 1227.60e6),
        code_sigma_m=0.5,
        carrier_sigma_m=0.003,
        require_los=not args.no_los,
        earth_margin_m=50e3,
        fov_half_angle_deg=None if args.no_fov else 80.0,
        dropout_prob=0.01,
        cycle_slip_gap_s=30.0,
        cycle_slip_prob_per_min=1e-4,
        sat_clock_bias_m=sat_clock_bias_m,
    )

    if args.stations_itrf:
        stations = load_stations_itrf_csv(args.stations_itrf, t0)
    elif args.stations:
        stations = load_stations_csv(args.stations)
    else:
        stations = default_ilrs_stations()

    slr = None
    if not args.skip_slr:
        slr_mask = slr_availability_mask(
            truth_states,
            t_grid,
            stations,
            min_elevation_deg=20.0,
            require_night=False,
            t0=t0,
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

    arcs = build_arcs(args.duration_s, args.dt_s, args.arc_s, args.overlap_s)
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

    results = run_pod_uq_multi_arc(
        prop_det,
        prop_stm,
        x0_truth,
        x0_guess,
        P0_guess,
        t_grid,
        env,
        arcs,
        gnss_sets=gnss_sets,
        slr_sets=None if all(s is None for s in slr_sets) else slr_sets,
        max_iter=6,
    )

    summaries = summarize_multi_arc(results, horizons_s=(6 * 3600.0, 24 * 3600.0))
    for i, summary in enumerate(summaries):
        print(f"Arc {i+1}: radial_rms_m={summary['radial_rms_m']:.3f} coverage={summary['coverage_3sigma']}")
        for h, sig, err in zip(summary["horizons_s"], summary["rtn_sigma"], summary["rtn_error"]):
            print(f"  t={h:.1f}s sigma={sig} err={err}")

    if args.antex:
        antex = parse_antex(args.antex)
        print(f"Loaded ANTEX entries: {len(antex)} (not applied to measurements)")


if __name__ == "__main__":
    main()
