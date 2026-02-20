#!/usr/bin/env python3
"""Run mission vs attitude case studies for UQ."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scripts.attitude_uq_modes import summarize_attitude_uq_modes
    from scripts.model_discrepancy import apply_density_model_discrepancy
    from scripts.payload_impact import compute_payload_impact_metrics, payload_metrics_enabled
    from scripts.summary_metadata import build_run_metadata
    from scripts.uq_parameter_channels import apply_uq_parameter_channels
except ModuleNotFoundError:
    from attitude_uq_modes import summarize_attitude_uq_modes
    from model_discrepancy import apply_density_model_discrepancy
    from payload_impact import compute_payload_impact_metrics, payload_metrics_enabled
    from summary_metadata import build_run_metadata
    from uq_parameter_channels import apply_uq_parameter_channels


def build_env(
    t_grid: np.ndarray,
    EnvInputs,
    density: float,
    temperature_K: float,
    particle_mass_kg: float,
    eta1_rad: np.ndarray | None = None,
    eta2_rad: np.ndarray | None = None,
) -> list:
    env = []
    for i, _ in enumerate(t_grid):
        e = EnvInputs()
        e.density = density
        e.temperature_K = temperature_K
        e.particles_mass_kg = particle_mass_kg
        e.wind_I = np.zeros(3)
        e.sun_position_I_m = np.zeros(3)
        e.moon_position_I_m = np.zeros(3)
        e.magnetic_field_I_T = np.zeros(3)
        e.tide_loading_accel_I_m_s2 = np.zeros(3)
        e.srp_scale = 1.0
        e.albedo_ir_scale = 1.0
        if eta1_rad is not None:
            e.eta1_rad = float(eta1_rad[i])
        if eta2_rad is not None:
            e.eta2_rad = float(eta2_rad[i])
        env.append(e)
    return env


def daily_mean(times: np.ndarray, values: np.ndarray):
    days = times.astype("datetime64[D]")
    unique_days, inv = np.unique(days, return_inverse=True)
    sums = np.zeros(len(unique_days), dtype=float)
    counts = np.zeros(len(unique_days), dtype=float)
    for idx, val in zip(inv, values):
        sums[idx] += val
        counts[idx] += 1.0
    means = sums / np.maximum(counts, 1.0)
    return unique_days.astype("datetime64[ns]"), means


def eci_to_geodetic_series(t_grid_s: np.ndarray, t0_utc: dt.datetime, r_eci: np.ndarray):
    from vleo_uq import ecef_to_geodetic
    from vleo_uq.pod_uq import eci_to_ecef
    from vleo_uq.env_sources import _gmst_rad

    lat_deg = np.zeros(len(t_grid_s), dtype=float)
    lon_deg = np.zeros(len(t_grid_s), dtype=float)
    alt_m = np.zeros(len(t_grid_s), dtype=float)
    t_grid_dt = []
    for i, t in enumerate(t_grid_s):
        ts = t0_utc + dt.timedelta(seconds=float(t))
        t_grid_dt.append(np.datetime64(ts))
        gmst = _gmst_rad(ts)
        r_ecef = eci_to_ecef(r_eci[i], gmst)
        lat_rad, lon_rad, alt = ecef_to_geodetic(r_ecef[0], r_ecef[1], r_ecef[2])
        lat_deg[i] = np.rad2deg(lat_rad)
        lon_deg[i] = np.rad2deg(lon_rad)
        alt_m[i] = alt
    return np.array(t_grid_dt, dtype="datetime64[ns]"), lat_deg, lon_deg, alt_m


def build_env_from_sources(
    t_grid_s: np.ndarray,
    t0_utc: dt.datetime,
    x0: np.ndarray,
    prop_det,
    EnvInputs,
    eta1_rad: np.ndarray | None,
    eta2_rad: np.ndarray | None,
    space_weather_source: str,
    omni_path: str,
    hwm14_lib: str,
    hwm14_data: str,
    interpolation: str,
):
    from vleo_uq import (
        DataProvenance,
        EnvSeriesBuilder,
        HWM14SharedLibBackend,
        HWM14WindModel,
        NRLMSIS21DensityModel,
        TimeSeries,
        WindOutputs,
        build_space_weather_dataset,
        parse_omni2_indices,
    )
    from vleo_uq.env_sources import DataSourceError

    env_guess = build_env(
        t_grid_s,
        EnvInputs,
        density=1e-12,
        temperature_K=1000.0,
        particle_mass_kg=28.0 * 1.6605390689252e-27,
        eta1_rad=eta1_rad,
        eta2_rad=eta2_rad,
    )
    det_guess = prop_det.propagate(x0, t_grid_s, env_guess)
    if not np.isfinite(det_guess).all():
        bad_rows = int(np.sum(~np.isfinite(det_guess).all(axis=1)))
        raise DataSourceError(
            f"deterministic pre-propagation produced non-finite states "
            f"(bad_rows={bad_rows}/{det_guess.shape[0]}). "
            "Use a finer scenario dt (ATT cases typically require <=5 s)."
        )
    t_grid_dt, lat_deg, lon_deg, alt_m = eci_to_geodetic_series(
        t_grid_s, t0_utc, det_guess[:, 0:3]
    )
    if (not np.isfinite(lat_deg).all()) or (not np.isfinite(lon_deg).all()) or (not np.isfinite(alt_m).all()):
        bad = int(np.sum(~(np.isfinite(lat_deg) & np.isfinite(lon_deg) & np.isfinite(alt_m))))
        raise DataSourceError(
            f"geodetic pre-processing produced non-finite samples (bad={bad}/{lat_deg.size})."
        )

    if space_weather_source != "omni":
        raise ValueError("only omni source is supported for now")

    prov = DataProvenance(
        name="NASA_OMNI",
        url="https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat",
        retrieved_utc="offline",
    )
    omni = parse_omni2_indices(
        omni_path,
        provenance_f107=prov,
        provenance_ap=prov,
        provenance_kp=prov,
        provenance_dst=prov,
        allow_missing=True,
    )
    f107_days, f107_means = daily_mean(omni["f107"].times, omni["f107"].values)
    day_start = f107_days.min()
    day_end = f107_days.max()
    full_days = np.arange(day_start, day_end + np.timedelta64(1, "D"), np.timedelta64(1, "D"))
    full_vals = np.interp(
        full_days.astype("datetime64[D]").astype("int64"),
        f107_days.astype("datetime64[D]").astype("int64"),
        f107_means,
    )
    f107_daily = TimeSeries(
        times=full_days.astype("datetime64[ns]"),
        values=full_vals,
        provenance=prov,
        units=omni["f107"].units,
    )

    space_weather = build_space_weather_dataset(
        t_grid=t_grid_dt,
        f107_series=f107_daily,
        f107a81_series=None,
        ap_series=omni["ap"],
        kp_series=omni["kp"],
        dst_series=omni["dst"],
        interpolation=interpolation,
        compute_f107a81_flag=True,
    )

    class _ZeroWindModel:
        def evaluate(
            self,
            t_grid,
            lat_deg,
            lon_deg,
            alt_m,
            indices,
        ):
            return WindOutputs(wind_I=np.zeros((len(t_grid), 3), dtype=float))

    strict_hwm14 = os.environ.get("VLEO_HWM14_STRICT", "0").strip() in {"1", "true", "TRUE", "yes", "YES"}
    try:
        backend = HWM14SharedLibBackend(
            lib_path=hwm14_lib,
            data_dir=hwm14_data,
            output_frame="eci",
        )
        wind_model = HWM14WindModel(backend)
    except (DataSourceError, OSError) as exc:
        if strict_hwm14:
            raise
        print(
            f"[env] warning: HWM14 unavailable ({exc}); falling back to zero wind.",
            flush=True,
        )
        wind_model = _ZeroWindModel()
    density_model = NRLMSIS21DensityModel()

    builder = EnvSeriesBuilder(
        density_model=density_model,
        wind_model=wind_model,
        space_weather=space_weather,
    )
    series = builder.build(
        t_grid=t_grid_dt,
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        alt_m=alt_m,
        r_eci=det_guess[:, 0:3],
        interpolation=interpolation,
    )

    if eta1_rad is not None or eta2_rad is not None:
        for i, entry in enumerate(series.env_inputs):
            if eta1_rad is not None:
                entry.eta1_rad = float(eta1_rad[i])
            if eta2_rad is not None:
                entry.eta2_rad = float(eta2_rad[i])

    return series


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


def wing_profile(t_grid: np.ndarray, amp_deg: float, period_s: float) -> tuple[np.ndarray, np.ndarray]:
    if amp_deg == 0.0:
        zeros = np.zeros_like(t_grid)
        return zeros, zeros
    amp = np.deg2rad(amp_deg)
    phase = 2.0 * np.pi * t_grid / period_s
    eta1 = amp * np.sin(phase)
    eta2 = amp * np.cos(phase)
    return eta1, eta2


def wing_profile_axis_constant(t_grid: np.ndarray, wing_deg: float, axis: str) -> tuple[np.ndarray, np.ndarray]:
    axis_key = str(axis).strip().lower()
    d = np.deg2rad(float(wing_deg))
    eta1 = np.zeros_like(t_grid, dtype=float)
    eta2 = np.zeros_like(t_grid, dtype=float)
    if axis_key == "roll":
        eta1[:] = d
        eta2[:] = -d
    elif axis_key == "pitch":
        eta1[:] = d
        eta2[:] = d
    elif axis_key == "yaw":
        eta1[:] = d
        eta2[:] = 0.0
    else:
        raise ValueError(f"unsupported maneuver axis '{axis}'")
    return eta1, eta2


def wing_profile_for_scenario(
    t_grid: np.ndarray,
    scenario: dict | None,
    default_amp_deg: float,
    default_period_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    def as_float(val) -> float | None:
        try:
            f = float(val)
        except Exception:
            return None
        return f if np.isfinite(f) else None

    if not isinstance(scenario, dict):
        return wing_profile(t_grid, default_amp_deg, default_period_s)

    eta1_const = as_float(scenario.get("wing_constant_eta1_deg"))
    eta2_const = as_float(scenario.get("wing_constant_eta2_deg"))
    if eta1_const is not None and eta2_const is not None:
        eta1 = np.full_like(t_grid, np.deg2rad(eta1_const), dtype=float)
        eta2 = np.full_like(t_grid, np.deg2rad(eta2_const), dtype=float)
        return eta1, eta2

    wing_const = as_float(scenario.get("wing_constant_deg"))
    axis = str(scenario.get("maneuver_axis", "")).strip().lower()
    if wing_const is not None and axis in {"roll", "pitch", "yaw"}:
        return wing_profile_axis_constant(t_grid, wing_const, axis)
    if wing_const is not None:
        eta = np.full_like(t_grid, np.deg2rad(wing_const), dtype=float)
        return eta, eta.copy()

    return wing_profile(t_grid, default_amp_deg, default_period_s)


def _normalize_vec(v: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if not np.isfinite(n) or n <= 0.0:
        out = np.array(fallback, dtype=float)
        nf = float(np.linalg.norm(out))
        if not np.isfinite(nf) or nf <= 0.0:
            return np.array([1.0, 0.0, 0.0], dtype=float)
        return out / nf
    return np.array(v, dtype=float) / n


def _quat_normalize_wxyz(q: np.ndarray) -> np.ndarray:
    q = np.array(q, dtype=float).reshape(4)
    n = float(np.linalg.norm(q))
    if not np.isfinite(n) or n <= 0.0:
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    return q / n


def _quat_mul_wxyz(q_a: np.ndarray, q_b: np.ndarray) -> np.ndarray:
    wa, xa, ya, za = _quat_normalize_wxyz(q_a)
    wb, xb, yb, zb = _quat_normalize_wxyz(q_b)
    return np.array(
        [
            wa * wb - xa * xb - ya * yb - za * zb,
            wa * xb + xa * wb + ya * zb - za * yb,
            wa * yb - xa * zb + ya * wb + za * xb,
            wa * zb + xa * yb - ya * xb + za * wb,
        ],
        dtype=float,
    )


def _quat_from_euler_zyx_deg(euler_zyx_deg: np.ndarray | list | tuple) -> np.ndarray:
    # Input order is [yaw_z, pitch_y, roll_x] in degrees.
    vals = np.array(euler_zyx_deg, dtype=float).reshape(3)
    yaw, pitch, roll = np.deg2rad(vals)
    cy = np.cos(0.5 * yaw)
    sy = np.sin(0.5 * yaw)
    cp = np.cos(0.5 * pitch)
    sp = np.sin(0.5 * pitch)
    cr = np.cos(0.5 * roll)
    sr = np.sin(0.5 * roll)
    q = np.array(
        [
            cr * cp * cy + sr * sp * sy,
            sr * cp * cy - cr * sp * sy,
            cr * sp * cy + sr * cp * sy,
            cr * cp * sy - sr * sp * cy,
        ],
        dtype=float,
    )
    return _quat_normalize_wxyz(q)


def _quat_from_dcm(R: np.ndarray) -> np.ndarray:
    R = np.array(R, dtype=float).reshape(3, 3)
    tr = float(np.trace(R))
    if tr > 0.0:
        s = np.sqrt(tr + 1.0) * 2.0
        w = 0.25 * s
        x = (R[2, 1] - R[1, 2]) / s
        y = (R[0, 2] - R[2, 0]) / s
        z = (R[1, 0] - R[0, 1]) / s
    elif (R[0, 0] > R[1, 1]) and (R[0, 0] > R[2, 2]):
        s = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2.0
        w = (R[2, 1] - R[1, 2]) / s
        x = 0.25 * s
        y = (R[0, 1] + R[1, 0]) / s
        z = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2.0
        w = (R[0, 2] - R[2, 0]) / s
        x = (R[0, 1] + R[1, 0]) / s
        y = 0.25 * s
        z = (R[1, 2] + R[2, 1]) / s
    else:
        s = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2.0
        w = (R[1, 0] - R[0, 1]) / s
        x = (R[0, 2] + R[2, 0]) / s
        y = (R[1, 2] + R[2, 1]) / s
        z = 0.25 * s
    return _quat_normalize_wxyz(np.array([w, x, y, z], dtype=float))


def _flow_x_alignment_quaternion(r0_eci_m: np.ndarray, v0_eci_m_s: np.ndarray) -> np.ndarray:
    # Define body +X along relative flow direction at t0: v_rel = wind - v ~= -v.
    bx_I = _normalize_vec(-np.array(v0_eci_m_s, dtype=float), np.array([1.0, 0.0, 0.0], dtype=float))
    h_I = np.cross(np.array(r0_eci_m, dtype=float), np.array(v0_eci_m_s, dtype=float))
    bz_I = _normalize_vec(h_I, np.array([0.0, 0.0, 1.0], dtype=float))
    by_I = _normalize_vec(np.cross(bz_I, bx_I), np.array([0.0, 1.0, 0.0], dtype=float))
    bz_I = _normalize_vec(np.cross(bx_I, by_I), bz_I)
    # R_BI maps inertial vectors to body vectors: v_B = R_BI * v_I.
    R_BI = np.vstack([bx_I, by_I, bz_I])
    return _quat_from_dcm(R_BI)


def scenario_initial_w0_rad_s(scenario: dict | None) -> np.ndarray | None:
    if not isinstance(scenario, dict):
        return None
    w_rad = scenario.get("initial_w_BI_B_rad_s")
    if isinstance(w_rad, (list, tuple)) and len(w_rad) == 3:
        try:
            return np.asarray(w_rad, dtype=float)
        except Exception:
            return None
    w_deg = scenario.get("initial_w_BI_B_deg_s")
    if isinstance(w_deg, (list, tuple)) and len(w_deg) == 3:
        try:
            return np.deg2rad(np.asarray(w_deg, dtype=float))
        except Exception:
            return None
    name = str(scenario.get("name", "")).lower()
    if "att_a1_detumble" in name:
        return np.deg2rad(np.array([8.0, -6.0, 10.0], dtype=float))
    return None


def save_json(path: Path, payload: dict) -> None:
    def to_jsonable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        return obj

    data = {k: to_jsonable(v) for k, v in payload.items()}
    path.write_text(json.dumps(data, indent=2))


def align_quaternions_to_reference(q_series: np.ndarray, ref_series: np.ndarray) -> np.ndarray:
    aligned = q_series.copy()
    for i in range(q_series.shape[0]):
        if np.dot(aligned[i], ref_series[i]) < 0.0:
            aligned[i] *= -1.0
    return aligned


def maybe_import_matplotlib():
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except Exception:
        return None
    return plt


def plot_case(
    outdir: Path,
    name: str,
    t_grid: np.ndarray,
    det_states: np.ndarray,
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    ut_mean: np.ndarray,
    stm_mean: np.ndarray,
    plot_all_states: bool,
    formation_det_sep: np.ndarray | None = None,
    formation_mc_sep: np.ndarray | None = None,
) -> None:
    plt = maybe_import_matplotlib()
    if plt is None:
        print("[case] matplotlib not available; skipping plots.")
        return

    case_dir = outdir / name
    case_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(det_states[:, 0], det_states[:, 1], det_states[:, 2], label="det")
    ax.plot(mc_mean[:, 0], mc_mean[:, 1], mc_mean[:, 2], label="mc_mean")
    ax.plot(ut_mean[:, 0], ut_mean[:, 1], ut_mean[:, 2], label="ut_mean")
    ax.plot(stm_mean[:, 0], stm_mean[:, 1], stm_mean[:, 2], label="stm_mean")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(case_dir / "trajectory_3d.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    labels = ["x", "y", "z"]
    for i, ax in enumerate(axes):
        ax.plot(t_grid, det_states[:, i], label="det")
        ax.plot(t_grid, mc_mean[:, i], label="mc_mean")
        ax.plot(t_grid, ut_mean[:, i], label="ut_mean")
        ax.plot(t_grid, stm_mean[:, i], label="stm_mean")
        ax.fill_between(
            t_grid,
            mc_mean[:, i] - mc_std[:, i],
            mc_mean[:, i] + mc_std[:, i],
            color="C1",
            alpha=0.2,
            label="mc_std" if i == 0 else None,
        )
        ax.set_ylabel(f"r_{labels[i]} [m]")
    axes[-1].set_xlabel("t [s]")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "position_components.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    labels = ["x", "y", "z"]
    for i, ax in enumerate(axes):
        ax.plot(t_grid, det_states[:, 3 + i], label="det")
        ax.plot(t_grid, mc_mean[:, 3 + i], label="mc_mean")
        ax.plot(t_grid, ut_mean[:, 3 + i], label="ut_mean")
        ax.plot(t_grid, stm_mean[:, 3 + i], label="stm_mean")
        ax.fill_between(
            t_grid,
            mc_mean[:, 3 + i] - mc_std[:, 3 + i],
            mc_mean[:, 3 + i] + mc_std[:, 3 + i],
            color="C1",
            alpha=0.2,
            label="mc_std" if i == 0 else None,
        )
        ax.set_ylabel(f"v_{labels[i]} [m/s]")
    axes[-1].set_xlabel("t [s]")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "velocity_components.png", dpi=150)
    plt.close(fig)

    q_ref = det_states[:, 6:10]
    mc_q = align_quaternions_to_reference(mc_mean[:, 6:10], q_ref)
    ut_q = align_quaternions_to_reference(ut_mean[:, 6:10], q_ref)
    stm_q = align_quaternions_to_reference(stm_mean[:, 6:10], q_ref)
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8, 8))
    labels = ["w", "x", "y", "z"]
    for i, ax in enumerate(axes):
        ax.plot(t_grid, q_ref[:, i], label="det")
        ax.plot(t_grid, mc_q[:, i], label="mc_mean")
        ax.plot(t_grid, ut_q[:, i], label="ut_mean")
        ax.plot(t_grid, stm_q[:, i], label="stm_mean")
        ax.set_ylabel(f"q_{labels[i]}")
    axes[-1].set_xlabel("t [s]")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "attitude_quaternion.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    labels = ["x", "y", "z"]
    for i, ax in enumerate(axes):
        ax.plot(t_grid, det_states[:, 10 + i], label="det")
        ax.plot(t_grid, mc_mean[:, 10 + i], label="mc_mean")
        ax.plot(t_grid, ut_mean[:, 10 + i], label="ut_mean")
        ax.plot(t_grid, stm_mean[:, 10 + i], label="stm_mean")
        ax.set_ylabel(f"w_{labels[i]} [rad/s]")
    axes[-1].set_xlabel("t [s]")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "attitude_rates.png", dpi=150)
    plt.close(fig)

    if formation_det_sep is not None and formation_mc_sep is not None:
        fig, ax = plt.subplots(figsize=(8, 4))
        for idx in range(formation_det_sep.shape[0]):
            ax.plot(t_grid, formation_det_sep[idx], label=f"det obj{idx+2}")
        for idx in range(formation_mc_sep.shape[0]):
            ax.plot(t_grid, formation_mc_sep[idx], linestyle="--", label=f"mc obj{idx+2}")
        ax.set_xlabel("t [s]")
        ax.set_ylabel("separation [m]")
        ax.legend()
        fig.tight_layout()
        fig.savefig(case_dir / "formation_separation.png", dpi=150)
        plt.close(fig)

    if plot_all_states:
        fig, axes = plt.subplots(6, 3, sharex=True, figsize=(10, 14))
        axes = axes.ravel()
        for i in range(18):
            ax = axes[i]
            ax.plot(t_grid, det_states[:, i], label="det")
            ax.plot(t_grid, mc_mean[:, i], label="mc_mean")
            ax.plot(t_grid, ut_mean[:, i], label="ut_mean")
            ax.plot(t_grid, stm_mean[:, i], label="stm_mean")
            ax.set_ylabel(f"x[{i}]")
        axes[-1].set_xlabel("t [s]")
        axes[0].legend()
        fig.tight_layout()
        fig.savefig(case_dir / "all_states.png", dpi=150)
        plt.close(fig)


def build_initial_state(
    mu: float,
    prop_state_size: int,
    r_offset: np.ndarray | None = None,
    w0: np.ndarray | None = None,
    orbit: dict | None = None,
    attitude_reference: str | None = None,
    initial_attitude_euler_deg_zyx: list[float] | tuple[float, float, float] | np.ndarray | None = None,
) -> np.ndarray:
    def rot1(angle: float) -> np.ndarray:
        c = np.cos(angle)
        s = np.sin(angle)
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])

    def rot3(angle: float) -> np.ndarray:
        c = np.cos(angle)
        s = np.sin(angle)
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

    def kepler_to_eci(
        mu_val: float,
        a_m: float,
        e: float,
        i_rad: float,
        raan_rad: float,
        argp_rad: float,
        M_rad: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        # Solve Kepler equation for eccentric anomaly E.
        E = float(M_rad)
        for _ in range(20):
            f = E - e * np.sin(E) - M_rad
            fp = 1.0 - e * np.cos(E)
            dE = -f / fp
            E += dE
            if abs(dE) < 1e-12:
                break
        cosE = np.cos(E)
        sinE = np.sin(E)
        fac = np.sqrt(max(0.0, 1.0 - e * e))
        # True anomaly and perifocal state.
        nu = np.arctan2(fac * sinE, cosE - e)
        p = a_m * (1.0 - e * e)
        r_pf = np.array(
            [
                p * np.cos(nu) / (1.0 + e * np.cos(nu)),
                p * np.sin(nu) / (1.0 + e * np.cos(nu)),
                0.0,
            ]
        )
        v_pf = np.array(
            [
                -np.sqrt(mu_val / p) * np.sin(nu),
                np.sqrt(mu_val / p) * (e + np.cos(nu)),
                0.0,
            ]
        )
        Q = rot3(raan_rad) @ rot1(i_rad) @ rot3(argp_rad)
        return Q @ r_pf, Q @ v_pf

    # Default: ISS-like near-circular orbit if nothing is supplied.
    default_orbit = {
        "type": "kepler",
        "alt_km": 420.0,
        "e": 0.0007,
        "i_deg": 51.6,
        "raan_deg": 0.0,
        "argp_deg": 0.0,
        "M_deg": 0.0,
    }
    orb = dict(default_orbit)
    if orbit:
        orb.update(orbit)

    if str(orb.get("type", "kepler")).lower() == "eci_state":
        r0 = np.array(orb.get("r_eci_m", [7000e3, 0.0, 0.0]), dtype=float)
        v0 = np.array(orb.get("v_eci_m_s", [0.0, np.sqrt(mu / np.linalg.norm(r0)), 0.0]), dtype=float)
    else:
        re_m = 6378137.0
        if orb.get("a_km") is not None:
            a_m = 1e3 * float(orb["a_km"])
        else:
            a_m = re_m + 1e3 * float(orb.get("alt_km", 420.0))
        e = float(orb.get("e", 0.0))
        i_rad = np.deg2rad(float(orb.get("i_deg", 51.6)))
        raan_rad = np.deg2rad(float(orb.get("raan_deg", 0.0)))
        argp_rad = np.deg2rad(float(orb.get("argp_deg", 0.0)))
        if orb.get("nu_deg") is not None:
            nu_rad = np.deg2rad(float(orb.get("nu_deg", 0.0)))
            # Convert nu->M for consistent element flow.
            E = 2.0 * np.arctan2(np.sqrt(1.0 - e) * np.sin(0.5 * nu_rad), np.sqrt(1.0 + e) * np.cos(0.5 * nu_rad))
            M_rad = E - e * np.sin(E)
        else:
            M_rad = np.deg2rad(float(orb.get("M_deg", 0.0)))
        r0, v0 = kepler_to_eci(mu, a_m, e, i_rad, raan_rad, argp_rad, M_rad)

    if r_offset is not None:
        r0 = r0 + np.array(r_offset, dtype=float)

    x0 = np.zeros(prop_state_size)
    x0[0:3] = r0
    x0[3:6] = v0
    q_BI = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    ref = str(attitude_reference or "").strip().lower()
    if ref in {"flow_x", "lvlh_flow_x", "att_flow_x"}:
        q_BI = _flow_x_alignment_quaternion(r0, v0)
    if initial_attitude_euler_deg_zyx is not None:
        vals = np.array(initial_attitude_euler_deg_zyx, dtype=float).reshape(-1)
        if vals.size == 3 and np.isfinite(vals).all():
            q_offset = _quat_from_euler_zyx_deg(vals)
            q_BI = _quat_mul_wxyz(q_offset, q_BI)
    x0[6:10] = _quat_normalize_wxyz(q_BI)
    if w0 is None:
        w0 = np.array([0.0, 0.0, 0.01])
    x0[10:13] = np.array(w0, dtype=float)
    return x0


def expand_covariance(P0: np.ndarray, target_size: int) -> np.ndarray:
    if P0.shape[0] != P0.shape[1]:
        raise ValueError("P0 must be square")
    if P0.shape[0] > target_size:
        raise ValueError("P0 larger than target state size")
    if P0.shape[0] == target_size:
        return P0
    out = np.zeros((target_size, target_size), dtype=float)
    out[: P0.shape[0], : P0.shape[1]] = P0
    return out


def sample_mc_initial_states(
    x0: np.ndarray,
    P0: np.ndarray,
    particles: int,
    rng: np.random.Generator,
    freeze_attitude: bool,
) -> np.ndarray:
    state_size = x0.shape[0]
    if particles <= 1:
        return x0[None, :].copy()
    P = expand_covariance(P0, state_size)
    P = 0.5 * (P + P.T)
    jitter = 0.0
    for _ in range(6):
        try:
            L = np.linalg.cholesky(P + jitter * np.eye(state_size))
            break
        except np.linalg.LinAlgError:
            diag_max = float(np.max(np.diag(P))) if P.size else 1.0
            jitter = (1e-12 if jitter == 0.0 else jitter * 10.0) * max(1.0, diag_max)
    else:
        raise ValueError("P0 is not positive definite")
    noise = rng.standard_normal((particles, state_size))
    X0 = x0 + noise @ L.T
    q_ref = x0[6:10].copy()
    if np.linalg.norm(q_ref) == 0.0:
        q_ref = np.array([1.0, 0.0, 0.0, 0.0])
    for i in range(particles):
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


def clone_propagator_config(PropagatorConfig, base_config, overrides: dict[str, Any], freeze_attitude: bool):
    cfg = PropagatorConfig()
    attrs = [
        "mu_earth_m3_s2",
        "use_j2_perturbation",
        "j2_earth",
        "use_j3_perturbation",
        "j3_earth",
        "use_j4_perturbation",
        "j4_earth",
        "gravity_fd_step_m",
        "earth_equatorial_radius_m",
        "use_sun_third_body",
        "use_moon_third_body",
        "mu_sun_m3_s2",
        "mu_moon_m3_s2",
        "sun_ephemeris_scale",
        "moon_ephemeris_scale",
        "use_srp_acceleration",
        "srp_cr",
        "srp_area_m2",
        "solar_pressure_1au_n_m2",
        "astronomical_unit_m",
        "use_albedo_ir_acceleration",
        "albedo_ir_cr",
        "albedo_ir_area_m2",
        "albedo_pressure_n_m2",
        "earth_ir_pressure_n_m2",
        "use_tide_loading_acceleration",
        "tide_loading_scale",
        "use_magnetic_torque",
        "residual_dipole_B_A_m2",
        "residual_dipole_scale",
        "magnetic_field_scale",
        "rtol",
        "atol",
        "min_step_s",
        "max_step_s",
        "step_safety",
        "min_step_factor",
        "max_step_factor",
        "max_steps",
        "rho_fast_tau_s",
        "rho_fast_sigma",
        "rho_bias_tau_s",
        "rho_bias_sigma",
        "wind_tau_s",
        "wind_sigma",
        "rng_seed",
        "verbose",
        "progress_stride",
        "particle_stride",
        "debug_state",
        "debug_stride",
        "freeze_attitude",
    ]
    for key in attrs:
        if hasattr(base_config, key):
            setattr(cfg, key, getattr(base_config, key))
    for key, value in (overrides or {}).items():
        if not hasattr(cfg, key):
            print(f"[case] warning: unknown propagator override '{key}' ignored")
            continue
        setattr(cfg, key, value)
    cfg.freeze_attitude = bool(freeze_attitude)
    return cfg


def parse_events_config(events_raw):
    if not events_raw:
        return []
    from vleo_uq import PropagationEvent

    events = []
    for idx, item in enumerate(events_raw):
        if not isinstance(item, dict):
            raise ValueError(f"event[{idx}] must be an object")
        t_s = float(item["t_s"])
        kind = str(item["kind"])
        params = dict(item.get("params", {}))
        events.append(PropagationEvent(t_s=t_s, kind=kind, params=params))
    return events


def apply_event_to_mean_cov(mean: np.ndarray, cov: np.ndarray, event) -> tuple[np.ndarray, np.ndarray]:
    from vleo_uq import apply_event_to_state

    mean_new = apply_event_to_state(mean, event)
    cov_new = np.array(cov, dtype=float, copy=True)
    if event.kind == "delta_v":
        sigma_dv = event.params.get("sigma_dv_m_s")
        if sigma_dv is not None:
            arr = np.array(sigma_dv, dtype=float)
            if arr.ndim == 0:
                arr = np.full(3, float(arr))
            if arr.shape[0] != 3:
                raise ValueError("sigma_dv_m_s must be scalar or length-3")
            for i in range(3):
                cov_new[3 + i, 3 + i] += float(arr[i]) ** 2
    elif event.kind == "attitude_reset":
        if bool(event.params.get("reset_covariance", True)):
            q_idx = [6, 7, 8, 9]
            cov_new[q_idx, :] = 0.0
            cov_new[:, q_idx] = 0.0
            cov_new[np.ix_(q_idx, q_idx)] = np.eye(4) * 1e-12
            if "w_BI_B" in event.params:
                w_idx = [10, 11, 12]
                cov_new[w_idx, :] = 0.0
                cov_new[:, w_idx] = 0.0
                cov_new[np.ix_(w_idx, w_idx)] = np.eye(3) * 1e-12
    return mean_new, cov_new


def propagate_all_with_events(
    prop_det,
    prop_mc,
    prop_ut,
    prop_stm,
    x0: np.ndarray,
    X0: np.ndarray,
    P0: np.ndarray,
    t_grid: np.ndarray,
    env: list,
    augment_process_noise: bool,
    events: list,
    event_time_tol_s: float,
    rng: np.random.Generator,
    events_sample_uncertainty: bool,
):
    from vleo_uq import (
        apply_event_to_state,
        event_indices,
        propagate_with_events,
        sample_event_realizations,
    )

    nt = len(t_grid)
    if events_sample_uncertainty:
        print("[case] events uncertainty sampling enabled for MC")
        mc_states = np.zeros((nt, X0.shape[0], X0.shape[1]), dtype=float)
        for p in range(X0.shape[0]):
            evp = sample_event_realizations(events, rng)
            tol = max(event_time_tol_s, 0.5 * float(np.min(np.diff(t_grid))) + 1e-9) if nt > 1 else event_time_tol_s
            mc_states[:, p, :] = propagate_with_events(
                prop_det, X0[p], t_grid, env, evp, time_tol_s=tol
            )
    else:
        ev_map = event_indices(t_grid, events, time_tol_s=event_time_tol_s)
        event_idxs = sorted(i for i in ev_map.keys() if 0 < i < nt)
        boundaries = [0] + event_idxs + [nt - 1]

        X_curr = X0.copy()
        if 0 in ev_map:
            for ev in ev_map[0]:
                for p in range(X_curr.shape[0]):
                    X_curr[p] = apply_event_to_state(X_curr[p], ev)
        mc_states = np.zeros((nt, X_curr.shape[0], X_curr.shape[1]), dtype=float)
        mc_states[0] = X_curr
        for bi in range(len(boundaries) - 1):
            s = boundaries[bi]
            e = boundaries[bi + 1]
            seg = prop_mc.propagate(X_curr, t_grid[s : e + 1], env[s : e + 1])
            if s == 0:
                mc_states[s : e + 1] = seg
            else:
                mc_states[s + 1 : e + 1] = seg[1:]
            X_curr = np.array(seg[-1], dtype=float, copy=True)
            if e in ev_map:
                for ev in ev_map[e]:
                    for p in range(X_curr.shape[0]):
                        X_curr[p] = apply_event_to_state(X_curr[p], ev)
                mc_states[e] = X_curr

    # Deterministic path always uses nominal events.
    det_states = propagate_with_events(
        prop_det, x0, t_grid, env, events, time_tol_s=event_time_tol_s
    )

    # UT/STM use deterministic event map at grid times.
    ev_map = event_indices(t_grid, events, time_tol_s=event_time_tol_s)
    event_idxs = sorted(i for i in ev_map.keys() if 0 < i < nt)
    boundaries = [0] + event_idxs + [nt - 1]

    n_state = x0.shape[0]
    ut_mean = np.zeros((nt, n_state), dtype=float)
    ut_cov = np.zeros((nt, n_state, n_state), dtype=float)
    stm_mean = np.zeros((nt, n_state), dtype=float)
    stm_cov = np.zeros((nt, n_state, n_state), dtype=float)

    ut_x = x0.copy()
    ut_P = P0.copy()
    stm_x = x0.copy()
    stm_P = P0.copy()
    if 0 in ev_map:
        for ev in ev_map[0]:
            ut_x, ut_P = apply_event_to_mean_cov(ut_x, ut_P, ev)
            stm_x, stm_P = apply_event_to_mean_cov(stm_x, stm_P, ev)
    ut_mean[0] = ut_x
    ut_cov[0] = ut_P
    stm_mean[0] = stm_x
    stm_cov[0] = stm_P

    for bi in range(len(boundaries) - 1):
        s = boundaries[bi]
        e = boundaries[bi + 1]
        seg_t = t_grid[s : e + 1]
        seg_env = env[s : e + 1]

        seg_ut_mean, seg_ut_cov = prop_ut.propagate(
            ut_x, ut_P, seg_t, seg_env, augment_process_noise=augment_process_noise
        )
        seg_stm_mean, seg_stm_cov = prop_stm.propagate(stm_x, stm_P, seg_t, seg_env)

        if s == 0:
            ut_mean[s : e + 1] = seg_ut_mean
            ut_cov[s : e + 1] = seg_ut_cov
            stm_mean[s : e + 1] = seg_stm_mean
            stm_cov[s : e + 1] = seg_stm_cov
        else:
            ut_mean[s + 1 : e + 1] = seg_ut_mean[1:]
            ut_cov[s + 1 : e + 1] = seg_ut_cov[1:]
            stm_mean[s + 1 : e + 1] = seg_stm_mean[1:]
            stm_cov[s + 1 : e + 1] = seg_stm_cov[1:]

        ut_x = np.array(seg_ut_mean[-1], dtype=float, copy=True)
        ut_P = np.array(seg_ut_cov[-1], dtype=float, copy=True)
        stm_x = np.array(seg_stm_mean[-1], dtype=float, copy=True)
        stm_P = np.array(seg_stm_cov[-1], dtype=float, copy=True)

        if e in ev_map:
            for ev in ev_map[e]:
                ut_x, ut_P = apply_event_to_mean_cov(ut_x, ut_P, ev)
                stm_x, stm_P = apply_event_to_mean_cov(stm_x, stm_P, ev)
            ut_mean[e] = ut_x
            ut_cov[e] = ut_P
            stm_mean[e] = stm_x
            stm_cov[e] = stm_P

    return det_states, mc_states, ut_mean, ut_cov, stm_mean, stm_cov


def run_case(
    name: str,
    outdir: Path,
    duration_s: float,
    dt_s: float,
    particles: int,
    config,
    props,
    EnvInputs,
    eta1_rad: np.ndarray | None,
    eta2_rad: np.ndarray | None,
    P0: np.ndarray,
    augment_process_noise: bool,
    seed: int,
    do_pod: bool = False,
    pod_arc_s: float = 1800.0,
    pod_overlap_s: float = 600.0,
    pod_skip_slr: bool = True,
    pod_sp3: str | None = None,
    pod_eop: str | None = None,
    pod_use_sat_clock: bool = False,
    pod_estimator: str = "batch",
    pod_enkf_members: int = 64,
    pod_enkf_seed: int | None = None,
    pod_enkf_inflation: float = 1.0,
    pod_enkf_use_carrier: bool = False,
    r_offset: np.ndarray | None = None,
    use_env_sources: bool = False,
    start_utc: dt.datetime | None = None,
    space_weather_source: str = "omni",
    omni_path: str = "data/space_weather/omni/omni2_all_years.dat",
    hwm14_lib: str = "data/hwm14/libhwm14.so",
    hwm14_data: str = "data/hwm14",
    interpolation: str = "nearest",
    events: list | None = None,
    event_time_tol_s: float | None = None,
    events_sample_uncertainty: bool = False,
    density_scale: float = 1.0,
    pod_disable_gnss: bool = False,
    pod_gnss_cycle_slip_gap_s: float | None = 30.0,
    pod_gnss_cycle_slip_prob_per_min: float = 1e-4,
    pod_gnss_code_outlier_prob: float = 0.0,
    pod_gnss_code_outlier_sigma_scale: float = 25.0,
    pod_gnss_carrier_outlier_prob: float = 0.0,
    pod_gnss_carrier_outlier_sigma_scale: float = 25.0,
    pod_slr_weather_on: bool = False,
    pod_slr_weather_clear_prob: float = 0.8,
    pod_slr_weather_p_stay_clear: float = 0.985,
    pod_slr_weather_p_stay_blocked: float = 0.93,
    pod_slr_weather_seed: int | None = None,
    pod_measurement_latency_s: float = 0.0,
    pod_measurement_latency_jitter_s: float = 0.0,
    pod_measurement_latency_seed: int | None = None,
    pod_ops_outage_on: bool = False,
    pod_ops_outage_rate_per_hour: float = 0.0,
    pod_ops_outage_mean_duration_s: float = 0.0,
    initial_w0_rad_s: np.ndarray | None = None,
    orbit: dict | None = None,
    attitude_reference: str | None = None,
    initial_attitude_euler_deg_zyx: list[float] | tuple[float, float, float] | np.ndarray | None = None,
    scenario_meta: dict | None = None,
) -> dict:
    prop_det, prop_mc, prop_ut, prop_stm = props
    t_grid = np.arange(0.0, duration_s + dt_s, dt_s)

    mu = config.mu_earth_m3_s2
    if getattr(config, "freeze_attitude", False):
        w0 = np.zeros(3)
    else:
        w0 = initial_w0_rad_s
    x0 = build_initial_state(
        mu,
        prop_det.state_size,
        r_offset=r_offset,
        w0=w0,
        orbit=orbit,
        attitude_reference=attitude_reference,
        initial_attitude_euler_deg_zyx=initial_attitude_euler_deg_zyx,
    )
    if use_env_sources:
        if start_utc is None:
            raise ValueError("start_utc must be provided when using env sources")
        env_series = build_env_from_sources(
            t_grid_s=t_grid,
            t0_utc=start_utc,
            x0=x0,
            prop_det=prop_det,
            EnvInputs=EnvInputs,
            eta1_rad=eta1_rad,
            eta2_rad=eta2_rad,
            space_weather_source=space_weather_source,
            omni_path=omni_path,
            hwm14_lib=hwm14_lib,
            hwm14_data=hwm14_data,
            interpolation=interpolation,
        )
        env = env_series.env_inputs
        if density_scale != 1.0:
            for e in env:
                e.density *= density_scale
    else:
        env = build_env(
            t_grid,
            EnvInputs,
            density=1e-12,
            temperature_K=1000.0,
            particle_mass_kg=28.0 * 1.6605390689252e-27,
            eta1_rad=eta1_rad,
            eta2_rad=eta2_rad,
        )
        if density_scale != 1.0:
            for e in env:
                e.density *= density_scale
    discrepancy_meta = apply_density_model_discrepancy(
        env,
        t_grid,
        scenario=scenario_meta if isinstance(scenario_meta, dict) else {"name": name},
        seed=int(seed),
    )
    rng = np.random.default_rng(seed)
    freeze_attitude = bool(getattr(config, "freeze_attitude", False))
    X0 = sample_mc_initial_states(x0, P0, particles, rng, freeze_attitude)
    parsed_events = parse_events_config(events)
    if parsed_events:
        tol = float(event_time_tol_s) if event_time_tol_s is not None else (0.5 * dt_s + 1e-9)
        print(f"[case] {name} events enabled count={len(parsed_events)} tol={tol:.3f}s")
        print(f"[case] {name} segmented MC/UT/STM/DET start")
        det_states, mc_states, ut_mean, ut_cov, stm_mean, stm_cov = propagate_all_with_events(
            prop_det=prop_det,
            prop_mc=prop_mc,
            prop_ut=prop_ut,
            prop_stm=prop_stm,
            x0=x0,
            X0=X0,
            P0=P0,
            t_grid=t_grid,
            env=env,
            augment_process_noise=augment_process_noise,
            events=parsed_events,
            event_time_tol_s=tol,
            rng=rng,
            events_sample_uncertainty=events_sample_uncertainty,
        )
        print(f"[case] {name} segmented MC/UT/STM/DET done")
    else:
        print(f"[case] {name} MC start particles={particles} steps={len(t_grid)}")
        mc_states = prop_mc.propagate(X0, t_grid, env)
        print(f"[case] {name} MC done")
        print(f"[case] {name} UT start")
        ut_mean, ut_cov = prop_ut.propagate(x0, P0, t_grid, env, augment_process_noise=augment_process_noise)
        print(f"[case] {name} UT done")
        print(f"[case] {name} STM start")
        stm_mean, stm_cov = prop_stm.propagate(x0, P0, t_grid, env)
        print(f"[case] {name} STM done")
        print(f"[case] {name} DET start")
        det_states = prop_det.propagate(x0, t_grid, env)
        print(f"[case] {name} DET done")

    mc_mean = mc_states.mean(axis=1)
    mc_std = mc_states.std(axis=1, ddof=1)
    mc_cov_final = np.cov(mc_states[-1], rowvar=False)

    ut_mean_err = float(np.linalg.norm(mc_mean[-1] - ut_mean[-1]))
    stm_mean_err = float(np.linalg.norm(mc_mean[-1] - stm_mean[-1]))
    ut_cov_err = float(np.linalg.norm(mc_cov_final - ut_cov[-1]) / np.linalg.norm(mc_cov_final))
    stm_cov_err = float(np.linalg.norm(mc_cov_final - stm_cov[-1]) / np.linalg.norm(mc_cov_final))

    case_dir = outdir / name
    case_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        case_dir / "mc_ut_stm.npz",
        t_grid=t_grid,
        det_states=det_states,
        mc_states=mc_states,
        mc_mean=mc_mean,
        mc_std=mc_std,
        mc_cov_final=mc_cov_final,
        ut_mean=ut_mean,
        ut_cov=ut_cov,
        stm_mean=stm_mean,
        stm_cov=stm_cov,
    )
    summary = {
        "case": name,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "particles": particles,
        "augment_process_noise": augment_process_noise,
        "density_scale": density_scale,
        "ut_mean_err": ut_mean_err,
        "stm_mean_err": stm_mean_err,
        "ut_cov_rel_err": ut_cov_err,
        "stm_cov_rel_err": stm_cov_err,
        "seed": seed,
        "events_count": len(parsed_events),
        "model_discrepancy": discrepancy_meta,
        "uq_parameter_draw": (
            dict((scenario_meta or {}).get("uq_parameter_draw", {}))
            if isinstance(scenario_meta, dict)
            else {}
        ),
        "metadata": build_run_metadata(
            scenario_meta if isinstance(scenario_meta, dict) else {"name": name},
            seed,
        ),
    }

    if payload_metrics_enabled(default=True):
        max_samples = int(os.environ.get("VLEO_PAYLOAD_MAX_SAMPLES", "25"))
        summary["payload_impact"] = compute_payload_impact_metrics(
            det_states=det_states,
            mc_states=mc_states,
            t_grid=t_grid,
            max_samples=max_samples,
        )
    if not freeze_attitude:
        summary["attitude_uq_modes"] = summarize_attitude_uq_modes(
            t_grid=t_grid,
            det_states=det_states,
            mc_states=mc_states,
            seed=int(seed),
            fast_tau_s=scenario_meta.get("att_uq_fast_tau_s", 200.0)
            if isinstance(scenario_meta, dict)
            else 200.0,
            fast_sigma_rad=scenario_meta.get("att_uq_fast_sigma_rad", 1e-4)
            if isinstance(scenario_meta, dict)
            else 1e-4,
        )

    if use_env_sources:
        save_json(
            case_dir / "space_weather.json",
            {
                "source": space_weather_source,
                "f107": env_series.space_weather["f107"],
                "f107a81": env_series.space_weather["f107a81"],
                "ap": env_series.space_weather["ap"],
                "kp": env_series.space_weather["kp"],
                "dst": env_series.space_weather["dst"],
            },
        )

    if do_pod:
        try:
            from vleo_uq import (
                default_ilrs_stations,
                parse_eop_all,
                parse_sp3,
                run_pod_uq_multi_arc,
                simulate_gnss_measurements,
                simulate_slr_measurements,
                slr_availability_mask,
                slice_gnss_measurements,
                slice_slr_measurements,
                summarize_multi_arc,
            )

            def build_arcs(total_s: float, dt_s: float, arc_s: float, overlap_s: float):
                arcs = []
                step = max(1, int(round((arc_s - overlap_s) / dt_s)))
                arc_len = max(2, int(round(arc_s / dt_s)))
                n = int(round(total_s / dt_s)) + 1
                start = 0
                while start < n - 1:
                    end = min(n, start + arc_len)
                    arcs.append((start, end))
                    start += step
                    if end == n:
                        break
                return arcs

            rng = np.random.default_rng(seed)
            latency_s = float(
                (scenario_meta or {}).get("pod_measurement_latency_s", pod_measurement_latency_s)
                if isinstance(scenario_meta, dict)
                else pod_measurement_latency_s
            )
            latency_jitter_s = float(
                (scenario_meta or {}).get("pod_measurement_latency_jitter_s", pod_measurement_latency_jitter_s)
                if isinstance(scenario_meta, dict)
                else pod_measurement_latency_jitter_s
            )
            latency_seed = int(
                (scenario_meta or {}).get(
                    "pod_measurement_latency_seed",
                    seed if pod_measurement_latency_seed is None else pod_measurement_latency_seed,
                )
                if isinstance(scenario_meta, dict)
                else (seed if pod_measurement_latency_seed is None else pod_measurement_latency_seed)
            )
            ops_outage_on = bool(
                (scenario_meta or {}).get("pod_ops_outage_on", pod_ops_outage_on)
                if isinstance(scenario_meta, dict)
                else pod_ops_outage_on
            )
            ops_outage_rate_per_hour = float(
                (scenario_meta or {}).get("pod_ops_outage_rate_per_hour", pod_ops_outage_rate_per_hour)
                if isinstance(scenario_meta, dict)
                else pod_ops_outage_rate_per_hour
            )
            ops_outage_mean_duration_s = float(
                (scenario_meta or {}).get(
                    "pod_ops_outage_mean_duration_s",
                    pod_ops_outage_mean_duration_s,
                )
                if isinstance(scenario_meta, dict)
                else pod_ops_outage_mean_duration_s
            )
            gnss_positions = None
            sat_clock_bias_m = None
            if (not pod_disable_gnss) and pod_sp3:
                sp3_path = Path(pod_sp3)
                if sp3_path.exists():
                    eph = parse_sp3(str(sp3_path), constellation_prefix=("G",))
                    t0 = eph.epochs[0]
                    t_grid_dt = [t0 + dt.timedelta(seconds=float(sec)) for sec in t_grid]
                    eop = None
                    if pod_eop:
                        eop_path = Path(pod_eop)
                        if eop_path.exists():
                            eop = parse_eop_all(str(eop_path))
                    gnss_positions = eph.sample(t_grid_dt, frame="eci", eop=eop)
                    if pod_use_sat_clock:
                        clocks = eph.sample_clock(t_grid_dt)
                        if clocks is not None:
                            sat_clock_bias_m = clocks * 299792458.0

            if (not pod_disable_gnss) and gnss_positions is None:
                gnss_positions = synthetic_gnss_positions(t_grid, config.mu_earth_m3_s2)
            gnss = None
            if not pod_disable_gnss:
                gnss = simulate_gnss_measurements(
                    det_states,
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
                    cycle_slip_gap_s=pod_gnss_cycle_slip_gap_s,
                    cycle_slip_prob_per_min=float(pod_gnss_cycle_slip_prob_per_min),
                    code_outlier_prob=float(pod_gnss_code_outlier_prob),
                    code_outlier_sigma_scale=float(pod_gnss_code_outlier_sigma_scale),
                    carrier_outlier_prob=float(pod_gnss_carrier_outlier_prob),
                    carrier_outlier_sigma_scale=float(pod_gnss_carrier_outlier_sigma_scale),
                    sat_clock_bias_m=sat_clock_bias_m,
                    ops_outage_on=ops_outage_on,
                    ops_outage_rate_per_hour=ops_outage_rate_per_hour,
                    ops_outage_mean_duration_s=ops_outage_mean_duration_s,
                )

            slr = None
            slr_mask = None
            if not pod_skip_slr:
                stations = default_ilrs_stations()
                weather_seed = pod_slr_weather_seed
                if weather_seed is None:
                    weather_seed = int(seed) + 7919
                slr_mask = slr_availability_mask(
                    det_states,
                    t_grid,
                    stations,
                    min_elevation_deg=20.0,
                    require_night=False,
                    t0=None,
                    weather_enabled=bool(pod_slr_weather_on),
                    weather_clear_prob=float(pod_slr_weather_clear_prob),
                    weather_p_stay_clear=float(pod_slr_weather_p_stay_clear),
                    weather_p_stay_blocked=float(pod_slr_weather_p_stay_blocked),
                    weather_seed=int(weather_seed),
                )
                slr = simulate_slr_measurements(
                    det_states,
                    t_grid,
                    stations,
                    sigma_m=0.01,
                    cadence_s=20.0,
                    availability_mask=slr_mask,
                    rng=rng,
                    ops_outage_on=ops_outage_on,
                    ops_outage_rate_per_hour=ops_outage_rate_per_hour,
                    ops_outage_mean_duration_s=ops_outage_mean_duration_s,
                )

            arcs = build_arcs(duration_s, dt_s, pod_arc_s, pod_overlap_s)
            gnss_sets = [] if gnss is not None else None
            slr_sets = []
            for start, end in arcs:
                if gnss is not None:
                    gnss_sets.append(slice_gnss_measurements(gnss, start, end))
                if slr is not None:
                    try:
                        slr_sets.append(slice_slr_measurements(slr, start, end))
                    except ValueError:
                        slr_sets.append(None)
                else:
                    slr_sets.append(None)

            if gnss is None and all(s is None for s in slr_sets):
                raise ValueError("POD requires at least one measurement source (GNSS and/or SLR)")

            x0_guess = x0.copy()
            x0_guess[0:3] += np.array([50.0, -20.0, 10.0])
            P0_guess = np.diag(np.ones(prop_stm.state_size))

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
                estimator=str(pod_estimator),
                enkf_members=int(pod_enkf_members),
                enkf_seed=int(seed if pod_enkf_seed is None else pod_enkf_seed),
                enkf_inflation=float(pod_enkf_inflation),
                enkf_use_carrier=bool(pod_enkf_use_carrier),
                measurement_latency_s=latency_s,
                measurement_latency_jitter_s=latency_jitter_s,
                measurement_latency_seed=latency_seed,
            )
            pod_summaries = summarize_multi_arc(pod_results, horizons_s=(1800.0, 3600.0))

            pod_dir = case_dir / "pod"
            pod_dir.mkdir(parents=True, exist_ok=True)
            for i, result in enumerate(pod_results):
                np.savez(
                    pod_dir / f"pod_arc_{i+1}.npz",
                    t_grid=result.t_grid,
                    rtn_error=result.rtn_error,
                    rtn_sigma=result.rtn_sigma,
                    radial_rms_m=result.radial_rms_m,
                )
                try:
                    import matplotlib.pyplot as plt

                    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
                    ax.plot(result.t_grid, result.rtn_error[:, 0], label="R error")
                    ax.plot(result.t_grid, result.rtn_error[:, 1], label="T error")
                    ax.plot(result.t_grid, result.rtn_error[:, 2], label="N error")
                    ax.plot(result.t_grid, 3.0 * result.rtn_sigma[:, 0], "--", label="3σ R")
                    ax.plot(result.t_grid, 3.0 * result.rtn_sigma[:, 1], "--", label="3σ T")
                    ax.plot(result.t_grid, 3.0 * result.rtn_sigma[:, 2], "--", label="3σ N")
                    ax.set_xlabel("t [s]")
                    ax.set_ylabel("RTN [m]")
                    ax.legend()
                    ax.grid(True, alpha=0.3)
                    fig.tight_layout()
                    fig.savefig(pod_dir / f"pod_arc_{i+1}_rtn.png", dpi=150)
                    plt.close(fig)
                except Exception:
                    pass
            save_json(pod_dir / "summary.json", {"pod_summaries": pod_summaries})
            summary["pod"] = {
                "arcs": len(pod_results),
                "arc_s": pod_arc_s,
                "overlap_s": pod_overlap_s,
                "estimator": str(pod_estimator),
                "enkf_members": int(pod_enkf_members),
                "enkf_inflation": float(pod_enkf_inflation),
                "enkf_use_carrier": bool(pod_enkf_use_carrier),
                "gnss_cycle_slip_gap_s": None if pod_gnss_cycle_slip_gap_s is None else float(pod_gnss_cycle_slip_gap_s),
                "gnss_cycle_slip_prob_per_min": float(pod_gnss_cycle_slip_prob_per_min),
                "gnss_code_outlier_prob": float(pod_gnss_code_outlier_prob),
                "gnss_code_outlier_sigma_scale": float(pod_gnss_code_outlier_sigma_scale),
                "gnss_carrier_outlier_prob": float(pod_gnss_carrier_outlier_prob),
                "gnss_carrier_outlier_sigma_scale": float(pod_gnss_carrier_outlier_sigma_scale),
                "slr_weather_on": bool(pod_slr_weather_on),
                "slr_weather_clear_prob": float(pod_slr_weather_clear_prob),
                "slr_weather_p_stay_clear": float(pod_slr_weather_p_stay_clear),
                "slr_weather_p_stay_blocked": float(pod_slr_weather_p_stay_blocked),
                "slr_weather_seed": None if pod_slr_weather_seed is None else int(pod_slr_weather_seed),
                "measurement_latency_s": float(latency_s),
                "measurement_latency_jitter_s": float(latency_jitter_s),
                "measurement_latency_seed": int(latency_seed),
                "ops_outage_on": bool(ops_outage_on),
                "ops_outage_rate_per_hour": float(ops_outage_rate_per_hour),
                "ops_outage_mean_duration_s": float(ops_outage_mean_duration_s),
            }
            if gnss is not None:
                summary["pod"]["gnss_meas_count"] = int(gnss.values.shape[0])
                if gnss.is_outlier is not None:
                    summary["pod"]["gnss_outlier_count"] = int(np.sum(gnss.is_outlier))
                if gnss.is_cycle_slip is not None:
                    summary["pod"]["gnss_cycle_slip_count"] = int(np.sum(gnss.is_cycle_slip))
                summary["pod"]["gnss_ambiguity_total"] = int(gnss.ambiguities_cycles.shape[0])
            if slr_mask is not None:
                summary["pod"]["slr_available_fraction"] = float(np.mean(slr_mask))
        except Exception as exc:
            summary["pod_error"] = str(exc)
    save_json(case_dir / "summary.json", summary)
    print(f"[case] {name} done -> {case_dir}")
    return {
        "t_grid": t_grid,
        "det_states": det_states,
        "mc_mean": mc_mean,
        "mc_std": mc_std,
        "ut_mean": ut_mean,
        "stm_mean": stm_mean,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Mission vs attitude case studies.")
    parser.add_argument("--case", choices=("mission", "attitude", "both"), default="both")
    parser.add_argument("--config", default=None, help="Optional JSON scenario config.")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--outdir", default=None, help="Optional output directory.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--augment_process_noise", action="store_true")
    parser.add_argument("--plot", action="store_true", help="Generate plots for each scenario.")
    parser.add_argument("--plot_all_states", action="store_true", help="Plot all 18 state components.")
    parser.add_argument("--mission_duration_s", type=float, default=3 * 3600.0)
    parser.add_argument("--mission_dt_s", type=float, default=60.0)
    parser.add_argument("--mission_particles", type=int, default=32)
    parser.add_argument("--attitude_duration_s", type=float, default=600.0)
    parser.add_argument("--attitude_dt_s", type=float, default=5.0)
    parser.add_argument("--attitude_particles", type=int, default=32)
    parser.add_argument("--wing_amp_deg", type=float, default=10.0)
    parser.add_argument("--wing_period_s", type=float, default=120.0)
    parser.add_argument("--ut_alpha", type=float, default=1.0)
    parser.add_argument("--ut_beta", type=float, default=2.0)
    parser.add_argument("--ut_kappa", type=float, default=0.0)
    parser.add_argument("--pod_uq", action="store_true", help="Run POD UQ for each scenario.")
    parser.add_argument("--pod_arc_s", type=float, default=1800.0)
    parser.add_argument("--pod_overlap_s", type=float, default=600.0)
    parser.add_argument("--pod_skip_slr", action="store_true", help="Skip SLR simulation for POD.")
    parser.add_argument("--pod_sp3", default=None, help="Optional SP3 file for GNSS ephemerides.")
    parser.add_argument("--pod_eop", default="data/eop/EOP-All.txt", help="Optional EOP file for SP3.")
    parser.add_argument("--pod_use_sat_clock", action="store_true", help="Use satellite clock bias from SP3.")
    parser.add_argument(
        "--pod_estimator",
        choices=("batch", "enkf"),
        default="batch",
        help="POD estimator backend.",
    )
    parser.add_argument("--pod_enkf_members", type=int, default=64, help="EnKF ensemble members for POD.")
    parser.add_argument("--pod_enkf_inflation", type=float, default=1.0, help="EnKF multiplicative inflation.")
    parser.add_argument(
        "--pod_enkf_use_carrier",
        action="store_true",
        help="Include carrier observations in POD EnKF (estimates ambiguities).",
    )
    parser.add_argument("--pod_gnss_cycle_slip_gap_s", type=float, default=30.0)
    parser.add_argument("--pod_gnss_cycle_slip_prob_per_min", type=float, default=1e-4)
    parser.add_argument("--pod_gnss_code_outlier_prob", type=float, default=0.0)
    parser.add_argument("--pod_gnss_code_outlier_sigma_scale", type=float, default=25.0)
    parser.add_argument("--pod_gnss_carrier_outlier_prob", type=float, default=0.0)
    parser.add_argument("--pod_gnss_carrier_outlier_sigma_scale", type=float, default=25.0)
    parser.add_argument("--pod_slr_weather_on", action="store_true")
    parser.add_argument("--pod_slr_weather_clear_prob", type=float, default=0.8)
    parser.add_argument("--pod_slr_weather_p_stay_clear", type=float, default=0.985)
    parser.add_argument("--pod_slr_weather_p_stay_blocked", type=float, default=0.93)
    parser.add_argument("--pod_slr_weather_seed", type=int, default=None)
    parser.add_argument("--use_env_sources", action="store_true", help="Use offline space weather + winds.")
    parser.add_argument("--start_utc", default=None, help="UTC start time, e.g. 2025-01-24T00:00:00.")
    parser.add_argument("--space_weather_source", choices=("omni",), default="omni")
    parser.add_argument("--omni_path", default="data/space_weather/omni/omni2_all_years.dat")
    parser.add_argument("--hwm14_lib", default="data/hwm14/libhwm14.so")
    parser.add_argument("--hwm14_data", default="data/hwm14")
    parser.add_argument("--env_interpolation", default="nearest")
    args = parser.parse_args()

    rank = int(os.environ.get("SLURM_PROCID", os.environ.get("OMPI_COMM_WORLD_RANK", "0")))
    world = int(os.environ.get("SLURM_NTASKS", os.environ.get("OMPI_COMM_WORLD_SIZE", "1")))
    if world > 1 and rank != 0:
        print(f"[case] rank {rank}/{world} idle (run on rank 0 only)")
        return

    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))

    start_utc = None
    if args.start_utc:
        start_utc = dt.datetime.fromisoformat(args.start_utc)

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

    geom = default_geometry()
    aero = AeroAdapter()
    aero.init(geom)
    default_vehicle = VehicleParams()

    def make_config(seed: int, freeze_attitude: bool) -> PropagatorConfig:
        cfg = PropagatorConfig()
        cfg.rtol = 1e-6
        cfg.atol = 1e-6
        cfg.max_step_s = 120.0
        cfg.rng_seed = seed
        cfg.freeze_attitude = freeze_attitude
        return cfg

    config_mission = make_config(args.seed, True)
    config_att = make_config(args.seed, False)

    def vehicle_from_scenario(scenario: dict | None):
        vehicle = VehicleParams()
        if not scenario:
            return vehicle
        mass_kg = scenario.get("spacecraft_mass_kg")
        if mass_kg is not None:
            try:
                mass = float(mass_kg)
                if np.isfinite(mass) and mass > 0.0:
                    vehicle.mass_kg = mass
            except Exception:
                pass
        inertia = scenario.get("spacecraft_inertia_kgm2")
        if isinstance(inertia, (list, tuple)):
            try:
                arr = np.asarray(inertia, dtype=float)
            except Exception:
                arr = np.array([], dtype=float)
            if arr.size == 3:
                vehicle.inertia_B = np.diag(arr)
            elif arr.size == 9:
                vehicle.inertia_B = arr.reshape(3, 3)
        return vehicle

    def build_props_for_config(cfg, ut_alpha: float, ut_beta: float, ut_kappa: float, vehicle=None):
        vehicle_use = vehicle if vehicle is not None else default_vehicle
        return (
            DeterministicPropagator(aero, vehicle_use, cfg),
            EnsemblePropagatorMC(aero, vehicle_use, cfg),
            SigmaPointPropagatorUT(aero, vehicle_use, cfg, ut_alpha, ut_beta, ut_kappa),
            StmPropagator(aero, vehicle_use, cfg),
        )

    props_mission = build_props_for_config(config_mission, args.ut_alpha, args.ut_beta, args.ut_kappa)
    props_att = build_props_for_config(config_att, args.ut_alpha, args.ut_beta, args.ut_kappa)

    outdir = (
        Path(args.outdir)
        if args.outdir
        else Path("results") / f"case_studies_{dt.datetime.now(dt.timezone.utc).strftime('%Y%m%d_%H%M%S')}"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    # P0 for mission: attitude-rate covariance set to zero.
    P0_mission = np.diag(
        [10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [0.0] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3
    )
    # P0 for attitude case: allow attitude-rate uncertainty.
    P0_att = np.diag(
        [10.0 ** 2] * 3 + [0.01 ** 2] * 3 + [1e-6] * 4 + [1e-4] * 3 + [0.05 ** 2] * 2 + [1.0 ** 2] * 3
    )

    def resolve_start_utc(value: str | None) -> dt.datetime | None:
        if value is None:
            return start_utc
        return dt.datetime.fromisoformat(value)

    if args.config:
        raw = json.loads(Path(args.config).read_text(encoding="utf-8"))
        scenarios = raw.get("scenarios", [])
        if not scenarios:
            raise ValueError("config must include scenarios")

        def scenario_cfg_and_props(scenario: dict, seed_value: int, freeze_attitude: bool):
            scenario_eff, _ = apply_uq_parameter_channels(scenario, seed=int(seed_value))
            base = config_att if freeze_attitude is False else config_mission
            overrides = scenario_eff.get("propagator_overrides", {})
            cfg = clone_propagator_config(
                PropagatorConfig=PropagatorConfig,
                base_config=base,
                overrides=overrides,
                freeze_attitude=freeze_attitude,
            )
            cfg.rng_seed = int(seed_value)
            ut_alpha = float(scenario_eff.get("ut_alpha", args.ut_alpha))
            ut_beta = float(scenario_eff.get("ut_beta", args.ut_beta))
            ut_kappa = float(scenario_eff.get("ut_kappa", args.ut_kappa))
            vehicle = vehicle_from_scenario(scenario_eff)
            return cfg, build_props_for_config(cfg, ut_alpha, ut_beta, ut_kappa, vehicle=vehicle), scenario_eff

        for scenario in scenarios:
            scenario_base = scenario
            name = scenario["name"]
            typ = scenario["type"]
            duration_s = float(scenario.get("duration_s", 3600.0))
            dt_s = float(scenario.get("dt_s", 60.0))
            particles = int(scenario.get("particles", 32))
            augment = bool(scenario.get("augment_process_noise", False))
            if typ == "mission":
                cfg_run, props_run, scenario = scenario_cfg_and_props(scenario_base, args.seed, True)
                res = run_case(
                    name,
                    outdir,
                    duration_s,
                    dt_s,
                    particles,
                    cfg_run,
                    props_run,
                    EnvInputs,
                    eta1_rad=None,
                    eta2_rad=None,
                    P0=P0_mission,
                    augment_process_noise=augment,
                    seed=args.seed,
                    do_pod=bool(scenario.get("pod_uq", args.pod_uq)),
                    pod_arc_s=float(scenario.get("pod_arc_s", args.pod_arc_s)),
                    pod_overlap_s=float(scenario.get("pod_overlap_s", args.pod_overlap_s)),
                    pod_skip_slr=bool(scenario.get("pod_skip_slr", args.pod_skip_slr)),
                    pod_sp3=scenario.get("pod_sp3", args.pod_sp3),
                    pod_eop=scenario.get("pod_eop", args.pod_eop),
                    pod_use_sat_clock=bool(scenario.get("pod_use_sat_clock", args.pod_use_sat_clock)),
                    pod_estimator=str(scenario.get("pod_estimator", args.pod_estimator)),
                    pod_enkf_members=int(scenario.get("pod_enkf_members", args.pod_enkf_members)),
                    pod_enkf_seed=scenario.get("pod_enkf_seed"),
                    pod_enkf_inflation=float(scenario.get("pod_enkf_inflation", args.pod_enkf_inflation)),
                    pod_enkf_use_carrier=bool(
                        scenario.get("pod_enkf_use_carrier", args.pod_enkf_use_carrier)
                    ),
                    use_env_sources=bool(scenario.get("use_env_sources", args.use_env_sources)),
                    start_utc=resolve_start_utc(scenario.get("start_utc")),
                    space_weather_source=scenario.get("space_weather_source", args.space_weather_source),
                    omni_path=scenario.get("omni_path", args.omni_path),
                    hwm14_lib=scenario.get("hwm14_lib", args.hwm14_lib),
                    hwm14_data=scenario.get("hwm14_data", args.hwm14_data),
                    interpolation=scenario.get("env_interpolation", args.env_interpolation),
                    events=scenario.get("events"),
                    event_time_tol_s=scenario.get("event_time_tol_s"),
                    events_sample_uncertainty=bool(scenario.get("events_sample_uncertainty", False)),
                    density_scale=float(scenario.get("density_scale", 1.0)),
                    pod_disable_gnss=bool(scenario.get("pod_disable_gnss", False)),
                    pod_gnss_cycle_slip_gap_s=scenario.get(
                        "pod_gnss_cycle_slip_gap_s", args.pod_gnss_cycle_slip_gap_s
                    ),
                    pod_gnss_cycle_slip_prob_per_min=float(
                        scenario.get("pod_gnss_cycle_slip_prob_per_min", args.pod_gnss_cycle_slip_prob_per_min)
                    ),
                    pod_gnss_code_outlier_prob=float(
                        scenario.get("pod_gnss_code_outlier_prob", args.pod_gnss_code_outlier_prob)
                    ),
                    pod_gnss_code_outlier_sigma_scale=float(
                        scenario.get(
                            "pod_gnss_code_outlier_sigma_scale",
                            args.pod_gnss_code_outlier_sigma_scale,
                        )
                    ),
                    pod_gnss_carrier_outlier_prob=float(
                        scenario.get("pod_gnss_carrier_outlier_prob", args.pod_gnss_carrier_outlier_prob)
                    ),
                    pod_gnss_carrier_outlier_sigma_scale=float(
                        scenario.get(
                            "pod_gnss_carrier_outlier_sigma_scale",
                            args.pod_gnss_carrier_outlier_sigma_scale,
                        )
                    ),
                    pod_slr_weather_on=bool(scenario.get("pod_slr_weather_on", args.pod_slr_weather_on)),
                    pod_slr_weather_clear_prob=float(
                        scenario.get("pod_slr_weather_clear_prob", args.pod_slr_weather_clear_prob)
                    ),
                    pod_slr_weather_p_stay_clear=float(
                        scenario.get("pod_slr_weather_p_stay_clear", args.pod_slr_weather_p_stay_clear)
                    ),
                    pod_slr_weather_p_stay_blocked=float(
                        scenario.get("pod_slr_weather_p_stay_blocked", args.pod_slr_weather_p_stay_blocked)
                    ),
                    pod_slr_weather_seed=scenario.get("pod_slr_weather_seed", args.pod_slr_weather_seed),
                    initial_w0_rad_s=scenario_initial_w0_rad_s(scenario),
                    orbit=scenario.get("orbit"),
                    attitude_reference=str(scenario.get("attitude_reference", "flow_x")),
                    initial_attitude_euler_deg_zyx=scenario.get("initial_attitude_euler_deg_zyx"),
                    scenario_meta=scenario,
                )
                if args.plot:
                    plot_case(
                        outdir,
                        name,
                        res["t_grid"],
                        res["det_states"],
                        res["mc_mean"],
                        res["mc_std"],
                        res["ut_mean"],
                        res["stm_mean"],
                        args.plot_all_states,
                    )
            elif typ == "formation":
                offsets = scenario.get("formation_offsets_m", [[0.0, 0.0, 0.0]])
                results = []
                for idx, offset in enumerate(offsets):
                    cfg_run, props_run, scenario = scenario_cfg_and_props(scenario_base, args.seed + idx, True)
                    res = run_case(
                        f"{name}/obj_{idx+1:02d}",
                        outdir,
                        duration_s,
                        dt_s,
                        particles,
                        cfg_run,
                        props_run,
                        EnvInputs,
                        eta1_rad=None,
                        eta2_rad=None,
                        P0=P0_mission,
                        augment_process_noise=augment,
                        seed=args.seed + idx,
                        r_offset=np.array(offset, dtype=float),
                        do_pod=bool(scenario.get("pod_uq", args.pod_uq)),
                        pod_arc_s=float(scenario.get("pod_arc_s", args.pod_arc_s)),
                        pod_overlap_s=float(scenario.get("pod_overlap_s", args.pod_overlap_s)),
                        pod_skip_slr=bool(scenario.get("pod_skip_slr", args.pod_skip_slr)),
                        pod_sp3=scenario.get("pod_sp3", args.pod_sp3),
                        pod_eop=scenario.get("pod_eop", args.pod_eop),
                        pod_use_sat_clock=bool(scenario.get("pod_use_sat_clock", args.pod_use_sat_clock)),
                        pod_estimator=str(scenario.get("pod_estimator", args.pod_estimator)),
                        pod_enkf_members=int(scenario.get("pod_enkf_members", args.pod_enkf_members)),
                        pod_enkf_seed=scenario.get("pod_enkf_seed"),
                        pod_enkf_inflation=float(scenario.get("pod_enkf_inflation", args.pod_enkf_inflation)),
                        pod_enkf_use_carrier=bool(
                            scenario.get("pod_enkf_use_carrier", args.pod_enkf_use_carrier)
                        ),
                        use_env_sources=bool(scenario.get("use_env_sources", args.use_env_sources)),
                        start_utc=resolve_start_utc(scenario.get("start_utc")),
                        space_weather_source=scenario.get("space_weather_source", args.space_weather_source),
                        omni_path=scenario.get("omni_path", args.omni_path),
                        hwm14_lib=scenario.get("hwm14_lib", args.hwm14_lib),
                        hwm14_data=scenario.get("hwm14_data", args.hwm14_data),
                        interpolation=scenario.get("env_interpolation", args.env_interpolation),
                        events=scenario.get("events"),
                        event_time_tol_s=scenario.get("event_time_tol_s"),
                        events_sample_uncertainty=bool(scenario.get("events_sample_uncertainty", False)),
                        density_scale=float(scenario.get("density_scale", 1.0)),
                        pod_disable_gnss=bool(scenario.get("pod_disable_gnss", False)),
                        pod_gnss_cycle_slip_gap_s=scenario.get(
                            "pod_gnss_cycle_slip_gap_s", args.pod_gnss_cycle_slip_gap_s
                        ),
                        pod_gnss_cycle_slip_prob_per_min=float(
                            scenario.get("pod_gnss_cycle_slip_prob_per_min", args.pod_gnss_cycle_slip_prob_per_min)
                        ),
                        pod_gnss_code_outlier_prob=float(
                            scenario.get("pod_gnss_code_outlier_prob", args.pod_gnss_code_outlier_prob)
                        ),
                        pod_gnss_code_outlier_sigma_scale=float(
                            scenario.get(
                                "pod_gnss_code_outlier_sigma_scale",
                                args.pod_gnss_code_outlier_sigma_scale,
                            )
                        ),
                        pod_gnss_carrier_outlier_prob=float(
                            scenario.get("pod_gnss_carrier_outlier_prob", args.pod_gnss_carrier_outlier_prob)
                        ),
                        pod_gnss_carrier_outlier_sigma_scale=float(
                            scenario.get(
                                "pod_gnss_carrier_outlier_sigma_scale",
                                args.pod_gnss_carrier_outlier_sigma_scale,
                            )
                        ),
                        pod_slr_weather_on=bool(scenario.get("pod_slr_weather_on", args.pod_slr_weather_on)),
                        pod_slr_weather_clear_prob=float(
                            scenario.get("pod_slr_weather_clear_prob", args.pod_slr_weather_clear_prob)
                        ),
                        pod_slr_weather_p_stay_clear=float(
                            scenario.get("pod_slr_weather_p_stay_clear", args.pod_slr_weather_p_stay_clear)
                        ),
                        pod_slr_weather_p_stay_blocked=float(
                            scenario.get("pod_slr_weather_p_stay_blocked", args.pod_slr_weather_p_stay_blocked)
                        ),
                        pod_slr_weather_seed=scenario.get("pod_slr_weather_seed", args.pod_slr_weather_seed),
                        initial_w0_rad_s=scenario_initial_w0_rad_s(scenario),
                        orbit=scenario.get("orbit"),
                        scenario_meta=scenario,
                    )
                    results.append(res)
                    if args.plot:
                        plot_case(
                            outdir,
                            f"{name}/obj_{idx+1:02d}",
                            res["t_grid"],
                            res["det_states"],
                            res["mc_mean"],
                            res["mc_std"],
                            res["ut_mean"],
                            res["stm_mean"],
                            args.plot_all_states,
                        )
                if results:
                    chief = results[0]
                    t_grid = chief["t_grid"]
                    det_sep = []
                    mc_sep = []
                    for follower in results[1:]:
                        det = np.linalg.norm(follower["det_states"][:, 0:3] - chief["det_states"][:, 0:3], axis=1)
                        mc = np.linalg.norm(follower["mc_mean"][:, 0:3] - chief["mc_mean"][:, 0:3], axis=1)
                        det_sep.append(det)
                        mc_sep.append(mc)
                    if det_sep:
                        det_sep = np.vstack(det_sep)
                        mc_sep = np.vstack(mc_sep)
                        np.savez(
                            outdir / name / "formation_metrics.npz",
                            t_grid=t_grid,
                            det_separation_m=det_sep,
                            mc_mean_separation_m=mc_sep,
                        )
                        if args.plot:
                            plot_case(
                                outdir,
                                name,
                                t_grid,
                                chief["det_states"],
                                chief["mc_mean"],
                                chief["mc_std"],
                                chief["ut_mean"],
                                chief["stm_mean"],
                                args.plot_all_states,
                                formation_det_sep=det_sep,
                                formation_mc_sep=mc_sep,
                            )
                        save_json(
                            outdir / name / "summary.json",
                            {
                                "case": name,
                                "type": "formation",
                                "num_objects": len(results),
                                "duration_s": duration_s,
                                "dt_s": dt_s,
                                "particles": particles,
                                "metadata": build_run_metadata(scenario_base, args.seed),
                            },
                        )
            elif typ == "attitude":
                t_grid = np.arange(0.0, duration_s + dt_s, dt_s)
                eta1_rad, eta2_rad = wing_profile_for_scenario(
                    t_grid,
                    scenario,
                    default_amp_deg=float(scenario.get("wing_amp_deg", args.wing_amp_deg)),
                    default_period_s=float(scenario.get("wing_period_s", args.wing_period_s)),
                )
                cfg_run, props_run, scenario = scenario_cfg_and_props(scenario_base, args.seed, False)
                res = run_case(
                    name,
                    outdir,
                    duration_s,
                    dt_s,
                    particles,
                    cfg_run,
                    props_run,
                    EnvInputs,
                    eta1_rad=eta1_rad,
                    eta2_rad=eta2_rad,
                    P0=P0_att,
                    augment_process_noise=augment,
                    seed=args.seed,
                    do_pod=bool(scenario.get("pod_uq", args.pod_uq)),
                    pod_arc_s=float(scenario.get("pod_arc_s", args.pod_arc_s)),
                    pod_overlap_s=float(scenario.get("pod_overlap_s", args.pod_overlap_s)),
                    pod_skip_slr=bool(scenario.get("pod_skip_slr", args.pod_skip_slr)),
                    pod_sp3=scenario.get("pod_sp3", args.pod_sp3),
                    pod_eop=scenario.get("pod_eop", args.pod_eop),
                    pod_use_sat_clock=bool(scenario.get("pod_use_sat_clock", args.pod_use_sat_clock)),
                    pod_estimator=str(scenario.get("pod_estimator", args.pod_estimator)),
                    pod_enkf_members=int(scenario.get("pod_enkf_members", args.pod_enkf_members)),
                    pod_enkf_seed=scenario.get("pod_enkf_seed"),
                    pod_enkf_inflation=float(scenario.get("pod_enkf_inflation", args.pod_enkf_inflation)),
                    pod_enkf_use_carrier=bool(
                        scenario.get("pod_enkf_use_carrier", args.pod_enkf_use_carrier)
                    ),
                    use_env_sources=bool(scenario.get("use_env_sources", args.use_env_sources)),
                    start_utc=resolve_start_utc(scenario.get("start_utc")),
                    space_weather_source=scenario.get("space_weather_source", args.space_weather_source),
                    omni_path=scenario.get("omni_path", args.omni_path),
                    hwm14_lib=scenario.get("hwm14_lib", args.hwm14_lib),
                    hwm14_data=scenario.get("hwm14_data", args.hwm14_data),
                    interpolation=scenario.get("env_interpolation", args.env_interpolation),
                    events=scenario.get("events"),
                    event_time_tol_s=scenario.get("event_time_tol_s"),
                    events_sample_uncertainty=bool(scenario.get("events_sample_uncertainty", False)),
                    density_scale=float(scenario.get("density_scale", 1.0)),
                    pod_disable_gnss=bool(scenario.get("pod_disable_gnss", False)),
                    pod_gnss_cycle_slip_gap_s=scenario.get(
                        "pod_gnss_cycle_slip_gap_s", args.pod_gnss_cycle_slip_gap_s
                    ),
                    pod_gnss_cycle_slip_prob_per_min=float(
                        scenario.get("pod_gnss_cycle_slip_prob_per_min", args.pod_gnss_cycle_slip_prob_per_min)
                    ),
                    pod_gnss_code_outlier_prob=float(
                        scenario.get("pod_gnss_code_outlier_prob", args.pod_gnss_code_outlier_prob)
                    ),
                    pod_gnss_code_outlier_sigma_scale=float(
                        scenario.get(
                            "pod_gnss_code_outlier_sigma_scale",
                            args.pod_gnss_code_outlier_sigma_scale,
                        )
                    ),
                    pod_gnss_carrier_outlier_prob=float(
                        scenario.get("pod_gnss_carrier_outlier_prob", args.pod_gnss_carrier_outlier_prob)
                    ),
                    pod_gnss_carrier_outlier_sigma_scale=float(
                        scenario.get(
                            "pod_gnss_carrier_outlier_sigma_scale",
                            args.pod_gnss_carrier_outlier_sigma_scale,
                        )
                    ),
                    pod_slr_weather_on=bool(scenario.get("pod_slr_weather_on", args.pod_slr_weather_on)),
                    pod_slr_weather_clear_prob=float(
                        scenario.get("pod_slr_weather_clear_prob", args.pod_slr_weather_clear_prob)
                    ),
                    pod_slr_weather_p_stay_clear=float(
                        scenario.get("pod_slr_weather_p_stay_clear", args.pod_slr_weather_p_stay_clear)
                    ),
                    pod_slr_weather_p_stay_blocked=float(
                        scenario.get("pod_slr_weather_p_stay_blocked", args.pod_slr_weather_p_stay_blocked)
                    ),
                    pod_slr_weather_seed=scenario.get("pod_slr_weather_seed", args.pod_slr_weather_seed),
                    initial_w0_rad_s=scenario_initial_w0_rad_s(scenario),
                    orbit=scenario.get("orbit"),
                    scenario_meta=scenario,
                )
                if args.plot:
                    plot_case(
                        outdir,
                        name,
                        res["t_grid"],
                        res["det_states"],
                        res["mc_mean"],
                        res["mc_std"],
                        res["ut_mean"],
                        res["stm_mean"],
                        args.plot_all_states,
                    )
            else:
                raise ValueError(f"unknown scenario type: {typ}")
        return

    if args.case in {"mission", "both"}:
        res = run_case(
            "mission",
            outdir,
            args.mission_duration_s,
            args.mission_dt_s,
            args.mission_particles,
            config_mission,
            props_mission,
            EnvInputs,
            eta1_rad=None,
            eta2_rad=None,
            P0=P0_mission,
            augment_process_noise=args.augment_process_noise,
            seed=args.seed,
            do_pod=args.pod_uq,
            pod_arc_s=args.pod_arc_s,
            pod_overlap_s=args.pod_overlap_s,
            pod_skip_slr=args.pod_skip_slr,
            pod_sp3=args.pod_sp3,
            pod_eop=args.pod_eop,
            pod_use_sat_clock=args.pod_use_sat_clock,
            pod_estimator=args.pod_estimator,
            pod_enkf_members=args.pod_enkf_members,
            pod_enkf_seed=None,
            pod_enkf_inflation=args.pod_enkf_inflation,
            pod_enkf_use_carrier=args.pod_enkf_use_carrier,
            pod_gnss_cycle_slip_gap_s=args.pod_gnss_cycle_slip_gap_s,
            pod_gnss_cycle_slip_prob_per_min=args.pod_gnss_cycle_slip_prob_per_min,
            pod_gnss_code_outlier_prob=args.pod_gnss_code_outlier_prob,
            pod_gnss_code_outlier_sigma_scale=args.pod_gnss_code_outlier_sigma_scale,
            pod_gnss_carrier_outlier_prob=args.pod_gnss_carrier_outlier_prob,
            pod_gnss_carrier_outlier_sigma_scale=args.pod_gnss_carrier_outlier_sigma_scale,
            pod_slr_weather_on=args.pod_slr_weather_on,
            pod_slr_weather_clear_prob=args.pod_slr_weather_clear_prob,
            pod_slr_weather_p_stay_clear=args.pod_slr_weather_p_stay_clear,
            pod_slr_weather_p_stay_blocked=args.pod_slr_weather_p_stay_blocked,
            pod_slr_weather_seed=args.pod_slr_weather_seed,
            use_env_sources=args.use_env_sources,
            start_utc=start_utc,
            space_weather_source=args.space_weather_source,
            omni_path=args.omni_path,
            hwm14_lib=args.hwm14_lib,
            hwm14_data=args.hwm14_data,
            interpolation=args.env_interpolation,
            scenario_meta={"name": "mission", "catalog_scenario_id": "MISSION"},
        )
        if args.plot:
            plot_case(
                outdir,
                "mission",
                res["t_grid"],
                res["det_states"],
                res["mc_mean"],
                res["mc_std"],
                res["ut_mean"],
                res["stm_mean"],
                args.plot_all_states,
            )

    if args.case in {"attitude", "both"}:
        t_grid = np.arange(0.0, args.attitude_duration_s + args.attitude_dt_s, args.attitude_dt_s)
        eta1_rad, eta2_rad = wing_profile(t_grid, args.wing_amp_deg, args.wing_period_s)
        res = run_case(
            "attitude",
            outdir,
            args.attitude_duration_s,
            args.attitude_dt_s,
            args.attitude_particles,
            config_att,
            props_att,
            EnvInputs,
            eta1_rad=eta1_rad,
            eta2_rad=eta2_rad,
            P0=P0_att,
            augment_process_noise=args.augment_process_noise,
            seed=args.seed,
            do_pod=args.pod_uq,
            pod_arc_s=args.pod_arc_s,
            pod_overlap_s=args.pod_overlap_s,
            pod_skip_slr=args.pod_skip_slr,
            pod_sp3=args.pod_sp3,
            pod_eop=args.pod_eop,
            pod_use_sat_clock=args.pod_use_sat_clock,
            pod_estimator=args.pod_estimator,
            pod_enkf_members=args.pod_enkf_members,
            pod_enkf_seed=None,
            pod_enkf_inflation=args.pod_enkf_inflation,
            pod_enkf_use_carrier=args.pod_enkf_use_carrier,
            pod_gnss_cycle_slip_gap_s=args.pod_gnss_cycle_slip_gap_s,
            pod_gnss_cycle_slip_prob_per_min=args.pod_gnss_cycle_slip_prob_per_min,
            pod_gnss_code_outlier_prob=args.pod_gnss_code_outlier_prob,
            pod_gnss_code_outlier_sigma_scale=args.pod_gnss_code_outlier_sigma_scale,
            pod_gnss_carrier_outlier_prob=args.pod_gnss_carrier_outlier_prob,
            pod_gnss_carrier_outlier_sigma_scale=args.pod_gnss_carrier_outlier_sigma_scale,
            pod_slr_weather_on=args.pod_slr_weather_on,
            pod_slr_weather_clear_prob=args.pod_slr_weather_clear_prob,
            pod_slr_weather_p_stay_clear=args.pod_slr_weather_p_stay_clear,
            pod_slr_weather_p_stay_blocked=args.pod_slr_weather_p_stay_blocked,
            pod_slr_weather_seed=args.pod_slr_weather_seed,
            use_env_sources=args.use_env_sources,
            start_utc=start_utc,
            space_weather_source=args.space_weather_source,
            omni_path=args.omni_path,
            hwm14_lib=args.hwm14_lib,
            hwm14_data=args.hwm14_data,
            interpolation=args.env_interpolation,
            attitude_reference="flow_x",
            scenario_meta={"name": "attitude", "catalog_scenario_id": "ATTITUDE"},
        )
        if args.plot:
            plot_case(
                outdir,
                "attitude",
                res["t_grid"],
                res["det_states"],
                res["mc_mean"],
                res["mc_std"],
                res["ut_mean"],
                res["stm_mean"],
                args.plot_all_states,
            )


if __name__ == "__main__":
    main()
