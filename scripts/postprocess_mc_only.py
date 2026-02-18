from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np

try:
    from scripts.summary_metadata import build_run_metadata
except ModuleNotFoundError:
    from summary_metadata import build_run_metadata


def save_json(path: Path, payload: dict) -> None:
    def to_jsonable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        return obj

    data = {k: to_jsonable(v) for k, v in payload.items()}
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def maybe_import_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    return plt


def quaternion_to_euler_zyx_deg(q_wxyz: np.ndarray) -> np.ndarray:
    q = np.array(q_wxyz, dtype=float).reshape(-1)
    if q.size != 4:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    n = np.linalg.norm(q)
    if not np.isfinite(n) or n <= 0.0:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    w, x, y, z = q / n
    # ZYX (yaw-pitch-roll) intrinsic convention.
    sinr_cosp = 2.0 * (w * x + y * z)
    cosr_cosp = 1.0 - 2.0 * (x * x + y * y)
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    sinp = 2.0 * (w * y - z * x)
    sinp = np.clip(sinp, -1.0, 1.0)
    pitch = np.arcsin(sinp)

    siny_cosp = 2.0 * (w * z + x * y)
    cosy_cosp = 1.0 - 2.0 * (y * y + z * z)
    yaw = np.arctan2(siny_cosp, cosy_cosp)
    return np.degrees(np.array([roll, pitch, yaw], dtype=float))


def plot_mc_case(
    case_dir: Path,
    t_grid: np.ndarray,
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    plot_all_states: bool,
    q_particles: np.ndarray | None = None,
) -> None:
    plt = maybe_import_matplotlib()
    if plt is None:
        print("[mc-only] matplotlib unavailable; skipping plots", flush=True)
        return

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    labels = ["x", "y", "z"]
    for i, ax in enumerate(axes):
        ax.plot(t_grid, mc_mean[:, i], label="mc_mean")
        ax.fill_between(
            t_grid,
            mc_mean[:, i] - mc_std[:, i],
            mc_mean[:, i] + mc_std[:, i],
            color="C0",
            alpha=0.2,
            label="mc_std" if i == 0 else None,
        )
        ax.set_ylabel(f"r_{labels[i]} [m]")
    axes[-1].set_xlabel("t [s]")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "mc_position_components.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    for i, ax in enumerate(axes):
        k = 3 + i
        ax.plot(t_grid, mc_mean[:, k], label="mc_mean")
        ax.fill_between(
            t_grid,
            mc_mean[:, k] - mc_std[:, k],
            mc_mean[:, k] + mc_std[:, k],
            color="C1",
            alpha=0.2,
            label="mc_std" if i == 0 else None,
        )
        ax.set_ylabel(f"v_{labels[i]} [m/s]")
    axes[-1].set_xlabel("t [s]")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "mc_velocity_components.png", dpi=150)
    plt.close(fig)

    r = mc_mean[:, 0:3]
    re_m = 6378137.0
    alt_m = np.linalg.norm(r, axis=1) - re_m
    safe_r = np.maximum(np.linalg.norm(r, axis=1), 1.0)
    radial_std_m = np.sqrt(np.sum((r / safe_r[:, None]) ** 2 * (mc_std[:, 0:3] ** 2), axis=1))

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_grid, alt_m, label="mc_mean altitude")
    ax.fill_between(
        t_grid,
        alt_m - radial_std_m,
        alt_m + radial_std_m,
        color="C2",
        alpha=0.2,
        label="approx. radial std",
    )
    ax.set_xlabel("t [s]")
    ax.set_ylabel("altitude [m]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(case_dir / "mc_altitude.png", dpi=150)
    plt.close(fig)

    if mc_mean.shape[1] >= 13:
        fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
        for i, ax in enumerate(axes):
            k = 10 + i
            mean_deg = np.degrees(mc_mean[:, k])
            std_deg = np.degrees(mc_std[:, k])
            ax.plot(t_grid, mean_deg, label="mc_mean")
            ax.fill_between(
                t_grid,
                mean_deg - std_deg,
                mean_deg + std_deg,
                color="C3",
                alpha=0.2,
                label="mc_std" if i == 0 else None,
            )
            ax.set_ylabel(f"w_{labels[i]} [deg/s]")
        axes[-1].set_xlabel("t [s]")
        axes[0].legend()
        fig.tight_layout()
        fig.savefig(case_dir / "mc_attitude_rates.png", dpi=150)
        plt.close(fig)

    if mc_mean.shape[1] >= 10:
        eul = np.array([quaternion_to_euler_zyx_deg(mc_mean[i, 6:10]) for i in range(mc_mean.shape[0])], dtype=float)
        eul_std = None
        mean_label = "mc_mean"
        std_label = "mc_std"
        if isinstance(q_particles, np.ndarray) and q_particles.ndim == 3 and q_particles.shape[0] == t_grid.size:
            particle_stats = euler_stats_from_particles_deg(q_particles)
            if particle_stats is not None:
                eul, eul_std = particle_stats
                mean_label = "particle_mean"
                std_label = "particle_std"
        fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
        eul_labels = ["roll", "pitch", "yaw"]
        for i, ax in enumerate(axes):
            ax.plot(t_grid, eul[:, i], label=mean_label)
            if eul_std is not None:
                ax.fill_between(
                    t_grid,
                    eul[:, i] - eul_std[:, i],
                    eul[:, i] + eul_std[:, i],
                    color="C4",
                    alpha=0.2,
                    label=std_label if i == 0 else None,
                )
            ax.set_ylabel(f"{eul_labels[i]} [deg]")
        axes[-1].set_xlabel("t [s]")
        axes[0].legend()
        fig.tight_layout()
        fig.savefig(case_dir / "mc_attitude_euler_deg.png", dpi=150)
        plt.close(fig)

    if plot_all_states:
        n_state = mc_mean.shape[1]
        n_cols = 3
        n_rows = int(np.ceil(n_state / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, figsize=(11, 3 * n_rows))
        axes = np.atleast_1d(axes).ravel()
        for i in range(n_state):
            ax = axes[i]
            ax.plot(t_grid, mc_mean[:, i], lw=1.0)
            ax.fill_between(
                t_grid,
                mc_mean[:, i] - mc_std[:, i],
                mc_mean[:, i] + mc_std[:, i],
                alpha=0.2,
            )
            ax.set_ylabel(f"x[{i}]")
        for i in range(n_state, len(axes)):
            axes[i].axis("off")
        axes[min(n_state - 1, len(axes) - 1)].set_xlabel("t [s]")
        fig.tight_layout()
        fig.savefig(case_dir / "mc_all_states.png", dpi=150)
        plt.close(fig)


def plot_formation_separation(case_dir: Path, t_grid: np.ndarray, separations_m: np.ndarray) -> None:
    plt = maybe_import_matplotlib()
    if plt is None:
        print("[mc-only] matplotlib unavailable; skipping formation separation plot", flush=True)
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    for i in range(separations_m.shape[0]):
        ax.plot(t_grid, separations_m[i], label=f"chief-follower {i+1}")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("separation [m]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(case_dir / "mc_formation_separation.png", dpi=150)
    plt.close(fig)


def plot_maneuver_torque(case_dir: Path, t_grid: np.ndarray, diag: dict) -> None:
    plt = maybe_import_matplotlib()
    if plt is None:
        return

    rep_tau = np.asarray(diag.get("rep_torque_B_Nm"))
    if rep_tau.ndim != 2 or rep_tau.shape[0] != t_grid.size or rep_tau.shape[1] != 3:
        return
    det_tau = np.asarray(diag.get("det_torque_B_Nm")) if "det_torque_B_Nm" in diag else None
    rep_norm = np.linalg.norm(rep_tau, axis=1)
    det_norm = np.linalg.norm(det_tau, axis=1) if isinstance(det_tau, np.ndarray) and det_tau.shape == rep_tau.shape else None

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(9, 10))
    labels = ["x", "y", "z"]
    for i in range(3):
        axes[i].plot(t_grid, rep_tau[:, i], label="MC representative")
        if isinstance(det_tau, np.ndarray) and det_tau.shape == rep_tau.shape:
            axes[i].plot(t_grid, det_tau[:, i], label="deterministic", alpha=0.8)
        axes[i].set_ylabel(f"tau_{labels[i]} [Nm]")
        axes[i].grid(True, alpha=0.3)
    axes[3].plot(t_grid, rep_norm, label="MC representative")
    if det_norm is not None:
        axes[3].plot(t_grid, det_norm, label="deterministic", alpha=0.8)
    axes[3].set_ylabel("|tau| [Nm]")
    axes[3].set_xlabel("t [s]")
    axes[3].grid(True, alpha=0.3)
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(case_dir / "maneuver_torque_timeseries.png", dpi=150)
    plt.close(fig)


def load_mc_stats(path: Path):
    with np.load(path) as data:
        t_grid = np.array(data["t_grid"], dtype=float)
        mc_mean = np.array(data["mc_mean"], dtype=float)
        mc_std = np.array(data["mc_std"], dtype=float)
        mc_cov_final = np.array(data["mc_cov_final"], dtype=float)
        n_total = int(data["n_total"]) if "n_total" in data else -1
    return t_grid, mc_mean, mc_std, mc_cov_final, n_total


def load_attitude_particles_from_shards(shard_dir: Path, t_grid_ref: np.ndarray) -> dict | None:
    paths = sorted(shard_dir.glob("ensemble_rank*.npz"))
    if not paths:
        return None

    q_parts = []
    w_parts = []
    n_ref = int(t_grid_ref.size)
    for path in paths:
        try:
            with np.load(path) as data:
                X = np.array(data["X"], dtype=float)
        except Exception:
            continue
        if X.ndim != 3 or X.shape[1] <= 0 or X.shape[2] < 13:
            continue
        n_use = min(n_ref, int(X.shape[0]))
        if n_use <= 0:
            continue
        q_parts.append(X[:n_use, :, 6:10])
        w_parts.append(X[:n_use, :, 10:13])

    if not q_parts:
        return None

    n_common = min(int(q.shape[0]) for q in q_parts)
    q_all = np.concatenate([q[:n_common] for q in q_parts], axis=1)
    w_all = np.concatenate([w[:n_common] for w in w_parts], axis=1)
    return {"t_grid": np.array(t_grid_ref[:n_common], dtype=float), "q": q_all, "w": w_all}


def load_maneuver_diagnostics(shard_dir: Path, t_grid_ref: np.ndarray) -> dict | None:
    paths = sorted(shard_dir.glob("maneuver_diag_rank*.npz"))
    if not paths:
        return None
    path = paths[0]
    try:
        with np.load(path) as data:
            out = {k: np.array(data[k]) for k in data.files}
    except Exception:
        return None

    nt = int(t_grid_ref.size)
    for key in (
        "rep_progress_axis_deg",
        "rep_eta1_cmd_rad",
        "rep_eta2_cmd_rad",
        "rep_eta_cmd_deg",
        "rep_eta_target_deg",
        "rep_eta_rate_deg_s",
        "rep_torque_B_Nm",
        "rep_phase_index",
        "det_progress_axis_deg",
        "det_eta1_cmd_rad",
        "det_eta2_cmd_rad",
        "det_eta_cmd_deg",
        "det_eta_target_deg",
        "det_eta_rate_deg_s",
        "det_torque_B_Nm",
        "det_phase_index",
    ):
        if key not in out:
            continue
        arr = np.asarray(out[key])
        if arr.ndim == 1:
            out[key] = arr[:nt]
        elif arr.ndim >= 2:
            out[key] = arr[:nt, ...]
    out["t_grid"] = np.array(t_grid_ref[: nt], dtype=float)
    return out


def _unwrap_axis_trace_from_quat_particles(q_particles: np.ndarray, axis: int) -> np.ndarray:
    nt, npart, _ = q_particles.shape
    trace = np.full((nt, npart), np.nan, dtype=float)
    for j in range(npart):
        e0 = quaternion_to_euler_zyx_deg(q_particles[0, j])
        prev = float(e0[axis])
        if not np.isfinite(prev):
            continue
        trace[0, j] = prev
        for i in range(1, nt):
            ei = quaternion_to_euler_zyx_deg(q_particles[i, j])
            curr = float(ei[axis])
            if not np.isfinite(curr):
                trace[i, j] = trace[i - 1, j]
                continue
            d = ((curr - prev + 180.0) % 360.0) - 180.0
            trace[i, j] = trace[i - 1, j] + d
            prev = curr
    return trace


def euler_stats_from_particles_deg(q_particles: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    if q_particles.ndim != 3 or q_particles.shape[2] < 4:
        return None
    nt = int(q_particles.shape[0])
    eul_trace = np.full((nt, q_particles.shape[1], 3), np.nan, dtype=float)
    for axis in range(3):
        eul_trace[:, :, axis] = _unwrap_axis_trace_from_quat_particles(q_particles, axis)

    finite = np.isfinite(eul_trace)
    counts = np.sum(finite, axis=1)
    sums = np.sum(np.where(finite, eul_trace, 0.0), axis=1)
    means = np.divide(sums, counts, out=np.full((nt, 3), np.nan, dtype=float), where=counts > 0)
    centered_sq = np.where(finite, (eul_trace - means[:, None, :]) ** 2, 0.0)
    var = np.divide(np.sum(centered_sq, axis=1), counts, out=np.full((nt, 3), np.nan, dtype=float), where=counts > 0)
    std = np.sqrt(np.maximum(var, 0.0))
    return means, std


def build_attitude_metric_overrides(
    *,
    requested_metrics: list[str],
    t_grid: np.ndarray,
    q_particles: np.ndarray | None,
    w_particles: np.ndarray | None,
) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if q_particles is None or q_particles.ndim != 3 or q_particles.shape[0] != t_grid.size:
        return out

    axis_idx = {"roll": 0, "pitch": 1, "yaw": 2}
    axis_trace_cache: dict[int, np.ndarray] = {}

    def axis_trace(axis_name: str) -> np.ndarray:
        idx = axis_idx[axis_name]
        if idx not in axis_trace_cache:
            axis_trace_cache[idx] = _unwrap_axis_trace_from_quat_particles(q_particles, idx)
        return axis_trace_cache[idx]

    for metric in requested_metrics:
        m = str(metric)
        m_time = re.fullmatch(r"time_to_(roll|pitch|yaw)_(\d+(?:\.\d+)?)_deg", m)
        m_over = re.fullmatch(r"overshoot_(roll|pitch|yaw)_(\d+(?:\.\d+)?)_deg", m)
        m_peak = re.fullmatch(r"peak_(roll|pitch|yaw)_rate_deg_s", m)
        m_rms_exc = re.fullmatch(r"rms_(roll|pitch|yaw)_excursion_deg", m)

        if m_time:
            axis_name = m_time.group(1)
            target_deg = float(m_time.group(2))
            tr = axis_trace(axis_name)
            progress = np.abs(tr - tr[0:1, :])
            vals = []
            for j in range(progress.shape[1]):
                idx = np.where(np.isfinite(progress[:, j]) & (progress[:, j] >= target_deg))[0]
                if idx.size > 0:
                    vals.append(float(t_grid[int(idx[0])]))
            if vals:
                out[m] = {
                    "value": float(np.nanmedian(np.array(vals, dtype=float))),
                    "status": "computed_from_particles",
                    "note": f"median over reached particles ({len(vals)}/{progress.shape[1]})",
                    "unit": "s",
                }
            else:
                out[m] = {"value": None, "status": "not_reached", "note": "no particle reached target", "unit": "s"}
            continue

        if m_over:
            axis_name = m_over.group(1)
            target_deg = float(m_over.group(2))
            tr = axis_trace(axis_name)
            progress = np.abs(tr - tr[0:1, :])
            max_prog = np.nanmax(progress, axis=0)
            over = np.maximum(0.0, max_prog - target_deg)
            over = over[np.isfinite(over)]
            if over.size > 0:
                out[m] = {
                    "value": float(np.nanmedian(over)),
                    "status": "computed_from_particles",
                    "note": "median overshoot across particles",
                    "unit": "deg",
                }
            else:
                out[m] = {"value": None, "status": "unavailable_in_mc_only"}
            continue

        if m_peak and w_particles is not None and w_particles.ndim == 3 and w_particles.shape[0] == t_grid.size:
            axis_name = m_peak.group(1)
            ax = axis_idx[axis_name]
            peak_per_particle = np.nanmax(np.abs(np.degrees(w_particles[:, :, ax])), axis=0)
            peak_per_particle = peak_per_particle[np.isfinite(peak_per_particle)]
            if peak_per_particle.size > 0:
                out[m] = {
                    "value": float(np.nanmedian(peak_per_particle)),
                    "status": "computed_from_particles",
                    "note": "median of per-particle peak rates",
                    "unit": "deg/s",
                }
            else:
                out[m] = {"value": None, "status": "unavailable_in_mc_only"}
            continue

        if m_rms_exc:
            axis_name = m_rms_exc.group(1)
            tr = axis_trace(axis_name)
            delta = tr - tr[0:1, :]
            rms = np.sqrt(np.nanmean(delta ** 2, axis=0))
            rms = rms[np.isfinite(rms)]
            if rms.size > 0:
                out[m] = {
                    "value": float(np.nanmedian(rms)),
                    "status": "computed_from_particles",
                    "note": "median RMS excursion across particles",
                    "unit": "deg",
                }
            else:
                out[m] = {"value": None, "status": "unavailable_in_mc_only"}

    return out


def prepare_mc_series_for_analysis(
    *,
    t_grid: np.ndarray,
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    mc_cov_final: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    reentry_alt_m = float(os.environ.get("VLEO_REENTRY_ALTITUDE_M", "120000.0"))
    info = {
        "raw_steps": int(t_grid.size),
        "raw_duration_s": float(t_grid[-1] - t_grid[0]) if t_grid.size > 1 else 0.0,
        "reentry_detected": False,
        "reentry_altitude_threshold_m": reentry_alt_m if reentry_alt_m > 0.0 else None,
        "reentry_time_s": None,
        "reentry_altitude_m": None,
        "analysis_stop_reason": "full_span",
    }
    if reentry_alt_m <= 0.0 or mc_mean.shape[1] < 3:
        info["analysis_steps"] = int(t_grid.size)
        info["analysis_duration_s"] = float(t_grid[-1] - t_grid[0]) if t_grid.size > 1 else 0.0
        return t_grid, mc_mean, mc_std, mc_cov_final, info

    r = mc_mean[:, 0:3]
    re_m = 6378137.0
    alt_m = np.linalg.norm(r, axis=1) - re_m
    idx = np.where(np.isfinite(alt_m) & (alt_m <= reentry_alt_m))[0]
    if idx.size == 0:
        info["analysis_steps"] = int(t_grid.size)
        info["analysis_duration_s"] = float(t_grid[-1] - t_grid[0]) if t_grid.size > 1 else 0.0
        return t_grid, mc_mean, mc_std, mc_cov_final, info

    i0 = int(idx[0])
    end = max(1, i0 + 1)
    t_eff = np.array(t_grid[:end], dtype=float)
    mean_eff = np.array(mc_mean[:end], dtype=float)
    std_eff = np.array(mc_std[:end], dtype=float)
    cov_eff = np.diag(np.square(std_eff[-1])) if std_eff.ndim == 2 and std_eff.shape[0] > 0 else mc_cov_final

    info["reentry_detected"] = True
    info["reentry_time_s"] = float(t_grid[i0])
    info["reentry_altitude_m"] = float(alt_m[i0])
    info["analysis_stop_reason"] = "reentry"
    info["analysis_steps"] = int(t_eff.size)
    info["analysis_duration_s"] = float(t_eff[-1] - t_eff[0]) if t_eff.size > 1 else 0.0
    return t_eff, mean_eff, std_eff, cov_eff, info


def summarize_mc_state(
    *,
    name: str,
    scenario_type: str,
    requested_metrics: list[str],
    t_grid: np.ndarray,
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    mc_cov_final: np.ndarray,
    n_total: int,
    analysis_info: dict | None = None,
) -> dict:
    r = mc_mean[:, 0:3]
    v = mc_mean[:, 3:6]
    re_m = 6378137.0
    alt_m = np.linalg.norm(r, axis=1) - re_m
    speed_mps = np.linalg.norm(v, axis=1)
    w_norm = np.linalg.norm(mc_mean[:, 10:13], axis=1) if mc_mean.shape[1] >= 13 else np.array([])

    summary = {
        "case": name,
        "mode": "mc_only",
        "type": scenario_type,
        "requested_catalog_metrics": requested_metrics,
        "duration_s": float(t_grid[-1] - t_grid[0]) if t_grid.size > 1 else 0.0,
        "steps": int(t_grid.size),
        "n_mc_effective": int(n_total),
        "final_altitude_m": float(alt_m[-1]),
        "min_altitude_m": float(np.nanmin(alt_m)),
        "max_altitude_m": float(np.nanmax(alt_m)),
        "final_speed_mps": float(speed_mps[-1]),
        "max_speed_mps": float(np.nanmax(speed_mps)),
        "final_pos_std_norm_m": float(np.linalg.norm(mc_std[-1, 0:3])),
        "final_vel_std_norm_mps": float(np.linalg.norm(mc_std[-1, 3:6])),
        "final_cov_trace": float(np.trace(mc_cov_final)) if mc_cov_final.ndim == 2 else float("nan"),
        "finite_mc_mean": bool(np.isfinite(mc_mean).all()),
        "finite_mc_std": bool(np.isfinite(mc_std).all()),
        "finite_mc_cov_final": bool(np.isfinite(mc_cov_final).all()),
    }
    if w_norm.size:
        summary["initial_body_rate_rad_s"] = [float(x) for x in mc_mean[0, 10:13]]
        summary["initial_body_rate_norm_rad_s"] = float(w_norm[0])
        summary["initial_body_rate_norm_deg_s"] = float(np.degrees(w_norm[0]))
        summary["final_rate_norm_rad_s"] = float(w_norm[-1])
        summary["max_rate_norm_rad_s"] = float(np.nanmax(w_norm))
    if analysis_info:
        summary.update(analysis_info)
    return summary


def quaternion_error_deg_from_initial(mc_mean: np.ndarray) -> np.ndarray:
    if mc_mean.shape[1] < 10:
        return np.zeros(mc_mean.shape[0], dtype=float)
    q = np.array(mc_mean[:, 6:10], dtype=float)
    q0 = np.array(q[0], dtype=float)
    n0 = np.linalg.norm(q0)
    if n0 == 0.0:
        q0 = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    else:
        q0 = q0 / n0
    out = np.zeros(q.shape[0], dtype=float)
    for i in range(q.shape[0]):
        qi = q[i]
        ni = np.linalg.norm(qi)
        if ni == 0.0:
            qi = q0
        else:
            qi = qi / ni
        dot = float(np.clip(np.abs(np.dot(qi, q0)), -1.0, 1.0))
        out[i] = np.degrees(2.0 * np.arccos(dot))
    return out


def build_catalog_core_outputs(
    *,
    requested_metrics: list[str],
    t_grid: np.ndarray,
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    separations_m: np.ndarray | None = None,
    q_particles: np.ndarray | None = None,
    w_particles: np.ndarray | None = None,
) -> dict:
    outputs: dict[str, dict] = {}
    duration_s = float(t_grid[-1] - t_grid[0]) if t_grid.size > 1 else 0.0
    eps = 1e-12

    r = mc_mean[:, 0:3] if mc_mean.shape[1] >= 3 else np.zeros((mc_mean.shape[0], 3), dtype=float)
    v = mc_mean[:, 3:6] if mc_mean.shape[1] >= 6 else np.zeros((mc_mean.shape[0], 3), dtype=float)
    re_m = 6378137.0
    alt_m = np.linalg.norm(r, axis=1) - re_m
    speed_mps = np.linalg.norm(v, axis=1)
    pos_std_norm = np.linalg.norm(mc_std[:, 0:3], axis=1) if mc_std.shape[1] >= 3 else np.zeros_like(t_grid)
    vel_std_norm = np.linalg.norm(mc_std[:, 3:6], axis=1) if mc_std.shape[1] >= 6 else np.zeros_like(t_grid)
    rate_norm = np.linalg.norm(mc_mean[:, 10:13], axis=1) if mc_mean.shape[1] >= 13 else np.zeros_like(t_grid)
    rate_norm_deg = np.degrees(rate_norm)
    qerr_deg = quaternion_error_deg_from_initial(mc_mean)
    euler_deg_unwrapped = None
    if mc_mean.shape[1] >= 10:
        euler_deg = np.array([quaternion_to_euler_zyx_deg(mc_mean[i, 6:10]) for i in range(mc_mean.shape[0])], dtype=float)
        euler_deg_unwrapped = np.degrees(np.unwrap(np.radians(euler_deg), axis=0))
    detumble_thresh_deg_s = float(np.clip(float(os.environ.get("VLEO_DETUMBLE_RATE_THRESH_DEG_S", "0.5")), 0.01, 30.0))

    sep = separations_m
    sep_primary = sep[0] if isinstance(sep, np.ndarray) and sep.ndim == 2 and sep.shape[0] > 0 else None
    sep_threshold_m = float(os.environ.get("VLEO_CLOSE_APPROACH_THRESHOLD_M", "100.0"))
    target_tol_m = float(os.environ.get("VLEO_RECONFIG_TOL_M", "10.0"))
    attitude_overrides = build_attitude_metric_overrides(
        requested_metrics=requested_metrics,
        t_grid=t_grid,
        q_particles=q_particles,
        w_particles=w_particles,
    )

    def put(metric: str, value, status: str, note: str = "", unit: str | None = None) -> None:
        rec = {"value": value, "status": status}
        if note:
            rec["note"] = note
        if unit:
            rec["unit"] = unit
        outputs[metric] = rec

    for metric in requested_metrics:
        m = str(metric)
        if m in attitude_overrides:
            outputs[m] = attitude_overrides[m]
            continue
        m_time = re.fullmatch(r"time_to_(roll|pitch|yaw)_(\d+(?:\.\d+)?)_deg", m)
        m_over = re.fullmatch(r"overshoot_(roll|pitch|yaw)_(\d+(?:\.\d+)?)_deg", m)
        m_peak = re.fullmatch(r"peak_(roll|pitch|yaw)_rate_deg_s", m)
        m_rms_exc = re.fullmatch(r"rms_(roll|pitch|yaw)_excursion_deg", m)

        if m_time:
            if euler_deg_unwrapped is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                axis_name = m_time.group(1)
                target_deg = float(m_time.group(2))
                axis = {"roll": 0, "pitch": 1, "yaw": 2}[axis_name]
                delta = np.abs(euler_deg_unwrapped[:, axis] - euler_deg_unwrapped[0, axis])
                idx = np.where(delta >= target_deg)[0]
                val = float(t_grid[int(idx[0])]) if idx.size > 0 else None
                put(m, val, "computed" if val is not None else "not_reached", "from unwrapped Euler ZYX", "s")
            continue
        if m_over:
            if euler_deg_unwrapped is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                axis_name = m_over.group(1)
                target_deg = float(m_over.group(2))
                axis = {"roll": 0, "pitch": 1, "yaw": 2}[axis_name]
                delta = np.abs(euler_deg_unwrapped[:, axis] - euler_deg_unwrapped[0, axis])
                val = float(max(0.0, np.nanmax(delta) - target_deg))
                put(m, val, "computed", "from unwrapped Euler ZYX", "deg")
            continue
        if m_peak:
            axis_name = m_peak.group(1)
            axis = {"roll": 0, "pitch": 1, "yaw": 2}[axis_name]
            if mc_mean.shape[1] >= 13:
                val = float(np.nanmax(np.abs(np.degrees(mc_mean[:, 10 + axis]))))
                put(m, val, "computed", unit="deg/s")
            else:
                put(m, None, "unavailable_in_mc_only")
            continue
        if m_rms_exc:
            if euler_deg_unwrapped is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                axis_name = m_rms_exc.group(1)
                axis = {"roll": 0, "pitch": 1, "yaw": 2}[axis_name]
                delta = euler_deg_unwrapped[:, axis] - euler_deg_unwrapped[0, axis]
                val = float(np.sqrt(np.nanmean(delta ** 2)))
                put(m, val, "computed", "from unwrapped Euler ZYX", "deg")
            continue

        if m == "time_to_detumble":
            idx = np.where(rate_norm_deg <= detumble_thresh_deg_s)[0]
            val = float(t_grid[int(idx[0])]) if idx.size > 0 else None
            put(m, val, "computed" if val is not None else "not_reached", f"threshold={detumble_thresh_deg_s:.3f} deg/s", "s")
        elif m == "max_body_rate":
            put(m, float(np.nanmax(rate_norm_deg)), "computed", unit="deg/s")
        elif m == "rms_body_rate":
            put(m, float(np.sqrt(np.nanmean(rate_norm_deg ** 2))), "computed", unit="deg/s")
        elif m in {"rms_pointing_error", "rms_yaw_error"}:
            put(m, float(np.sqrt(np.nanmean(qerr_deg ** 2))), "proxy", "quat angle vs initial attitude", "deg")
        elif m in {"p95_pointing_error", "p95_yaw_error"}:
            put(m, float(np.nanpercentile(qerr_deg, 95.0)), "proxy", "quat angle vs initial attitude", "deg")
        elif m == "max_pointing_error":
            put(m, float(np.nanmax(qerr_deg)), "proxy", "quat angle vs initial attitude", "deg")
        elif m == "delta_altitude":
            put(m, float(alt_m[-1] - alt_m[0]), "computed", unit="m")
        elif m == "mean_da_dt":
            val = float((alt_m[-1] - alt_m[0]) / max(duration_s, eps))
            put(m, val, "computed", unit="m/s")
        elif m == "pos_rms_eci":
            put(m, float(np.sqrt(np.nanmean(pos_std_norm ** 2))), "proxy", "from MC spread (std norm)", "m")
        elif m == "vel_rms_eci":
            put(m, float(np.sqrt(np.nanmean(vel_std_norm ** 2))), "proxy", "from MC spread (std norm)", "m/s")
        elif m in {"along_track_sigma", "cross_track_sigma", "radial_sigma"}:
            axis = {"radial_sigma": 0, "along_track_sigma": 1, "cross_track_sigma": 2}[m]
            val = float(mc_std[-1, axis]) if mc_std.shape[1] > axis else None
            put(m, val, "proxy", "ECI component std proxy at final epoch", "m")
        elif m == "expected_pos_rms":
            put(m, float(np.nanmean(pos_std_norm)), "proxy", "mean position std norm", "m")
        elif m == "percentile_pos_rms":
            put(m, float(np.nanpercentile(pos_std_norm, 95.0)), "proxy", "p95 position std norm", "m")
        elif m in {"sep_error_rms", "lvhl_rms_error"}:
            if sep_primary is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                err = sep_primary - sep_primary[0]
                put(m, float(np.sqrt(np.nanmean(err ** 2))), "proxy", "chief-follower-1 separation error", "m")
        elif m == "sep_error_p95":
            if sep_primary is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                err = np.abs(sep_primary - sep_primary[0])
                put(m, float(np.nanpercentile(err, 95.0)), "proxy", "chief-follower-1 separation error", "m")
        elif m in {"min_inter_sat_range", "min_range_distribution", "collision_margin_min"}:
            if sep is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                put(m, float(np.nanmin(sep)), "computed", unit="m")
        elif m in {"reconfig_time", "time_to_target_sep"}:
            if sep_primary is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                target_sep = None
                target_note = "target=initial_separation"
                if m == "time_to_target_sep":
                    env_target = os.environ.get("VLEO_TARGET_SEP_M", "").strip()
                    if env_target:
                        try:
                            target_sep = float(env_target)
                            target_note = f"target={target_sep:.3f} m (env)"
                        except ValueError:
                            target_sep = None
                if target_sep is None:
                    target_sep = float(sep_primary[0])

                err = np.abs(sep_primary - target_sep)
                within = err <= target_tol_m
                if within.size == 0:
                    val = None
                elif target_note == "target=initial_separation":
                    # Reconfiguration time is meaningful only after leaving the tolerance tube.
                    out_idx = np.where(~within)[0]
                    if out_idx.size == 0:
                        val = 0.0
                        put(
                            m,
                            val,
                            "proxy",
                            f"{target_note}; tol={target_tol_m:.3f} m; always within tolerance",
                            "s",
                        )
                        continue
                    reacq_idx = np.where(within & (np.arange(within.size) > out_idx[0]))[0]
                    val = float(t_grid[int(reacq_idx[0])]) if reacq_idx.size > 0 else None
                else:
                    idx = np.where(within & (np.arange(within.size) > 0))[0]
                    val = float(t_grid[int(idx[0])]) if idx.size > 0 else None

                put(
                    m,
                    val,
                    "proxy" if val is not None else "not_reached",
                    f"{target_note}; tol={target_tol_m:.3f} m",
                    "s",
                )
        elif m in {"close_approach_count"}:
            if sep is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                below = sep < sep_threshold_m
                events = int(np.sum(np.diff(below.astype(int), axis=1) == 1))
                put(m, events, "proxy", f"threshold={sep_threshold_m:.3f} m")
        elif m in {
            "p_min_range_below_threshold",
            "box_violation_prob",
            "constraint_violation_prob",
            "formation_violation_prob",
            "estimator_divergence_prob",
            "overshoot_prob",
        }:
            if sep is None:
                put(m, None, "unavailable_in_mc_only")
            else:
                val = 1.0 if np.nanmin(sep) < sep_threshold_m else 0.0
                put(m, float(val), "proxy", "mean-trajectory indicator only (not Monte Carlo event probability)")
        else:
            put(m, None, "unavailable_in_mc_only")

    return outputs


def _parse_reached_fraction(note: str | None) -> float | None:
    if not note:
        return None
    m = re.search(r"\((\d+)\s*/\s*(\d+)\)", str(note))
    if not m:
        return None
    den = int(m.group(2))
    if den <= 0:
        return None
    return float(int(m.group(1)) / den)


def build_maneuver_actuator_stats(
    *,
    diag: dict,
    t_grid: np.ndarray,
    core_outputs: dict,
) -> dict:
    axis_raw = diag.get("axis")
    axis = str(np.asarray(axis_raw).item()) if axis_raw is not None else ""
    series_prefix = "det" if f"det_phase_index" in diag else "rep"
    phase_idx = np.asarray(diag.get(f"{series_prefix}_phase_index"))
    eta_cmd_deg = np.asarray(diag.get(f"{series_prefix}_eta_cmd_deg"))
    eta_target_deg = np.asarray(diag.get(f"{series_prefix}_eta_target_deg"))
    eta_rate_deg_s = np.asarray(diag.get(f"{series_prefix}_eta_rate_deg_s"))
    progress_axis_deg = np.asarray(diag.get(f"{series_prefix}_progress_axis_deg"))
    eta_max_deg = float(np.asarray(diag.get("eta_max_deg")).item()) if "eta_max_deg" in diag else None
    eta_rate_max_deg_s = float(np.asarray(diag.get("eta_rate_max_deg_s")).item()) if "eta_rate_max_deg_s" in diag else None
    eta_accel_max_deg_s2 = float(np.asarray(diag.get("eta_accel_max_deg_s2")).item()) if "eta_accel_max_deg_s2" in diag else None
    half_turn_deg = float(np.asarray(diag.get("half_turn_deg")).item()) if "half_turn_deg" in diag else None
    full_turn_deg = float(np.asarray(diag.get("full_turn_deg")).item()) if "full_turn_deg" in diag else None

    stats: dict[str, object] = {
        "source": series_prefix,
        "axis": axis,
        "control_switch_count": int(np.sum(np.diff(phase_idx) != 0)) if phase_idx.ndim == 1 and phase_idx.size > 1 else None,
        "command_max_abs_deg": float(np.nanmax(np.abs(eta_cmd_deg))) if eta_cmd_deg.ndim == 1 and eta_cmd_deg.size > 0 else None,
        "target_max_abs_deg": float(np.nanmax(np.abs(eta_target_deg))) if eta_target_deg.ndim == 1 and eta_target_deg.size > 0 else None,
        "command_rate_max_abs_deg_s": float(np.nanmax(np.abs(eta_rate_deg_s))) if eta_rate_deg_s.ndim == 1 and eta_rate_deg_s.size > 0 else None,
        "progress_max_abs_deg": float(np.nanmax(np.abs(progress_axis_deg))) if progress_axis_deg.ndim == 1 and progress_axis_deg.size > 0 else None,
    }

    if eta_cmd_deg.ndim == 1 and eta_cmd_deg.size > 1:
        dcmd = np.diff(eta_cmd_deg)
        dt = np.diff(t_grid)
        dcmd_rate = np.divide(dcmd, dt, out=np.full_like(dcmd, np.nan), where=dt > 0.0)
        stats["command_jump_max_abs_deg"] = float(np.nanmax(np.abs(dcmd)))
        stats["command_jump_rate_max_abs_deg_s"] = float(np.nanmax(np.abs(dcmd_rate)))

    if eta_max_deg is not None and eta_cmd_deg.ndim == 1 and eta_cmd_deg.size > 0:
        sat_mask = np.abs(eta_cmd_deg) >= max(0.0, eta_max_deg - 1e-9)
        stats["eta_max_deg"] = float(eta_max_deg)
        stats["saturation_time_fraction"] = float(np.mean(sat_mask))
    if eta_rate_max_deg_s is not None:
        stats["eta_rate_max_deg_s"] = float(eta_rate_max_deg_s)
    if eta_accel_max_deg_s2 is not None:
        stats["eta_accel_max_deg_s2"] = float(eta_accel_max_deg_s2)

    if progress_axis_deg.ndim == 1 and progress_axis_deg.size > 0:
        prog_abs = np.abs(progress_axis_deg)
        if half_turn_deg is not None:
            stats["half_turn_reached"] = bool(np.any(prog_abs >= half_turn_deg))
        if full_turn_deg is not None:
            reached_full = np.where(prog_abs >= full_turn_deg)[0]
            stats["full_turn_reached"] = bool(reached_full.size > 0)
            if reached_full.size > 0:
                stats["time_to_full_turn_s"] = float(t_grid[int(reached_full[0])])

    metric_180 = core_outputs.get(f"time_to_{axis}_180_deg", {}) if axis else {}
    frac_180 = _parse_reached_fraction(metric_180.get("note"))
    if frac_180 is not None:
        stats["reached_fraction_180"] = frac_180

    return stats


def write_mc_case_outputs(
    *,
    case_dir: Path,
    summary_name: str,
    scenario_type: str,
    requested_metrics: list[str],
    t_grid: np.ndarray,
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    mc_cov_final: np.ndarray,
    n_total: int,
    do_plot: bool,
    plot_all_states: bool,
    q_particles: np.ndarray | None = None,
    w_particles: np.ndarray | None = None,
    maneuver_diag: dict | None = None,
    scenario: dict | None = None,
    seed: int | None = None,
) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    t_eff, mean_eff, std_eff, cov_eff, analysis_info = prepare_mc_series_for_analysis(
        t_grid=t_grid,
        mc_mean=mc_mean,
        mc_std=mc_std,
        mc_cov_final=mc_cov_final,
    )
    case_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        case_dir / "mc_only.npz",
        t_grid=t_eff,
        mc_mean=mean_eff,
        mc_std=std_eff,
        mc_cov_final=cov_eff,
        n_total=n_total,
        raw_steps=analysis_info.get("raw_steps", int(t_grid.size)),
        analysis_steps=analysis_info.get("analysis_steps", int(t_eff.size)),
    )
    summary = summarize_mc_state(
        name=summary_name,
        scenario_type=scenario_type,
        requested_metrics=requested_metrics,
        t_grid=t_eff,
        mc_mean=mean_eff,
        mc_std=std_eff,
        mc_cov_final=cov_eff,
        n_total=n_total,
        analysis_info=analysis_info,
    )
    q_eff = None
    w_eff = None
    if isinstance(q_particles, np.ndarray) and q_particles.ndim == 3 and q_particles.shape[0] >= t_eff.size:
        q_eff = q_particles[: t_eff.size]
    if isinstance(w_particles, np.ndarray) and w_particles.ndim == 3 and w_particles.shape[0] >= t_eff.size:
        w_eff = w_particles[: t_eff.size]

    summary["catalog_core_outputs"] = build_catalog_core_outputs(
        requested_metrics=requested_metrics,
        t_grid=t_eff,
        mc_mean=mean_eff,
        mc_std=std_eff,
        q_particles=q_eff,
        w_particles=w_eff,
    )
    if isinstance(q_eff, np.ndarray):
        summary["attitude_metric_source"] = "per_particle_quaternions"
        summary["attitude_particles_used"] = int(q_eff.shape[1])
    if isinstance(maneuver_diag, dict):
        diag_save = {}
        for key, value in maneuver_diag.items():
            arr = np.asarray(value)
            if arr.ndim >= 1 and arr.shape[0] >= t_eff.size:
                diag_save[key] = arr[: t_eff.size]
            else:
                diag_save[key] = arr
        np.savez_compressed(case_dir / "maneuver_torque_diagnostics.npz", **diag_save)
        rep_tau = np.asarray(diag_save.get("rep_torque_B_Nm"))
        det_tau = np.asarray(diag_save.get("det_torque_B_Nm")) if "det_torque_B_Nm" in diag_save else None
        rep_norm = np.linalg.norm(rep_tau, axis=1) if rep_tau.ndim == 2 and rep_tau.shape[1] == 3 else np.array([])
        det_norm = np.linalg.norm(det_tau, axis=1) if isinstance(det_tau, np.ndarray) and det_tau.shape == rep_tau.shape else np.array([])
        summary["maneuver_torque_diagnostics"] = {
            "available": True,
            "max_rep_torque_norm_Nm": float(np.nanmax(rep_norm)) if rep_norm.size > 0 else None,
            "max_det_torque_norm_Nm": float(np.nanmax(det_norm)) if det_norm.size > 0 else None,
            "source": "runtime diagnostic (Euler rigid-body torque reconstruction)",
        }
        summary["maneuver_actuator_stats"] = build_maneuver_actuator_stats(
            diag=diag_save,
            t_grid=t_eff,
            core_outputs=summary.get("catalog_core_outputs", {}),
        )
    summary["metadata"] = build_run_metadata(
        scenario if isinstance(scenario, dict) else {"name": summary_name},
        seed,
    )
    save_json(case_dir / "summary.json", summary)
    if do_plot:
        plot_mc_case(case_dir, t_eff, mean_eff, std_eff, plot_all_states, q_particles=q_eff)
        if isinstance(maneuver_diag, dict):
            plot_maneuver_torque(case_dir, t_eff, maneuver_diag)
    return summary, t_eff, mean_eff, std_eff


def main() -> None:
    parser = argparse.ArgumentParser(description="Postprocess MC-only outputs for one scenario.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--plot_all_states", action="store_true")
    args = parser.parse_args()

    raw = json.loads(Path(args.config).read_text(encoding="utf-8"))
    scenarios = raw.get("scenarios", [])
    scenario = next((s for s in scenarios if s.get("name") == args.scenario), None)
    if scenario is None:
        raise ValueError(f"scenario '{args.scenario}' not found in config")

    name = str(scenario["name"])
    scenario_type = str(scenario.get("type", "unknown"))
    requested_metrics = [str(m) for m in scenario.get("catalog_output_metrics", [])]
    outdir = Path(args.outdir)
    case_dir = outdir / name
    case_dir.mkdir(parents=True, exist_ok=True)

    if scenario_type == "formation":
        offsets = scenario.get("formation_offsets_m", [[0.0, 0.0, 0.0]])
        if not isinstance(offsets, list) or not offsets:
            offsets = [[0.0, 0.0, 0.0]]
        object_summaries = []
        object_means = []
        object_stds = []
        object_t_grids = []

        for idx in range(len(offsets)):
            obj = f"obj_{idx+1:02d}"
            stats_path = outdir / name / f"{obj}_mc_stats.npz"
            if not stats_path.exists():
                raise FileNotFoundError(stats_path)
            t_grid, mc_mean, mc_std, mc_cov_final, n_total = load_mc_stats(stats_path)
            obj_dir = case_dir / obj
            obj_summary, t_eff, mean_eff, std_eff = write_mc_case_outputs(
                case_dir=obj_dir,
                summary_name=f"{name}/{obj}",
                scenario_type=scenario_type,
                requested_metrics=requested_metrics,
                t_grid=t_grid,
                mc_mean=mc_mean,
                mc_std=mc_std,
                mc_cov_final=mc_cov_final,
                n_total=n_total,
                do_plot=args.plot,
                plot_all_states=args.plot_all_states,
                scenario=scenario,
                seed=args.seed,
            )
            object_summaries.append(obj_summary)
            object_t_grids.append(t_eff)
            object_means.append(mean_eff)
            object_stds.append(std_eff)

        scenario_summary = {
            "case": name,
            "mode": "mc_only",
            "type": scenario_type,
            "requested_catalog_metrics": requested_metrics,
            "num_objects": int(len(object_summaries)),
            "object_summaries": object_summaries,
        }

        sep_arr = None
        t_ref = object_t_grids[0] if object_t_grids else np.array([0.0], dtype=float)
        if object_t_grids:
            n_common = min(int(t.size) for t in object_t_grids)
            n_common = max(1, n_common)
            t_ref = object_t_grids[0][:n_common]
            object_means = [m[:n_common] for m in object_means]
            object_stds = [s[:n_common] for s in object_stds]
            scenario_summary["analysis_common_steps"] = int(n_common)
            scenario_summary["analysis_common_duration_s"] = float(t_ref[-1] - t_ref[0]) if t_ref.size > 1 else 0.0

        if len(object_means) > 1:
            chief = object_means[0][:, 0:3]
            sep = []
            for follower in object_means[1:]:
                sep.append(np.linalg.norm(follower[:, 0:3] - chief, axis=1))
            sep_arr = np.vstack(sep)
            np.savez(case_dir / "formation_mc_metrics.npz", t_grid=t_ref, mc_mean_separation_m=sep_arr)
            scenario_summary["formation_metrics"] = {
                "separation_final_m": [float(x) for x in sep_arr[:, -1]],
                "separation_min_m": [float(x) for x in np.nanmin(sep_arr, axis=1)],
                "separation_max_m": [float(x) for x in np.nanmax(sep_arr, axis=1)],
            }
            if args.plot:
                plot_formation_separation(case_dir, t_ref, sep_arr)

        scenario_summary["catalog_core_outputs"] = build_catalog_core_outputs(
            requested_metrics=requested_metrics,
            t_grid=t_ref,
            mc_mean=object_means[0] if object_means else np.zeros((1, 18), dtype=float),
            mc_std=object_stds[0] if object_stds else np.zeros((1, 18), dtype=float),
            separations_m=sep_arr,
        )
        scenario_summary["metadata"] = build_run_metadata(scenario, args.seed)

        save_json(case_dir / "summary.json", scenario_summary)
        print(f"[mc-only] wrote {case_dir}/summary.json", flush=True)
        return

    mc_stats_path = outdir / f"{name}_mc_stats.npz"
    if not mc_stats_path.exists():
        raise FileNotFoundError(mc_stats_path)

    t_grid, mc_mean, mc_std, mc_cov_final, n_total = load_mc_stats(mc_stats_path)
    shard_dir = outdir / "mc_shards" / name
    need_attitude_particles = bool(
        scenario_type == "attitude"
        or any(
            re.fullmatch(r"(time_to|overshoot|peak|rms)_(roll|pitch|yaw).*", m)
            for m in requested_metrics
        )
    )
    particle_series = load_attitude_particles_from_shards(shard_dir, t_grid) if need_attitude_particles else None
    maneuver_diag = load_maneuver_diagnostics(shard_dir, t_grid) if scenario_type == "attitude" else None
    write_mc_case_outputs(
        case_dir=case_dir,
        summary_name=name,
        scenario_type=scenario_type,
        requested_metrics=requested_metrics,
        t_grid=t_grid,
        mc_mean=mc_mean,
        mc_std=mc_std,
        mc_cov_final=mc_cov_final,
        n_total=n_total,
        do_plot=args.plot,
        plot_all_states=args.plot_all_states,
        q_particles=particle_series["q"] if isinstance(particle_series, dict) else None,
        w_particles=particle_series["w"] if isinstance(particle_series, dict) else None,
        maneuver_diag=maneuver_diag,
        scenario=scenario,
        seed=args.seed,
    )
    print(f"[mc-only] wrote {case_dir}/summary.json", flush=True)


if __name__ == "__main__":
    main()
