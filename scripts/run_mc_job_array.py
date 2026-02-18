from __future__ import annotations

import argparse
import importlib
import os
import time
from pathlib import Path
from typing import Callable, Dict, Tuple

import numpy as np


def parse_scenario(spec: str) -> Callable[[], Dict[str, object]]:
    if ":" not in spec:
        raise ValueError("scenario must be in module:function form")
    module_name, func_name = spec.split(":", 1)
    module = importlib.import_module(module_name)
    func = getattr(module, func_name)
    if not callable(func):
        raise ValueError("scenario is not callable")
    return func


def shard_from_env() -> Tuple[int, int]:
    if "SLURM_ARRAY_TASK_ID" in os.environ and "SLURM_ARRAY_TASK_COUNT" in os.environ:
        return int(os.environ["SLURM_ARRAY_TASK_ID"]), int(os.environ["SLURM_ARRAY_TASK_COUNT"])
    if "SLURM_PROCID" in os.environ and "SLURM_NTASKS" in os.environ:
        return int(os.environ["SLURM_PROCID"]), int(os.environ["SLURM_NTASKS"])
    if "OMPI_COMM_WORLD_RANK" in os.environ and "OMPI_COMM_WORLD_SIZE" in os.environ:
        return int(os.environ["OMPI_COMM_WORLD_RANK"]), int(os.environ["OMPI_COMM_WORLD_SIZE"])
    return -1, -1


def propagate_chunked(
    propagator,
    x0: np.ndarray,
    t_grid: np.ndarray,
    env,
    chunk_steps: int,
    tag: str,
) -> np.ndarray:
    nt = len(t_grid)
    if nt == 0:
        raise ValueError("t_grid must have at least 1 entry")
    if nt == 1:
        out = np.zeros((1, x0.shape[0], x0.shape[1]), dtype=float)
        out[0] = x0
        return out

    out = np.zeros((nt, x0.shape[0], x0.shape[1]), dtype=float)
    x_curr = x0
    i0 = 0
    while i0 < nt - 1:
        i1 = min(nt - 1, i0 + chunk_steps)
        t_chunk = t_grid[i0 : i1 + 1] - t_grid[i0]
        env_chunk = env[i0 : i1 + 1]
        t_chunk_start = time.time()
        out_chunk = propagator.propagate(x_curr, t_chunk, env_chunk)
        out[i0 : i1 + 1] = out_chunk
        x_curr = out_chunk[-1]
        elapsed = time.time() - t_chunk_start
        print(f"{tag} step {i1+1}/{nt} chunk_s={elapsed:.1f}")
        i0 = i1
    return out


def _as_finite_float(value) -> float | None:
    try:
        out = float(value)
    except Exception:
        return None
    if not np.isfinite(out):
        return None
    return out


def _axis_index(axis: str) -> int:
    key = str(axis).strip().lower()
    return {"roll": 0, "pitch": 1, "yaw": 2}[key]


def _maneuver_control_mode(scenario_cfg: dict) -> str:
    explicit = str(scenario_cfg.get("maneuver_control_mode", "")).strip().lower()
    if explicit in {"point", "rate", "open_loop"}:
        return explicit
    toggles = scenario_cfg.get("catalog_control_toggles")
    if isinstance(toggles, (list, tuple, set)):
        keys = {str(t).strip().upper() for t in toggles}
        if "CTL_ATT_AERO_RATE" in keys:
            return "rate"
        if "CTL_ATT_AERO_POINT" in keys:
            return "point"
    return "open_loop"


def _quat_to_euler_zyx_deg(q_wxyz: np.ndarray) -> np.ndarray:
    q = np.asarray(q_wxyz, dtype=float).reshape(-1)
    if q.size != 4:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    n = np.linalg.norm(q)
    if not np.isfinite(n) or n <= 0.0:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    w, x, y, z = q / n
    sinr_cosp = 2.0 * (w * x + y * z)
    cosr_cosp = 1.0 - 2.0 * (x * x + y * y)
    roll = np.arctan2(sinr_cosp, cosr_cosp)
    sinp = np.clip(2.0 * (w * y - z * x), -1.0, 1.0)
    pitch = np.arcsin(sinp)
    siny_cosp = 2.0 * (w * z + x * y)
    cosy_cosp = 1.0 - 2.0 * (y * y + z * z)
    yaw = np.arctan2(siny_cosp, cosy_cosp)
    return np.degrees(np.array([roll, pitch, yaw], dtype=float))


def _wrapped_delta_deg(curr: float, prev: float) -> float:
    d = curr - prev
    return ((d + 180.0) % 360.0) - 180.0


def _axis_eta_cmd(axis: str, wing_deg: float) -> tuple[float, float]:
    axis_key = str(axis).strip().lower()
    d = np.deg2rad(float(wing_deg))
    if axis_key == "roll":
        return float(d), float(-d)
    if axis_key == "pitch":
        return float(d), float(d)
    if axis_key == "yaw":
        return float(d), 0.0
    raise ValueError(f"unsupported maneuver axis '{axis}'")


def _inertia_matrix_from_scenario(scenario_cfg: dict) -> np.ndarray:
    inertia = scenario_cfg.get("spacecraft_inertia_kgm2")
    if isinstance(inertia, (list, tuple)):
        arr = np.asarray(inertia, dtype=float).reshape(-1)
        if arr.size == 3 and np.isfinite(arr).all():
            return np.diag(arr)
        if arr.size == 9 and np.isfinite(arr).all():
            return arr.reshape(3, 3)
    # Baseline fallback used across ATT scenarios in the catalog.
    return np.diag(np.array([0.15, 0.12, 0.20], dtype=float))


def _extract_maneuver_runtime_config(scenario_cfg: dict) -> dict | None:
    def cfg_float(*keys: str) -> float | None:
        for key in keys:
            v = _as_finite_float(scenario_cfg.get(key))
            if v is not None:
                return float(v)
        return None

    def cfg_or_env(default: float, env_key: str, *keys: str) -> float:
        v = cfg_float(*keys)
        if v is not None:
            return float(v)
        e = _as_finite_float(os.environ.get(env_key))
        if e is not None:
            return float(e)
        return float(default)

    axis = str(scenario_cfg.get("maneuver_axis", "")).strip().lower()
    wing_const = _as_finite_float(scenario_cfg.get("wing_constant_deg"))
    if axis not in {"roll", "pitch", "yaw"} or wing_const is None or wing_const == 0.0:
        return None
    half_turn = _as_finite_float(scenario_cfg.get("maneuver_half_turn_deg"))
    full_turn = _as_finite_float(scenario_cfg.get("maneuver_full_turn_deg"))
    if half_turn is None:
        half_turn = 180.0
    if full_turn is None:
        full_turn = 360.0
    half_turn = max(0.0, float(half_turn))
    full_turn = max(float(full_turn), half_turn + 1e-6)
    eta_max_deg = cfg_or_env(110.0, "VLEO_ETA_MAX_DEG", "eta_max_deg", "eta_limit_deg")
    eta_limit_rad = cfg_float("eta_limit_rad")
    if eta_limit_rad is not None:
        eta_max_deg = float(np.degrees(abs(eta_limit_rad)))
    eta_rate_max_deg_s = cfg_or_env(5.0, "VLEO_ETA_RATE_MAX_DEG_S", "eta_rate_max_deg_s")
    eta_accel_max_deg_s2 = cfg_or_env(
        1.0,
        "VLEO_ETA_ACCEL_MAX_DEG_S2",
        "eta_accel_max_deg_s2",
        "eta_accel_max_deg_s2_1deg",
    )
    control_mode = _maneuver_control_mode(scenario_cfg)
    point_kp_cmd_per_deg = cfg_or_env(
        0.5,
        "VLEO_ATT_POINT_KP_CMD_PER_DEG",
        "att_point_kp_cmd_per_deg",
        "maneuver_point_kp_cmd_per_deg",
    )
    point_kd_cmd_per_deg_s = cfg_or_env(
        2.0,
        "VLEO_ATT_POINT_KD_CMD_PER_DEG_S",
        "att_point_kd_cmd_per_deg_s",
        "maneuver_point_kd_cmd_per_deg_s",
    )
    rate_target_deg_s = cfg_or_env(
        1.0,
        "VLEO_ATT_RATE_TARGET_DEG_S",
        "att_rate_target_deg_s",
        "maneuver_rate_target_deg_s",
    )
    rate_kp_cmd_per_deg_s = cfg_or_env(
        12.0,
        "VLEO_ATT_RATE_KP_CMD_PER_DEG_S",
        "att_rate_kp_cmd_per_deg_s",
        "maneuver_rate_kp_cmd_per_deg_s",
    )
    return {
        "axis": axis,
        "control_mode": control_mode,
        "wing_constant_deg": float(wing_const),
        "half_turn_deg": half_turn,
        "full_turn_deg": full_turn,
        "eta_max_deg": max(0.0, float(abs(eta_max_deg))),
        "eta_rate_max_deg_s": max(0.0, float(abs(eta_rate_max_deg_s))),
        "eta_accel_max_deg_s2": max(0.0, float(abs(eta_accel_max_deg_s2))),
        "point_kp_cmd_per_deg": max(0.0, float(abs(point_kp_cmd_per_deg))),
        "point_kd_cmd_per_deg_s": max(0.0, float(abs(point_kd_cmd_per_deg_s))),
        "rate_target_deg_s": max(0.0, float(abs(rate_target_deg_s))),
        "rate_kp_cmd_per_deg_s": max(0.0, float(abs(rate_kp_cmd_per_deg_s))),
    }


def _propagate_member_with_axis_runtime(
    det_propagator,
    x0: np.ndarray,
    t_grid: np.ndarray,
    env,
    *,
    axis: str,
    control_mode: str,
    wing_constant_deg: float,
    half_turn_deg: float,
    full_turn_deg: float,
    eta_max_deg: float,
    eta_rate_max_deg_s: float,
    eta_accel_max_deg_s2: float,
    point_kp_cmd_per_deg: float,
    point_kd_cmd_per_deg_s: float,
    rate_target_deg_s: float,
    rate_kp_cmd_per_deg_s: float,
    inertia_B: np.ndarray,
) -> tuple[np.ndarray, dict]:
    nt = int(t_grid.size)
    nx = int(x0.size)
    if nt < 2:
        out = np.zeros((nt, nx), dtype=float)
        if nt == 1:
            out[0] = x0
        diag = {
            "axis": axis,
            "control_mode": str(control_mode),
            "half_turn_deg": float(half_turn_deg),
            "full_turn_deg": float(full_turn_deg),
            "progress_axis_deg": np.zeros(nt, dtype=float),
            "eta1_cmd_rad": np.zeros(nt, dtype=float),
            "eta2_cmd_rad": np.zeros(nt, dtype=float),
            "eta_cmd_deg": np.zeros(nt, dtype=float),
            "eta_target_deg": np.zeros(nt, dtype=float),
            "eta_rate_deg_s": np.zeros(nt, dtype=float),
            "axis_rate_deg_s": np.zeros(nt, dtype=float),
            "torque_B_Nm": np.full((nt, 3), np.nan, dtype=float),
            "phase_index": np.zeros(nt, dtype=int),
            "eta_max_deg": float(eta_max_deg),
            "eta_rate_max_deg_s": float(eta_rate_max_deg_s),
            "eta_accel_max_deg_s2": float(eta_accel_max_deg_s2),
            "point_kp_cmd_per_deg": float(point_kp_cmd_per_deg),
            "point_kd_cmd_per_deg_s": float(point_kd_cmd_per_deg_s),
            "rate_target_deg_s": float(rate_target_deg_s),
            "rate_kp_cmd_per_deg_s": float(rate_kp_cmd_per_deg_s),
        }
        return out, diag

    out = np.zeros((nt, nx), dtype=float)
    out[0] = x0

    control_mode = str(control_mode).strip().lower()
    if control_mode not in {"point", "rate", "open_loop"}:
        control_mode = "open_loop"
    axis_idx = _axis_index(axis)
    wing_sign = 1.0 if wing_constant_deg >= 0.0 else -1.0
    eta_max_deg = max(0.0, float(eta_max_deg))
    eta_rate_max_deg_s = max(0.0, float(eta_rate_max_deg_s))
    eta_accel_max_deg_s2 = max(0.0, float(eta_accel_max_deg_s2))
    point_kp_cmd_per_deg = max(0.0, float(point_kp_cmd_per_deg))
    point_kd_cmd_per_deg_s = max(0.0, float(point_kd_cmd_per_deg_s))
    rate_target_deg_s = max(0.0, float(rate_target_deg_s))
    rate_kp_cmd_per_deg_s = max(0.0, float(rate_kp_cmd_per_deg_s))
    wing_mag = min(abs(float(wing_constant_deg)), eta_max_deg if eta_max_deg > 0.0 else abs(float(wing_constant_deg)))

    progress_deg = np.zeros(nt, dtype=float)
    eta1_cmd = np.zeros(nt, dtype=float)
    eta2_cmd = np.zeros(nt, dtype=float)
    eta_cmd_deg = np.zeros(nt, dtype=float)
    eta_target_deg = np.zeros(nt, dtype=float)
    eta_rate_deg_s = np.zeros(nt, dtype=float)
    axis_rate_deg_s = np.zeros(nt, dtype=float)
    torque_B = np.full((nt, 3), np.nan, dtype=float)
    phase_index = np.zeros(nt, dtype=int)

    prev_angle = float(_quat_to_euler_zyx_deg(out[0, 6:10])[axis_idx])
    cmd_deg_state = 0.0
    cmd_rate_deg_s_state = 0.0
    phase = 0

    def _target_cmd_deg(prog_axis_deg: float, axis_rate_cmd_deg_s: float, phase_idx: int) -> float:
        signed_prog = wing_sign * prog_axis_deg
        signed_rate = wing_sign * axis_rate_cmd_deg_s
        if control_mode == "point":
            if phase_idx == 0:
                point_target_deg = half_turn_deg
            else:
                point_target_deg = full_turn_deg
            err_point_deg = point_target_deg - signed_prog
            cmd_signed_deg = point_kp_cmd_per_deg * err_point_deg - point_kd_cmd_per_deg_s * signed_rate
            return float(wing_sign * np.clip(cmd_signed_deg, -wing_mag, wing_mag))
        if control_mode == "rate":
            target_rate_deg_s = rate_target_deg_s if phase_idx == 0 else 0.0
            rate_err_deg_s = target_rate_deg_s - signed_rate
            cmd_signed_deg = rate_kp_cmd_per_deg_s * rate_err_deg_s
            return float(wing_sign * np.clip(cmd_signed_deg, -wing_mag, wing_mag))
        if phase_idx == 0:
            return float(wing_sign * wing_mag)
        if phase_idx == 1:
            return float(-wing_sign * wing_mag)
        return 0.0

    for i in range(nt - 1):
        dt_step = float(t_grid[i + 1] - t_grid[i])
        if dt_step <= 0.0:
            raise RuntimeError("t_grid must be strictly increasing")

        curr_axis_rate_deg_s = 0.0
        if out.shape[1] >= 13:
            curr_axis_rate_deg_s = float(np.degrees(out[i, 10 + axis_idx]))
        axis_rate_deg_s[i] = curr_axis_rate_deg_s
        target_deg = _target_cmd_deg(progress_deg[i], curr_axis_rate_deg_s, phase)

        cmd_deg_prev = cmd_deg_state
        if eta_rate_max_deg_s > 0.0 and eta_accel_max_deg_s2 > 0.0:
            err_deg = target_deg - cmd_deg_state
            des_rate_deg_s = np.clip(err_deg / dt_step, -eta_rate_max_deg_s, eta_rate_max_deg_s)
            d_rate_deg_s = des_rate_deg_s - cmd_rate_deg_s_state
            max_d_rate_deg_s = eta_accel_max_deg_s2 * dt_step
            d_rate_deg_s = float(np.clip(d_rate_deg_s, -max_d_rate_deg_s, max_d_rate_deg_s))
            cmd_rate_deg_s_state = float(np.clip(cmd_rate_deg_s_state + d_rate_deg_s, -eta_rate_max_deg_s, eta_rate_max_deg_s))
            cmd_deg_state = float(np.clip(cmd_deg_state + cmd_rate_deg_s_state * dt_step, -eta_max_deg, eta_max_deg))
            # Avoid numerical overshoot around target crossings.
            if (target_deg - cmd_deg_prev) * (target_deg - cmd_deg_state) <= 0.0:
                cmd_deg_state = float(target_deg)
                cmd_rate_deg_s_state = 0.0
        else:
            cmd_deg_state = float(np.clip(target_deg, -eta_max_deg, eta_max_deg))
            cmd_rate_deg_s_state = 0.0

        eta1, eta2 = _axis_eta_cmd(axis, cmd_deg_state)
        eta1_cmd[i] = eta1
        eta2_cmd[i] = eta2
        eta_cmd_deg[i] = cmd_deg_state
        eta_target_deg[i] = target_deg
        eta_rate_deg_s[i] = cmd_rate_deg_s_state
        phase_index[i] = phase

        env[i].eta1_rad = float(eta1)
        env[i].eta2_rad = float(eta2)
        env[i + 1].eta1_rad = float(eta1)
        env[i + 1].eta2_rad = float(eta2)

        out_step = det_propagator.propagate(out[i], np.array([0.0, dt_step], dtype=float), env[i : i + 2])
        out[i + 1] = out_step[-1]

        axis_angle = float(_quat_to_euler_zyx_deg(out[i + 1, 6:10])[axis_idx])
        if np.isfinite(axis_angle) and np.isfinite(prev_angle):
            progress_deg[i + 1] = progress_deg[i] + _wrapped_delta_deg(axis_angle, prev_angle)
            prev_angle = axis_angle
        else:
            progress_deg[i + 1] = progress_deg[i]

        if out.shape[1] >= 13:
            w_prev = out[i, 10:13]
            w_next = out[i + 1, 10:13]
            wdot = (w_next - w_prev) / dt_step
            torque_B[i + 1] = inertia_B @ wdot + np.cross(w_prev, inertia_B @ w_prev)

        signed_prog = wing_sign * progress_deg[i + 1]
        if control_mode == "point":
            if phase == 0 and signed_prog >= half_turn_deg:
                phase = 1
            if phase == 1 and signed_prog >= full_turn_deg:
                phase = 2
        elif control_mode == "rate":
            if phase == 0 and signed_prog >= full_turn_deg:
                phase = 1
        else:
            prog_abs = abs(progress_deg[i + 1])
            if phase == 0 and prog_abs >= half_turn_deg:
                phase = 1
            if phase == 1 and prog_abs >= full_turn_deg:
                phase = 2

    eta1_cmd[-1] = eta1_cmd[-2]
    eta2_cmd[-1] = eta2_cmd[-2]
    eta_cmd_deg[-1] = eta_cmd_deg[-2]
    eta_target_deg[-1] = eta_target_deg[-2]
    eta_rate_deg_s[-1] = eta_rate_deg_s[-2]
    axis_rate_deg_s[-1] = axis_rate_deg_s[-2]
    phase_index[-1] = phase

    diag = {
        "axis": axis,
        "control_mode": str(control_mode),
        "half_turn_deg": float(half_turn_deg),
        "full_turn_deg": float(full_turn_deg),
        "progress_axis_deg": progress_deg,
        "eta1_cmd_rad": eta1_cmd,
        "eta2_cmd_rad": eta2_cmd,
        "eta_cmd_deg": eta_cmd_deg,
        "eta_target_deg": eta_target_deg,
        "eta_rate_deg_s": eta_rate_deg_s,
        "axis_rate_deg_s": axis_rate_deg_s,
        "torque_B_Nm": torque_B,
        "phase_index": phase_index,
        "eta_max_deg": float(eta_max_deg),
        "eta_rate_max_deg_s": float(eta_rate_max_deg_s),
        "eta_accel_max_deg_s2": float(eta_accel_max_deg_s2),
        "point_kp_cmd_per_deg": float(point_kp_cmd_per_deg),
        "point_kd_cmd_per_deg_s": float(point_kd_cmd_per_deg_s),
        "rate_target_deg_s": float(rate_target_deg_s),
        "rate_kp_cmd_per_deg_s": float(rate_kp_cmd_per_deg_s),
    }
    return out, diag


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MC ensemble shard for job arrays.")
    parser.add_argument("--scenario", required=True, help="module:function returning scenario dict")
    parser.add_argument("--outdir", default="results/mc_shards", help="output directory")
    parser.add_argument("--shard", type=int, default=None, help="shard index (0-based)")
    parser.add_argument("--shards", type=int, default=None, help="total shards")
    parser.add_argument("--chunk_steps", type=int, default=0, help="chunk size in steps for progress")
    args = parser.parse_args()

    shard, shards = shard_from_env()
    if args.shard is not None:
        shard = args.shard
    if args.shards is not None:
        shards = args.shards
    if shard < 0 or shards <= 0:
        raise RuntimeError("provide shard/shards or set SLURM_ARRAY_TASK_ID/COUNT")

    scenario = parse_scenario(args.scenario)()
    for key in ("propagator", "t_grid", "env", "X0"):
        if key not in scenario:
            raise RuntimeError(f"scenario missing key '{key}'")

    propagator = scenario["propagator"]
    det_propagator = scenario.get("det_propagator")
    t_grid = np.asarray(scenario["t_grid"], dtype=float)
    env = scenario["env"]
    X0 = np.asarray(scenario["X0"], dtype=float)
    scenario_cfg = scenario.get("scenario_config") if isinstance(scenario.get("scenario_config"), dict) else {}
    maneuver_runtime_cfg = _extract_maneuver_runtime_config(scenario_cfg)
    inertia_B = _inertia_matrix_from_scenario(scenario_cfg) if maneuver_runtime_cfg is not None else None
    safe_mode = os.environ.get("VLEO_MC_SAFE_MODE", "0").strip().lower() in {"1", "true", "yes", "on"}

    n_total = X0.shape[0]
    start = (n_total * shard) // shards
    end = (n_total * (shard + 1)) // shards
    if start >= end:
        print(f"[mc] rank {shard}/{shards} empty slice (members={n_total}), writing empty shard")
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        out_path = outdir / f"ensemble_rank{shard:04d}.npz"
        np.savez_compressed(
            out_path,
            X=np.zeros((len(t_grid), 0, X0.shape[1]), dtype=float),
            t_grid=t_grid,
            start=start,
            end=end,
            shard=shard,
            shards=shards,
            failed=0,
            failed_members=np.array([], dtype=int),
        )
        print(f"Saved empty shard {shard}/{shards} to {out_path}")
        return

    X_slice = X0[start:end]
    t0 = time.time()
    print(
        f"[mc] rank {shard}/{shards} start={start} end={end} "
        f"members={end-start} steps={len(t_grid)}"
    )

    chunk_steps = int(args.chunk_steps)
    failed_members: list[int] = []
    representative_diag: dict | None = None
    representative_local_idx: int | None = None
    representative_global_idx: int | None = None
    deterministic_diag: dict | None = None
    if safe_mode and maneuver_runtime_cfg is not None and det_propagator is not None and inertia_B is not None:
        x0_nominal = np.asarray(scenario.get("x0_nominal", X0[0]), dtype=float).reshape(-1)
        try:
            _, deterministic_diag = _propagate_member_with_axis_runtime(
                det_propagator,
                x0_nominal,
                t_grid,
                env,
                axis=str(maneuver_runtime_cfg["axis"]),
                control_mode=str(maneuver_runtime_cfg["control_mode"]),
                wing_constant_deg=float(maneuver_runtime_cfg["wing_constant_deg"]),
                half_turn_deg=float(maneuver_runtime_cfg["half_turn_deg"]),
                full_turn_deg=float(maneuver_runtime_cfg["full_turn_deg"]),
                eta_max_deg=float(maneuver_runtime_cfg["eta_max_deg"]),
                eta_rate_max_deg_s=float(maneuver_runtime_cfg["eta_rate_max_deg_s"]),
                eta_accel_max_deg_s2=float(maneuver_runtime_cfg["eta_accel_max_deg_s2"]),
                point_kp_cmd_per_deg=float(maneuver_runtime_cfg["point_kp_cmd_per_deg"]),
                point_kd_cmd_per_deg_s=float(maneuver_runtime_cfg["point_kd_cmd_per_deg_s"]),
                rate_target_deg_s=float(maneuver_runtime_cfg["rate_target_deg_s"]),
                rate_kp_cmd_per_deg_s=float(maneuver_runtime_cfg["rate_kp_cmd_per_deg_s"]),
                inertia_B=inertia_B,
            )
        except Exception as det_exc:
            print(
                f"[mc] rank {shard}/{shards} deterministic maneuver diagnostic failed "
                f"({type(det_exc).__name__}: {det_exc})"
            )
    if safe_mode:
        if det_propagator is None:
            raise RuntimeError("VLEO_MC_SAFE_MODE=1 requires scenario key 'det_propagator'")
        print(f"[mc] rank {shard}/{shards} safe-mode enabled (deterministic per-member propagation)")
        kept = []
        for j in range(X_slice.shape[0]):
            try:
                xj0 = X_slice[j]
                if maneuver_runtime_cfg is not None and inertia_B is not None:
                    out_j, diag_j = _propagate_member_with_axis_runtime(
                        det_propagator,
                        xj0,
                        t_grid,
                        env,
                        axis=str(maneuver_runtime_cfg["axis"]),
                        control_mode=str(maneuver_runtime_cfg["control_mode"]),
                        wing_constant_deg=float(maneuver_runtime_cfg["wing_constant_deg"]),
                        half_turn_deg=float(maneuver_runtime_cfg["half_turn_deg"]),
                        full_turn_deg=float(maneuver_runtime_cfg["full_turn_deg"]),
                        eta_max_deg=float(maneuver_runtime_cfg["eta_max_deg"]),
                        eta_rate_max_deg_s=float(maneuver_runtime_cfg["eta_rate_max_deg_s"]),
                        eta_accel_max_deg_s2=float(maneuver_runtime_cfg["eta_accel_max_deg_s2"]),
                        point_kp_cmd_per_deg=float(maneuver_runtime_cfg["point_kp_cmd_per_deg"]),
                        point_kd_cmd_per_deg_s=float(maneuver_runtime_cfg["point_kd_cmd_per_deg_s"]),
                        rate_target_deg_s=float(maneuver_runtime_cfg["rate_target_deg_s"]),
                        rate_kp_cmd_per_deg_s=float(maneuver_runtime_cfg["rate_kp_cmd_per_deg_s"]),
                        inertia_B=inertia_B,
                    )
                    if representative_diag is None:
                        representative_diag = diag_j
                        representative_local_idx = int(j)
                        representative_global_idx = int(start + j)
                else:
                    out_j = det_propagator.propagate(xj0, t_grid, env)
                kept.append(out_j)
            except Exception as member_exc:
                failed_members.append(j)
                print(f"[mc] rank {shard}/{shards} dropped member local={j} ({type(member_exc).__name__}: {member_exc})")
        if kept:
            out = np.stack(kept, axis=1)
        else:
            out = np.zeros((len(t_grid), 0, X_slice.shape[1]), dtype=float)
    else:
        try:
            if chunk_steps > 0:
                out = propagate_chunked(
                    propagator,
                    X_slice,
                    t_grid,
                    env,
                    chunk_steps,
                    f"[mc] rank {shard}/{shards}",
                )
            else:
                out = propagator.propagate(X_slice, t_grid, env)
        except Exception as exc:
            # Keep the shard alive even if a subset of particles diverges.
            print(f"[mc] rank {shard}/{shards} bulk propagation failed ({type(exc).__name__}: {exc}), retrying member-wise")
            kept = []
            for j in range(X_slice.shape[0]):
                try:
                    Xj = X_slice[j : j + 1]
                    if chunk_steps > 0:
                        out_j = propagate_chunked(
                            propagator,
                            Xj,
                            t_grid,
                            env,
                            chunk_steps,
                            f"[mc] rank {shard}/{shards} member {j}",
                        )
                    else:
                        out_j = propagator.propagate(Xj, t_grid, env)
                    kept.append(out_j[:, 0, :])
                except Exception as member_exc:
                    failed_members.append(j)
                    print(f"[mc] rank {shard}/{shards} dropped member local={j} ({type(member_exc).__name__}: {member_exc})")
            if kept:
                out = np.stack(kept, axis=1)
            else:
                out = np.zeros((len(t_grid), 0, X_slice.shape[1]), dtype=float)

    elapsed = time.time() - t0
    print(f"[mc] rank {shard}/{shards} done in {elapsed:.1f}s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"ensemble_rank{shard:04d}.npz"
    np.savez_compressed(
        out_path,
        X=out,
        t_grid=t_grid,
        start=start,
        end=end,
        shard=shard,
        shards=shards,
        failed=int(len(failed_members) > 0),
        failed_members=np.array(failed_members, dtype=int),
    )
    print(f"Saved shard {shard}/{shards} to {out_path}")

    if representative_diag is not None:
        diag_path = outdir / f"maneuver_diag_rank{shard:04d}.npz"
        payload = {
            "t_grid": t_grid,
            "representative_local_idx": int(representative_local_idx if representative_local_idx is not None else -1),
            "representative_global_idx": int(representative_global_idx if representative_global_idx is not None else -1),
            "axis": np.array(str(representative_diag["axis"])),
            "control_mode": np.array(str(representative_diag.get("control_mode", "open_loop"))),
            "half_turn_deg": float(representative_diag["half_turn_deg"]),
            "full_turn_deg": float(representative_diag["full_turn_deg"]),
            "rep_progress_axis_deg": np.asarray(representative_diag["progress_axis_deg"], dtype=float),
            "rep_eta1_cmd_rad": np.asarray(representative_diag["eta1_cmd_rad"], dtype=float),
            "rep_eta2_cmd_rad": np.asarray(representative_diag["eta2_cmd_rad"], dtype=float),
            "rep_eta_cmd_deg": np.asarray(representative_diag["eta_cmd_deg"], dtype=float),
            "rep_eta_target_deg": np.asarray(representative_diag["eta_target_deg"], dtype=float),
            "rep_eta_rate_deg_s": np.asarray(representative_diag["eta_rate_deg_s"], dtype=float),
            "rep_axis_rate_deg_s": np.asarray(representative_diag["axis_rate_deg_s"], dtype=float),
            "rep_torque_B_Nm": np.asarray(representative_diag["torque_B_Nm"], dtype=float),
            "rep_phase_index": np.asarray(representative_diag["phase_index"], dtype=int),
            "eta_max_deg": float(representative_diag["eta_max_deg"]),
            "eta_rate_max_deg_s": float(representative_diag["eta_rate_max_deg_s"]),
            "eta_accel_max_deg_s2": float(representative_diag["eta_accel_max_deg_s2"]),
            "point_kp_cmd_per_deg": float(representative_diag["point_kp_cmd_per_deg"]),
            "point_kd_cmd_per_deg_s": float(representative_diag["point_kd_cmd_per_deg_s"]),
            "rate_target_deg_s": float(representative_diag["rate_target_deg_s"]),
            "rate_kp_cmd_per_deg_s": float(representative_diag["rate_kp_cmd_per_deg_s"]),
        }
        if deterministic_diag is not None:
            payload.update(
                {
                    "det_progress_axis_deg": np.asarray(deterministic_diag["progress_axis_deg"], dtype=float),
                    "det_eta1_cmd_rad": np.asarray(deterministic_diag["eta1_cmd_rad"], dtype=float),
                    "det_eta2_cmd_rad": np.asarray(deterministic_diag["eta2_cmd_rad"], dtype=float),
                    "det_eta_cmd_deg": np.asarray(deterministic_diag["eta_cmd_deg"], dtype=float),
                    "det_eta_target_deg": np.asarray(deterministic_diag["eta_target_deg"], dtype=float),
                    "det_eta_rate_deg_s": np.asarray(deterministic_diag["eta_rate_deg_s"], dtype=float),
                    "det_axis_rate_deg_s": np.asarray(deterministic_diag["axis_rate_deg_s"], dtype=float),
                    "det_torque_B_Nm": np.asarray(deterministic_diag["torque_B_Nm"], dtype=float),
                    "det_phase_index": np.asarray(deterministic_diag["phase_index"], dtype=int),
                }
            )
        np.savez_compressed(diag_path, **payload)
        print(f"Saved maneuver diagnostics for rank {shard}/{shards} to {diag_path}")


if __name__ == "__main__":
    main()
