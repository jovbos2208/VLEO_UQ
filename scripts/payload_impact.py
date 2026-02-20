from __future__ import annotations

from typing import Any

import numpy as np


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)):
        return bool(int(value) != 0)
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    return bool(default)


def _quat_to_dcm_wxyz(q_wxyz: np.ndarray) -> np.ndarray:
    q = np.asarray(q_wxyz, dtype=float).reshape(4)
    n = float(np.linalg.norm(q))
    if not np.isfinite(n) or n <= 0.0:
        return np.eye(3)
    w, x, y, z = q / n
    return np.array(
        [
            [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
            [2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x)],
            [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), 1.0 - 2.0 * (x * x + y * y)],
        ],
        dtype=float,
    )


def _ground_intersect_spherical(r_sc_eci: np.ndarray, los_eci: np.ndarray, earth_radius_m: float) -> np.ndarray:
    r = np.asarray(r_sc_eci, dtype=float).reshape(3)
    l = np.asarray(los_eci, dtype=float).reshape(3)
    l_norm = float(np.linalg.norm(l))
    if not np.isfinite(l_norm) or l_norm <= 0.0:
        return np.full(3, np.nan, dtype=float)
    l = l / l_norm
    b = float(np.dot(r, l))
    c = float(np.dot(r, r) - earth_radius_m * earth_radius_m)
    disc = b * b - c
    if disc < 0.0:
        return np.full(3, np.nan, dtype=float)
    s = -b - float(np.sqrt(disc))
    if s <= 0.0:
        s = -b + float(np.sqrt(disc))
    if s <= 0.0:
        return np.full(3, np.nan, dtype=float)
    return r + s * l


def _great_circle_distance_m(p_a: np.ndarray, p_b: np.ndarray, earth_radius_m: float) -> float:
    a = np.asarray(p_a, dtype=float).reshape(3)
    b = np.asarray(p_b, dtype=float).reshape(3)
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if not np.isfinite(na) or not np.isfinite(nb) or na <= 0.0 or nb <= 0.0:
        return float("nan")
    ca = np.clip(float(np.dot(a, b) / (na * nb)), -1.0, 1.0)
    return float(earth_radius_m * np.arccos(ca))


def _sample_indices(n: int, max_samples: int) -> np.ndarray:
    n = int(n)
    max_samples = max(1, int(max_samples))
    if n <= max_samples:
        return np.arange(n, dtype=int)
    idx = np.linspace(0, n - 1, max_samples, dtype=int)
    return np.unique(idx)


def compute_payload_impact_metrics(
    *,
    det_states: np.ndarray,
    mc_states: np.ndarray,
    t_grid: np.ndarray,
    boresight_body: np.ndarray | None = None,
    earth_radius_m: float = 6378137.0,
    max_samples: int = 25,
) -> dict[str, Any]:
    det = np.asarray(det_states, dtype=float)
    mc = np.asarray(mc_states, dtype=float)
    t = np.asarray(t_grid, dtype=float).reshape(-1)
    if det.ndim != 2 or mc.ndim != 3:
        return {"available": False, "reason": "invalid_shapes"}
    if det.shape[0] != mc.shape[0] or det.shape[0] != t.size:
        return {"available": False, "reason": "time_shape_mismatch"}
    if det.shape[1] < 10 or mc.shape[2] < 10:
        return {"available": False, "reason": "state_size_lt_10"}
    if boresight_body is None:
        boresight_body = np.array([0.0, 0.0, 1.0], dtype=float)
    b_body = np.asarray(boresight_body, dtype=float).reshape(3)
    b_norm = float(np.linalg.norm(b_body))
    if not np.isfinite(b_norm) or b_norm <= 0.0:
        b_body = np.array([0.0, 0.0, 1.0], dtype=float)
    else:
        b_body = b_body / b_norm

    idxs = _sample_indices(t.size, max_samples)
    rows: list[dict[str, Any]] = []
    p95_vals: list[float] = []
    for i in idxs:
        R_det = _quat_to_dcm_wxyz(det[i, 6:10])
        los_det = R_det.T @ b_body
        gp_det = _ground_intersect_spherical(det[i, 0:3], los_det, earth_radius_m)

        errs = []
        for j in range(mc.shape[1]):
            R_mc = _quat_to_dcm_wxyz(mc[i, j, 6:10])
            los_mc = R_mc.T @ b_body
            gp_mc = _ground_intersect_spherical(mc[i, j, 0:3], los_mc, earth_radius_m)
            if np.any(~np.isfinite(gp_mc)) or np.any(~np.isfinite(gp_det)):
                continue
            d = _great_circle_distance_m(gp_mc, gp_det, earth_radius_m)
            if np.isfinite(d):
                errs.append(float(d))
        errs_arr = np.asarray(errs, dtype=float)
        if errs_arr.size == 0:
            rows.append({"time_s": float(t[i]), "valid_count": 0})
            continue
        p95 = float(np.percentile(errs_arr, 95.0))
        p95_vals.append(p95)
        rows.append(
            {
                "time_s": float(t[i]),
                "valid_count": int(errs_arr.size),
                "mean_m": float(np.mean(errs_arr)),
                "p50_m": float(np.percentile(errs_arr, 50.0)),
                "p95_m": p95,
                "max_m": float(np.max(errs_arr)),
            }
        )

    valid_rows = [r for r in rows if int(r.get("valid_count", 0)) > 0]
    out = {
        "available": bool(len(valid_rows) > 0),
        "mode": "ensemble_ground_intersection",
        "sample_count": int(len(rows)),
        "valid_sample_count": int(len(valid_rows)),
        "valid_sample_fraction": float(len(valid_rows) / max(1, len(rows))),
        "boresight_body": [float(x) for x in b_body],
        "samples": rows,
    }
    if valid_rows:
        out["final_p95_geolocation_error_m"] = float(valid_rows[-1]["p95_m"])
        out["mean_p95_geolocation_error_m"] = float(np.mean(p95_vals)) if p95_vals else None
        out["max_p95_geolocation_error_m"] = float(np.max(p95_vals)) if p95_vals else None
    return out


def compute_payload_impact_proxy_from_stats(
    *,
    mc_mean: np.ndarray,
    mc_std: np.ndarray,
    earth_radius_m: float = 6378137.0,
) -> dict[str, Any]:
    mean = np.asarray(mc_mean, dtype=float)
    std = np.asarray(mc_std, dtype=float)
    if mean.ndim != 2 or std.ndim != 2 or mean.shape != std.shape:
        return {"available": False, "reason": "invalid_shapes"}
    if mean.shape[1] < 10:
        return {"available": False, "reason": "state_size_lt_10"}
    pos_sigma = float(np.linalg.norm(std[-1, 0:3])) if std.shape[1] >= 3 else float("nan")
    qvec_sigma = float(np.linalg.norm(std[-1, 7:10])) if std.shape[1] >= 10 else float("nan")
    theta_sigma = 2.0 * qvec_sigma if np.isfinite(qvec_sigma) else float("nan")
    geo_sigma = float(np.hypot(pos_sigma, earth_radius_m * theta_sigma)) if np.isfinite(theta_sigma) else float("nan")
    return {
        "available": bool(np.isfinite(geo_sigma)),
        "mode": "small_angle_proxy_from_final_std",
        "final_position_sigma_norm_m": pos_sigma,
        "final_attitude_sigma_proxy_rad": theta_sigma,
        "final_geolocation_sigma_proxy_m": geo_sigma,
        "final_p95_geolocation_error_m": float(1.96 * geo_sigma) if np.isfinite(geo_sigma) else None,
    }


def payload_metrics_enabled(default: bool = True) -> bool:
    import os

    return _as_bool(os.environ.get("VLEO_PAYLOAD_METRICS_ON"), default)
