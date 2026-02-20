"""POD UQ harness: truth propagation, measurement simulation, batch OD, and metrics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence
import datetime as dt

import numpy as np

from ._vleo_uq import DeterministicPropagator, StmPropagator
from .env_sources import _gmst_rad


@dataclass
class PodUqResult:
    x0_truth: np.ndarray
    x0_est: np.ndarray
    P0_post: np.ndarray
    truth_states: np.ndarray
    est_states: np.ndarray
    measurements: np.ndarray
    meas_indices: np.ndarray
    rtn_sigma: np.ndarray
    rtn_error: np.ndarray
    radial_rms_m: float
    coverage_3sigma: dict
    coverage_sigma: dict
    t_grid: Optional[np.ndarray] = None


def normalize_quaternion_state(x: np.ndarray) -> None:
    q = x[6:10]
    norm = np.linalg.norm(q)
    if norm == 0.0:
        raise ValueError("quaternion must be non-zero")
    x[6:10] = q / norm


def build_measurement_indices(t_grid: np.ndarray, cadence_s: Optional[float]) -> np.ndarray:
    if cadence_s is None:
        return np.arange(len(t_grid), dtype=int)
    if cadence_s <= 0.0:
        raise ValueError("cadence_s must be positive")
    indices = [0]
    next_t = t_grid[0] + cadence_s
    for i in range(1, len(t_grid)):
        if t_grid[i] >= next_t - 1e-9:
            indices.append(i)
            next_t += cadence_s
    return np.array(indices, dtype=int)


def simulate_position_measurements(
    states: np.ndarray,
    sigma_m: Sequence[float] | float,
    rng: np.random.Generator,
    indices: Optional[Iterable[int]] = None,
) -> np.ndarray:
    if indices is None:
        indices = range(states.shape[0])
    indices = np.array(list(indices), dtype=int)
    base = states[indices, 0:3]
    if np.isscalar(sigma_m):
        noise = rng.normal(0.0, float(sigma_m), size=base.shape)
    else:
        sigma = np.array(sigma_m, dtype=float).reshape(1, 3)
        noise = rng.normal(0.0, 1.0, size=base.shape) * sigma
    return base + noise


def rtn_basis(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    r_norm = np.linalg.norm(r)
    if r_norm == 0.0:
        raise ValueError("position magnitude must be non-zero")
    r_hat = r / r_norm
    h = np.cross(r, v)
    h_norm = np.linalg.norm(h)
    if h_norm == 0.0:
        raise ValueError("orbit normal magnitude must be non-zero")
    n_hat = h / h_norm
    t_hat = np.cross(n_hat, r_hat)
    return np.column_stack((r_hat, t_hat, n_hat))


def compute_rtn_errors(truth_states: np.ndarray, est_states: np.ndarray) -> np.ndarray:
    errors = est_states[:, 0:3] - truth_states[:, 0:3]
    rtn_err = np.zeros_like(errors)
    for i in range(errors.shape[0]):
        C = rtn_basis(truth_states[i, 0:3], truth_states[i, 3:6])
        rtn_err[i] = C.T @ errors[i]
    return rtn_err


def compute_rtn_sigmas(P0_post: np.ndarray, stm: np.ndarray, truth_states: np.ndarray) -> np.ndarray:
    nt = truth_states.shape[0]
    sigmas = np.zeros((nt, 3))
    for i in range(nt):
        Phi = stm[i]
        P_t = Phi @ P0_post @ Phi.T
        P_pos = P_t[0:3, 0:3]
        C = rtn_basis(truth_states[i, 0:3], truth_states[i, 3:6])
        P_rtn = C.T @ P_pos @ C
        sigmas[i] = np.sqrt(np.diag(P_rtn))
    return sigmas


def compute_coverage_sigma_levels(
    rtn_err: np.ndarray,
    rtn_sigma: np.ndarray,
    levels: Sequence[float] = (1.0, 2.0, 3.0),
) -> dict:
    coverage = {}
    labels = ["R", "T", "N"]
    for level in levels:
        if level <= 0.0:
            raise ValueError("sigma coverage levels must be positive")
        key = f"{int(level) if float(level).is_integer() else level}sigma"
        level_cov = {}
        for idx, label in enumerate(labels):
            within = np.abs(rtn_err[:, idx]) <= float(level) * rtn_sigma[:, idx]
            level_cov[label] = float(np.mean(within))
        coverage[key] = level_cov
    return coverage


def compute_coverage_3sigma(rtn_err: np.ndarray, rtn_sigma: np.ndarray) -> dict:
    return compute_coverage_sigma_levels(rtn_err, rtn_sigma, levels=(3.0,))["3sigma"]


def _to_datetimes(t_grid: Sequence, t0: Optional[dt.datetime]) -> list[dt.datetime]:
    if isinstance(t_grid, np.ndarray) and np.issubdtype(t_grid.dtype, np.datetime64):
        return [dt.datetime.utcfromtimestamp(t.astype("datetime64[ns]").astype(int) / 1e9) for t in t_grid]
    if isinstance(t_grid, (list, tuple)) and t_grid and isinstance(t_grid[0], dt.datetime):
        return list(t_grid)
    if t0 is None:
        raise ValueError("t0 is required when t_grid is not datetime-like")
    return [t0 + dt.timedelta(seconds=float(sec)) for sec in t_grid]


def batch_od_position(
    prop_stm: StmPropagator,
    x0_guess: np.ndarray,
    P0_guess: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence,
    meas_indices: np.ndarray,
    measurements: np.ndarray,
    sigma_m: Sequence[float] | float,
    max_iter: int = 6,
    tol: float = 1e-6,
    damping: float = 1e-8,
) -> tuple[np.ndarray, np.ndarray]:
    if np.isscalar(sigma_m):
        sigma_vec = np.array([float(sigma_m)] * 3)
    else:
        sigma_vec = np.array(sigma_m, dtype=float)
    if np.any(sigma_vec <= 0.0):
        raise ValueError("measurement sigma must be positive")
    W = np.diag(1.0 / (sigma_vec * sigma_vec))

    x_est = x0_guess.copy()
    normalize_quaternion_state(x_est)

    for _ in range(max_iter):
        mean, _cov, stm = prop_stm.propagate(x_est, P0_guess, t_grid, env, return_stm=True)

        N = np.zeros((prop_stm.state_size, prop_stm.state_size))
        b = np.zeros(prop_stm.state_size)

        for k, idx in enumerate(meas_indices):
            y_pred = mean[idx, 0:3]
            r = measurements[k] - y_pred
            A = stm[idx][0:3, :]
            N += A.T @ W @ A
            b += A.T @ W @ r

        N += damping * np.eye(prop_stm.state_size)

        try:
            dx = np.linalg.solve(N, b)
        except np.linalg.LinAlgError:
            dx = np.linalg.lstsq(N, b, rcond=None)[0]

        x_est += dx
        normalize_quaternion_state(x_est)

        if np.linalg.norm(dx) < tol:
            break

    try:
        P0_post = np.linalg.inv(N)
    except np.linalg.LinAlgError:
        P0_post = np.linalg.pinv(N)

    return x_est, P0_post


def run_pod_uq(
    prop_det: DeterministicPropagator,
    prop_stm: StmPropagator,
    x0_truth: np.ndarray,
    x0_guess: np.ndarray,
    P0_guess: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence,
    meas_sigma_m: Sequence[float] | float,
    meas_cadence_s: Optional[float] = None,
    measurement_latency_s: float = 0.0,
    measurement_latency_jitter_s: float = 0.0,
    measurement_latency_seed: int = 0,
    rng: Optional[np.random.Generator] = None,
    max_iter: int = 6,
) -> PodUqResult:
    if rng is None:
        rng = np.random.default_rng(0)

    truth_states = prop_det.propagate(x0_truth, t_grid, env)
    meas_indices = build_measurement_indices(t_grid, meas_cadence_s)
    rng_lat = np.random.default_rng(int(measurement_latency_seed))
    meas_model_indices = apply_measurement_latency(
        t_grid,
        meas_indices,
        latency_s=float(measurement_latency_s),
        jitter_s=float(measurement_latency_jitter_s),
        rng=rng_lat,
    )
    measurements = simulate_position_measurements(truth_states, meas_sigma_m, rng, meas_model_indices)

    x0_est, P0_post = batch_od_position(
        prop_stm,
        x0_guess,
        P0_guess,
        t_grid,
        env,
        meas_model_indices,
        measurements,
        meas_sigma_m,
        max_iter=max_iter,
    )

    est_states = prop_det.propagate(x0_est, t_grid, env)
    _, _, stm = prop_stm.propagate(x0_est, P0_post, t_grid, env, return_stm=True)

    rtn_error = compute_rtn_errors(truth_states, est_states)
    rtn_sigma = compute_rtn_sigmas(P0_post, stm, truth_states)
    radial_rms_m = float(np.sqrt(np.mean(rtn_error[:, 0] ** 2)))
    coverage_3sigma = compute_coverage_3sigma(rtn_error, rtn_sigma)
    coverage_sigma = compute_coverage_sigma_levels(rtn_error, rtn_sigma)

    return PodUqResult(
        x0_truth=x0_truth,
        x0_est=x0_est,
        P0_post=P0_post,
        truth_states=truth_states,
        est_states=est_states,
        measurements=measurements,
        meas_indices=meas_indices,
        rtn_sigma=rtn_sigma,
        rtn_error=rtn_error,
        radial_rms_m=radial_rms_m,
        coverage_3sigma=coverage_3sigma,
        coverage_sigma=coverage_sigma,
        t_grid=np.array(t_grid, dtype=float),
    )


MEAS_CODE = 1
MEAS_CARRIER = 2
MEAS_SLR = 3


@dataclass
class GroundStation:
    name: str
    lat_rad: float
    lon_rad: float
    alt_m: float = 0.0


@dataclass
class GnssMeasurements:
    values: np.ndarray
    types: np.ndarray
    time_indices: np.ndarray
    sat_indices: np.ndarray
    freq_indices: np.ndarray
    sigmas: np.ndarray
    sat_positions_eci: np.ndarray
    wavelengths_m: np.ndarray
    frequencies_hz: np.ndarray
    clock_bias_truth: np.ndarray
    clock_drift_truth: np.ndarray
    tropo_truth: np.ndarray
    ambiguity_ids: np.ndarray
    ambiguities_cycles: np.ndarray
    clock_bias_sigma_rw_m: float
    clock_drift_sigma_rw_m_s: float
    tropo_sigma_rw_m: float
    sat_clock_bias_m: Optional[np.ndarray] = None
    sat_pco_eci_m: Optional[np.ndarray] = None
    is_outlier: Optional[np.ndarray] = None
    is_cycle_slip: Optional[np.ndarray] = None


@dataclass
class SlrMeasurements:
    values: np.ndarray
    time_indices: np.ndarray
    station_indices: np.ndarray
    sigmas: np.ndarray
    stations: Sequence[GroundStation]
    station_positions_eci: np.ndarray
    range_bias_m: np.ndarray


@dataclass
class PodUqMeasurementResult:
    x0_truth: np.ndarray
    x0_est: np.ndarray
    P0_post: np.ndarray
    truth_states: np.ndarray
    est_states: np.ndarray
    measurements: dict
    rtn_sigma: np.ndarray
    rtn_error: np.ndarray
    radial_rms_m: float
    coverage_3sigma: dict
    coverage_sigma: dict
    nuisance: dict
    t_grid: Optional[np.ndarray] = None


def geodetic_to_ecef(lat_rad: float, lon_rad: float, alt_m: float) -> np.ndarray:
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)
    sin_lon = np.sin(lon_rad)
    cos_lon = np.cos(lon_rad)
    N = a / np.sqrt(1.0 - e2 * sin_lat * sin_lat)
    x = (N + alt_m) * cos_lat * cos_lon
    y = (N + alt_m) * cos_lat * sin_lon
    z = (N * (1.0 - e2) + alt_m) * sin_lat
    return np.array([x, y, z])


def ecef_to_geodetic(x_m: float, y_m: float, z_m: float) -> tuple[float, float, float]:
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    lon = np.arctan2(y_m, x_m)
    p = np.sqrt(x_m * x_m + y_m * y_m)
    if p == 0.0:
        lat = np.sign(z_m) * (np.pi / 2.0)
        alt = np.abs(z_m) - a * (1.0 - f)
        return lat, lon, alt
    lat = np.arctan2(z_m, p * (1.0 - e2))
    for _ in range(6):
        sin_lat = np.sin(lat)
        N = a / np.sqrt(1.0 - e2 * sin_lat * sin_lat)
        alt = p / np.cos(lat) - N
        lat = np.arctan2(z_m, p * (1.0 - e2 * (N / (N + alt))))
    sin_lat = np.sin(lat)
    N = a / np.sqrt(1.0 - e2 * sin_lat * sin_lat)
    alt = p / np.cos(lat) - N
    return lat, lon, alt


def rotation_z(theta_rad: float) -> np.ndarray:
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def ecef_to_eci(vec_ecef: np.ndarray, theta_rad: float) -> np.ndarray:
    return rotation_z(theta_rad) @ vec_ecef


def eci_to_ecef(vec_eci: np.ndarray, theta_rad: float) -> np.ndarray:
    return rotation_z(theta_rad).T @ vec_eci


def ecef_to_enu(vec_ecef: np.ndarray, lat_rad: float, lon_rad: float) -> np.ndarray:
    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)
    sin_lon = np.sin(lon_rad)
    cos_lon = np.cos(lon_rad)
    rot = np.array(
        [
            [-sin_lon, cos_lon, 0.0],
            [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
            [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
        ]
    )
    return rot @ vec_ecef


def build_time_mask(t_grid: np.ndarray, gaps: Optional[Sequence[tuple[float, float]]]) -> np.ndarray:
    mask = np.ones(len(t_grid), dtype=bool)
    if gaps:
        for start, end in gaps:
            mask &= ~((t_grid >= start) & (t_grid <= end))
    return mask


def apply_measurement_latency(
    t_grid: np.ndarray,
    indices: np.ndarray,
    latency_s: float = 0.0,
    jitter_s: float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    idx = np.array(indices, dtype=int).reshape(-1)
    if idx.size == 0:
        return idx
    if latency_s <= 0.0 and jitter_s <= 0.0:
        return idx
    if rng is None:
        rng = np.random.default_rng(0)
    t = np.asarray(t_grid, dtype=float).reshape(-1)
    if t.size == 0:
        return np.zeros_like(idx)
    jitter = np.zeros(idx.size, dtype=float)
    if jitter_s > 0.0:
        jitter = rng.normal(0.0, float(jitter_s), size=idx.size)
    delayed_t = t[idx] - float(latency_s) - jitter
    out = np.searchsorted(t, delayed_t, side="right") - 1
    out = np.clip(out, 0, t.size - 1)
    return out.astype(int)


def sample_operational_outage_mask(
    t_grid: np.ndarray,
    *,
    enabled: bool = False,
    start_rate_per_hour: float = 0.0,
    mean_duration_s: float = 0.0,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    t = np.asarray(t_grid, dtype=float).reshape(-1)
    n = t.size
    mask = np.ones(n, dtype=bool)
    if (not enabled) or n == 0 or start_rate_per_hour <= 0.0 or mean_duration_s <= 0.0:
        return mask
    if rng is None:
        rng = np.random.default_rng(0)

    blocked_until = -np.inf
    for i in range(n):
        ti = float(t[i])
        if ti < blocked_until:
            mask[i] = False
            continue
        if i == 0:
            dt = max(1.0, float(np.median(np.diff(t))) if n > 1 else 1.0)
        else:
            dt = max(0.0, float(t[i] - t[i - 1]))
        p_start = 1.0 - float(np.exp(-start_rate_per_hour * dt / 3600.0))
        if rng.random() < p_start:
            dur = float(rng.exponential(mean_duration_s))
            blocked_until = ti + max(1.0, dur)
            mask[i] = False
    return mask


def line_of_sight(r_sc: np.ndarray, r_tx: np.ndarray, r_earth_m: float = 6378137.0) -> bool:
    d = r_tx - r_sc
    denom = np.dot(d, d)
    if denom == 0.0:
        return False
    t = -np.dot(r_sc, d) / denom
    if t <= 0.0:
        return np.dot(r_sc, r_sc) > r_earth_m * r_earth_m
    if t >= 1.0:
        return np.dot(r_tx, r_tx) > r_earth_m * r_earth_m
    closest = r_sc + t * d
    return np.dot(closest, closest) > r_earth_m * r_earth_m


def simulate_gnss_measurements(
    truth_states: np.ndarray,
    t_grid: np.ndarray,
    sat_positions_eci: np.ndarray,
    code_sigma_m: Sequence[float] | float = 0.5,
    carrier_sigma_m: Sequence[float] | float = 0.003,
    frequencies_hz: Optional[Sequence[float]] = None,
    wavelengths_m: Optional[Sequence[float]] = None,
    clock_bias_sigma_rw_m: float = 0.5,
    clock_drift_sigma_rw_m_s: float = 1e-4,
    tropo_sigma_rw_m: float = 0.0,
    ambiguity_cycles: Optional[np.ndarray] = None,
    rng: Optional[np.random.Generator] = None,
    cadence_s: Optional[float] = 1.0,
    availability_mask: Optional[np.ndarray] = None,
    gaps: Optional[Sequence[tuple[float, float]]] = None,
    require_los: bool = True,
    earth_radius_m: float = 6378137.0,
    earth_margin_m: float = 50e3,
    fov_half_angle_deg: Optional[float] = 80.0,
    boresight_mode: str = "zenith",
    dropout_prob: float = 0.01,
    cycle_slip_gap_s: Optional[float] = 30.0,
    cycle_slip_prob_per_min: float = 1e-4,
    code_outlier_prob: float = 0.0,
    code_outlier_sigma_scale: float = 25.0,
    carrier_outlier_prob: float = 0.0,
    carrier_outlier_sigma_scale: float = 25.0,
    include_code: bool = True,
    include_carrier: bool = True,
    sat_clock_bias_m: Optional[np.ndarray] = None,
    sat_pco_eci_m: Optional[np.ndarray] = None,
    ops_outage_on: bool = False,
    ops_outage_rate_per_hour: float = 0.0,
    ops_outage_mean_duration_s: float = 0.0,
) -> GnssMeasurements:
    if rng is None:
        rng = np.random.default_rng(0)
    if sat_positions_eci.shape[0] != truth_states.shape[0]:
        raise ValueError("sat_positions_eci must align with t_grid length")
    n_sats = sat_positions_eci.shape[1]
    if n_sats == 0:
        raise ValueError("sat_positions_eci must include at least one satellite")

    c_m_s = 299792458.0
    if wavelengths_m is None:
        if frequencies_hz is None:
            frequencies_hz = (1575.42e6, 1227.60e6)
        wavelengths = np.array([c_m_s / f for f in frequencies_hz], dtype=float)
    else:
        wavelengths = np.array(wavelengths_m, dtype=float)
        if frequencies_hz is None:
            frequencies_hz = tuple(c_m_s / wavelengths)
    frequencies = np.array(frequencies_hz, dtype=float)
    if wavelengths.shape[0] != frequencies.shape[0]:
        raise ValueError("frequencies_hz and wavelengths_m must align")
    n_freq = wavelengths.shape[0]

    def expand_sigma(val: Sequence[float] | float, name: str) -> np.ndarray:
        if np.isscalar(val):
            return np.full(n_freq, float(val))
        arr = np.array(val, dtype=float)
        if arr.shape[0] != n_freq:
            raise ValueError(f"{name} must match number of frequencies")
        return arr

    code_sigmas = expand_sigma(code_sigma_m, "code_sigma_m")
    carrier_sigmas = expand_sigma(carrier_sigma_m, "carrier_sigma_m")

    if ambiguity_cycles is not None:
        ambiguity_cycles = np.array(ambiguity_cycles, dtype=int)
        if ambiguity_cycles.shape not in {(n_sats,), (n_sats, n_freq)}:
            raise ValueError("ambiguity_cycles must match (n_sats,) or (n_sats, n_freq)")

    clock_bias = np.zeros(len(t_grid))
    clock_drift = np.zeros(len(t_grid))
    tropo = np.zeros(len(t_grid))
    for i in range(1, len(t_grid)):
        dt = t_grid[i] - t_grid[i - 1]
        if dt <= 0.0:
            raise ValueError("t_grid must be strictly increasing")
        clock_drift[i] = clock_drift[i - 1] + clock_drift_sigma_rw_m_s * np.sqrt(dt) * rng.normal()
        clock_bias[i] = (
            clock_bias[i - 1]
            + clock_drift[i - 1] * dt
            + clock_bias_sigma_rw_m * np.sqrt(dt) * rng.normal()
        )
        tropo[i] = tropo[i - 1] + tropo_sigma_rw_m * np.sqrt(dt) * rng.normal()

    time_mask = build_time_mask(t_grid, gaps)
    time_indices = build_measurement_indices(t_grid, cadence_s)
    time_indices = time_indices[time_mask[time_indices]]
    ops_mask = sample_operational_outage_mask(
        t_grid,
        enabled=bool(ops_outage_on),
        start_rate_per_hour=float(max(0.0, ops_outage_rate_per_hour)),
        mean_duration_s=float(max(0.0, ops_outage_mean_duration_s)),
        rng=rng,
    )
    time_indices = time_indices[ops_mask[time_indices]]

    values = []
    types = []
    time_out = []
    sat_out = []
    freq_out = []
    sigmas = []
    amb_ids_out = []
    outlier_flags = []
    cycle_slip_flags = []

    if availability_mask is not None:
        availability_mask = np.array(availability_mask, dtype=bool)
        if availability_mask.shape not in {(len(t_grid),), (len(t_grid), n_sats)}:
            raise ValueError("availability_mask must be (nt,) or (nt, n_sats)")

    if sat_clock_bias_m is not None:
        sat_clock_bias_m = np.array(sat_clock_bias_m, dtype=float)
        if sat_clock_bias_m.shape != (len(t_grid), n_sats):
            raise ValueError("sat_clock_bias_m must be (nt, n_sats)")
    if sat_pco_eci_m is not None:
        sat_pco_eci_m = np.array(sat_pco_eci_m, dtype=float)
        if sat_pco_eci_m.shape != (len(t_grid), n_sats, 3):
            raise ValueError("sat_pco_eci_m must be (nt, n_sats, 3)")

    if fov_half_angle_deg is not None:
        fov_half_angle_rad = np.deg2rad(fov_half_angle_deg)
    else:
        fov_half_angle_rad = None

    def boresight_vector(r_sc: np.ndarray) -> Optional[np.ndarray]:
        if fov_half_angle_rad is None:
            return None
        if boresight_mode == "zenith":
            norm = np.linalg.norm(r_sc)
            if norm == 0.0:
                return None
            return r_sc / norm
        if boresight_mode == "none":
            return None
        raise ValueError("boresight_mode must be 'zenith' or 'none'")

    amb_ids = np.full((n_sats, n_freq), -1, dtype=int)
    amb_cycles = []
    last_obs_t = np.full((n_sats, n_freq), np.nan)

    def new_ambiguity(sat_idx: int, freq_idx: int) -> None:
        amb_ids[sat_idx, freq_idx] = len(amb_cycles)
        amb_cycles.append(rng.integers(10_000, 200_000))

    for s in range(n_sats):
        for f in range(n_freq):
            if ambiguity_cycles is not None:
                if ambiguity_cycles.ndim == 1:
                    value = ambiguity_cycles[s]
                else:
                    value = ambiguity_cycles[s, f]
                amb_ids[s, f] = len(amb_cycles)
                amb_cycles.append(int(value))
            else:
                new_ambiguity(s, f)

    for ti in time_indices:
        r_sc = truth_states[ti, 0:3]
        boresight = boresight_vector(r_sc)
        for s in range(n_sats):
            if availability_mask is not None:
                if availability_mask.ndim == 1 and not availability_mask[ti]:
                    continue
                if availability_mask.ndim == 2 and not availability_mask[ti, s]:
                    continue
            r_sat = sat_positions_eci[ti, s]
            if sat_pco_eci_m is not None:
                r_sat = r_sat + sat_pco_eci_m[ti, s]
            if require_los and not line_of_sight(r_sc, r_sat, earth_radius_m + earth_margin_m):
                continue
            if dropout_prob > 0.0 and rng.random() < dropout_prob:
                continue
            rho = np.linalg.norm(r_sc - r_sat)
            bias = clock_bias[ti] + tropo[ti]
            if sat_clock_bias_m is not None:
                bias += sat_clock_bias_m[ti, s]
            los = r_sat - r_sc
            los_norm = np.linalg.norm(los)
            if los_norm == 0.0:
                continue
            los_hat = los / los_norm

            if boresight is not None and fov_half_angle_rad is not None:
                angle = np.arccos(np.clip(np.dot(los_hat, boresight), -1.0, 1.0))
                if angle > fov_half_angle_rad:
                    continue

            for f in range(n_freq):
                if cycle_slip_gap_s is not None:
                    last_t = last_obs_t[s, f]
                    if np.isfinite(last_t) and (t_grid[ti] - last_t) > cycle_slip_gap_s:
                        new_ambiguity(s, f)
                        slipped = True
                    else:
                        slipped = False
                else:
                    slipped = False
                if cycle_slip_prob_per_min > 0.0:
                    slip_prob = cycle_slip_prob_per_min * (
                        0.0
                        if np.isnan(last_obs_t[s, f])
                        else (t_grid[ti] - last_obs_t[s, f]) / 60.0
                    )
                    slip_prob = float(np.clip(slip_prob, 0.0, 1.0))
                    if slip_prob > 0.0 and rng.random() < slip_prob:
                        new_ambiguity(s, f)
                        slipped = True

                if include_code:
                    code_noise = code_sigmas[f] * rng.normal()
                    code_is_outlier = False
                    if code_outlier_prob > 0.0 and rng.random() < code_outlier_prob:
                        code_noise += code_sigmas[f] * code_outlier_sigma_scale * rng.normal()
                        code_is_outlier = True
                    values.append(rho + bias + code_noise)
                    types.append(MEAS_CODE)
                    time_out.append(ti)
                    sat_out.append(s)
                    freq_out.append(f)
                    sigmas.append(code_sigmas[f])
                    amb_ids_out.append(-1)
                    outlier_flags.append(code_is_outlier)
                    cycle_slip_flags.append(False)

                if include_carrier:
                    amb_idx = amb_ids[s, f]
                    amb_cycles_val = amb_cycles[amb_idx]
                    carrier_noise = carrier_sigmas[f] * rng.normal()
                    carrier_is_outlier = False
                    if carrier_outlier_prob > 0.0 and rng.random() < carrier_outlier_prob:
                        carrier_noise += carrier_sigmas[f] * carrier_outlier_sigma_scale * rng.normal()
                        carrier_is_outlier = True
                    values.append(
                        rho + bias + amb_cycles_val * wavelengths[f] + carrier_noise
                    )
                    types.append(MEAS_CARRIER)
                    time_out.append(ti)
                    sat_out.append(s)
                    freq_out.append(f)
                    sigmas.append(carrier_sigmas[f])
                    amb_ids_out.append(amb_idx)
                    outlier_flags.append(carrier_is_outlier)
                    cycle_slip_flags.append(slipped)

                last_obs_t[s, f] = t_grid[ti]

    return GnssMeasurements(
        values=np.array(values, dtype=float),
        types=np.array(types, dtype=int),
        time_indices=np.array(time_out, dtype=int),
        sat_indices=np.array(sat_out, dtype=int),
        freq_indices=np.array(freq_out, dtype=int),
        sigmas=np.array(sigmas, dtype=float),
        sat_positions_eci=sat_positions_eci,
        wavelengths_m=wavelengths,
        frequencies_hz=frequencies,
        clock_bias_truth=clock_bias,
        clock_drift_truth=clock_drift,
        tropo_truth=tropo,
        ambiguity_ids=np.array(amb_ids_out, dtype=int),
        ambiguities_cycles=np.array(amb_cycles, dtype=int),
        clock_bias_sigma_rw_m=clock_bias_sigma_rw_m,
        clock_drift_sigma_rw_m_s=clock_drift_sigma_rw_m_s,
        tropo_sigma_rw_m=tropo_sigma_rw_m,
        sat_clock_bias_m=sat_clock_bias_m,
        sat_pco_eci_m=sat_pco_eci_m,
        is_outlier=np.array(outlier_flags, dtype=bool),
        is_cycle_slip=np.array(cycle_slip_flags, dtype=bool),
    )


def build_station_positions_eci(
    stations: Sequence[GroundStation],
    t_grid: np.ndarray,
    theta0_rad: float = 0.0,
    omega_earth_rad_s: float = 7.2921150e-5,
) -> np.ndarray:
    ecef = np.array([geodetic_to_ecef(s.lat_rad, s.lon_rad, s.alt_m) for s in stations])
    out = np.zeros((len(t_grid), len(stations), 3))
    t0 = t_grid[0]
    for i, t in enumerate(t_grid):
        theta = theta0_rad + omega_earth_rad_s * (t - t0)
        rot = rotation_z(theta)
        out[i] = (rot @ ecef.T).T
    return out


def elevation_mask(
    truth_states: np.ndarray,
    t_grid: np.ndarray,
    stations: Sequence[GroundStation],
    min_elevation_deg: float = 10.0,
    theta0_rad: float = 0.0,
    omega_earth_rad_s: float = 7.2921150e-5,
) -> np.ndarray:
    mask = np.zeros((len(t_grid), len(stations)), dtype=bool)
    t0 = t_grid[0]
    for i, t in enumerate(t_grid):
        theta = theta0_rad + omega_earth_rad_s * (t - t0)
        for j, st in enumerate(stations):
            r_ecef = geodetic_to_ecef(st.lat_rad, st.lon_rad, st.alt_m)
            r_sc_ecef = eci_to_ecef(truth_states[i, 0:3], theta)
            rho = r_sc_ecef - r_ecef
            enu = ecef_to_enu(rho, st.lat_rad, st.lon_rad)
            elev = np.arcsin(enu[2] / np.linalg.norm(enu))
            mask[i, j] = elev >= np.deg2rad(min_elevation_deg)
    return mask


def sun_position_eci(ts: dt.datetime) -> np.ndarray:
    jd = (ts - dt.datetime(2000, 1, 1, 12)).total_seconds() / 86400.0 + 2451545.0
    T = (jd - 2451545.0) / 36525.0
    L = np.deg2rad((280.460 + 36000.770 * T) % 360.0)
    g = np.deg2rad((357.528 + 35999.050 * T) % 360.0)
    lam = L + np.deg2rad(1.915) * np.sin(g) + np.deg2rad(0.020) * np.sin(2.0 * g)
    eps = np.deg2rad(23.4393 - 0.0130 * T)
    r_au = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2.0 * g)
    r_m = r_au * 149597870700.0
    x = r_m * np.cos(lam)
    y = r_m * np.cos(eps) * np.sin(lam)
    z = r_m * np.sin(eps) * np.sin(lam)
    return np.array([x, y, z])


def station_night_mask(
    t_grid: Sequence,
    stations: Sequence[GroundStation],
    t0: Optional[dt.datetime],
    twilight_deg: float = -6.0,
) -> np.ndarray:
    times = _to_datetimes(t_grid, t0)
    mask = np.zeros((len(times), len(stations)), dtype=bool)
    for i, ts in enumerate(times):
        gmst = _gmst_rad(ts)
        rot = rotation_z(gmst)
        sun_ecef = rot.T @ sun_position_eci(ts)
        for j, st in enumerate(stations):
            r_ecef = geodetic_to_ecef(st.lat_rad, st.lon_rad, st.alt_m)
            rho = sun_ecef - r_ecef
            enu = ecef_to_enu(rho, st.lat_rad, st.lon_rad)
            elev = np.arcsin(enu[2] / np.linalg.norm(enu))
            mask[i, j] = elev <= np.deg2rad(twilight_deg)
    return mask


def station_weather_mask(
    nt: int,
    n_stations: int,
    weather_clear_prob: float = 0.8,
    weather_p_stay_clear: float = 0.985,
    weather_p_stay_blocked: float = 0.93,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    if nt <= 0 or n_stations <= 0:
        raise ValueError("nt and n_stations must be positive")
    if rng is None:
        rng = np.random.default_rng(0)
    clear_prob = float(np.clip(weather_clear_prob, 0.0, 1.0))
    p_stay_clear = float(np.clip(weather_p_stay_clear, 0.0, 1.0))
    p_stay_blocked = float(np.clip(weather_p_stay_blocked, 0.0, 1.0))

    mask = np.zeros((nt, n_stations), dtype=bool)
    states = rng.random(n_stations) < clear_prob
    mask[0, :] = states
    for i in range(1, nt):
        r = rng.random(n_stations)
        stay_clear = states & (r < p_stay_clear)
        recover = (~states) & (r >= p_stay_blocked)
        states = stay_clear | recover
        mask[i, :] = states
    return mask


def slr_availability_mask(
    truth_states: np.ndarray,
    t_grid: np.ndarray,
    stations: Sequence[GroundStation],
    min_elevation_deg: float = 20.0,
    require_night: bool = True,
    twilight_deg: float = -6.0,
    t0: Optional[dt.datetime] = None,
    theta0_rad: float = 0.0,
    omega_earth_rad_s: float = 7.2921150e-5,
    weather_enabled: bool = False,
    weather_clear_prob: float = 0.8,
    weather_p_stay_clear: float = 0.985,
    weather_p_stay_blocked: float = 0.93,
    weather_seed: Optional[int] = None,
    weather_rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    mask = elevation_mask(
        truth_states,
        t_grid,
        stations,
        min_elevation_deg=min_elevation_deg,
        theta0_rad=theta0_rad,
        omega_earth_rad_s=omega_earth_rad_s,
    )
    if require_night:
        night = station_night_mask(t_grid, stations, t0, twilight_deg=twilight_deg)
        mask &= night
    if weather_enabled:
        if weather_rng is None:
            weather_rng = np.random.default_rng(0 if weather_seed is None else int(weather_seed))
        weather = station_weather_mask(
            nt=len(t_grid),
            n_stations=len(stations),
            weather_clear_prob=weather_clear_prob,
            weather_p_stay_clear=weather_p_stay_clear,
            weather_p_stay_blocked=weather_p_stay_blocked,
            rng=weather_rng,
        )
        mask &= weather
    return mask


def simulate_slr_measurements(
    truth_states: np.ndarray,
    t_grid: np.ndarray,
    stations: Sequence[GroundStation],
    sigma_m: float = 0.01,
    rng: Optional[np.random.Generator] = None,
    cadence_s: Optional[float] = 20.0,
    availability_mask: Optional[np.ndarray] = None,
    gaps: Optional[Sequence[tuple[float, float]]] = None,
    min_elevation_deg: float = 20.0,
    theta0_rad: float = 0.0,
    omega_earth_rad_s: float = 7.2921150e-5,
    range_bias_m: Optional[Sequence[float]] = None,
    range_bias_sigma_m: float = 0.005,
    ops_outage_on: bool = False,
    ops_outage_rate_per_hour: float = 0.0,
    ops_outage_mean_duration_s: float = 0.0,
) -> SlrMeasurements:
    if rng is None:
        rng = np.random.default_rng(0)
    if range_bias_m is None:
        if range_bias_sigma_m > 0.0:
            range_bias_m = rng.normal(0.0, range_bias_sigma_m, size=len(stations))
        else:
            range_bias_m = np.zeros(len(stations))
    range_bias_m = np.array(range_bias_m, dtype=float)
    if range_bias_m.shape[0] != len(stations):
        raise ValueError("range_bias_m must match station count")

    time_mask = build_time_mask(t_grid, gaps)
    time_indices = build_measurement_indices(t_grid, cadence_s)
    time_indices = time_indices[time_mask[time_indices]]
    ops_mask = sample_operational_outage_mask(
        t_grid,
        enabled=bool(ops_outage_on),
        start_rate_per_hour=float(max(0.0, ops_outage_rate_per_hour)),
        mean_duration_s=float(max(0.0, ops_outage_mean_duration_s)),
        rng=rng,
    )
    time_indices = time_indices[ops_mask[time_indices]]

    if availability_mask is None:
        availability_mask = elevation_mask(
            truth_states,
            t_grid,
            stations,
            min_elevation_deg=min_elevation_deg,
            theta0_rad=theta0_rad,
            omega_earth_rad_s=omega_earth_rad_s,
        )
    else:
        availability_mask = np.array(availability_mask, dtype=bool)
        if availability_mask.shape != (len(t_grid), len(stations)):
            raise ValueError("availability_mask must be (nt, n_stations)")

    station_positions_eci = build_station_positions_eci(
        stations, t_grid, theta0_rad=theta0_rad, omega_earth_rad_s=omega_earth_rad_s
    )

    values = []
    time_out = []
    station_out = []
    sigmas = []

    for ti in time_indices:
        r_sc = truth_states[ti, 0:3]
        for s in range(len(stations)):
            if not availability_mask[ti, s]:
                continue
            r_st = station_positions_eci[ti, s]
            rho = np.linalg.norm(r_sc - r_st)
            values.append(rho + range_bias_m[s] + sigma_m * rng.normal())
            time_out.append(ti)
            station_out.append(s)
            sigmas.append(sigma_m)

    return SlrMeasurements(
        values=np.array(values, dtype=float),
        time_indices=np.array(time_out, dtype=int),
        station_indices=np.array(station_out, dtype=int),
        sigmas=np.array(sigmas, dtype=float),
        stations=stations,
        station_positions_eci=station_positions_eci,
        range_bias_m=range_bias_m,
    )


def batch_od_measurements(
    prop_stm: StmPropagator,
    x0_guess: np.ndarray,
    P0_guess: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence,
    gnss: Optional[GnssMeasurements] = None,
    slr: Optional[SlrMeasurements] = None,
    estimate_clock_bias: bool = True,
    estimate_clock_drift: bool = False,
    estimate_tropo: bool = False,
    estimate_ambiguity: bool = True,
    estimate_slr_bias: bool = False,
    use_prior: bool = True,
    measurement_latency_s: float = 0.0,
    measurement_latency_jitter_s: float = 0.0,
    measurement_latency_seed: int = 0,
    max_iter: int = 6,
    tol: float = 1e-6,
    damping: float = 1e-8,
) -> tuple[np.ndarray, np.ndarray, dict]:
    n_state = prop_stm.state_size
    rng_latency = np.random.default_rng(int(measurement_latency_seed))
    gnss_time_effective = None
    slr_time_effective = None
    if gnss is not None:
        gnss_time_effective = apply_measurement_latency(
            t_grid,
            gnss.time_indices,
            latency_s=float(measurement_latency_s),
            jitter_s=float(measurement_latency_jitter_s),
            rng=rng_latency,
        )
    if slr is not None:
        slr_time_effective = apply_measurement_latency(
            t_grid,
            slr.time_indices,
            latency_s=float(measurement_latency_s),
            jitter_s=float(measurement_latency_jitter_s),
            rng=rng_latency,
        )

    offsets = {"state": 0}
    offset = n_state

    if gnss is not None and estimate_clock_bias:
        gnss_times = np.unique(gnss_time_effective)
        clock_bias_map = {t: i for i, t in enumerate(gnss_times)}
        offsets["clock_bias"] = offset
        offset += len(gnss_times)
    else:
        gnss_times = np.array([], dtype=int)
        clock_bias_map = {}

    if gnss is not None and estimate_clock_drift:
        offsets["clock_drift"] = offset
        offset += len(gnss_times)
    else:
        offsets["clock_drift"] = None

    if gnss is not None and estimate_tropo:
        offsets["tropo"] = offset
        offset += len(gnss_times)
    else:
        offsets["tropo"] = None

    has_carrier = gnss is not None and np.any(gnss.types == MEAS_CARRIER)
    if gnss is not None and estimate_ambiguity and has_carrier:
        offsets["ambiguity"] = offset
        offset += gnss.ambiguities_cycles.shape[0]
    else:
        offsets["ambiguity"] = None

    if slr is not None and estimate_slr_bias:
        offsets["slr_bias"] = offset
        offset += len(slr.stations)
    else:
        offsets["slr_bias"] = None

    n_params = offset

    p = np.zeros(n_params)
    p[:n_state] = x0_guess
    normalize_quaternion_state(p[:n_state])

    for _ in range(max_iter):
        mean, _cov, stm = prop_stm.propagate(p[:n_state], P0_guess, t_grid, env, return_stm=True)

        rows = []
        resids = []
        sigmas = []

        if gnss is not None:
            clock_bias_idx = None
            if offsets.get("clock_bias") is not None:
                clock_bias_idx = np.array(
                    [offsets["clock_bias"] + clock_bias_map[t] for t in gnss_time_effective],
                    dtype=int,
                )
            tropo_idx = None
            if offsets.get("tropo") is not None:
                tropo_idx = np.array(
                    [offsets["tropo"] + clock_bias_map[t] for t in gnss_time_effective],
                    dtype=int,
                )
            amb_idx = None
            if offsets.get("ambiguity") is not None:
                amb_idx = offsets["ambiguity"] + gnss.ambiguity_ids

            for k in range(len(gnss.values)):
                ti = int(gnss_time_effective[k])
                si = gnss.sat_indices[k]
                fi = gnss.freq_indices[k]
                r_sc = mean[ti, 0:3]
                r_sat = gnss.sat_positions_eci[ti, si]
                rho_vec = r_sc - r_sat
                rho = np.linalg.norm(rho_vec)
                rho_hat = rho_vec / rho

                pred = rho
                if clock_bias_idx is not None:
                    pred += p[clock_bias_idx[k]]
                if tropo_idx is not None:
                    pred += p[tropo_idx[k]]
                if gnss.types[k] == MEAS_CARRIER and amb_idx is not None:
                    pred += p[amb_idx[k]] * gnss.wavelengths_m[fi]

                r = gnss.values[k] - pred
                row = np.zeros(n_params)
                row[:n_state] = rho_hat @ stm[ti][0:3, :]
                if clock_bias_idx is not None:
                    row[clock_bias_idx[k]] = 1.0
                if tropo_idx is not None:
                    row[tropo_idx[k]] = 1.0
                if gnss.types[k] == MEAS_CARRIER and amb_idx is not None:
                    row[amb_idx[k]] = gnss.wavelengths_m[fi]

                rows.append(row)
                resids.append(r)
                sigmas.append(gnss.sigmas[k])

            if offsets.get("clock_bias") is not None and gnss.clock_bias_sigma_rw_m > 0.0:
                for i in range(len(gnss_times) - 1):
                    dt = t_grid[gnss_times[i + 1]] - t_grid[gnss_times[i]]
                    if dt <= 0.0:
                        continue
                    idx_prev = offsets["clock_bias"] + i
                    idx_next = offsets["clock_bias"] + i + 1
                    row = np.zeros(n_params)
                    row[idx_next] = 1.0
                    row[idx_prev] = -1.0
                    r = -(row @ p)
                    rows.append(row)
                    resids.append(r)
                    sigmas.append(gnss.clock_bias_sigma_rw_m * np.sqrt(dt))

            if offsets.get("clock_drift") is not None and gnss.clock_drift_sigma_rw_m_s > 0.0:
                for i in range(len(gnss_times) - 1):
                    dt = t_grid[gnss_times[i + 1]] - t_grid[gnss_times[i]]
                    if dt <= 0.0:
                        continue
                    idx_prev = offsets["clock_drift"] + i
                    idx_next = offsets["clock_drift"] + i + 1
                    row = np.zeros(n_params)
                    row[idx_next] = 1.0
                    row[idx_prev] = -1.0
                    r = -(row @ p)
                    rows.append(row)
                    resids.append(r)
                    sigmas.append(gnss.clock_drift_sigma_rw_m_s * np.sqrt(dt))

                if offsets.get("clock_bias") is not None:
                    for i in range(len(gnss_times) - 1):
                        dt = t_grid[gnss_times[i + 1]] - t_grid[gnss_times[i]]
                        if dt <= 0.0:
                            continue
                        idx_bias_next = offsets["clock_bias"] + i + 1
                        idx_bias_prev = offsets["clock_bias"] + i
                        idx_drift = offsets["clock_drift"] + i
                        row = np.zeros(n_params)
                        row[idx_bias_next] = 1.0
                        row[idx_bias_prev] = -1.0
                        row[idx_drift] = -dt
                        r = -(row @ p)
                        rows.append(row)
                        resids.append(r)
                        sigmas.append(gnss.clock_bias_sigma_rw_m * np.sqrt(dt))

            if offsets.get("tropo") is not None and gnss.tropo_sigma_rw_m > 0.0:
                for i in range(len(gnss_times) - 1):
                    dt = t_grid[gnss_times[i + 1]] - t_grid[gnss_times[i]]
                    if dt <= 0.0:
                        continue
                    idx_prev = offsets["tropo"] + i
                    idx_next = offsets["tropo"] + i + 1
                    row = np.zeros(n_params)
                    row[idx_next] = 1.0
                    row[idx_prev] = -1.0
                    r = -(row @ p)
                    rows.append(row)
                    resids.append(r)
                    sigmas.append(gnss.tropo_sigma_rw_m * np.sqrt(dt))

        if slr is not None:
            for k in range(len(slr.values)):
                ti = int(slr_time_effective[k])
                si = slr.station_indices[k]
                r_sc = mean[ti, 0:3]
                r_st = slr.station_positions_eci[ti, si]
                rho_vec = r_sc - r_st
                rho = np.linalg.norm(rho_vec)
                rho_hat = rho_vec / rho

                pred = rho
                if offsets.get("slr_bias") is not None:
                    pred += p[offsets["slr_bias"] + si]

                r = slr.values[k] - pred
                row = np.zeros(n_params)
                row[:n_state] = rho_hat @ stm[ti][0:3, :]
                if offsets.get("slr_bias") is not None:
                    row[offsets["slr_bias"] + si] = 1.0
                rows.append(row)
                resids.append(r)
                sigmas.append(slr.sigmas[k])

        if not rows:
            raise ValueError("no measurements provided for OD")

        A = np.vstack(rows)
        r = np.array(resids)
        w = 1.0 / (np.array(sigmas) ** 2)
        Aw = A * w[:, None]
        N = Aw.T @ A
        b = Aw.T @ r
        if use_prior:
            try:
                P0_inv = np.linalg.inv(P0_guess)
            except np.linalg.LinAlgError:
                P0_inv = np.linalg.pinv(P0_guess)
            N[:n_state, :n_state] += P0_inv
            b[:n_state] += P0_inv @ (x0_guess - p[:n_state])
        N += damping * np.eye(n_params)

        try:
            dx = np.linalg.solve(N, b)
        except np.linalg.LinAlgError:
            dx = np.linalg.lstsq(N, b, rcond=None)[0]

        p += dx
        normalize_quaternion_state(p[:n_state])

        if np.linalg.norm(dx) < tol:
            break

    try:
        P_post = np.linalg.inv(N)
    except np.linalg.LinAlgError:
        P_post = np.linalg.pinv(N)

    nuisance = {}
    nuisance["measurement_latency_s"] = float(measurement_latency_s)
    nuisance["measurement_latency_jitter_s"] = float(measurement_latency_jitter_s)
    nuisance["measurement_latency_seed"] = int(measurement_latency_seed)
    if gnss_time_effective is not None:
        nuisance["gnss_time_effective"] = np.array(gnss_time_effective, dtype=int)
    if slr_time_effective is not None:
        nuisance["slr_time_effective"] = np.array(slr_time_effective, dtype=int)
    if offsets.get("clock_bias") is not None:
        nuisance["clock_bias"] = p[offsets["clock_bias"] : offsets["clock_bias"] + len(gnss_times)]
    if offsets.get("clock_drift") is not None:
        nuisance["clock_drift"] = p[offsets["clock_drift"] : offsets["clock_drift"] + len(gnss_times)]
    if offsets.get("tropo") is not None:
        nuisance["tropo"] = p[offsets["tropo"] : offsets["tropo"] + len(gnss_times)]
    if offsets.get("ambiguity") is not None:
        n_amb = gnss.ambiguities_cycles.shape[0]
        nuisance["ambiguity"] = p[offsets["ambiguity"] : offsets["ambiguity"] + n_amb]
    if offsets.get("slr_bias") is not None:
        nuisance["slr_bias"] = p[offsets["slr_bias"] : offsets["slr_bias"] + len(slr.stations)]

    return p[:n_state], P_post[:n_state, :n_state], nuisance


def _sample_ensemble(
    x0: np.ndarray,
    P0: np.ndarray,
    members: int,
    rng: np.random.Generator,
) -> np.ndarray:
    if members < 1:
        raise ValueError("members must be >= 1")
    if members == 1:
        out = x0[None, :].copy()
        normalize_quaternion_state(out[0])
        return out
    P = 0.5 * (P0 + P0.T)
    jitter = 0.0
    for _ in range(6):
        try:
            L = np.linalg.cholesky(P + jitter * np.eye(P.shape[0]))
            break
        except np.linalg.LinAlgError:
            diag_max = float(np.max(np.diag(P))) if P.size else 1.0
            jitter = (1e-12 if jitter == 0.0 else jitter * 10.0) * max(1.0, diag_max)
    else:
        evals, evecs = np.linalg.eigh(P)
        eps = 1e-12 * max(1.0, float(np.max(evals)))
        evals = np.maximum(evals, eps)
        L = evecs @ np.diag(np.sqrt(evals))
    out = x0 + rng.standard_normal((members, x0.size)) @ L.T
    q_ref = x0[6:10].copy()
    if np.linalg.norm(q_ref) == 0.0:
        q_ref = np.array([1.0, 0.0, 0.0, 0.0])
    for i in range(out.shape[0]):
        q = out[i, 6:10]
        n = np.linalg.norm(q)
        if n == 0.0:
            q = q_ref.copy()
        else:
            q = q / n
        if np.dot(q, q_ref) < 0.0:
            q = -q
        out[i, 6:10] = q
    return out


def _ensemble_mean_state(X: np.ndarray) -> np.ndarray:
    mean = np.mean(X, axis=0)
    q_ref = X[0, 6:10]
    q = X[:, 6:10].copy()
    for i in range(q.shape[0]):
        if np.dot(q[i], q_ref) < 0.0:
            q[i] = -q[i]
    q_mean = np.mean(q, axis=0)
    nq = np.linalg.norm(q_mean)
    if nq == 0.0:
        q_mean = np.array([1.0, 0.0, 0.0, 0.0])
    else:
        q_mean = q_mean / nq
    mean[6:10] = q_mean
    return mean


def _state_covariance(X: np.ndarray) -> np.ndarray:
    if X.shape[0] <= 1:
        return np.zeros((X.shape[1], X.shape[1]), dtype=float)
    return np.cov(X, rowvar=False)


def _enkf_scalar_update(
    A: np.ndarray,
    y_pred: np.ndarray,
    z: float,
    sigma: float,
    rng: np.random.Generator,
    inflation: float,
) -> np.ndarray:
    if sigma <= 0.0:
        raise ValueError("measurement sigma must be positive")
    if A.shape[0] <= 1:
        return A
    if inflation > 0.0 and inflation != 1.0:
        mean_a = np.mean(A, axis=0)
        A = mean_a + inflation * (A - mean_a)
    xa = A - np.mean(A, axis=0)
    ya = y_pred - float(np.mean(y_pred))
    pyy = float((ya @ ya) / (A.shape[0] - 1) + sigma * sigma)
    if not np.isfinite(pyy) or pyy <= 0.0:
        return A
    pxy = (xa.T @ ya) / (A.shape[0] - 1)
    K = pxy / pyy
    perturb = rng.normal(0.0, sigma, size=A.shape[0])
    innovation = float(z) + perturb - y_pred
    return A + np.outer(innovation, K)


def run_pod_uq_measurements_enkf(
    prop_det: DeterministicPropagator,
    x0_truth: np.ndarray,
    x0_guess: np.ndarray,
    P0_guess: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence,
    gnss: Optional[GnssMeasurements] = None,
    slr: Optional[SlrMeasurements] = None,
    members: int = 64,
    seed: int = 0,
    estimate_clock_bias: bool = True,
    estimate_clock_drift: bool = False,
    estimate_tropo: bool = False,
    estimate_ambiguity: bool = False,
    estimate_slr_bias: bool = False,
    use_carrier: bool = False,
    inflation: float = 1.0,
    measurement_latency_s: float = 0.0,
    measurement_latency_jitter_s: float = 0.0,
    measurement_latency_seed: int = 0,
) -> PodUqMeasurementResult:
    if gnss is None and slr is None:
        raise ValueError("no measurements provided for OD")
    if members < 2:
        raise ValueError("EnKF requires at least two members")
    n_state = prop_det.state_size
    if P0_guess.shape != (n_state, n_state):
        raise ValueError("P0_guess must match propagator state size")

    truth_states = prop_det.propagate(x0_truth, t_grid, env)
    rng = np.random.default_rng(seed)
    rng_latency = np.random.default_rng(int(measurement_latency_seed))
    X = _sample_ensemble(x0_guess, P0_guess, members, rng)
    gnss_time_effective = None
    slr_time_effective = None
    if gnss is not None:
        gnss_time_effective = apply_measurement_latency(
            t_grid,
            gnss.time_indices,
            latency_s=float(measurement_latency_s),
            jitter_s=float(measurement_latency_jitter_s),
            rng=rng_latency,
        )
    if slr is not None:
        slr_time_effective = apply_measurement_latency(
            t_grid,
            slr.time_indices,
            latency_s=float(measurement_latency_s),
            jitter_s=float(measurement_latency_jitter_s),
            rng=rng_latency,
        )

    offsets = {"state": 0}
    d = n_state
    has_carrier = gnss is not None and np.any(gnss.types == MEAS_CARRIER)
    if gnss is not None and estimate_clock_bias:
        offsets["clock_bias"] = d
        d += 1
    else:
        offsets["clock_bias"] = None
    if gnss is not None and estimate_clock_drift:
        offsets["clock_drift"] = d
        d += 1
    else:
        offsets["clock_drift"] = None
    if gnss is not None and estimate_tropo:
        offsets["tropo"] = d
        d += 1
    else:
        offsets["tropo"] = None
    if gnss is not None and estimate_ambiguity and has_carrier:
        offsets["ambiguity"] = d
        d += int(gnss.ambiguities_cycles.shape[0])
    else:
        offsets["ambiguity"] = None
    if slr is not None and estimate_slr_bias:
        offsets["slr_bias"] = d
        d += len(slr.stations)
    else:
        offsets["slr_bias"] = None

    A = np.zeros((members, d), dtype=float)
    A[:, :n_state] = X
    if offsets["ambiguity"] is not None:
        n_amb = int(gnss.ambiguities_cycles.shape[0])
        idx0 = int(offsets["ambiguity"])
        # Start near zero with broad spread; carrier data will contract these states.
        A[:, idx0 : idx0 + n_amb] = rng.normal(0.0, 1e5, size=(members, n_amb))
    if offsets["slr_bias"] is not None:
        idx0 = int(offsets["slr_bias"])
        A[:, idx0 : idx0 + len(slr.stations)] = rng.normal(0.0, 0.1, size=(members, len(slr.stations)))

    gnss_by_time: dict[int, list[int]] = {}
    if gnss is not None:
        for k in range(len(gnss.values)):
            ti = int(gnss_time_effective[k])
            gnss_by_time.setdefault(ti, []).append(k)
    slr_by_time: dict[int, list[int]] = {}
    if slr is not None:
        for k in range(len(slr.values)):
            ti = int(slr_time_effective[k])
            slr_by_time.setdefault(ti, []).append(k)

    est_states = np.zeros((len(t_grid), n_state), dtype=float)
    state_cov = np.zeros((len(t_grid), n_state, n_state), dtype=float)

    for i in range(len(t_grid)):
        if i > 0:
            dt_s = float(t_grid[i] - t_grid[i - 1])
            if dt_s <= 0.0:
                raise ValueError("t_grid must be strictly increasing")
            seg_t = np.array([t_grid[i - 1], t_grid[i]], dtype=float)
            seg_env = [env[i - 1], env[i]]
            for p in range(members):
                x_prev = A[p, :n_state]
                x_next = prop_det.propagate(x_prev, seg_t, seg_env)[-1]
                normalize_quaternion_state(x_next)
                A[p, :n_state] = x_next

            # Propagate nuisance states as random walks.
            if gnss is not None and offsets["clock_drift"] is not None and gnss.clock_drift_sigma_rw_m_s > 0.0:
                idx = int(offsets["clock_drift"])
                A[:, idx] += gnss.clock_drift_sigma_rw_m_s * np.sqrt(dt_s) * rng.normal(size=members)
            if gnss is not None and offsets["clock_bias"] is not None and gnss.clock_bias_sigma_rw_m > 0.0:
                idx_b = int(offsets["clock_bias"])
                if offsets["clock_drift"] is not None:
                    idx_d = int(offsets["clock_drift"])
                    A[:, idx_b] += A[:, idx_d] * dt_s
                A[:, idx_b] += gnss.clock_bias_sigma_rw_m * np.sqrt(dt_s) * rng.normal(size=members)
            if gnss is not None and offsets["tropo"] is not None and gnss.tropo_sigma_rw_m > 0.0:
                idx_t = int(offsets["tropo"])
                A[:, idx_t] += gnss.tropo_sigma_rw_m * np.sqrt(dt_s) * rng.normal(size=members)

        # GNSS updates at epoch i.
        if gnss is not None and i in gnss_by_time:
            for k in gnss_by_time[i]:
                typ = int(gnss.types[k])
                if typ == MEAS_CARRIER and not use_carrier:
                    continue
                if typ == MEAS_CARRIER and offsets["ambiguity"] is None:
                    continue
                si = int(gnss.sat_indices[k])
                fi = int(gnss.freq_indices[k])
                r_sat = gnss.sat_positions_eci[i, si]
                if gnss.sat_pco_eci_m is not None:
                    r_sat = r_sat + gnss.sat_pco_eci_m[i, si]
                rho_vec = A[:, 0:3] - r_sat.reshape(1, 3)
                y_pred = np.linalg.norm(rho_vec, axis=1)
                if offsets["clock_bias"] is not None:
                    y_pred = y_pred + A[:, int(offsets["clock_bias"])]
                if offsets["tropo"] is not None:
                    y_pred = y_pred + A[:, int(offsets["tropo"])]
                if gnss.sat_clock_bias_m is not None:
                    y_pred = y_pred + float(gnss.sat_clock_bias_m[i, si])
                if typ == MEAS_CARRIER:
                    amb_id = int(gnss.ambiguity_ids[k])
                    if amb_id >= 0:
                        idx_a = int(offsets["ambiguity"]) + amb_id
                        y_pred = y_pred + A[:, idx_a] * float(gnss.wavelengths_m[fi])
                A = _enkf_scalar_update(
                    A,
                    y_pred=y_pred,
                    z=float(gnss.values[k]),
                    sigma=float(gnss.sigmas[k]),
                    rng=rng,
                    inflation=inflation,
                )
                for p in range(members):
                    normalize_quaternion_state(A[p, :n_state])

        # SLR updates at epoch i.
        if slr is not None and i in slr_by_time:
            for k in slr_by_time[i]:
                si = int(slr.station_indices[k])
                r_st = slr.station_positions_eci[i, si]
                rho_vec = A[:, 0:3] - r_st.reshape(1, 3)
                y_pred = np.linalg.norm(rho_vec, axis=1)
                if offsets["slr_bias"] is not None:
                    y_pred = y_pred + A[:, int(offsets["slr_bias"]) + si]
                A = _enkf_scalar_update(
                    A,
                    y_pred=y_pred,
                    z=float(slr.values[k]),
                    sigma=float(slr.sigmas[k]),
                    rng=rng,
                    inflation=inflation,
                )
                for p in range(members):
                    normalize_quaternion_state(A[p, :n_state])

        est_states[i] = _ensemble_mean_state(A[:, :n_state])
        state_cov[i] = _state_covariance(A[:, :n_state])

    x0_est = est_states[0].copy()
    P0_post = state_cov[0].copy()

    rtn_error = compute_rtn_errors(truth_states, est_states)
    rtn_sigma = np.zeros((len(t_grid), 3), dtype=float)
    for i in range(len(t_grid)):
        P_pos = state_cov[i, 0:3, 0:3]
        C = rtn_basis(truth_states[i, 0:3], truth_states[i, 3:6])
        P_rtn = C.T @ P_pos @ C
        rtn_sigma[i] = np.sqrt(np.maximum(np.diag(P_rtn), 0.0))
    radial_rms_m = float(np.sqrt(np.mean(rtn_error[:, 0] ** 2)))
    coverage_3sigma = compute_coverage_3sigma(rtn_error, rtn_sigma)
    coverage_sigma = compute_coverage_sigma_levels(rtn_error, rtn_sigma)

    measurements = {}
    if gnss is not None:
        measurements["gnss"] = gnss
    if slr is not None:
        measurements["slr"] = slr

    nuisance: dict[str, np.ndarray | float] = {
        "state_cov_final": state_cov[-1],
    }
    if offsets["clock_bias"] is not None:
        nuisance["clock_bias"] = float(np.mean(A[:, int(offsets["clock_bias"])]))
    if offsets["clock_drift"] is not None:
        nuisance["clock_drift"] = float(np.mean(A[:, int(offsets["clock_drift"])]))
    if offsets["tropo"] is not None:
        nuisance["tropo"] = float(np.mean(A[:, int(offsets["tropo"])]))
    if offsets["ambiguity"] is not None:
        n_amb = int(gnss.ambiguities_cycles.shape[0])
        idx0 = int(offsets["ambiguity"])
        nuisance["ambiguity"] = np.mean(A[:, idx0 : idx0 + n_amb], axis=0)
    if offsets["slr_bias"] is not None:
        idx0 = int(offsets["slr_bias"])
        nuisance["slr_bias"] = np.mean(A[:, idx0 : idx0 + len(slr.stations)], axis=0)

    return PodUqMeasurementResult(
        x0_truth=x0_truth,
        x0_est=x0_est,
        P0_post=P0_post,
        truth_states=truth_states,
        est_states=est_states,
        measurements=measurements,
        rtn_sigma=rtn_sigma,
        rtn_error=rtn_error,
        radial_rms_m=radial_rms_m,
        coverage_3sigma=coverage_3sigma,
        coverage_sigma=coverage_sigma,
        nuisance=nuisance,
        t_grid=np.array(t_grid, dtype=float),
    )


def run_pod_uq_measurements(
    prop_det: DeterministicPropagator,
    prop_stm: StmPropagator,
    x0_truth: np.ndarray,
    x0_guess: np.ndarray,
    P0_guess: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence,
    gnss: Optional[GnssMeasurements] = None,
    slr: Optional[SlrMeasurements] = None,
    max_iter: int = 6,
    estimate_clock_bias: bool = True,
    estimate_clock_drift: bool = False,
    estimate_tropo: bool = False,
    estimate_ambiguity: bool = True,
    estimate_slr_bias: bool = False,
    use_prior: bool = True,
    measurement_latency_s: float = 0.0,
    measurement_latency_jitter_s: float = 0.0,
    measurement_latency_seed: int = 0,
) -> PodUqMeasurementResult:
    truth_states = prop_det.propagate(x0_truth, t_grid, env)

    x0_est, P0_post, nuisance = batch_od_measurements(
        prop_stm,
        x0_guess,
        P0_guess,
        t_grid,
        env,
        gnss=gnss,
        slr=slr,
        max_iter=max_iter,
        estimate_clock_bias=estimate_clock_bias,
        estimate_clock_drift=estimate_clock_drift,
        estimate_tropo=estimate_tropo,
        estimate_ambiguity=estimate_ambiguity,
        estimate_slr_bias=estimate_slr_bias,
        use_prior=use_prior,
        measurement_latency_s=measurement_latency_s,
        measurement_latency_jitter_s=measurement_latency_jitter_s,
        measurement_latency_seed=measurement_latency_seed,
    )

    est_states = prop_det.propagate(x0_est, t_grid, env)
    _, _, stm = prop_stm.propagate(x0_est, P0_post, t_grid, env, return_stm=True)

    rtn_error = compute_rtn_errors(truth_states, est_states)
    rtn_sigma = compute_rtn_sigmas(P0_post, stm, truth_states)
    radial_rms_m = float(np.sqrt(np.mean(rtn_error[:, 0] ** 2)))
    coverage_3sigma = compute_coverage_3sigma(rtn_error, rtn_sigma)
    coverage_sigma = compute_coverage_sigma_levels(rtn_error, rtn_sigma)

    measurements = {}
    if gnss is not None:
        measurements["gnss"] = gnss
    if slr is not None:
        measurements["slr"] = slr

    return PodUqMeasurementResult(
        x0_truth=x0_truth,
        x0_est=x0_est,
        P0_post=P0_post,
        truth_states=truth_states,
        est_states=est_states,
        measurements=measurements,
        rtn_sigma=rtn_sigma,
        rtn_error=rtn_error,
        radial_rms_m=radial_rms_m,
        coverage_3sigma=coverage_3sigma,
        coverage_sigma=coverage_sigma,
        nuisance=nuisance,
        t_grid=np.array(t_grid, dtype=float),
    )


def run_pod_uq_multi_arc(
    prop_det: DeterministicPropagator,
    prop_stm: StmPropagator,
    x0_truth: np.ndarray,
    x0_guess: np.ndarray,
    P0_guess: np.ndarray,
    t_grid: np.ndarray,
    env: Sequence,
    arcs: Sequence[tuple[int, int]],
    gnss_sets: Optional[Sequence[Optional[GnssMeasurements]]] = None,
    slr_sets: Optional[Sequence[Optional[SlrMeasurements]]] = None,
    max_iter: int = 6,
    couple_arcs: bool = True,
    carry_covariance: bool = True,
    estimator: str = "batch",
    enkf_members: int = 64,
    enkf_seed: int = 0,
    enkf_inflation: float = 1.0,
    enkf_use_carrier: bool = False,
    measurement_latency_s: float = 0.0,
    measurement_latency_jitter_s: float = 0.0,
    measurement_latency_seed: int = 0,
) -> list[PodUqMeasurementResult]:
    truth_full = prop_det.propagate(x0_truth, t_grid, env)
    guess_full = prop_det.propagate(x0_guess, t_grid, env)
    if len(arcs) == 0:
        return []
    if couple_arcs:
        for i in range(1, len(arcs)):
            if arcs[i][0] < arcs[i - 1][0]:
                raise ValueError("arcs must be sorted by start index when couple_arcs=True")
    if gnss_sets is not None and len(gnss_sets) != len(arcs):
        raise ValueError("gnss_sets length must match arcs")
    if slr_sets is not None and len(slr_sets) != len(arcs):
        raise ValueError("slr_sets length must match arcs")

    estimator = estimator.lower().strip()
    if estimator not in {"batch", "enkf"}:
        raise ValueError("estimator must be 'batch' or 'enkf'")

    prev_start = None
    prev_x0_est = None
    prev_P0_post = None

    results = []
    for i, (start, end) in enumerate(arcs):
        if end <= start:
            raise ValueError("arc end must be greater than start")
        t_arc = t_grid[start:end]
        env_arc = env[start:end]
        if gnss_sets is not None:
            gnss = gnss_sets[i]
        else:
            gnss = None
        if slr_sets is not None:
            slr = slr_sets[i]
        else:
            slr = None
        x0_truth_arc = truth_full[start].copy()
        x0_guess_arc = guess_full[start].copy()
        P0_arc = P0_guess
        if couple_arcs and i > 0 and prev_x0_est is not None and prev_P0_post is not None and prev_start is not None:
            link_start = prev_start
            link_end = start
            if link_end < link_start:
                raise ValueError("arc start indices must be non-decreasing for coupled mode")
            if link_end == link_start:
                x0_guess_arc = prev_x0_est.copy()
                if carry_covariance:
                    P0_arc = prev_P0_post.copy()
            else:
                t_link = t_grid[link_start : link_end + 1]
                env_link = env[link_start : link_end + 1]
                mean_link, cov_link = prop_stm.propagate(
                    prev_x0_est, prev_P0_post, t_link, env_link, return_stm=False
                )
                x0_guess_arc = mean_link[-1].copy()
                if carry_covariance:
                    P0_arc = cov_link[-1].copy()
        if estimator == "enkf":
            result = run_pod_uq_measurements_enkf(
                prop_det=prop_det,
                x0_truth=x0_truth_arc,
                x0_guess=x0_guess_arc,
                P0_guess=P0_arc,
                t_grid=t_arc,
                env=env_arc,
                gnss=gnss,
                slr=slr,
                members=enkf_members,
                seed=enkf_seed + i,
                estimate_clock_bias=True,
                estimate_clock_drift=False,
                estimate_tropo=False,
                estimate_ambiguity=enkf_use_carrier,
                estimate_slr_bias=False,
                use_carrier=enkf_use_carrier,
                inflation=enkf_inflation,
                measurement_latency_s=measurement_latency_s,
                measurement_latency_jitter_s=measurement_latency_jitter_s,
                measurement_latency_seed=measurement_latency_seed + i,
            )
            # Keep the same arc-coupling convention as batch OD (estimate at arc start),
            # which also behaves well when arcs overlap.
            prev_start = start
            prev_x0_est = result.x0_est.copy()
            prev_P0_post = result.P0_post.copy()
        else:
            result = run_pod_uq_measurements(
                prop_det,
                prop_stm,
                x0_truth_arc,
                x0_guess_arc,
                P0_arc,
                t_arc,
                env_arc,
                gnss=gnss,
                slr=slr,
                max_iter=max_iter,
                measurement_latency_s=measurement_latency_s,
                measurement_latency_jitter_s=measurement_latency_jitter_s,
                measurement_latency_seed=measurement_latency_seed + i,
            )
            prev_start = start
            prev_x0_est = result.x0_est.copy()
            prev_P0_post = result.P0_post.copy()
        results.append(result)
    return results


def default_ilrs_stations() -> list[GroundStation]:
    return [
        GroundStation("ZIML", np.deg2rad(46.8772), np.deg2rad(7.4652), 951.2),
        GroundStation("YARL", np.deg2rad(-29.0464), np.deg2rad(115.3467), 244.0),
        GroundStation("HA4T", np.deg2rad(20.706489), np.deg2rad(-156.256921), 3056.272),
        GroundStation("GODL", np.deg2rad(39.0206), np.deg2rad(-76.82770), 19.184),
        GroundStation("YEBL", np.deg2rad(40.5245), np.deg2rad(-3.0905), 974.335),
    ]


def load_stations_csv(path: str) -> list[GroundStation]:
    stations = []
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        for line in handle:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = [p.strip() for p in raw.split(",")]
            if header is None and any(p.lower() in {"name", "lat_deg", "lon_deg", "lat_rad"} for p in parts):
                header = [p.lower() for p in parts]
                continue
            if header is None:
                header = ["name", "lat_deg", "lon_deg", "alt_m"]
            row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
            name = row.get("name", f"STA{len(stations)+1}")
            if "lat_rad" in row or "lon_rad" in row:
                lat = float(row.get("lat_rad", "0"))
                lon = float(row.get("lon_rad", "0"))
            else:
                lat = np.deg2rad(float(row.get("lat_deg", "0")))
                lon = np.deg2rad(float(row.get("lon_deg", "0")))
            alt = float(row.get("alt_m", "0") or 0.0)
            stations.append(GroundStation(name, lat, lon, alt))
    if not stations:
        raise ValueError("no stations parsed from CSV")
    return stations


def _decimal_year_to_datetime(year_dec: float) -> dt.datetime:
    year = int(year_dec)
    rem = year_dec - year
    start = dt.datetime(year, 1, 1)
    end = dt.datetime(year + 1, 1, 1)
    return start + (end - start) * rem


def load_stations_itrf_csv(path: str, t0: dt.datetime) -> list[GroundStation]:
    stations = []
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        for line in handle:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = [p.strip() for p in raw.split(",")]
            if header is None and any(p.lower() in {"x_m", "y_m", "z_m", "vx_m_s"} for p in parts):
                header = [p.lower() for p in parts]
                continue
            if header is None:
                header = ["name", "x_m", "y_m", "z_m", "vx_m_s", "vy_m_s", "vz_m_s", "epoch_year"]
            row = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
            name = row.get("name", f"STA{len(stations)+1}")
            x = float(row.get("x_m", "0"))
            y = float(row.get("y_m", "0"))
            z = float(row.get("z_m", "0"))
            vx = float(row.get("vx_m_s", "0") or 0.0)
            vy = float(row.get("vy_m_s", "0") or 0.0)
            vz = float(row.get("vz_m_s", "0") or 0.0)
            epoch_raw = row.get("epoch_year", "") or row.get("epoch", "")
            if epoch_raw:
                epoch = _decimal_year_to_datetime(float(epoch_raw))
            else:
                epoch = t0
            dt_s = (t0 - epoch).total_seconds()
            x_t = x + vx * dt_s
            y_t = y + vy * dt_s
            z_t = z + vz * dt_s
            lat, lon, alt = ecef_to_geodetic(x_t, y_t, z_t)
            stations.append(GroundStation(name, lat, lon, alt))
    if not stations:
        raise ValueError("no ITRF stations parsed from CSV")
    return stations


def _residual_stats(values_m: np.ndarray, sigmas_m: np.ndarray) -> dict | None:
    vals = np.asarray(values_m, dtype=float).reshape(-1)
    sig = np.asarray(sigmas_m, dtype=float).reshape(-1)
    if vals.size == 0 or sig.size == 0 or vals.size != sig.size:
        return None
    valid = np.isfinite(vals) & np.isfinite(sig) & (sig > 0.0)
    if not np.any(valid):
        return None
    v = vals[valid]
    s = sig[valid]
    n = v.size
    vn = v / s
    std_m = float(np.std(v, ddof=1)) if n > 1 else 0.0
    std_n = float(np.std(vn, ddof=1)) if n > 1 else 0.0
    return {
        "count": int(n),
        "mean_m": float(np.mean(v)),
        "std_m": std_m,
        "rms_m": float(np.sqrt(np.mean(v * v))),
        "mean_norm": float(np.mean(vn)),
        "std_norm": std_n,
        "rms_norm": float(np.sqrt(np.mean(vn * vn))),
        "p95_abs_norm": float(np.percentile(np.abs(vn), 95.0)),
    }


def _time_mapped_nuisance(value, ti: int, time_map: dict[int, int]) -> float:
    if value is None:
        return 0.0
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size == 0:
        return 0.0
    if arr.size == 1:
        return float(arr[0])
    idx = time_map.get(int(ti))
    if idx is not None and 0 <= idx < arr.size:
        return float(arr[idx])
    idx_fallback = int(np.clip(int(ti), 0, arr.size - 1))
    return float(arr[idx_fallback])


def measurement_residual_statistics(result: PodUqMeasurementResult) -> dict:
    out: dict[str, dict] = {}
    if not isinstance(result.measurements, dict):
        return out

    nuisance = result.nuisance if isinstance(result.nuisance, dict) else {}
    est_states = np.asarray(result.est_states, dtype=float)

    gnss = result.measurements.get("gnss")
    if isinstance(gnss, GnssMeasurements) and gnss.values.size > 0:
        gnss_time_map = {int(t): i for i, t in enumerate(np.unique(gnss.time_indices))}
        cb = nuisance.get("clock_bias")
        tr = nuisance.get("tropo")
        amb = np.asarray(nuisance.get("ambiguity", []), dtype=float).reshape(-1)

        gnss_res = np.zeros(gnss.values.size, dtype=float)
        for k in range(gnss.values.size):
            ti = int(gnss.time_indices[k])
            si = int(gnss.sat_indices[k])
            fi = int(gnss.freq_indices[k])
            r_sc = est_states[ti, 0:3]
            r_sat = np.asarray(gnss.sat_positions_eci[ti, si], dtype=float)
            if gnss.sat_pco_eci_m is not None:
                r_sat = r_sat + np.asarray(gnss.sat_pco_eci_m[ti, si], dtype=float)
            pred = float(np.linalg.norm(r_sc - r_sat))
            pred += _time_mapped_nuisance(cb, ti, gnss_time_map)
            pred += _time_mapped_nuisance(tr, ti, gnss_time_map)
            if gnss.sat_clock_bias_m is not None:
                pred += float(gnss.sat_clock_bias_m[ti, si])
            if int(gnss.types[k]) == MEAS_CARRIER:
                amb_id = int(gnss.ambiguity_ids[k])
                if amb_id >= 0 and amb_id < amb.size:
                    pred += float(amb[amb_id] * gnss.wavelengths_m[fi])
            gnss_res[k] = float(gnss.values[k] - pred)

        stats_all = _residual_stats(gnss_res, gnss.sigmas)
        if stats_all is not None:
            out["gnss_all"] = stats_all
        code_mask = gnss.types == MEAS_CODE
        if np.any(code_mask):
            stats_code = _residual_stats(gnss_res[code_mask], gnss.sigmas[code_mask])
            if stats_code is not None:
                out["gnss_code"] = stats_code
        carrier_mask = gnss.types == MEAS_CARRIER
        if np.any(carrier_mask):
            stats_carrier = _residual_stats(gnss_res[carrier_mask], gnss.sigmas[carrier_mask])
            if stats_carrier is not None:
                out["gnss_carrier"] = stats_carrier

    slr = result.measurements.get("slr")
    if isinstance(slr, SlrMeasurements) and slr.values.size > 0:
        slr_bias = np.asarray(nuisance.get("slr_bias", []), dtype=float).reshape(-1)
        slr_res = np.zeros(slr.values.size, dtype=float)
        for k in range(slr.values.size):
            ti = int(slr.time_indices[k])
            si = int(slr.station_indices[k])
            r_sc = est_states[ti, 0:3]
            r_st = np.asarray(slr.station_positions_eci[ti, si], dtype=float)
            pred = float(np.linalg.norm(r_sc - r_st))
            if slr_bias.size == 1:
                pred += float(slr_bias[0])
            elif slr_bias.size > si:
                pred += float(slr_bias[si])
            slr_res[k] = float(slr.values[k] - pred)
        stats_slr = _residual_stats(slr_res, slr.sigmas)
        if stats_slr is not None:
            out["slr"] = stats_slr

    return out


def summarize_arc_metrics(
    result: PodUqMeasurementResult,
    horizons_s: Sequence[float] = (6 * 3600.0, 24 * 3600.0),
    t_grid: Optional[np.ndarray] = None,
) -> dict:
    if t_grid is None:
        if result.t_grid is None:
            raise ValueError("t_grid is required for arc summary")
        t_grid = result.t_grid
    t_grid = np.array(t_grid, dtype=float)
    summaries = {
        "radial_rms_m": result.radial_rms_m,
        "coverage_3sigma": result.coverage_3sigma,
        "coverage_sigma": result.coverage_sigma,
        "measurement_residuals": measurement_residual_statistics(result),
        "horizons_s": [],
        "rtn_sigma": [],
        "rtn_error": [],
    }
    for h in horizons_s:
        if h < t_grid[0] or h > t_grid[-1]:
            continue
        idx = int(np.argmin(np.abs(t_grid - h)))
        summaries["horizons_s"].append(float(t_grid[idx]))
        summaries["rtn_sigma"].append(result.rtn_sigma[idx].tolist())
        summaries["rtn_error"].append(result.rtn_error[idx].tolist())
    return summaries


def summarize_multi_arc(
    results: Sequence[PodUqMeasurementResult],
    horizons_s: Sequence[float] = (6 * 3600.0, 24 * 3600.0),
) -> list[dict]:
    return [summarize_arc_metrics(r, horizons_s=horizons_s) for r in results]


def slice_gnss_measurements(
    gnss: GnssMeasurements,
    start_idx: int,
    end_idx: int,
) -> GnssMeasurements:
    mask = (gnss.time_indices >= start_idx) & (gnss.time_indices < end_idx)
    if not np.any(mask):
        raise ValueError("no GNSS measurements in arc")
    time_indices = gnss.time_indices[mask] - start_idx
    sat_positions = gnss.sat_positions_eci[start_idx:end_idx]
    clock_bias = gnss.clock_bias_truth[start_idx:end_idx]
    clock_drift = gnss.clock_drift_truth[start_idx:end_idx]
    tropo = gnss.tropo_truth[start_idx:end_idx]
    sat_clock_bias_m = None
    if gnss.sat_clock_bias_m is not None:
        sat_clock_bias_m = gnss.sat_clock_bias_m[start_idx:end_idx]
    sat_pco_eci_m = None
    if gnss.sat_pco_eci_m is not None:
        sat_pco_eci_m = gnss.sat_pco_eci_m[start_idx:end_idx]
    # Compact ambiguity IDs to the subset actually used in this arc.
    amb_ids = gnss.ambiguity_ids[mask].copy()
    amb_cycles = gnss.ambiguities_cycles
    valid = amb_ids >= 0
    if np.any(valid):
        used = np.unique(amb_ids[valid])
        remap = {int(old): i for i, old in enumerate(used.tolist())}
        for old, new in remap.items():
            amb_ids[amb_ids == old] = new
        amb_cycles = gnss.ambiguities_cycles[used]
    else:
        amb_cycles = np.zeros((0,), dtype=int)

    return GnssMeasurements(
        values=gnss.values[mask],
        types=gnss.types[mask],
        time_indices=time_indices,
        sat_indices=gnss.sat_indices[mask],
        freq_indices=gnss.freq_indices[mask],
        sigmas=gnss.sigmas[mask],
        sat_positions_eci=sat_positions,
        wavelengths_m=gnss.wavelengths_m,
        frequencies_hz=gnss.frequencies_hz,
        clock_bias_truth=clock_bias,
        clock_drift_truth=clock_drift,
        tropo_truth=tropo,
        ambiguity_ids=amb_ids,
        ambiguities_cycles=amb_cycles,
        clock_bias_sigma_rw_m=gnss.clock_bias_sigma_rw_m,
        clock_drift_sigma_rw_m_s=gnss.clock_drift_sigma_rw_m_s,
        tropo_sigma_rw_m=gnss.tropo_sigma_rw_m,
        sat_clock_bias_m=sat_clock_bias_m,
        sat_pco_eci_m=sat_pco_eci_m,
    )


def slice_slr_measurements(
    slr: SlrMeasurements,
    start_idx: int,
    end_idx: int,
) -> SlrMeasurements:
    mask = (slr.time_indices >= start_idx) & (slr.time_indices < end_idx)
    if not np.any(mask):
        raise ValueError("no SLR measurements in arc")
    time_indices = slr.time_indices[mask] - start_idx
    station_positions = slr.station_positions_eci[start_idx:end_idx]
    return SlrMeasurements(
        values=slr.values[mask],
        time_indices=time_indices,
        station_indices=slr.station_indices[mask],
        sigmas=slr.sigmas[mask],
        stations=slr.stations,
        station_positions_eci=station_positions,
        range_bias_m=slr.range_bias_m,
    )
