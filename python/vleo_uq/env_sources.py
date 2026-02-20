import csv
import ctypes
import datetime as dt
import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np
import json
import re


class DataSourceError(RuntimeError):
    pass


@dataclass(frozen=True)
class DataProvenance:
    name: str
    url: str
    retrieved_utc: str
    notes: str = ""

    def validate(self) -> None:
        if not self.name or not self.url or not self.retrieved_utc:
            raise DataSourceError("provenance must include name, url, and retrieved_utc")


@dataclass
class TimeSeries:
    times: np.ndarray
    values: np.ndarray
    provenance: DataProvenance
    units: str = ""

    def validate(self) -> None:
        if len(self.times) != len(self.values):
            raise DataSourceError("times and values length mismatch")
        if len(self.times) == 0:
            raise DataSourceError("time series is empty")
        if np.any(~np.isfinite(self.values)):
            raise DataSourceError("time series contains non-finite values")
        self.provenance.validate()


def _parse_time(value: str) -> np.datetime64:
    if value.endswith("Z"):
        value = value[:-1]
    try:
        if len(value) == 7 and value[4] == "-":
            value = value + "-01"
        elif len(value) == 4 and value.isdigit():
            value = value + "-01-01"
        parsed = dt.datetime.fromisoformat(value)
    except ValueError as exc:
        raise DataSourceError(f"invalid time format: {value}") from exc
    if parsed.tzinfo is not None:
        parsed = parsed.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return np.datetime64(parsed)


def _to_datetime64_array(times: Sequence) -> np.ndarray:
    if isinstance(times, np.ndarray) and np.issubdtype(times.dtype, np.datetime64):
        return times
    return np.array([_parse_time(str(t)) for t in times], dtype="datetime64[ns]")


def _interpolate_series(series: TimeSeries, t_grid: Sequence, method: str) -> np.ndarray:
    series.validate()
    if method not in {"linear", "nearest"}:
        raise DataSourceError("interpolation method must be 'linear' or 'nearest'")

    times = series.times.astype("datetime64[s]").astype("int64")
    target = _to_datetime64_array(t_grid).astype("datetime64[s]").astype("int64")

    if method == "linear":
        return np.interp(target, times, series.values)

    idx = np.searchsorted(times, target, side="left")
    idx = np.clip(idx, 0, len(times) - 1)
    left = np.clip(idx - 1, 0, len(times) - 1)
    choose_left = (target - times[left]) <= (times[idx] - target)
    return np.where(choose_left, series.values[left], series.values[idx])


def _datetime64_to_datetime(ts: np.datetime64) -> dt.datetime:
    seconds = ts.astype("datetime64[ns]").astype("int64") / 1e9
    return dt.datetime.utcfromtimestamp(seconds)


def _julian_date(timestamp: dt.datetime) -> float:
    year = timestamp.year
    month = timestamp.month
    day = timestamp.day + (
        timestamp.hour + (timestamp.minute + timestamp.second / 60.0) / 60.0
    ) / 24.0
    if month <= 2:
        year -= 1
        month += 12
    a = year // 100
    b = 2 - a + a // 4
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b - 1524.5
    return jd


def _gmst_rad(timestamp: dt.datetime) -> float:
    jd = _julian_date(timestamp)
    t = (jd - 2451545.0) / 36525.0
    gmst_sec = (
        67310.54841
        + (876600.0 * 3600 + 8640184.812866) * t
        + 0.093104 * t * t
        - 6.2e-6 * t * t * t
    )
    gmst_sec = gmst_sec % 86400.0
    return gmst_sec * (2.0 * math.pi / 86400.0)


def _enu_to_ecef(north: float, east: float, lat_rad: float, lon_rad: float) -> np.ndarray:
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)
    east_vec = np.array([-sin_lon, cos_lon, 0.0], dtype=float)
    north_vec = np.array([-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat], dtype=float)
    return east * east_vec + north * north_vec


def _ecef_to_eci(vec_ecef: np.ndarray, gmst_rad: float) -> np.ndarray:
    cos_g = math.cos(gmst_rad)
    sin_g = math.sin(gmst_rad)
    rot = np.array(
        [
            [cos_g, -sin_g, 0.0],
            [sin_g, cos_g, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    return rot @ vec_ecef


def _sun_position_eci_m(timestamp: dt.datetime) -> np.ndarray:
    jd = _julian_date(timestamp)
    t = (jd - 2451545.0) / 36525.0
    mean_long = math.radians((280.460 + 36000.770 * t) % 360.0)
    mean_anomaly = math.radians((357.528 + 35999.050 * t) % 360.0)
    lambda_ecl = (
        mean_long
        + math.radians(1.915) * math.sin(mean_anomaly)
        + math.radians(0.020) * math.sin(2.0 * mean_anomaly)
    )
    eps = math.radians(23.4393 - 0.0130 * t)
    r_au = 1.00014 - 0.01671 * math.cos(mean_anomaly) - 0.00014 * math.cos(2.0 * mean_anomaly)
    r_m = r_au * 149597870700.0
    return np.array(
        [
            r_m * math.cos(lambda_ecl),
            r_m * math.cos(eps) * math.sin(lambda_ecl),
            r_m * math.sin(eps) * math.sin(lambda_ecl),
        ],
        dtype=float,
    )


def _moon_position_eci_m(timestamp: dt.datetime) -> np.ndarray:
    jd = _julian_date(timestamp)
    d = jd - 2451545.0

    l0 = math.radians((218.316 + 13.176396 * d) % 360.0)
    mm = math.radians((134.963 + 13.064993 * d) % 360.0)
    ms = math.radians((357.529 + 0.98560028 * d) % 360.0)
    darg = math.radians((297.850 + 12.190749 * d) % 360.0)
    farg = math.radians((93.272 + 13.229350 * d) % 360.0)

    lon = (
        l0
        + math.radians(6.289) * math.sin(mm)
        + math.radians(1.274) * math.sin(2.0 * darg - mm)
        + math.radians(0.658) * math.sin(2.0 * darg)
        + math.radians(0.214) * math.sin(2.0 * mm)
        - math.radians(0.186) * math.sin(ms)
    )
    lat = (
        math.radians(5.128) * math.sin(farg)
        + math.radians(0.280) * math.sin(mm + farg)
        + math.radians(0.277) * math.sin(mm - farg)
        + math.radians(0.173) * math.sin(2.0 * darg - farg)
    )
    dist_m = (
        385001000.0
        - 20905000.0 * math.cos(mm)
        - 3699000.0 * math.cos(2.0 * darg - mm)
        - 2956000.0 * math.cos(2.0 * darg)
        - 570000.0 * math.cos(2.0 * mm)
    )

    eps = math.radians(23.4393 - 3.563e-7 * d)
    x_ecl = dist_m * math.cos(lat) * math.cos(lon)
    y_ecl = dist_m * math.cos(lat) * math.sin(lon)
    z_ecl = dist_m * math.sin(lat)
    return np.array(
        [
            x_ecl,
            y_ecl * math.cos(eps) - z_ecl * math.sin(eps),
            y_ecl * math.sin(eps) + z_ecl * math.cos(eps),
        ],
        dtype=float,
    )


def _earth_dipole_field_eci_t(r_eci_m: np.ndarray) -> np.ndarray:
    r = np.asarray(r_eci_m, dtype=float).reshape(3)
    r_norm = float(np.linalg.norm(r))
    if (not np.isfinite(r_norm)) or r_norm <= 0.0:
        return np.zeros(3, dtype=float)
    r_hat = r / r_norm
    # First-order Earth dipole approximation (aligned with +Z in ECI for simplicity).
    m_vec = np.array([0.0, 0.0, 7.94e22], dtype=float)  # [A m^2]
    mu0_over_4pi = 1.0e-7
    return (mu0_over_4pi / (r_norm ** 3)) * (3.0 * np.dot(m_vec, r_hat) * r_hat - m_vec)


def _tide_loading_accel_eci_m_s2(
    r_eci_m: np.ndarray,
    sun_eci_m: np.ndarray,
    moon_eci_m: np.ndarray,
) -> np.ndarray:
    r = np.asarray(r_eci_m, dtype=float).reshape(3)
    rn = float(np.linalg.norm(r))
    if (not np.isfinite(rn)) or rn <= 0.0:
        return np.zeros(3, dtype=float)
    rhat = r / rn

    sun = np.asarray(sun_eci_m, dtype=float).reshape(3)
    moon = np.asarray(moon_eci_m, dtype=float).reshape(3)
    sun_n = float(np.linalg.norm(sun))
    moon_n = float(np.linalg.norm(moon))
    if sun_n <= 0.0 or moon_n <= 0.0:
        return np.zeros(3, dtype=float)
    sun_hat = sun / sun_n
    moon_hat = moon / moon_n

    # Surrogate Earth tide/loading acceleration model:
    # a_tide ~ A * (3 (rhat·dhat) dhat - rhat), summed for Sun and Moon.
    amp_sun = 1.5e-8
    amp_moon = 2.5e-8
    sun_term = amp_sun * (3.0 * float(np.dot(rhat, sun_hat)) * sun_hat - rhat)
    moon_term = amp_moon * (3.0 * float(np.dot(rhat, moon_hat)) * moon_hat - rhat)
    return sun_term + moon_term


def _load_swpc_json(path: str) -> Tuple[str, list]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list) or len(data) == 0:
        raise DataSourceError("SWPC JSON file must contain a non-empty list")
    return path, data


def parse_swpc_json_series(
    path: str,
    value_key: str,
    provenance: DataProvenance,
    time_key: Optional[str] = None,
) -> TimeSeries:
    _, data = _load_swpc_json(path)

    if isinstance(data[0], list):
        header = [str(h) for h in data[0]]
        rows = data[1:]
        if time_key is None:
            for candidate in ("time_tag", "time-tag", "time"):
                if candidate in header:
                    time_key = candidate
                    break
        if time_key is None or time_key not in header:
            raise DataSourceError("time column not found in SWPC JSON table")
        if value_key not in header:
            raise DataSourceError(f"value column '{value_key}' not found in SWPC JSON table")
        t_idx = header.index(time_key)
        v_idx = header.index(value_key)
        times = _to_datetime64_array([row[t_idx] for row in rows])
        values = np.array([float(row[v_idx]) for row in rows], dtype=float)
        return TimeSeries(times=times, values=values, provenance=provenance, units="")

    if isinstance(data[0], dict):
        if time_key is None:
            for candidate in ("time_tag", "time-tag", "time"):
                if candidate in data[0]:
                    time_key = candidate
                    break
        if time_key is None:
            raise DataSourceError("time key not found in SWPC JSON objects")
        if value_key not in data[0]:
            raise DataSourceError(f"value key '{value_key}' not found in SWPC JSON objects")
        times = _to_datetime64_array([row[time_key] for row in data])
        values = np.array([float(row[value_key]) for row in data], dtype=float)
        return TimeSeries(times=times, values=values, provenance=provenance, units="")

    raise DataSourceError("SWPC JSON format must be a list of lists or list of dicts")


def parse_swpc_json_series_multi(
    path: str,
    value_keys: Dict[str, DataProvenance],
    time_key: Optional[str] = None,
) -> Dict[str, TimeSeries]:
    _, data = _load_swpc_json(path)
    series = {}

    if isinstance(data[0], list):
        header = [str(h) for h in data[0]]
        rows = data[1:]
        if time_key is None:
            for candidate in ("time_tag", "time-tag", "time"):
                if candidate in header:
                    time_key = candidate
                    break
        if time_key is None or time_key not in header:
            raise DataSourceError("time column not found in SWPC JSON table")
        t_idx = header.index(time_key)
        times = _to_datetime64_array([row[t_idx] for row in rows])
        for key, prov in value_keys.items():
            if key not in header:
                raise DataSourceError(f"value column '{key}' not found in SWPC JSON table")
            v_idx = header.index(key)
            values = np.array([float(row[v_idx]) for row in rows], dtype=float)
            series[key] = TimeSeries(times=times, values=values, provenance=prov, units="")
        return series

    if isinstance(data[0], dict):
        if time_key is None:
            for candidate in ("time_tag", "time-tag", "time"):
                if candidate in data[0]:
                    time_key = candidate
                    break
        if time_key is None:
            raise DataSourceError("time key not found in SWPC JSON objects")
        times = _to_datetime64_array([row[time_key] for row in data])
        for key, prov in value_keys.items():
            if key not in data[0]:
                raise DataSourceError(f"value key '{key}' not found in SWPC JSON objects")
            values = np.array([float(row[key]) for row in data], dtype=float)
            series[key] = TimeSeries(times=times, values=values, provenance=prov, units="")
        return series

    raise DataSourceError("SWPC JSON format must be a list of lists or list of dicts")


def parse_swpc_45day_ap_f107(
    path: str,
    provenance_f107: DataProvenance,
    provenance_ap: DataProvenance,
) -> Tuple[TimeSeries, TimeSeries, Optional[TimeSeries]]:
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            tokens = stripped.split()
            if len(tokens) < 4:
                continue
            if tokens[0].isdigit() and len(tokens[0]) == 8:
                year = int(tokens[0][:4])
                month = int(tokens[0][4:6])
                day = int(tokens[0][6:8])
                offset = 1
            elif (
                len(tokens[0]) == 4
                and tokens[0].isdigit()
                and tokens[1].isdigit()
                and tokens[2].isdigit()
            ):
                year = int(tokens[0])
                month = int(tokens[1])
                day = int(tokens[2])
                offset = 3
            else:
                continue
            try:
                nums = [float(t) for t in tokens[offset:]]
            except ValueError:
                continue
            if len(nums) < 2:
                raise DataSourceError(
                    "SWPC 45-day file must include at least Ap and F10.7 columns"
                )
            rows.append((dt.date(year, month, day), nums))

    if not rows:
        raise DataSourceError("SWPC 45-day file contained no parseable data")

    times = _to_datetime64_array([r[0].isoformat() for r in rows])
    ap = np.array([r[1][0] for r in rows], dtype=float)
    f107 = np.array([r[1][1] for r in rows], dtype=float)
    f107a = None
    if any(len(r[1]) >= 3 for r in rows):
        f107a = np.array([r[1][2] if len(r[1]) >= 3 else np.nan for r in rows], dtype=float)
        if np.any(~np.isfinite(f107a)):
            raise DataSourceError("SWPC 45-day file missing f107a for some rows")

    f107_series = TimeSeries(times=times, values=f107, provenance=provenance_f107, units="sfu")
    ap_series = TimeSeries(times=times, values=ap, provenance=provenance_ap, units="nT")
    f107a_series = None
    if f107a is not None:
        f107a_series = TimeSeries(
            times=times, values=f107a, provenance=provenance_f107, units="sfu"
        )

    return f107_series, ap_series, f107a_series


def parse_swpc_f107_csv(
    path: str,
    provenance_f107: DataProvenance,
    provenance_f107a81: Optional[DataProvenance] = None,
) -> Tuple[TimeSeries, Optional[TimeSeries]]:
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise DataSourceError("SWPC CSV missing header")
        if "time" not in reader.fieldnames or "f107" not in reader.fieldnames:
            raise DataSourceError("SWPC CSV requires 'time' and 'f107' columns")
        for row in reader:
            rows.append(row)
    if not rows:
        raise DataSourceError("SWPC CSV contained no data rows")

    times = _to_datetime64_array([r["time"] for r in rows])
    f107 = np.array([float(r["f107"]) for r in rows], dtype=float)
    f107_series = TimeSeries(times=times, values=f107, provenance=provenance_f107, units="sfu")

    f107a_series = None
    if "f107a81" in rows[0]:
        if provenance_f107a81 is None:
            raise DataSourceError("f107a81 column present but no provenance provided")
        f107a = np.array([float(r["f107a81"]) for r in rows], dtype=float)
        f107a_series = TimeSeries(
            times=times, values=f107a, provenance=provenance_f107a81, units="sfu"
        )

    return f107_series, f107a_series


def parse_gfz_kp_ap(
    path: str,
    provenance_kp: DataProvenance,
    provenance_ap: DataProvenance,
) -> Tuple[TimeSeries, TimeSeries]:
    def parse_kp(token: str) -> float:
        token = token.strip()
        if not token:
            raise DataSourceError("empty Kp token")
        if token[-1] in {"+", "-", "o"} and token[0].isdigit():
            base = float(token[0])
            if token[-1] == "+":
                return base + 0.3
            if token[-1] == "-":
                return base - 0.3
            return base
        return float(token)

    times = []
    kp_vals = []
    ap_vals = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if len(tokens) < 3 + 16:
                continue
            if not (tokens[0].isdigit() and tokens[1].isdigit() and tokens[2].isdigit()):
                continue
            year = int(tokens[0])
            month = int(tokens[1])
            day = int(tokens[2])
            kp_tokens = tokens[3:11]
            ap_tokens = tokens[11:19]
            if len(kp_tokens) != 8 or len(ap_tokens) != 8:
                raise DataSourceError("GFZ Kp/ap file must include 8 Kp and 8 ap values per day")

            base = dt.datetime(year, month, day, 0, 0, 0)
            for i in range(8):
                times.append(base + dt.timedelta(hours=3 * i))
                kp_vals.append(parse_kp(kp_tokens[i]))
                ap_vals.append(float(ap_tokens[i]))

    if not times:
        raise DataSourceError("GFZ Kp/ap file contained no parseable data")

    times_arr = _to_datetime64_array([t.isoformat() for t in times])
    kp_series = TimeSeries(times=times_arr, values=np.array(kp_vals), provenance=provenance_kp, units="Kp")
    ap_series = TimeSeries(times=times_arr, values=np.array(ap_vals), provenance=provenance_ap, units="nT")
    return kp_series, ap_series


def parse_kyoto_dst_hourly(path: str, provenance_dst: DataProvenance) -> TimeSeries:
    times = []
    values = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if len(tokens) < 27:
                continue
            if not (tokens[0].isdigit() and tokens[1].isdigit() and tokens[2].isdigit()):
                continue
            year = int(tokens[0])
            month = int(tokens[1])
            day = int(tokens[2])
            vals = tokens[3:27]
            if len(vals) != 24:
                raise DataSourceError("Kyoto Dst file must include 24 hourly values per day")
            base = dt.datetime(year, month, day, 0, 0, 0)
            for i, token in enumerate(vals):
                times.append(base + dt.timedelta(hours=i))
                values.append(float(token))

    if not times:
        raise DataSourceError("Kyoto Dst file contained no parseable data")

    return TimeSeries(
        times=_to_datetime64_array([t.isoformat() for t in times]),
        values=np.array(values),
        provenance=provenance_dst,
        units="nT",
    )


def parse_celestrak_sw_all(
    path: str,
    block: str,
    flux_variant: str,
    f107a_variant: str,
    ap_mode: str,
    kp_mode: str,
    provenance_f107: DataProvenance,
    provenance_f107a81: DataProvenance,
    provenance_ap: DataProvenance,
    provenance_kp: DataProvenance,
) -> Dict[str, TimeSeries]:
    block = block.upper()
    if block not in {"OBSERVED", "DAILY_PREDICTED", "MONTHLY_PREDICTED"}:
        raise DataSourceError("block must be OBSERVED, DAILY_PREDICTED, or MONTHLY_PREDICTED")
    if flux_variant not in {"observed", "adjusted"}:
        raise DataSourceError("flux_variant must be 'observed' or 'adjusted'")
    if f107a_variant not in {"ctr81", "lst81"}:
        raise DataSourceError("f107a_variant must be 'ctr81' or 'lst81'")
    if ap_mode not in {"daily", "3hour"}:
        raise DataSourceError("ap_mode must be 'daily' or '3hour'")
    if kp_mode not in {"daily", "3hour"}:
        raise DataSourceError("kp_mode must be 'daily' or '3hour'")

    data_lines = []
    in_block = False
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("BEGIN "):
                in_block = stripped.split("BEGIN", 1)[1].strip().upper() == block
                continue
            if stripped.startswith("END "):
                in_block = False
                continue
            if not in_block or not stripped or stripped.startswith("#"):
                continue
            data_lines.append(stripped)

    if not data_lines:
        raise DataSourceError(f"no data lines found for block {block}")

    times = []
    f107_vals = []
    f107a_vals = []
    ap_vals = []
    kp_vals = []

    for line in data_lines:
        parts = line.split()
        if len(parts) < 33:
            raise DataSourceError("CelesTrak SW-All line has unexpected column count")
        year = int(parts[0])
        month = int(parts[1])
        day = int(parts[2])
        times.append(dt.date(year, month, day))

        kp_3hour = [float(p) / 10.0 for p in parts[5:13]]
        ap_3hour = [float(p) for p in parts[14:22]]
        ap_daily = float(parts[22])

        adj_f107 = float(parts[26])
        adj_ctr81 = float(parts[28])
        adj_lst81 = float(parts[29])
        obs_f107 = float(parts[30])
        obs_ctr81 = float(parts[31])
        obs_lst81 = float(parts[32])

        if flux_variant == "adjusted":
            f107_vals.append(adj_f107)
        else:
            f107_vals.append(obs_f107)

        if f107a_variant == "ctr81":
            f107a_vals.append(adj_ctr81 if flux_variant == "adjusted" else obs_ctr81)
        else:
            f107a_vals.append(adj_lst81 if flux_variant == "adjusted" else obs_lst81)

        if ap_mode == "daily":
            ap_vals.append(ap_daily)
        else:
            ap_vals.extend(ap_3hour)

        if kp_mode == "daily":
            kp_vals.append(float(parts[13]) / 10.0)
        else:
            kp_vals.extend(kp_3hour)

    if ap_mode == "3hour" or kp_mode == "3hour":
        expanded_times = []
        for t in times:
            base = dt.datetime(t.year, t.month, t.day, 0, 0, 0)
            for i in range(8):
                expanded_times.append(base + dt.timedelta(hours=3 * i))
        if ap_mode == "3hour" and len(ap_vals) != len(expanded_times):
            raise DataSourceError("CelesTrak 3-hour Ap length mismatch")
        if kp_mode == "3hour" and len(kp_vals) != len(expanded_times):
            raise DataSourceError("CelesTrak 3-hour Kp length mismatch")
        ap_times = expanded_times if ap_mode == "3hour" else times
        kp_times = expanded_times if kp_mode == "3hour" else times
    else:
        ap_times = times
        kp_times = times

    return {
        "f107": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times]),
            values=np.array(f107_vals),
            provenance=provenance_f107,
            units="sfu",
        ),
        "f107a81": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times]),
            values=np.array(f107a_vals),
            provenance=provenance_f107a81,
            units="sfu",
        ),
        "ap": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in ap_times]),
            values=np.array(ap_vals),
            provenance=provenance_ap,
            units="nT",
        ),
        "kp": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in kp_times]),
            values=np.array(kp_vals),
            provenance=provenance_kp,
            units="Kp",
        ),
    }


def parse_silso_sn_daily(
    path: str,
    provenance_sn: DataProvenance,
    allow_missing: bool = False,
) -> TimeSeries:
    times = []
    values = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 5:
                continue
            year = int(parts[0])
            month = int(parts[1])
            day = int(parts[2])
            sn = float(parts[4])
            if not allow_missing and sn < 0.0:
                raise DataSourceError("SILSO SN contains missing values; set allow_missing=True to accept")
            times.append(dt.date(year, month, day))
            values.append(sn)

    if not times:
        raise DataSourceError("SILSO SN file contained no parseable data")

    return TimeSeries(
        times=_to_datetime64_array([t.isoformat() for t in times]),
        values=np.array(values),
        provenance=provenance_sn,
        units="SN",
    )


def parse_swpc_3day_forecast_kp(
    path: str,
    provenance_kp: DataProvenance,
) -> TimeSeries:
    issued_year = None
    header_dates = []
    kp_dates = []
    kp_values = []
    in_table = False
    header_seen = False

    month_tokens = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith(":Issued:"):
                tokens = stripped.split()
                # Format: :Issued: 2026 Jan 21 1230 UTC
                for token in tokens:
                    if token.isdigit() and len(token) == 4:
                        issued_year = int(token)
                        break
                continue

            if stripped.startswith("NOAA Kp index breakdown"):
                in_table = True
                continue

            if in_table:
                if not stripped:
                    if header_seen:
                        break
                    continue
                if stripped.split() and stripped.split()[0] in month_tokens:
                    header = stripped.split()
                    if len(header) < 6:
                        raise DataSourceError("Kp forecast header is incomplete")
                    header_dates = [(header[i], header[i + 1]) for i in range(0, len(header), 2)]
                    header_seen = True
                    continue

                if stripped[0].isdigit() and "-" in stripped:
                    columns = re.split(r"\\s{2,}", stripped)
                    if len(columns) < 4:
                        continue
                    if issued_year is None:
                        raise DataSourceError("Issued year not found in 3-day forecast file")
                    time_range = columns[0]
                    values = []
                    for col in columns[1:4]:
                        match = re.search(r"[-+]?\\d+(?:\\.\\d+)?", col)
                        if not match:
                            raise DataSourceError("Kp forecast value is not numeric")
                        values.append(match.group(0))
                    if len(header_dates) < 3:
                        raise DataSourceError("Kp forecast header not found")
                    for (month_name, day_str), value_str in zip(header_dates[:3], values[:3]):
                        month_num = dt.datetime.strptime(month_name, "%b").month
                        day_num = int(day_str)
                        start_hour = int(time_range.split("-")[0][:2])
                        ts = dt.datetime(issued_year, month_num, day_num, start_hour, 0, 0)
                        kp_values.append(float(value_str))
                        kp_dates.append(ts)

    if not kp_values:
        raise DataSourceError("no Kp forecast values found in 3-day forecast file")

    return TimeSeries(
        times=_to_datetime64_array([t.isoformat() for t in kp_dates]),
        values=np.array(kp_values),
        provenance=provenance_kp,
        units="Kp",
    )


def parse_swpc_3day_geomag_forecast(
    path: str,
    provenance_kp: DataProvenance,
    provenance_ap: DataProvenance,
) -> Dict[str, TimeSeries]:
    issued_year = None
    kp_dates = []
    kp_values = []
    ap_dates = []
    ap_values = []
    in_kp_table = False
    header_dates = []
    header_seen = False

    month_tokens = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith(":Issued:"):
                tokens = stripped.split()
                for token in tokens:
                    if token.isdigit() and len(token) == 4:
                        issued_year = int(token)
                        break
                continue

            if stripped.startswith("Predicted Ap"):
                # Format: Predicted Ap 21 Jan-23 Jan 042-018-010
                parts = stripped.split()
                if issued_year is None:
                    raise DataSourceError("Issued year not found in geomag forecast file")
                try:
                    start_day = int(parts[2])
                    start_month = dt.datetime.strptime(parts[3], "%b").month
                    ap_vals = parts[-1].split("-")
                    for i, ap_str in enumerate(ap_vals):
                        ap_dates.append(dt.date(issued_year, start_month, start_day + i))
                        ap_values.append(float(ap_str))
                except (ValueError, IndexError) as exc:
                    raise DataSourceError("failed to parse predicted Ap line") from exc
                continue

            if stripped.startswith("NOAA Kp index forecast"):
                in_kp_table = True
                continue

            if in_kp_table:
                if not stripped:
                    if header_seen:
                        break
                    continue
                if stripped.split() and stripped.split()[0] in month_tokens:
                    header = stripped.split()
                    if len(header) < 6:
                        raise DataSourceError("Kp forecast header is incomplete")
                    header_dates = [(header[i], header[i + 1]) for i in range(0, len(header), 2)]
                    header_seen = True
                    continue
                if stripped[0].isdigit() and "-" in stripped:
                    parts = stripped.split()
                    if len(parts) < 4:
                        continue
                    if issued_year is None:
                        raise DataSourceError("Issued year not found in geomag forecast file")
                    time_range = parts[0]
                    values = parts[1:4]
                    if len(header_dates) < 3:
                        raise DataSourceError("Kp forecast header not found")
                    for (month_name, day_str), value_str in zip(header_dates[:3], values[:3]):
                        month_num = dt.datetime.strptime(month_name, "%b").month
                        day_num = int(day_str)
                        start_hour = int(time_range.split("-")[0][:2])
                        ts = dt.datetime(issued_year, month_num, day_num, start_hour, 0, 0)
                        kp_dates.append(ts)
                        kp_values.append(float(value_str))

    if not kp_values:
        raise DataSourceError("no Kp forecast values found in geomag forecast file")
    if not ap_values:
        raise DataSourceError("no Ap forecast values found in geomag forecast file")

    return {
        "kp": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in kp_dates]),
            values=np.array(kp_values),
            provenance=provenance_kp,
            units="Kp",
        ),
        "ap": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in ap_dates]),
            values=np.array(ap_values),
            provenance=provenance_ap,
            units="nT",
        ),
    }


def parse_omni2_indices(
    path: str,
    provenance_f107: DataProvenance,
    provenance_ap: DataProvenance,
    provenance_kp: DataProvenance,
    provenance_dst: DataProvenance,
    provenance_sunspot: Optional[DataProvenance] = None,
    allow_missing: bool = False,
) -> Dict[str, TimeSeries]:
    times_f107 = []
    values_f107 = []
    times_ap = []
    values_ap = []
    times_kp = []
    values_kp = []
    times_dst = []
    values_dst = []
    times_sn = []
    values_sn = []

    def append_value(times: list, values: list, tstamp: dt.datetime, value: float, valid: bool, name: str):
        if not valid:
            if allow_missing:
                return
            raise DataSourceError(f"OMNI2 contains missing {name} value")
        times.append(tstamp)
        values.append(value)

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if len(tokens) < 51:
                continue
            if not (tokens[0].isdigit() and tokens[1].isdigit() and tokens[2].isdigit()):
                continue
            year = int(tokens[0])
            doy = int(tokens[1])
            hour = int(tokens[2])
            tstamp = dt.datetime(year, 1, 1, 0, 0, 0) + dt.timedelta(days=doy - 1, hours=hour)

            kp_raw = int(tokens[38])
            kp_valid = kp_raw != 99
            append_value(times_kp, values_kp, tstamp, kp_raw / 10.0, kp_valid, "Kp")

            sn_raw = int(tokens[39])
            sn_valid = sn_raw != 999
            if provenance_sunspot is not None:
                append_value(times_sn, values_sn, tstamp, float(sn_raw), sn_valid, "sunspot number")

            dst_raw = int(tokens[40])
            dst_valid = dst_raw != 99999
            append_value(times_dst, values_dst, tstamp, float(dst_raw), dst_valid, "Dst")

            ap_raw = int(tokens[49])
            ap_valid = ap_raw != 999
            append_value(times_ap, values_ap, tstamp, float(ap_raw), ap_valid, "Ap")

            f107_raw = float(tokens[50])
            f107_valid = f107_raw < 999.0 and np.isfinite(f107_raw)
            append_value(times_f107, values_f107, tstamp, f107_raw, f107_valid, "F10.7")

    if not times_f107 or not times_ap or not times_kp or not times_dst:
        raise DataSourceError("OMNI2 file contained no parseable index data")

    series = {
        "f107": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times_f107]),
            values=np.array(values_f107, dtype=float),
            provenance=provenance_f107,
            units="sfu",
        ),
        "ap": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times_ap]),
            values=np.array(values_ap, dtype=float),
            provenance=provenance_ap,
            units="nT",
        ),
        "kp": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times_kp]),
            values=np.array(values_kp, dtype=float),
            provenance=provenance_kp,
            units="Kp",
        ),
        "dst": TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times_dst]),
            values=np.array(values_dst, dtype=float),
            provenance=provenance_dst,
            units="nT",
        ),
    }

    if provenance_sunspot is not None and times_sn:
        series["sunspot"] = TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times_sn]),
            values=np.array(values_sn, dtype=float),
            provenance=provenance_sunspot,
            units="SN",
        )

    return series


def parse_jb2008_indices(
    path: str,
    provenance: Dict[str, DataProvenance],
    allow_missing: bool = False,
) -> Dict[str, TimeSeries]:
    fields = ("f10", "f10b", "s10", "s10b", "m10", "m10b", "y10", "y10b", "dstdtc")
    if not set(fields).issubset(provenance.keys()):
        missing = sorted(set(fields) - set(provenance.keys()))
        raise DataSourceError(f"JB2008 provenance missing fields: {missing}")

    times = []
    values = {key: [] for key in fields}
    expect_indices = False
    current_time: Optional[dt.datetime] = None

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if len(tokens) >= 5 and tokens[0].isdigit() and len(tokens[0]) == 4 and tokens[1].isdigit():
                year = int(tokens[0])
                doy = int(tokens[1])
                hour = int(tokens[2])
                minute = int(tokens[3])
                second = float(tokens[4])
                current_time = dt.datetime(year, 1, 1, 0, 0, 0) + dt.timedelta(
                    days=doy - 1, hours=hour, minutes=minute, seconds=second
                )
                expect_indices = True
                continue

            if expect_indices and current_time is not None:
                parts = stripped.split()
                if len(parts) < 10:
                    raise DataSourceError("JB2008 index line is malformed")
                try:
                    nums = [float(v) for v in parts[:10]]
                except ValueError as exc:
                    raise DataSourceError("JB2008 index values are not numeric") from exc
                index_vals = nums[1:10]
                if any(not np.isfinite(v) for v in index_vals):
                    if allow_missing:
                        expect_indices = False
                        current_time = None
                        continue
                    raise DataSourceError("JB2008 index line contains non-finite values")
                times.append(current_time)
                for key, val in zip(fields, index_vals):
                    values[key].append(val)
                expect_indices = False
                current_time = None

    if not times:
        raise DataSourceError("JB2008 file contained no index records")

    series = {}
    for key in fields:
        unit = "nT" if key == "dstdtc" else "sfu"
        series[key] = TimeSeries(
            times=_to_datetime64_array([t.isoformat() for t in times]),
            values=np.array(values[key], dtype=float),
            provenance=provenance[key],
            units=unit,
        )
    return series


def compute_f107a81(times: Sequence, f107: np.ndarray) -> np.ndarray:
    t = _to_datetime64_array(times).astype("datetime64[D]").astype("int64")
    if len(t) != len(f107):
        raise DataSourceError("f107 length mismatch")
    if len(t) < 81:
        raise DataSourceError("need at least 81 days of f107 to compute f107a81")
    if not np.all(np.diff(t) == 1):
        raise DataSourceError("f107a81 computation requires daily cadence")

    window = np.ones(81, dtype=float) / 81.0
    padded = np.pad(f107, (40, 40), mode="edge")
    return np.convolve(padded, window, mode="valid")


def build_space_weather_dataset(
    t_grid: Sequence,
    f107_series: TimeSeries,
    ap_series: TimeSeries,
    kp_series: TimeSeries,
    dst_series: TimeSeries,
    f107a81_series: Optional[TimeSeries] = None,
    interpolation: Optional[str] = None,
    compute_f107a81_flag: bool = False,
) -> "SpaceWeatherDataset":
    if interpolation is None:
        raise DataSourceError("interpolation method must be provided explicitly")

    f107_series.validate()
    ap_series.validate()
    kp_series.validate()
    dst_series.validate()

    if f107a81_series is None and not compute_f107a81_flag:
        raise DataSourceError("f107a81 series is required or compute_f107a81_flag must be True")

    if f107a81_series is None and compute_f107a81_flag:
        f107a81_vals = compute_f107a81(f107_series.times, f107_series.values)
        f107a81_series = TimeSeries(
            times=f107_series.times,
            values=f107a81_vals,
            provenance=f107_series.provenance,
            units=f107_series.units,
        )

    f107a81_series.validate()

    dataset = SpaceWeatherDataset(
        times=_to_datetime64_array(t_grid),
        f107=_interpolate_series(f107_series, t_grid, interpolation),
        f107a81=_interpolate_series(f107a81_series, t_grid, interpolation),
        ap=_interpolate_series(ap_series, t_grid, interpolation),
        kp=_interpolate_series(kp_series, t_grid, interpolation),
        dst=_interpolate_series(dst_series, t_grid, interpolation),
        provenance={
            "f107": f107_series.provenance,
            "f107a81": f107a81_series.provenance,
            "ap": ap_series.provenance,
            "kp": kp_series.provenance,
            "dst": dst_series.provenance,
        },
    )
    dataset.validate()
    dataset.validate_recommended_sources()
    return dataset

@dataclass
class SpaceWeatherDataset:
    times: np.ndarray
    f107: np.ndarray
    f107a81: np.ndarray
    ap: np.ndarray
    kp: np.ndarray
    dst: np.ndarray
    provenance: Dict[str, DataProvenance]

    def validate(self) -> None:
        n = len(self.times)
        for name, arr in (
            ("f107", self.f107),
            ("f107a81", self.f107a81),
            ("ap", self.ap),
            ("kp", self.kp),
            ("dst", self.dst),
        ):
            if len(arr) != n:
                raise DataSourceError(f"{name} length does not match times")
            if np.any(~np.isfinite(arr)):
                raise DataSourceError(f"{name} contains non-finite values")
            prov = self.provenance.get(name)
            if prov is None:
                raise DataSourceError(f"missing provenance for {name}")
            prov.validate()

        if n == 0:
            raise DataSourceError("space weather dataset is empty")

    def validate_recommended_sources(self) -> None:
        allowed = {
            "f107": ("NOAA_SWPC", "CELESTRAK", "NASA_OMNI"),
            "f107a81": ("NOAA_SWPC", "CELESTRAK", "NASA_OMNI"),
            "ap": ("GFZ", "CELESTRAK", "NOAA_SWPC", "NASA_OMNI"),
            "kp": ("GFZ", "CELESTRAK", "NOAA_SWPC", "NASA_OMNI"),
            "dst": ("KYOTO", "CELESTRAK", "NOAA_NCEI", "NASA_OMNI"),
        }
        for key, options in allowed.items():
            prov = self.provenance.get(key)
            if prov is None:
                continue
            if not any(opt in prov.name for opt in options):
                raise DataSourceError(
                    f"{key} provenance '{prov.name}' is not from recommended sources {options}"
                )
    @classmethod
    def from_csv(cls, path: str, provenance: Dict[str, DataProvenance], time_column: str = "time"):
        required = {"f107", "f107a81", "ap", "kp", "dst"}
        rows = []
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            if time_column not in reader.fieldnames:
                raise DataSourceError(f"missing time column '{time_column}'")
            if not required.issubset(set(reader.fieldnames or [])):
                missing = sorted(required - set(reader.fieldnames or []))
                raise DataSourceError(f"missing required columns: {missing}")
            for row in reader:
                rows.append(row)

        if not rows:
            raise DataSourceError("no rows in space weather CSV")

        times = _to_datetime64_array([r[time_column] for r in rows])
        return cls(
            times=times,
            f107=np.array([float(r["f107"]) for r in rows]),
            f107a81=np.array([float(r["f107a81"]) for r in rows]),
            ap=np.array([float(r["ap"]) for r in rows]),
            kp=np.array([float(r["kp"]) for r in rows]),
            dst=np.array([float(r["dst"]) for r in rows]),
            provenance=provenance,
        )

    def interpolate(self, t_grid: Sequence, method: Optional[str]) -> Dict[str, np.ndarray]:
        if method not in {"linear", "nearest"}:
            raise DataSourceError("interpolation method must be 'linear' or 'nearest'")

        times = self.times.astype("datetime64[s]").astype("int64")
        target = _to_datetime64_array(t_grid).astype("datetime64[s]").astype("int64")

        def interp(values: np.ndarray) -> np.ndarray:
            if method == "linear":
                return np.interp(target, times, values)
            idx = np.searchsorted(times, target, side="left")
            idx = np.clip(idx, 0, len(times) - 1)
            left = np.clip(idx - 1, 0, len(times) - 1)
            choose_left = (target - times[left]) <= (times[idx] - target)
            return np.where(choose_left, values[left], values[idx])

        return {
            "f107": interp(self.f107),
            "f107a81": interp(self.f107a81),
            "ap": interp(self.ap),
            "kp": interp(self.kp),
            "dst": interp(self.dst),
        }


@dataclass
class DensityOutputs:
    density_kg_m3: np.ndarray
    temperature_K: np.ndarray
    particles_mass_kg: np.ndarray


@dataclass
class WindOutputs:
    wind_I: np.ndarray


class DensityModel:
    def evaluate(
        self,
        t_grid: Sequence,
        lat_deg: np.ndarray,
        lon_deg: np.ndarray,
        alt_m: np.ndarray,
        indices: Dict[str, np.ndarray],
    ) -> DensityOutputs:
        raise NotImplementedError


class WindModel:
    def evaluate(
        self,
        t_grid: Sequence,
        lat_deg: np.ndarray,
        lon_deg: np.ndarray,
        alt_m: np.ndarray,
        indices: Dict[str, np.ndarray],
    ) -> WindOutputs:
        raise NotImplementedError


class NRLMSIS21DensityModel(DensityModel):
    def __init__(self, *, version: float = 2.1) -> None:
        try:
            import pymsis.msis as msis
        except ImportError as exc:
            raise DataSourceError("pymsis is required for NRLMSIS 2.1 evaluation") from exc
        self._msis = msis
        self._version = version

    def evaluate(
        self,
        t_grid: Sequence,
        lat_deg: np.ndarray,
        lon_deg: np.ndarray,
        alt_m: np.ndarray,
        indices: Dict[str, np.ndarray],
    ) -> DensityOutputs:
        for key in ("f107", "f107a81", "ap"):
            if key not in indices:
                raise DataSourceError(f"missing required index '{key}' for NRLMSIS")

        dates = _to_datetime64_array(t_grid)
        lats = np.asarray(lat_deg, dtype=float)
        lons = np.asarray(lon_deg, dtype=float)
        alts_km = np.asarray(alt_m, dtype=float) / 1000.0

        if not (len(dates) == len(lats) == len(lons) == len(alts_km)):
            raise DataSourceError("NRLMSIS inputs must have matching lengths")

        ap_vals = np.asarray(indices["ap"], dtype=float)
        if ap_vals.ndim == 1:
            ap_vals = np.tile(ap_vals.reshape(-1, 1), (1, 7))
        if ap_vals.shape[0] != len(dates):
            raise DataSourceError("NRLMSIS ap series length mismatch")

        out = self._msis.calculate(
            dates,
            lons,
            lats,
            alts_km,
            f107s=np.asarray(indices["f107"], dtype=float),
            f107as=np.asarray(indices["f107a81"], dtype=float),
            aps=ap_vals,
            version=self._version,
        )

        if out.shape[-1] < int(self._msis.Variable.TEMPERATURE) + 1:
            raise DataSourceError("NRLMSIS output shape is unexpected")

        rho = out[:, self._msis.Variable.MASS_DENSITY]
        temperature = out[:, self._msis.Variable.TEMPERATURE]
        species = out[
            :,
            [
                self._msis.Variable.N2,
                self._msis.Variable.O2,
                self._msis.Variable.O,
                self._msis.Variable.HE,
                self._msis.Variable.H,
                self._msis.Variable.AR,
                self._msis.Variable.N,
                self._msis.Variable.ANOMALOUS_O,
                self._msis.Variable.NO,
            ],
        ]
        n_total = np.sum(species, axis=1)
        if np.any(n_total <= 0.0):
            raise DataSourceError("NRLMSIS returned non-positive total number density")
        mean_mass = rho / n_total

        return DensityOutputs(
            density_kg_m3=np.asarray(rho, dtype=float),
            temperature_K=np.asarray(temperature, dtype=float),
            particles_mass_kg=np.asarray(mean_mass, dtype=float),
        )


class DTM2013DensityModel(DensityModel):
    def __init__(self) -> None:
        raise DataSourceError(
            "DTM2013 evaluation requires a licensed implementation; provide a project-specific adapter."
        )


class JB2008DensityModel(DensityModel):
    def __init__(self) -> None:
        raise DataSourceError(
            "JB2008 evaluation requires an official implementation; provide a project-specific adapter."
        )


class HWM14SharedLibBackend:
    def __init__(
        self,
        lib_path: Optional[str] = None,
        data_dir: Optional[str] = None,
        symbol: Optional[str] = None,
        output_frame: str = "eci",
    ) -> None:
        if lib_path is None:
            lib_path = os.environ.get("HWM14_LIB") or os.environ.get("VLEO_HWM14_LIB")
        if lib_path is None:
            default_path = os.path.join("data", "hwm14", "libhwm14.so")
            if os.path.isfile(default_path):
                lib_path = default_path
        if lib_path is None or not os.path.isfile(lib_path):
            raise DataSourceError("HWM14 shared library not found; set HWM14_LIB or pass lib_path")

        if data_dir is None:
            candidate = os.path.join("data", "hwm14")
            if os.path.isdir(candidate):
                data_dir = candidate
        if data_dir is not None:
            required = {"dwm07b104i.dat", "gd2qd.dat", "hwm123114.bin"}
            present = set(os.listdir(data_dir)) if os.path.isdir(data_dir) else set()
            missing = sorted(required - present)
            if missing:
                raise DataSourceError(f"HWM14 data directory missing files: {missing}")
            os.environ.setdefault("HWMPATH", os.path.abspath(data_dir))
            os.environ.setdefault("HWM14_DATA", os.path.abspath(data_dir))

        self._output_frame = output_frame
        if output_frame not in {"eci", "ecef"}:
            raise DataSourceError("HWM14 output_frame must be 'eci' or 'ecef'")

        try:
            lib = ctypes.CDLL(lib_path)
        except OSError as exc:
            raise DataSourceError(
                f"failed to load HWM14 shared library '{lib_path}': {exc}"
            ) from exc
        candidates = (
            [symbol]
            if symbol is not None
            else ["__hwm14_module_MOD_hwm14", "hwm14_", "hwm14", "HWM14"]
        )
        func = None
        for name in candidates:
            if name is None:
                continue
            try:
                func = getattr(lib, name)
                break
            except AttributeError:
                continue
        if func is None:
            raise DataSourceError("HWM14 symbol not found in shared library")
        func.argtypes = [
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
        ]
        self._func = func

    def evaluate(
        self,
        t_grid: Sequence,
        lat_deg: np.ndarray,
        lon_deg: np.ndarray,
        alt_m: np.ndarray,
        indices: Dict[str, np.ndarray],
    ) -> WindOutputs:
        if "ap" not in indices:
            raise DataSourceError("HWM14 requires ap index input")

        times = _to_datetime64_array(t_grid)
        lats = np.asarray(lat_deg, dtype=float)
        lons = np.asarray(lon_deg, dtype=float)
        alts_km = np.asarray(alt_m, dtype=float) / 1000.0
        ap_vals = np.asarray(indices["ap"], dtype=float)

        if not (len(times) == len(lats) == len(lons) == len(alts_km) == len(ap_vals)):
            raise DataSourceError("HWM14 inputs must have matching lengths")

        wind_out = np.zeros((len(times), 3), dtype=float)
        for i, ts in enumerate(times):
            timestamp = _datetime64_to_datetime(ts)
            doy = timestamp.timetuple().tm_yday
            iyd = (timestamp.year % 100) * 1000 + doy
            sec = (
                timestamp.hour * 3600.0
                + timestamp.minute * 60.0
                + timestamp.second
                + timestamp.microsecond / 1e6
            )

            iyd_c = ctypes.c_int(iyd)
            sec_c = ctypes.c_float(sec)
            alt_c = ctypes.c_float(float(alts_km[i]))
            lat_c = ctypes.c_float(float(lats[i]))
            lon_c = ctypes.c_float(float(lons[i]))
            stl_c = ctypes.c_float(0.0)
            f107a_c = ctypes.c_float(0.0)
            f107_c = ctypes.c_float(0.0)
            ap = (ctypes.c_float * 2)(0.0, float(ap_vals[i]))
            w = (ctypes.c_float * 2)()

            self._func(
                ctypes.byref(iyd_c),
                ctypes.byref(sec_c),
                ctypes.byref(alt_c),
                ctypes.byref(lat_c),
                ctypes.byref(lon_c),
                ctypes.byref(stl_c),
                ctypes.byref(f107a_c),
                ctypes.byref(f107_c),
                ap,
                w,
            )

            north = float(w[0])
            east = float(w[1])
            lat_rad = math.radians(lats[i])
            lon_rad = math.radians(lons[i])
            wind_ecef = _enu_to_ecef(north, east, lat_rad, lon_rad)
            if self._output_frame == "eci":
                gmst = _gmst_rad(timestamp)
                wind_out[i] = _ecef_to_eci(wind_ecef, gmst)
            else:
                wind_out[i] = wind_ecef

        return WindOutputs(wind_I=wind_out)


class HWM14WindModel(WindModel):
    def __init__(self, backend) -> None:
        if backend is None or not hasattr(backend, "evaluate"):
            raise DataSourceError("HWM14 backend with evaluate() is required")
        self._backend = backend

    def evaluate(
        self,
        t_grid: Sequence,
        lat_deg: np.ndarray,
        lon_deg: np.ndarray,
        alt_m: np.ndarray,
        indices: Dict[str, np.ndarray],
    ) -> WindOutputs:
        out = self._backend.evaluate(
            t_grid=t_grid,
            lat_deg=np.asarray(lat_deg, dtype=float),
            lon_deg=np.asarray(lon_deg, dtype=float),
            alt_m=np.asarray(alt_m, dtype=float),
            indices=indices,
        )
        if not isinstance(out, WindOutputs):
            raise DataSourceError("HWM14 backend must return WindOutputs")
        return out


class SrpEclipseModel:
    def evaluate(self, t_grid: Sequence, r_eci: np.ndarray) -> Dict[str, np.ndarray]:
        raise NotImplementedError


@dataclass
class EnvSeries:
    env_inputs: Sequence
    space_weather: Dict[str, np.ndarray]
    srp_scale: Optional[np.ndarray]
    in_eclipse: Optional[np.ndarray]


class EnvSeriesBuilder:
    def __init__(
        self,
        density_model: DensityModel,
        wind_model: WindModel,
        space_weather: SpaceWeatherDataset,
        srp_eclipse_model: Optional[SrpEclipseModel] = None,
    ) -> None:
        if density_model is None or wind_model is None:
            raise DataSourceError("density_model and wind_model are required")
        if space_weather is None:
            raise DataSourceError("space_weather dataset is required")
        space_weather.validate()
        space_weather.validate_recommended_sources()
        self._density_model = density_model
        self._wind_model = wind_model
        self._space_weather = space_weather
        self._srp_eclipse_model = srp_eclipse_model

    def build(
        self,
        t_grid: Sequence,
        lat_deg: np.ndarray,
        lon_deg: np.ndarray,
        alt_m: np.ndarray,
        r_eci: Optional[np.ndarray] = None,
        interpolation: Optional[str] = None,
    ) -> EnvSeries:
        from vleo_uq import EnvInputs

        if interpolation is None:
            raise DataSourceError("interpolation method must be provided explicitly")
        indices = self._space_weather.interpolate(t_grid, interpolation)
        dens = self._density_model.evaluate(t_grid, lat_deg, lon_deg, alt_m, indices)
        wind = self._wind_model.evaluate(t_grid, lat_deg, lon_deg, alt_m, indices)

        times = _to_datetime64_array(t_grid)
        n = len(times)
        if len(dens.density_kg_m3) != n or len(wind.wind_I) != n:
            raise DataSourceError("model outputs must match t_grid length")

        srp_scale = None
        in_eclipse = None
        srp_scale_arr = np.ones(n, dtype=float)
        if self._srp_eclipse_model is not None:
            if r_eci is None:
                raise DataSourceError("r_eci is required for SRP/eclipse evaluation")
            srp = self._srp_eclipse_model.evaluate(t_grid, r_eci)
            srp_scale = srp.get("srp_scale")
            in_eclipse = srp.get("in_eclipse")
            if srp_scale is None or in_eclipse is None:
                raise DataSourceError("srp_eclipse_model must return srp_scale and in_eclipse")
            srp_scale_arr = np.asarray(srp_scale, dtype=float).reshape(-1)
            if srp_scale_arr.size != n:
                raise DataSourceError("srp_scale length must match t_grid")
            srp_scale_arr = np.clip(srp_scale_arr, 0.0, 1.5)

        r_eci_arr = None
        if r_eci is not None:
            r_eci_arr = np.asarray(r_eci, dtype=float)
            if r_eci_arr.shape != (n, 3):
                raise DataSourceError("r_eci must have shape (len(t_grid), 3)")

        env = []
        for i in range(n):
            ts = _datetime64_to_datetime(times[i])
            sun_i = _sun_position_eci_m(ts)
            moon_i = _moon_position_eci_m(ts)
            e = EnvInputs()
            e.density = float(dens.density_kg_m3[i])
            e.temperature_K = float(dens.temperature_K[i])
            e.particles_mass_kg = float(dens.particles_mass_kg[i])
            e.wind_I = wind.wind_I[i]
            e.sun_position_I_m = sun_i
            e.moon_position_I_m = moon_i
            e.srp_scale = float(srp_scale_arr[i])
            e.albedo_ir_scale = 1.0
            if r_eci_arr is not None:
                e.magnetic_field_I_T = _earth_dipole_field_eci_t(r_eci_arr[i])
                e.tide_loading_accel_I_m_s2 = _tide_loading_accel_eci_m_s2(r_eci_arr[i], sun_i, moon_i)
            env.append(e)

        return EnvSeries(env_inputs=env, space_weather=indices, srp_scale=srp_scale, in_eclipse=in_eclipse)
