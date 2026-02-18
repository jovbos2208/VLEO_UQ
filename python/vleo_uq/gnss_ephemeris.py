"""SP3 ephemeris parsing and interpolation helpers (offline GNSS support)."""

from __future__ import annotations

from dataclasses import dataclass
import datetime as dt
from typing import Iterable, Optional, Sequence

import numpy as np

from .env_sources import _gmst_rad


def _parse_epoch(line: str) -> dt.datetime:
    parts = line[1:].strip().split()
    if len(parts) < 6:
        raise ValueError("invalid SP3 epoch line")
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])
    hour = int(parts[3])
    minute = int(parts[4])
    sec = float(parts[5])
    sec_int = int(sec)
    micro = int(round((sec - sec_int) * 1e6))
    return dt.datetime(year, month, day, hour, minute, sec_int, micro)


def _parse_position(line: str) -> tuple[str, np.ndarray, float]:
    sat_id = line[1:4].strip()
    parts = line[4:].strip().split()
    if len(parts) < 3:
        raise ValueError("invalid SP3 position line")
    pos_km = np.array([float(parts[0]), float(parts[1]), float(parts[2])], dtype=float)
    clock_s = np.nan
    if len(parts) >= 4:
        try:
            clock_s = float(parts[3]) * 1e-6
        except ValueError:
            clock_s = np.nan
    return sat_id, pos_km * 1000.0, clock_s


def _rotation_z(theta_rad: float) -> np.ndarray:
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rotation_x(theta_rad: float) -> np.ndarray:
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _rotation_y(theta_rad: float) -> np.ndarray:
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _julian_date(ts: dt.datetime) -> float:
    year = ts.year
    month = ts.month
    day = ts.day
    hour = ts.hour + ts.minute / 60.0 + ts.second / 3600.0 + ts.microsecond / 3.6e9
    if month <= 2:
        year -= 1
        month += 12
    A = int(year / 100)
    B = 2 - A + int(A / 4)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    return jd + hour / 24.0


def _gmst_rad_ut1(ts_utc: dt.datetime, ut1_utc_s: float) -> float:
    ts_ut1 = ts_utc + dt.timedelta(seconds=float(ut1_utc_s))
    jd_ut1 = _julian_date(ts_ut1)
    T = (jd_ut1 - 2451545.0) / 36525.0
    gmst_sec = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * T
        + 0.093104 * T * T
        - 6.2e-6 * T * T * T
    )
    gmst_sec = gmst_sec % 86400.0
    return gmst_sec * (2.0 * np.pi / 86400.0)


@dataclass
class EopSeries:
    epochs: Sequence[dt.datetime]
    x_rad: np.ndarray
    y_rad: np.ndarray
    ut1_utc_s: np.ndarray

    def _times_sec(self, t0: Optional[dt.datetime] = None) -> np.ndarray:
        if t0 is None:
            t0 = self.epochs[0]
        return np.array([(t - t0).total_seconds() for t in self.epochs], dtype=float)

    def interpolate(self, ts: dt.datetime) -> tuple[float, float, float]:
        t0 = self.epochs[0]
        t_ref = self._times_sec(t0)
        t_query = (ts - t0).total_seconds()
        x = np.interp(t_query, t_ref, self.x_rad)
        y = np.interp(t_query, t_ref, self.y_rad)
        ut1 = np.interp(t_query, t_ref, self.ut1_utc_s)
        return float(x), float(y), float(ut1)


def parse_eop_all(path: str) -> EopSeries:
    epochs = []
    x_rad = []
    y_rad = []
    ut1_utc_s = []
    arcsec_to_rad = np.pi / (180.0 * 3600.0)
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#") or line.startswith("BEGIN") or line.startswith("END"):
                continue
            parts = line.strip().split()
            if len(parts) < 11:
                continue
            try:
                year = int(parts[0])
                month = int(parts[1])
                day = int(parts[2])
                x_arcsec = float(parts[4])
                y_arcsec = float(parts[5])
                ut1 = float(parts[6])
            except ValueError:
                continue
            epochs.append(dt.datetime(year, month, day))
            x_rad.append(x_arcsec * arcsec_to_rad)
            y_rad.append(y_arcsec * arcsec_to_rad)
            ut1_utc_s.append(ut1)
    if not epochs:
        raise ValueError("no EOP records parsed")
    return EopSeries(
        epochs=epochs,
        x_rad=np.array(x_rad, dtype=float),
        y_rad=np.array(y_rad, dtype=float),
        ut1_utc_s=np.array(ut1_utc_s, dtype=float),
    )


@dataclass
class AntennaPhaseCenter:
    sat_id: str
    pco_by_freq_m: dict[str, np.ndarray]


def parse_antex(path: str) -> dict[str, AntennaPhaseCenter]:
    entries: dict[str, AntennaPhaseCenter] = {}
    current_id = None
    current_pco: dict[str, np.ndarray] = {}
    current_freq = None

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if "START OF ANTENNA" in line:
                current_id = None
                current_pco = {}
                current_freq = None
                continue
            if "END OF ANTENNA" in line:
                if current_id:
                    entries[current_id] = AntennaPhaseCenter(current_id, current_pco)
                current_id = None
                current_pco = {}
                current_freq = None
                continue
            if "TYPE / SERIAL NO" in line:
                parts = line[:60].split()
                sat = None
                for token in parts:
                    if token and token[0] in {"G", "E", "R", "C", "J", "I", "S"} and token[1:].isdigit():
                        sat = token
                        break
                current_id = sat
                continue
            if "START OF FREQUENCY" in line:
                current_freq = line[:6].strip()
                continue
            if "END OF FREQUENCY" in line:
                current_freq = None
                continue
            if "NORTH / EAST / UP" in line and current_id and current_freq:
                parts = line[:60].split()
                if len(parts) >= 3:
                    pco_mm = np.array([float(parts[0]), float(parts[1]), float(parts[2])], dtype=float)
                    current_pco[current_freq] = pco_mm * 1e-3

    return entries

@dataclass
class Sp3Ephemeris:
    epochs: Sequence[dt.datetime]
    sat_ids: Sequence[str]
    positions_ecef_m: np.ndarray
    clocks_s: Optional[np.ndarray] = None

    def _times_sec(self, t0: Optional[dt.datetime] = None) -> np.ndarray:
        if t0 is None:
            t0 = self.epochs[0]
        return np.array([(t - t0).total_seconds() for t in self.epochs], dtype=float)

    def sample(
        self,
        t_grid: Sequence,
        t0: Optional[dt.datetime] = None,
        frame: str = "eci",
        eop: Optional[EopSeries] = None,
    ) -> np.ndarray:
        if frame not in {"ecef", "eci"}:
            raise ValueError("frame must be 'ecef' or 'eci'")

        if isinstance(t_grid, np.ndarray) and np.issubdtype(t_grid.dtype, np.datetime64):
            t_grid_dt = [dt.datetime.utcfromtimestamp(t.astype("datetime64[ns]").astype(int) / 1e9)
                         for t in t_grid]
        elif isinstance(t_grid, (list, tuple)) and t_grid and isinstance(t_grid[0], dt.datetime):
            t_grid_dt = list(t_grid)
        else:
            if t0 is None:
                t0 = self.epochs[0]
            t_grid_dt = None
            t_query = np.array(t_grid, dtype=float)

        if t_grid_dt is not None:
            if t0 is None:
                t0 = self.epochs[0]
            t_query = np.array([(t - t0).total_seconds() for t in t_grid_dt], dtype=float)

        t_ref = self._times_sec(t0)
        nt = len(t_query)
        ns = len(self.sat_ids)
        out = np.full((nt, ns, 3), np.nan, dtype=float)

        for s in range(ns):
            for k in range(3):
                series = self.positions_ecef_m[:, s, k]
                valid = np.isfinite(series)
                if np.count_nonzero(valid) < 2:
                    continue
                out[:, s, k] = np.interp(t_query, t_ref[valid], series[valid])

        if frame == "ecef":
            return out

        if t_grid_dt is None:
            t_grid_dt = [t0 + dt.timedelta(seconds=float(sec)) for sec in t_query]

        for i, ts in enumerate(t_grid_dt):
            if eop is None:
                gmst = _gmst_rad(ts)
                rot = _rotation_z(gmst)
                out[i] = (rot @ out[i].T).T
            else:
                x_pole, y_pole, ut1_utc = eop.interpolate(ts)
                gmst = _gmst_rad_ut1(ts, ut1_utc)
                rot = _rotation_z(gmst) @ _rotation_y(-x_pole) @ _rotation_x(-y_pole)
                out[i] = (rot @ out[i].T).T
        return out

    def sample_clock(
        self,
        t_grid: Sequence,
        t0: Optional[dt.datetime] = None,
    ) -> Optional[np.ndarray]:
        if self.clocks_s is None:
            return None
        if isinstance(t_grid, np.ndarray) and np.issubdtype(t_grid.dtype, np.datetime64):
            t_grid_dt = [dt.datetime.utcfromtimestamp(t.astype("datetime64[ns]").astype(int) / 1e9)
                         for t in t_grid]
        elif isinstance(t_grid, (list, tuple)) and t_grid and isinstance(t_grid[0], dt.datetime):
            t_grid_dt = list(t_grid)
        else:
            if t0 is None:
                t0 = self.epochs[0]
            t_grid_dt = None
            t_query = np.array(t_grid, dtype=float)

        if t_grid_dt is not None:
            if t0 is None:
                t0 = self.epochs[0]
            t_query = np.array([(t - t0).total_seconds() for t in t_grid_dt], dtype=float)

        t_ref = self._times_sec(t0)
        nt = len(t_query)
        ns = len(self.sat_ids)
        out = np.full((nt, ns), np.nan, dtype=float)
        for s in range(ns):
            series = self.clocks_s[:, s]
            valid = np.isfinite(series)
            if np.count_nonzero(valid) < 2:
                continue
            out[:, s] = np.interp(t_query, t_ref[valid], series[valid])
        return out


def parse_sp3(
    path: str,
    constellation_prefix: Optional[Iterable[str]] = None,
) -> Sp3Ephemeris:
    epochs: list[dt.datetime] = []
    epoch_maps: list[dict[str, np.ndarray]] = []
    clock_maps: list[dict[str, float]] = []
    sat_ids: list[str] = []
    sat_set = set()

    prefixes = None
    if constellation_prefix is not None:
        prefixes = tuple(constellation_prefix)

    current: dict[str, np.ndarray] = {}
    current_clock: dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith("*"):
                if epochs:
                    epoch_maps.append(current)
                    clock_maps.append(current_clock)
                    current = {}
                    current_clock = {}
                epochs.append(_parse_epoch(line))
            elif line.startswith("P"):
                sat_id, pos, clock_s = _parse_position(line)
                if prefixes is not None and not sat_id.startswith(prefixes):
                    continue
                current[sat_id] = pos
                current_clock[sat_id] = clock_s
                if sat_id not in sat_set:
                    sat_set.add(sat_id)
                    sat_ids.append(sat_id)

    if epochs:
        epoch_maps.append(current)
        clock_maps.append(current_clock)

    if not epochs:
        raise ValueError("no epochs parsed from SP3")
    if not sat_ids:
        raise ValueError("no satellite positions parsed from SP3")

    positions = np.full((len(epochs), len(sat_ids), 3), np.nan, dtype=float)
    clocks = np.full((len(epochs), len(sat_ids)), np.nan, dtype=float)
    id_to_idx = {sid: i for i, sid in enumerate(sat_ids)}
    for ei, mapping in enumerate(epoch_maps):
        for sid, pos in mapping.items():
            positions[ei, id_to_idx[sid]] = pos
        for sid, clk in clock_maps[ei].items():
            clocks[ei, id_to_idx[sid]] = clk

    return Sp3Ephemeris(
        epochs=epochs,
        sat_ids=sat_ids,
        positions_ecef_m=positions,
        clocks_s=clocks,
    )
