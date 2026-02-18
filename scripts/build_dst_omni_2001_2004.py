from __future__ import annotations

import csv
import datetime as dt
from pathlib import Path


def parse_omni_dst(path: Path):
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 41:
                continue
            if not (parts[0].isdigit() and parts[1].isdigit() and parts[2].isdigit()):
                continue
            year = int(parts[0])
            doy = int(parts[1])
            hour = int(parts[2])
            dst_raw = int(parts[40])
            if dst_raw == 99999:
                continue
            timestamp = dt.datetime(year, 1, 1, 0, 0, 0) + dt.timedelta(
                days=doy - 1, hours=hour
            )
            yield timestamp, float(dst_raw)


def main() -> None:
    omni_path = Path("data/space_weather/omni/omni2_all_years.dat")
    out_path = Path("data/space_weather/dst_omni/dst_2001_2004_hourly.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    start = dt.datetime(2001, 1, 1)
    end = dt.datetime(2004, 12, 31, 23, 0, 0)

    rows = []
    for ts, dst in parse_omni_dst(omni_path):
        if start <= ts <= end:
            rows.append((ts.isoformat(), dst))

    if not rows:
        raise RuntimeError("no OMNI Dst samples found for 2001-2004")

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "dst"])
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()
