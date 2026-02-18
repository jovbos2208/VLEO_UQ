from __future__ import annotations

import numpy as np

from vleo_uq import (
    DataProvenance,
    HWM14SharedLibBackend,
    HWM14WindModel,
    NRLMSIS21DensityModel,
    TimeSeries,
    build_space_weather_dataset,
    parse_omni2_indices,
)


def daily_mean(times: np.ndarray, values: np.ndarray) -> TimeSeries:
    days = times.astype("datetime64[D]")
    unique_days, inv = np.unique(days, return_inverse=True)
    sums = np.zeros(len(unique_days), dtype=float)
    counts = np.zeros(len(unique_days), dtype=float)
    for idx, val in zip(inv, values):
        sums[idx] += val
        counts[idx] += 1.0
    means = sums / np.maximum(counts, 1.0)
    return TimeSeries(
        times=unique_days.astype("datetime64[ns]"),
        values=means,
        provenance=DataProvenance(
            name="NASA_OMNI",
            url="https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat",
            retrieved_utc="offline",
        ),
        units="sfu",
    )


def main() -> None:
    omni_path = "data/space_weather/omni/omni2_all_years.dat"
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

    # F10.7 and F10.7a81 should be daily; OMNI provides hourly, so we average per day.
    f107_daily = daily_mean(omni["f107"].times, omni["f107"].values)

    # Build a time grid for one day (hourly) as an example.
    t_grid = omni["f107"].times[:24]

    space_weather = build_space_weather_dataset(
        t_grid=t_grid,
        f107_series=f107_daily,
        f107a81_series=None,
        ap_series=omni["ap"],
        kp_series=omni["kp"],
        dst_series=omni["dst"],
        interpolation="nearest",
        compute_f107a81_flag=True,
    )

    backend = HWM14SharedLibBackend(
        lib_path="data/hwm14/libhwm14.so",
        data_dir="data/hwm14",
        output_frame="eci",
    )
    wind_model = HWM14WindModel(backend)
    density_model = NRLMSIS21DensityModel()

    # Placeholder geometry inputs for demonstration only.
    lat_deg = np.zeros(len(t_grid))
    lon_deg = np.zeros(len(t_grid))
    alt_m = np.full(len(t_grid), 400e3)

    indices = space_weather.interpolate(t_grid, "nearest")
    wind = wind_model.evaluate(t_grid, lat_deg, lon_deg, alt_m, indices)
    dens = density_model.evaluate(t_grid, lat_deg, lon_deg, alt_m, indices)

    print("wind_I shape:", wind.wind_I.shape)
    print("density shape:", dens.density_kg_m3.shape)


if __name__ == "__main__":
    main()
