Offline Usage Notes

This project is designed to run without internet access once data and model
artifacts are staged into the repository.

Required staged assets (current layout)
- HWM14 shared library and data: `data/hwm14/libhwm14.so`, `data/hwm14/dwm07b104i.dat`,
  `data/hwm14/gd2qd.dat`, `data/hwm14/hwm123114.bin`.
- GNSS ephemerides for POD: `data/gnss/*.SP3` (decompressed).
- Earth orientation: `data/eop/EOP-All.txt`.
- Space weather indices:
  - OMNI2: `data/space_weather/omni/omni2_all_years.dat`
  - GFZ Kp/ap: `data/space_weather/gfz/Kp_ap_Ap_SN_F107_since_1932.txt`
  - Dst: `data/space_weather/dst/dstYYMM.for.request` (as needed)
  - SWPC forecast/flux: `data/space_weather/swpc_forecast/45-day-ap-forecast.txt`
  - CelesTrak SW-All: `data/space_weather/SW-All.txt`
  - JB2008 proxies: `data/space_weather/jb2008/JB2008_AUTO_OUTPUT.DAT`

Build HWM14 offline
1) Clone `jacobwilliams/HWM14` once on a machine with internet.
2) Use `scripts/build_hwm14.sh /path/to/HWM14` to compile and stage files into `data/hwm14/`.

Offline Python installation
1) Run `scripts/build_wheelhouse.sh` on a machine with internet.
2) Copy the resulting `wheelhouse/` directory to the cluster.
3) Run `scripts/install_offline.sh /path/to/wheelhouse` on the cluster.

No network access should be required at runtime if the above files are present.

Remaining manual/licensed items
- DTM2013 and JB2008 model executables/bindings (licensed).
- HASDM density database (restricted).
- Optional: EGM2008 gravity and NAIF planetary SPKs (e.g., de440) for SRP/eclipse geometry.
 
Notes
- Offline space-weather + wind coupling requires `pymsis` (NRLMSIS 2.1) and the HWM14 library.

Multi-node MC (MPI/Slurm)
- Use `scripts/run_mc_job_array.py` with MPI ranks or Slurm tasks.
- Helper launcher: `scripts/run_mc_mpi.sh`.
- Scenario builders live in `scripts/mc_scenarios.py` (mission/formation/attitude).

Multi-node UT (MPI/Slurm)
- Use `scripts/run_ut_job_array.py` with MPI ranks or Slurm tasks.
- Helper launcher: `scripts/run_ut_mpi.sh`.
- Scenario builders live in `scripts/ut_scenarios.py` (mission/formation/attitude).
