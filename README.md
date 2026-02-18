# VLEO_UQ

Uncertainty quantification framework for very low Earth orbit (VLEO) dynamics, attitude coupling, and POD-oriented analysis.

## Scope

This repository combines:
- C++ dynamics/force-moment core (`cpp/`, `force_moment_module_v4/`)
- Python orchestration and analysis (`python/`, `scripts/`)
- scenario definitions (`configs/`, `scenario_catalog.md`)

Key workflows include deterministic + stochastic propagation, MC/UT comparisons, POD UQ post-processing, and validation gates.

## Repository Layout

- `cpp/`: native core and bindings support code
- `python/vleo_uq/`: Python package
- `scripts/`: run/merge/postprocess/validation utilities
- `configs/`: case-study configuration files
- `tests/`: unit/regression tests
- `results/`: generated run outputs

`data/`, `results/`, `wheelhouse/`, and `slurm_scripts/` are local runtime artifacts and are ignored by git in this repository setup.

## Prerequisites

- Python `>=3.10`
- CMake `>=3.21`
- A C++20-capable compiler
- (cluster runs) SLURM environment and required modules

## Local Setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -e .
```

## Run Locally

Single entrypoint:

```bash
bash run_scenarios.sh 1000 4 configs/case_studies.json configs/case_studies_1000.json
```

Direct call:

```bash
python scripts/run_case_studies.py \
  --config configs/case_studies_1000.json \
  --threads 4 \
  --plot \
  --pod_uq
```

Outputs are written under `results/`.

## Required Data (not versioned)

Create this local layout:

```text
data/
  hwm14/
    libhwm14.so
    dwm07b104i.dat
    gd2qd.dat
    hwm123114.bin
  space_weather/omni/
    omni2_all_years.dat
  eop/
    EOP-All.txt
  gnss/
    ... SP3 files ...
```

1. OMNI2 space weather:

```bash
mkdir -p data/space_weather/omni
wget -O data/space_weather/omni/omni2_all_years.dat \
  https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat
```

2. HWM14 shared library + data:

```bash
mkdir -p data/hwm14
bash scripts/build_hwm14.sh /path/to/HWM14 data/hwm14
```

3. EOP file for SP3 interpolation:

```bash
mkdir -p data/eop
wget -O data/eop/EOP-All.txt https://celestrak.org/SpaceData/EOP-All.txt
```

4. GNSS SP3 orbit files:

```bash
bash scripts/download_sp3_weeks.sh
find data/gnss -type f \( -name '*.gz' -o -name '*.Z' \) -exec gunzip -f {} \;
```

At minimum, make sure one POD SP3 file exists at the path you pass to `--pod_sp3`
(for example `data/gnss/quiet/IGS0OPSFIN_20250240000_01D_15M_ORB.SP3`).

To run with these defaults, set:

```bash
export VLEO_OMNI_PATH=data/space_weather/omni/omni2_all_years.dat
export VLEO_HWM14_LIB=data/hwm14/libhwm14.so
export VLEO_HWM14_DATA=data/hwm14
```

## Validation Gates

Current run pipeline includes:
- POD harness validation (`scripts/validate_pod_uq_harness.py`)
- measurement realism gate (`scripts/validate_measurement_realism_gate.py`)

## Documentation

- `Plan.md`
- `Plan_progression.md`
- `docs/combined_implementation_check_plan.md`
