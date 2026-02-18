# Decision Log

This file freezes baseline technical decisions required by `docs/combined_implementation_check_plan.md`.

## 2026-02-18: Sprint 0 baseline freeze

### DG-1: 18-state definition

Status: accepted

Decision:

- Use the current propagated state layout from C++ bindings and propagator implementation:
  - `0:3`   `r_I` [m], ECI position
  - `3:6`   `v_I` [m/s], ECI velocity
  - `6:10`  `q_BI` [unitless], quaternion `[w,x,y,z]` mapping inertial -> body
  - `10:13` `w_BI_B` [rad/s], body angular rate
  - `13`    `log_rho_fast` [log scale], fast OU density term
  - `14`    `log_rho_bias` [log scale], slow OU/bias density term
  - `15:18` `wind_bias_I` [m/s], inertial wind bias OU term

Evidence:

- `cpp/include/vleo/propagator.hpp` (`kBaseStateSize=13`, `kLatentSize=5`, `kStateSize=18`).
- `cpp/bindings/bindings.cpp` state layout string in propagator docstring.
- `cpp/src/propagator.cpp` latent index definitions and state unpacking.

### DG-2: Baseline GSI model

Status: accepted

Decision:

- Baseline aerodynamic/GSI model for runs is the currently implemented accommodation-coefficient formulation with `temperature_ratio_method=1` default.
- `GSI_C1_CLL` remains catalog-level only until a CLL model is explicitly implemented in `force_moment_module_v4` and wired into runtime toggles.

Evidence:

- `force_moment_module_v4/src/aerodynamics.cpp` contains accommodation-driven aerodynamic model with selectable `temperature_ratio_method` and no CLL branch.
- `force_moment_module_v4/include/env_config.h` and `cpp/include/vleo/propagator.hpp` default `temperature_ratio_method=1`.
- `scenario_catalog.md` marks CLL as "if implemented".

### DG-3: Primary OD baseline

Status: accepted

Decision:

- Primary POD OD estimator baseline is `batch` (batch least-squares per arc).
- EnKF remains an optional comparison mode when explicitly selected.

Evidence:

- `scripts/catalog_to_case_config.py` maps estimator to `enkf` only when `EST_ENKF` is present; otherwise defaults to `batch`.
- `scripts/run_case_studies.py`, `scripts/postprocess_scenario.py`, and `scripts/postprocess_uq_pod.py` default `pod_estimator` to `batch`.
- `python/vleo_uq/pod_uq.py` implements batch OD path (`batch_od_measurements`, `run_pod_uq_measurements`).

### DG-4: Validation truth definition

Status: accepted

Decision:

- Baseline truth for POD/UQ validation is internally generated dynamics truth from deterministic propagation of `x0_truth`, with synthetic measurement generation on top.
- Validation metrics are state-space and residual-based against that truth baseline.
- Accelerometer-assisted truth and DSMC-subset truth are deferred extension tracks.

Evidence:

- `python/vleo_uq/pod_uq.py` uses `truth_states = prop_det.propagate(x0_truth, t_grid, env)` and evaluates RTN errors/coverage against that trajectory.
- `Plan_progression.md` lists accelerometer-assisted/DSMC truth as open options.

## Change control

- Any change to these four decisions must add a new dated entry in this file and update affected configs/docs.
