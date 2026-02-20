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

## 2026-02-19: Mathematical-description scope freeze (MD-01)

Status: accepted

Decision:

- For all `Partial`/`Missing` items in `docs/mathematical_description_implementation_checklist.md`, each item is now explicitly classified as:
  - `Implement in framework`, or
  - `Defer with explicit thesis-scope rationale`.

Immediate implementation decisions (completed):

- `MD-02`: composition uncertainty channel implemented as a stochastic OU/RW surrogate on `particles_mass_kg` in `scripts/model_discrepancy.py`.
- `MD-03`: storm regime-switch/jump disturbance channel implemented as a Poisson-start piecewise jump process in `scripts/model_discrepancy.py`.

Deferred-to-next-phase decisions:

- Higher-fidelity force-model additions (`J2+`, third-body, SRP acceleration coupling, magnetic torque) remain planned under `MD-04..MD-07`.
- Rarefied-regime uncertainty formalization and continuous-thrust noise modeling remain planned under `MD-08..MD-09`.
- POD latency model and full magnetometer integration remain planned under `MD-10..MD-11`.
- AO/plasma/structural/fuel/software-human channels remain planned-or-deferred under `MD-12..MD-14` depending on thesis scope.

Evidence:

- `docs/mathematical_description_execution_roadmap.md`
- `docs/mathematical_description_implementation_checklist.md`
- `scripts/model_discrepancy.py`

## 2026-02-20: Core force-model closure pass (MD-04..MD-07)

Status: accepted

Decision:

- Implement core force-model completeness in propagated dynamics for thesis-phase runs:
  - Zonal gravity expanded to optional `J2/J3/J4`.
  - Third-body acceleration supports Sun/Moon with ephemeris-scale uncertainty controls.
  - SRP acceleration added with `srp_scale` coupling from environment time series.
  - Magnetic residual dipole torque added with simple dipole-field surrogate and scale channels.

Scope notes:

- SRP is now propagated; albedo/IR remains out of scope for this pass.
- Third-body ephemeris uncertainty is represented through scale channels (not full covariance ephemeris model).

Evidence:

- `cpp/include/vleo/propagator.hpp`
- `cpp/src/propagator.cpp`
- `cpp/bindings/bindings.cpp`
- `python/vleo_uq/env_sources.py`
- `scripts/uq_parameter_channels.py`
- `tests/test_propagator_regression.py`

## 2026-02-20: Estimation and surrogate-channel closure pass (MD-10..MD-14)

Status: accepted

Decision:

- Complete the remaining main implementation gaps with explicit surrogate channels and tests:
  - POD latency model implemented as delayed measurement-time remapping with configurable jitter and seed.
  - Operational outage process implemented for GNSS/SLR measurement availability to represent software/human disruptions.
  - Magnetometer measurement model integrated across full-state MEKF/UKF/ENKF.
  - AO erosion surrogate implemented as stochastic modulation of aerodynamic accommodation angles.
  - Plasma/charging surrogate implemented as stochastic scaling of magnetic-field input.
  - Structural flex surrogate implemented as stochastic scaling on `srp_scale` and aero-angle coupling.
  - Fuel gauging uncertainty represented as an explicit mass-uncertainty surrogate channel.

Scope notes:

- Earth tides/loading and albedo/IR-specific SRP extensions remain outside this closure pass.
- Surrogate channels are designed for campaign-scale UQ realism, not high-fidelity multiphysics replacement.

Evidence:

- `python/vleo_uq/pod_uq.py`
- `python/vleo_uq/attitude_filter.py`
- `python/vleo_uq/__init__.py`
- `scripts/model_discrepancy.py`
- `scripts/uq_parameter_channels.py`
- `docs/model_discrepancy.md`
- `tests/test_pod_measurements.py`
- `tests/test_attitude_filter_fullstate.py`
- `tests/test_model_discrepancy.py`
- `tests/test_uq_parameter_channels.py`

## 2026-02-20: Final closure pass for remaining physics/events gaps

Status: accepted

Decision:

- Close the last remaining checklist gaps with explicit runtime implementations:
  - Add albedo/IR acceleration channel in core propagator.
  - Add tide/loading acceleration channel in core propagator.
  - Add explicit rarefied-regime switching surrogate channel in discrepancy model.
  - Add continuous finite-burn event model (acceleration-duration driven) with stochastic sampling.

Scope notes:

- New channels are intentionally lightweight and campaign-friendly.
- Focus is implementation completeness; high-fidelity calibration remains separate validation work.

Evidence:

- `cpp/include/vleo/propagator.hpp`
- `cpp/src/propagator.cpp`
- `cpp/bindings/bindings.cpp`
- `python/vleo_uq/events.py`
- `python/vleo_uq/env_sources.py`
- `scripts/model_discrepancy.py`
- `tests/test_propagator_regression.py`
- `tests/test_events.py`
- `tests/test_model_discrepancy.py`
