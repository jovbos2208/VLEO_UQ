# Mathematical Description Implementation Checklist

Source: `Protocoll/4S Paper/Mathematical_description.md`
Purpose: track whether each mathematical uncertainty model described there is implemented in this framework.
Execution roadmap: `docs/mathematical_description_execution_roadmap.md`

Legend:
- `Check`: `[x]` implemented, `[ ]` not fully implemented
- `Status`: `Implemented`, `Partial`, `Missing`

## A) Stochastic Building Blocks

| Check | Item | Status | Implementation Evidence | Gap / Note |
|---|---|---|---|---|
| [x] | White Gaussian measurement noise | Implemented | `python/vleo_uq/pod_uq.py` (`simulate_gnss_measurements`, `simulate_slr_measurements`), `python/vleo_uq/attitude_filter.py` | Used in GNSS/SLR/gyro/star-tracker simulation |
| [x] | Random walk (RW) bias states | Implemented | `python/vleo_uq/pod_uq.py` (clock/tropo RW), `cpp/src/propagator.cpp` (`tau<=0` latent update) | RW appears in both measurement and process models |
| [x] | First-order Gauss-Markov / OU | Implemented | `cpp/src/propagator.cpp` (latent OU update), `scripts/model_discrepancy.py` | Core latent uncertainty model |
| [x] | Discrete AR(1) form of GM(1) | Implemented | `cpp/src/propagator.cpp` (phi = exp(-dt/tau)), `scripts/model_discrepancy.py` | Equivalent discrete OU step used |
| [x] | Multiplicative/lognormal scaling (e.g., density) | Implemented | `cpp/src/propagator.cpp` (`rho_eff = rho*exp(log terms)`), `scripts/uq_parameter_channels.py`, `scripts/model_discrepancy.py` | Positive-definite scaling implemented |
| [x] | Heavy-tailed process option | Implemented | `scripts/model_discrepancy.py` (`_student_t_unit_variance`) | Implemented for model-discrepancy channel |

## B) VLEO Uncertainty Sources (from Section "Mathematical and quantitative models...")

| Check | Item | Status | Implementation Evidence | Gap / Note |
|---|---|---|---|---|
| [x] | Neutral thermospheric density uncertainty | Implemented | `python/vleo_uq/env_sources.py` (`NRLMSIS21DensityModel`), `cpp/src/propagator.cpp` (log-density latent states), `scripts/model_discrepancy.py` | Includes deterministic model + stochastic correction |
| [x] | Thermospheric wind uncertainty | Implemented | `python/vleo_uq/env_sources.py` (`HWM14WindModel`), `cpp/src/propagator.cpp` (wind bias latent states) | Wind bias OU channel present |
| [x] | Composition variability as explicit stochastic state | Implemented | `python/vleo_uq/env_sources.py`, `scripts/model_discrepancy.py` (composition OU/RW surrogate on `particles_mass_kg`) | Implemented as stochastic surrogate channel (no new propagated state dimension) |
| [x] | GSI / aerodynamic parameter uncertainty | Implemented | `force_moment_module_v4/*` (energy accommodation and aero model), `cpp/include/vleo/propagator.hpp` (`eta1_rad`, `eta2_rad`, `temperature_ratio_method`) | Aerodynamic parameterization present |
| [x] | Rarefied regime / Knudsen switching uncertainty | Implemented | `scripts/model_discrepancy.py` (`rarefied_regime_switch_*`), `force_moment_module_v4/*`, `tests/test_model_discrepancy.py` | Explicit stochastic regime-switch surrogate implemented via `temperature_ratio_method` switching |
| [x] | Atomic oxygen flux and erosion | Implemented | `scripts/model_discrepancy.py` (`ao_erosion_*` surrogate channel), `tests/test_model_discrepancy.py` | Implemented as stochastic surrogate on aerodynamic accommodation angles (`eta1_rad`, `eta2_rad`) |
| [x] | Plasma interactions and charging | Implemented | `scripts/model_discrepancy.py` (`plasma_charging_*` surrogate channel), `tests/test_model_discrepancy.py` | Implemented as stochastic surrogate on magnetic-field input scaling |
| [x] | Stochastic aerodynamic force/torque residual channel | Implemented | `cpp/src/propagator.cpp` latent OU processes, `scripts/model_discrepancy.py` | Implemented via latent and discrepancy channels |
| [x] | Space-weather index uncertainty drivers (F10.7/Kp/Ap/Dst ingestion) | Implemented | `python/vleo_uq/env_sources.py` (`parse_omni2_indices`, forecast parsers), `scripts/run_case_studies.py` | Ingestion/interpolation implemented |
| [x] | Explicit storm jump/regime-switch model | Implemented | `scripts/model_discrepancy.py` (`storm_jump_*` Poisson-start piecewise jump regime) | Piecewise storm-jump disturbance channel implemented |
| [x] | Higher-order geopotential coefficient uncertainty | Implemented | `cpp/src/propagator.cpp`, `cpp/include/vleo/propagator.hpp`, `cpp/bindings/bindings.cpp`, `scripts/uq_parameter_channels.py` | J2/J3/J4 zonal gravity options with coefficient-uncertainty parameter channels |
| [x] | Earth tides/loading uncertainty | Implemented | `cpp/src/propagator.cpp`, `cpp/include/vleo/propagator.hpp`, `cpp/bindings/bindings.cpp`, `python/vleo_uq/env_sources.py`, `tests/test_propagator_regression.py` | Tide/loading acceleration channel propagated with scale uncertainty path |
| [x] | Third-body ephemeris uncertainty | Implemented | `cpp/src/propagator.cpp`, `cpp/include/vleo/propagator.hpp`, `cpp/bindings/bindings.cpp`, `python/vleo_uq/env_sources.py`, `scripts/uq_parameter_channels.py` | Sun/Moon third-body force with ephemeris-scale uncertainty channels (`sun_ephemeris_scale`, `moon_ephemeris_scale`) |
| [x] | SRP / albedo / IR force uncertainty | Implemented | `cpp/src/propagator.cpp`, `cpp/include/vleo/propagator.hpp`, `cpp/bindings/bindings.cpp`, `scripts/uq_parameter_channels.py`, `tests/test_propagator_regression.py` | SRP + albedo/IR acceleration channels with uncertainty parameterization implemented |
| [x] | Magnetic field + residual dipole uncertainty | Implemented | `cpp/src/propagator.cpp`, `cpp/include/vleo/propagator.hpp`, `python/vleo_uq/env_sources.py`, `cpp/bindings/bindings.cpp`, `scripts/uq_parameter_channels.py` | Dipole-field surrogate + residual dipole torque and uncertainty scale channels implemented |
| [x] | Structural flexing / thermal deformation | Implemented | `scripts/model_discrepancy.py` (`structural_flex_*` surrogate channel), `tests/test_model_discrepancy.py` | Implemented as stochastic surrogate on `srp_scale` and accommodation-angle coupling |
| [x] | Mass and inertia uncertainty | Implemented | `scripts/uq_parameter_channels.py` (`uq_mass_*`, `uq_inertia_rel_sigma`) | Parameter-randomization channel exists |
| [x] | Fuel gauging uncertainty | Implemented | `scripts/uq_parameter_channels.py` (`uq_fuel_gauging_sigma_kg`, `uq_fuel_gauging_rel_sigma`), `tests/test_uq_parameter_channels.py` | Implemented as mass-uncertainty surrogate (no explicit propellant tank state) |
| [x] | Continuous thrust magnitude/direction/noise model | Implemented | `python/vleo_uq/events.py` (`finite_burn` event with acceleration/duration and stochastic sampling), `tests/test_events.py` | Continuous finite-burn represented via per-step pulse expansion and stochastic burn-parameter sampling |

## C) Measurement, Estimation, and Operations Uncertainty

| Check | Item | Status | Implementation Evidence | Gap / Note |
|---|---|---|---|---|
| [x] | GNSS code/phase noise, multipath-like outliers, cycle slips | Implemented | `python/vleo_uq/pod_uq.py` (`simulate_gnss_measurements`) | Includes dropout, slips, outlier inflation |
| [x] | SLR measurement noise and station bias | Implemented | `python/vleo_uq/pod_uq.py` (`simulate_slr_measurements`, `estimate_slr_bias`) | Includes weather-driven availability mask |
| [x] | Star tracker + IMU + magnetometer full sensor suite | Implemented | `python/vleo_uq/attitude_filter.py`, `tests/test_attitude_filter_fullstate.py` | MEKF/UKF/ENKF full-state paths support magnetometer measurements |
| [x] | Timing/clock error models | Implemented | `python/vleo_uq/pod_uq.py` (clock bias/drift RW states) | Implemented in truth + OD nuisance estimation |
| [x] | Parameter correlation/identifiability handling in OD | Implemented | `python/vleo_uq/pod_uq.py` (batch WLS/MAP with nuisance states, damping/priors) | Correlations handled via normal equations/prior terms |
| [x] | Operational latency model | Implemented | `python/vleo_uq/pod_uq.py` (`apply_measurement_latency`, latency-aware OD paths), `tests/test_pod_measurements.py` | Explicit delayed measurement mapping with configurable latency/jitter/seed |
| [x] | Maneuver execution uncertainty | Implemented | `python/vleo_uq/events.py`, `scripts/run_case_studies.py` | `delta_v` timing/scale perturbations supported |
| [x] | Software/human factor statistical model | Implemented | `python/vleo_uq/pod_uq.py` (`sample_operational_outage_mask`, GNSS/SLR outage controls), `tests/test_pod_measurements.py` | Implemented as Poisson-start operational outage process affecting measurement availability |
| [x] | Generic measurement bias + stochastic noise abstraction | Implemented | `python/vleo_uq/pod_uq.py` (clock/tropo/ambiguity/slr bias states) | Cross-cutting bias/noise model present |

## D) Worked Example Readiness (from "Worked quantitative examples")

| Check | Item | Status | Implementation Evidence | Gap / Note |
|---|---|---|---|---|
| [x] | Density scale to drag acceleration example can be reproduced | Implemented | `cpp/src/propagator.cpp`, `scripts/model_discrepancy.py`, `scripts/uq_parameter_channels.py` | All required terms available |
| [x] | Random-walk drag-like coefficient behavior reproducible (via stochastic channels) | Implemented | `scripts/model_discrepancy.py`, `cpp/src/propagator.cpp` latent channels | Direct `C_D` state not explicit, but equivalent multiplicative channel exists |
| [x] | GNSS carrier noise to position covariance study | Implemented | `python/vleo_uq/pod_uq.py` + `scripts/postprocess_uq_pod.py` | UQ/POD pipeline computes coverage and RTN metrics |
| [x] | Clock bias to range error and RW/GM modeling | Implemented | `python/vleo_uq/pod_uq.py` (clock state evolution and estimation) | Fully represented |

## E) Quick Completion Metrics

- Implemented: 38
- Partial: 0
- Missing: 0

## F) Priority Closure Order (recommended)

1. All checklist items are now implemented in-framework.
2. Remaining work is calibration/validation depth (not implementation completeness).
