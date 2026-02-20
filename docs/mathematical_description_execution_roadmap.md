# Mathematical Description Implementation Roadmap

Source checklist: `docs/mathematical_description_implementation_checklist.md`
Source theory: `Protocoll/4S Paper/Mathematical_description.md`

## Goal

Close all `Partial` and `Missing` checklist items by either:
- implementing them in the framework, or
- explicitly deferring them with thesis-scope justification.

## Priorities

- `P0` required for strong paper-grade orbit/attitude UQ claims.
- `P1` important extensions for robustness and operational realism.
- `P2` extended research channels; implement only if in scope, otherwise defer with rationale.

## Progress (2026-02-20)

| ID | Progress |
|---|---|
| `MD-01` | Done (scope freeze entry added to decision log) |
| `MD-02` | Done (composition discrepancy channel implemented in `scripts/model_discrepancy.py`) |
| `MD-03` | Done (storm jump/regime channel implemented in `scripts/model_discrepancy.py`) |
| `MD-04` | Done (optional J2 geopotential term added in propagator + regression test) |
| `MD-05` | Done (optional Sun/Moon third-body acceleration added in propagator + regression test) |
| `MD-06` | Done (SRP acceleration coupling added in core propagator with `srp_scale`) |
| `MD-07` | Done (magnetic residual dipole torque channel added in propagated dynamics) |
| `MD-08` | Done (explicit rarefied-regime switch surrogate implemented in discrepancy pipeline) |
| `MD-09` | Done (continuous finite-burn event model implemented via pulse expansion + stochastic burn sampling) |
| `MD-10` | Done (explicit latency/jitter measurement mapping wired through POD OD paths) |
| `MD-11` | Done (magnetometer measurement model integrated in MEKF/UKF/ENKF full-state filters) |
| `MD-12` | Done (AO erosion surrogate channel added to model discrepancy pipeline) |
| `MD-13` | Done (plasma/charging surrogate channel added to model discrepancy pipeline) |
| `MD-14` | Done (structural/fuel/software-human surrogate channels added: structural flex, fuel gauging mass surrogate, operational outage model) |

## Issue List

| ID | Priority | Work Item | Target Status | Main Code Targets | Validation Artifact | Effort |
|---|---|---|---|---|---|---|
| MD-01 | P0 | Scope freeze for all `Partial`/`Missing` entries | Implement or explicitly defer each item | `docs/decision_log.md`, `docs/mathematical_description_implementation_checklist.md` | Updated decision log + checklist | S |
| MD-02 | P0 | Add composition uncertainty channel (state or surrogate) | Implemented | `python/vleo_uq/env_sources.py`, `scripts/model_discrepancy.py`, `scripts/uq_parameter_channels.py`, `scripts/run_case_studies.py` | Composition toggle test + sensitivity report | M |
| MD-03 | P0 | Add storm regime-switch/jump model for thermosphere disturbance | Implemented | `scripts/model_discrepancy.py`, `python/vleo_uq/env_sources.py`, `scripts/run_case_studies.py` | Storm jump validation plot/report | M |
| MD-04 | P0 | Add at least J2 geopotential option in propagator | Implemented | `cpp/include/vleo/propagator.hpp`, `cpp/src/propagator.cpp`, `cpp/bindings/bindings.cpp` | New propagator regression test (J2 on/off) | M |
| MD-05 | P0 | Add third-body acceleration (Sun/Moon) with toggles | Implemented | `cpp/include/vleo/propagator.hpp`, `cpp/src/propagator.cpp`, `cpp/bindings/bindings.cpp`, `python/vleo_uq/env_sources.py` | Deterministic delta check + tolerance sweep update | M |
| MD-06 | P0 | Add SRP acceleration coupling using existing eclipse/srp scale | Implemented | `cpp/include/vleo/propagator.hpp`, `cpp/src/propagator.cpp`, `scripts/run_case_studies.py`, `python/vleo_uq/env_sources.py` | SRP on/off scenario comparison | M |
| MD-07 | P1 | Add magnetic residual dipole disturbance torque model | Implemented | `cpp/include/vleo/propagator.hpp`, `cpp/src/propagator.cpp`, `scripts/uq_parameter_channels.py` | Attitude disturbance test case | M |
| MD-08 | P1 | Add rarefied-regime uncertainty representation (Knudsen/model-form blend) | Implemented | `scripts/model_discrepancy.py`, `force_moment_module_v4/src/aerodynamics.cpp`, `tests/test_model_discrepancy.py` | Coefficient uncertainty sweep report | M |
| MD-09 | P1 | Add continuous finite-burn thrust/noise model (beyond impulsive delta-v) | Implemented | `python/vleo_uq/events.py`, `scripts/run_case_studies.py`, `tests/test_events.py` | Maneuver realism gate extension | M |
| MD-10 | P1 | Add explicit delayed-measurement/latency model in POD | Implemented | `python/vleo_uq/pod_uq.py`, `scripts/postprocess_uq_pod.py` | Latency stress scenario report | M |
| MD-11 | P1 | Complete sensor suite extension with magnetometer model | Implemented | `python/vleo_uq/attitude_filter.py`, `python/vleo_uq/__init__.py`, `scripts/example_attitude_filter*.py` | Attitude filter consistency test | S |
| MD-12 | P2 | AO erosion surrogate channel affecting aero coefficients | Implemented | `scripts/model_discrepancy.py`, `docs/model_discrepancy.md`, `tests/test_model_discrepancy.py` | AO surrogate sensitivity run | M |
| MD-13 | P2 | Plasma/charging surrogate disturbance channel | Implemented | `scripts/model_discrepancy.py`, `docs/model_discrepancy.md`, `tests/test_model_discrepancy.py` | Charging disturbance scenario report | M |
| MD-14 | P2 | Structural/fuel/software-human uncertainty channels (surrogate or defer) | Implemented | `scripts/uq_parameter_channels.py`, `scripts/model_discrepancy.py`, `python/vleo_uq/pod_uq.py`, `tests/test_uq_parameter_channels.py`, `tests/test_pod_measurements.py` | Surrogate-channel regression tests | M |

## Execution Order

1. `MD-01` (scope freeze).
2. `MD-02` + `MD-03` (thermosphere uncertainty completeness).
3. `MD-04` + `MD-05` + `MD-06` (core force-model completeness).
4. `MD-07` + `MD-08` + `MD-09` (attitude and aero robustness extensions).
5. `MD-10` + `MD-11` (measurement/estimation realism extensions).
6. `MD-12` + `MD-13` + `MD-14` (extended channels or formal deferral).

## Suggested Sprinting

- Sprint A (P0-1): `MD-01..MD-03`
- Sprint B (P0-2): `MD-04..MD-06`
- Sprint C (P1): `MD-07..MD-11`
- Sprint D (P2): `MD-12..MD-14`

## Definition of Done

- Each issue in `MD-01..MD-14` is implemented with tests/docs in current framework state.
- `docs/mathematical_description_implementation_checklist.md` has no `Partial` or `Missing` entries.
- Validation scripts/reports can now focus on calibration and campaign evidence generation.
