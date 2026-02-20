# Combined Implementation Check Plan

Scope: unify `Plan.md` and `Plan_progression.md` into one executable plan that verifies every requested capability, document, and validation gate.

## 1) Control Objectives

- [x] Build and expose the hybrid C++/Python stack with reproducible high-fidelity propagation.
- [x] Verify deterministic, MC, UT, and STM behavior against explicit acceptance criteria.
- [x] Verify uncertainty source coverage (initial state, parameters, latent processes, discrete events, model discrepancy).
- [x] Verify POD UQ, attitude UQ, payload impact, and sensitivity workflows.
- [x] Verify maneuver realism status points `P1..P8`.
- [x] Verify documentation, CI/HPC metadata, and scenario campaign traceability.

## 2) Decision Gates (must be closed first)

From `Plan_progression.md` section 8:

- [x] DG-1: Freeze exact 18-state definition (ordering, units, frames).
- [x] DG-2: Freeze baseline GSI model (Maxwell vs CLL vs Tuttas generalized).
- [x] DG-3: Freeze OD baseline estimator (EKF vs EnKF vs batch LS).
- [x] DG-4: Freeze validation truth definition (residual-only vs accelerometer-assisted vs DSMC subset).

Exit criterion: one signed decision record in `docs/decision_log.md` and referenced in all configs.

## 3) Traceability Matrix (what must be checked)

### 3.1 Architecture, ABI, and Build

- [x] A-1 (`Plan.md` A): C++/Python boundary follows performance split (no Python callback in propagation inner loop).
- [x] A-2 (`Plan.md` B): force/torque adapter ABI implemented (`EnvInputs`, `VehicleParams`, concrete `AeroAdapter` boundary).
- [x] A-3 (`Plan.md` B): thread safety guaranteed (immutable/shared-safe model or per-thread clone).
- [x] A-4 (`Plan.md` C, Ticket 1): pybind11 + CMake + scikit-build-core build/import path works (`pip install -e .`, import extension).
- [x] A-5 (`Plan_progression` 1.1): architecture truth page exists and reflects current code path.

### 3.2 Dynamics and Propagation Methods

- [x] D-1 (`Plan.md` D1, Ticket 2): deterministic SP truth baseline with high-order integrator and event handling.
- [x] D-2 (`Plan.md` D2, Ticket 3): MC ensemble propagator supports contiguous batch input/output and parallel execution.
- [x] D-3 (`Plan.md` D2, Ticket 5): UT sigma-point propagator implemented and compared to MC.
- [x] D-4 (`Plan.md` D2, Ticket 6): STM/covariance propagation available and consistency-checked.
- [x] D-5 (`Plan.md` H): integrator tolerance convergence demonstrated.

### 3.3 Uncertainty Source Coverage

- [x] UQ-1 (`Plan.md` E1): uncertain initial state sampling implemented and reproducible.
- [x] UQ-2 (`Plan.md` E2/E3): constant and slowly varying parameter uncertainty channels implemented.
- [x] UQ-3 (`Plan.md` E4, Ticket 4): latent OU/random-walk processes integrated with exact/validated discretization.
- [x] UQ-4 (`Plan.md` E5): discrete random events modelled (maneuver errors, mode transition timing).
- [x] UQ-5 (`Plan_progression` 3.2): explicit model discrepancy channel implemented with runtime toggle.

### 3.4 Measurement, POD, Attitude, Payload

- [x] MP-1 (`Plan_progression` 2.1): publication-grade measurement-model docs for GNSS/SLR/accelerometer are complete.
- [x] MP-2 (`Plan_progression` 2.1): cycle-slip/outlier module active for C9 stress tests.
- [x] MP-3 (`Plan_progression` 2.1): SLR station weather availability model active for C8.
- [x] MP-4 (`Plan.md` F1, Ticket 7): POD UQ harness runs truth -> measurements -> OD -> posterior metrics.
- [x] MP-5 (`Plan.md` F2): attitude UQ chain supports both fast and full-fidelity modes.
- [x] MP-6 (`Plan.md` F3, Ticket 8): payload impact metrics computed from orbit+attitude uncertainty.

### 3.5 Maneuver Realism (`P1..P8`)

- [x] MR-1 (`P1`): actuator limits implemented.
- [x] MR-2 (`P2`): finite-rate/finite-accel command evolution implemented.
- [x] MR-3 (`P3`): closed-loop `CTL_ATT_AERO_RATE` / `CTL_ATT_AERO_POINT` completed.
- [x] MR-4 (`P4`): maneuver diagnostics in postprocess completed.
- [x] MR-5 (`P5`): recommended wing-sweep range enforced as default policy (if approved).
- [x] MR-6 (`P6`): static authority map gate runs before time-domain sweeps.
- [x] MR-7 (`P7`): maneuver-study initial rate covariance tightening active.
- [x] MR-8 (`P8`): formal pass/fail validation gates enforced in pipeline.

### 3.6 Sensitivity, Validation, CI/HPC

- [x] SV-1 (`Plan.md` G, Ticket 9): Morris screening + Sobol workflow implemented.
- [x] SV-2 (`Plan_progression` 3.1): grouped then expanded variance attribution reported for key QoIs.
- [x] SV-3 (`Plan.md` H): MC convergence, coverage calibration, and ablation checks are automated.
- [x] SV-4 (`Plan_progression` 4.1): numerics documentation includes integrator and stochastic discretization checks.
- [x] SV-5 (`Plan_progression` 7): CI small deterministic run, unified output schema, metadata enforcement.

### 3.7 Scenario Campaign Progression

- [x] SC-1 (`Plan_progression` Phase A): run ATT_A1, ATT_A2, OD_C1/C2/C3 with DET/MC/UT as specified.
- [x] SC-2 (`Plan_progression` Phase B): run ATT_A7, OD_C9, OD_C8, FORM_B6 stress cases.
- [x] SC-3 (`Plan_progression` Phase C): run FORM_B2/B3/B5/B7 mission-level formation studies.
- [x] SC-4 (`Plan_progression` 5): experiment matrix file tracks scenario x toggles x estimator x UQ x seeds.
  Evidence run: `results/run_campaign_quick_20260219d` via `scripts/run_phase_campaign.py`.

## 4) Execution Sequence (implementation order)

### Sprint 0 (Week 1): Freeze definitions and audit baseline

- [x] Close DG-1..DG-4.
- [x] Refresh `docs/architecture.md` and add missing state/coupling references.
- [x] Record current status of `P1..P8` with evidence links.

### Sprint 1 (Week 1-2): Deterministic and build reliability

- [x] Complete deterministic regression tests (propagation-only, aero on/off, measurement sanity).
- [x] Validate build/import path and adapter interface checks.
- [x] Add CI deterministic smoke run and metadata assertions.

### Sprint 2 (Week 2-4): Stochastic core and propagation parity

- [x] Finalize OU exact discretization tests and seeded reproducibility.
- [x] Run MC/UT/STM parity suite on controlled cases; compute coverage vs MC benchmark.
- [x] Add integrator tolerance and MC convergence studies with archived plots/tables.

### Sprint 3 (Week 4-6): Measurement/POD and maneuver realism gaps

- [x] Complete `MR-3`, `MR-6`, `MR-8` (and `MR-5` if policy approved).
- [x] Implement cycle-slip/outlier and weather-limited SLR modules.
- [x] Run POD UQ harness and publish RTN sigma + RMS + coverage outputs.

### Sprint 4 (Week 6-8): Sensitivity and campaign runs

- [x] Implement Sobol/Morris postprocess and grouped attribution outputs.
- [x] Run Phase A/B/C scenario matrix with standardized summaries.
- [x] Add model discrepancy toggle study and calibration evidence.

## 5) Validation Gates (hard pass/fail)

- [x] VG-1 Numerics: tolerance sweep converges (state deltas below threshold at final and key horizons).
- [x] VG-2 Stochastic processes: OU mean/variance/autocorrelation match theory within tolerance.
- [x] VG-3 Method agreement: UT/STM covariance realism (1-sigma/2-sigma/3-sigma coverage) acceptable against MC.
- [x] VG-4 Maneuver realism: saturation, rate, reachability, and authority-map gate all pass before batch sweeps.
- [x] VG-5 Measurement realism: GNSS/SLR/accelerometer residual behavior consistent with configured noise/bias models.
- [x] VG-6 Reproducibility: same seed gives same statistics independent of thread count.

## 6) Required Deliverables (to close both plans)

Docs:

- [x] `docs/architecture.md` (updated)
- [x] `docs/state_and_coupling.md`
- [x] `docs/maneuver_realism_status.md`
- [x] `docs/measurement_models.md`
- [x] `docs/numerics.md`
- [x] `docs/model_discrepancy.md`
- [x] `docs/experiment_matrix.md`
- [x] `docs/decision_log.md`

Code/tests/pipeline:

- [x] deterministic regression tests for propagator+aero
- [x] OU exact discretization validation test
- [x] cycle-slip/outlier injection (C9)
- [x] SLR weather availability model (C8)
- [x] Sobol/Morris postprocessing script(s)
- [x] CI deterministic smoke run + metadata enforcement

Outputs:

- [x] standardized per-run JSON summary with hash/scenario/seed/toggles
- [x] MC convergence report
- [x] coverage report (MC vs UT vs STM)
- [x] maneuver realism report (`P1..P8` status + pass/fail)
- [x] Phase A/B/C scenario campaign report

## 7) Done Definition

Both source plans are considered checked when:

- [x] every checklist item in sections 3-6 is either complete or explicitly waived by decision record,
- [x] all validation gates `VG-1..VG-6` pass,
- [x] scenario campaign outputs and documents are reproducible from repository scripts.
