# Further Progression Plan — VLEO Aerodynamics & UQ Framework (Roadmap)

This document defines the **next development and research steps** for the hybrid C++/Python VLEO simulation + UQ framework, based on the current Codex overview:

- **C++ core**: `propagator.hpp/.cpp`, `aero_adapter.cpp`  
  - 18-state coupled propagation (orbit + attitude + latent atmosphere states)  
  - mesh-based aero force/torque with shadowing and hinge rotations  
  - deterministic / MC / UT / STM support  
- **Scenario layer**: `scenario_catalog.md`, `catalog_to_case_config.py`  
  - catalog → runnable JSON configs, ATT/FORM/OD mapping  
- **Orchestration & HPC**: SLURM scripts + job-array runners  
  - DET/STM/MC/UT parallel runs, shard merging, postprocessing  
- **Environment + POD + sensors**: Python package (`run_case_studies.py`, etc.)  
  - space-weather ingestion, winds, GNSS/SLR simulation, POD UQ, attitude filtering (MEKF/UKF/EnKF)

---

# Maneuver Realism Execution Status (Current Sprint)

This section tracks execution of the ATT maneuver realism plan (wing-angle evolution and physical plausibility).

- [x] `P1` Actuator realism limits in runtime maneuver control  
  Implemented in `scripts/run_mc_job_array.py` with:
  - `eta_max_deg = 110`
  - `eta_rate_max_deg_s = 5`
  - `eta_accel_max_deg_s2 = 1`
  Limits are configurable per-scenario and via env vars.

- [x] `P2` Replace instantaneous sign flips with finite-rate actuator evolution  
  Implemented as rate/acceleration-limited command tracking toward phase targets in `scripts/run_mc_job_array.py`.

- [x] `P3` Improve control semantics beyond pure waveform fallback  
  Implemented closed-loop maneuver semantics in `scripts/run_mc_job_array.py` for:
  - `CTL_ATT_AERO_RATE` (axis-rate tracking command),
  - `CTL_ATT_AERO_POINT` (angle/rate feedback command).
  Actuator rate/accel limits remain enforced in the same runtime path.

- [x] `P4` Add maneuver robustness diagnostics  
  Added actuator/maneuver stats in `scripts/postprocess_mc_only.py`, including:
  control switch count, max command jumps/rates, saturation fraction, full-turn reach indicators, and reached fraction parsing.

- [x] `P5` Wing-sweep physical-domain guidance  
  Added explicit `in_recommended_range` (`|wing_deg| <= 60`) flag in `scripts/summarize_att_wing_sweep.py`.  
  Default sweep-policy enforcement is now active in `scripts/build_att_wing_sweep_config.py` and `ATT_baseline.slurm` with explicit opt-in override.

- [x] `P6` Static authority map gate before time-domain sweeps  
  Implemented in `scripts/validate_maneuver_realism_gate.py` (`--stage pre`) and integrated in `run_all_A_scenarios.slurm`.

- [x] `P7` Tighten maneuver-study initial attitude-rate spread  
  Implemented scenario-aware covariance tightening in:
  - `scripts/mc_scenarios.py`
  - `scripts/ut_scenarios.py`
  with `initial_rate_sigma_deg_s` support (default `0.5 deg/s` for maneuver scenarios).

- [x] `P8` Formal validation gates (pass/fail criteria) in pipeline  
  Implemented in `scripts/validate_maneuver_realism_gate.py` (`--stage post`) and integrated in `run_all_A_scenarios.slurm` as hard-stop checks.

---

## 1) Immediate priorities (1–2 weeks): Make the framework auditable

### 1.1 Create a single “architecture truth” page
**Deliverable:** `docs/architecture.md` with:
- module boundaries (C++ core vs Python orchestration)
- ownership map (which file implements what)
- run pipeline (config → sampling → propagation → sensors → estimation → outputs)

**Action items**
- Document the 18-state definition explicitly (state ordering, units, frames).
- Add a diagram of the call chain from Python to C++ (adapter boundary).

### 1.2 Add a deterministic regression test suite
**Goal:** ensure changes do not silently break dynamics/aero/measurement models.

**Deliverables**
- `tests/test_propagator_regression.*`
- fixed-seed test cases:
  - propagation-only (no measurements)
  - propagation + aero (shadowing on/off)
  - GNSS/SLR measurement simulation sanity checks

**Acceptance criteria**
- deterministic trajectory matches baseline within tolerance (pos/vel/attitude).
- aero force/torque magnitudes remain within physically plausible bands.

---

## 2) Core model completeness (2–6 weeks): Measurement models + state coupling

### 2.1 Measurement models: make them publication-grade
You use **raw GNSS code + carrier**. This must be explicit in theory and code.

**Deliverables**
- `docs/measurement_models.md` with:
  - GNSS pseudorange + carrier phase equations, biases, clock model, ambiguity states
  - cycle-slip/outlier model (for stress scenario C9)
  - SLR range equation + station bias model + visibility masks
  - accelerometer measurement equation + bias/random walk model (scenario C3)

**Implementation checks**
- Confirm which terms are modeled vs neglected:
  - ionosphere/troposphere handling (even if simplified in VLEO)
  - relativistic/light-time corrections (likely negligible but state explicitly)
  - antenna phase center offsets and lever arms

### 2.2 Explicit coupling of latent atmosphere states into aero
Codify exactly how the latent states enter the force model:
- log-density factor: ρ = ρ_model · exp(δρ)
- wind bias: w = w_model + δw

**Deliverable**
- `docs/state_and_coupling.md`: state vector + dynamics + where each state enters forces/torques and measurement models.

---

## 3) UQ methodology strengthening (1–2 months): beyond MC/UT runs

### 3.1 Variance decomposition / sensitivity analysis
Thesis-quality UQ needs attribution:
- “How much uncertainty comes from density vs GSI vs attitude vs measurement bias?”

**Deliverables**
- Global sensitivity (Sobol indices) for QoIs:
  - drag acceleration a_D, lifetime/decay, pointing error, formation separation
- Local sensitivity maps using STM/Jacobians already supported

**Implementation path**
- Add `postprocess_sobol.py`
- Use existing MC samples; compute indices on QoIs.

### 3.2 Explicit model-form uncertainty channel
Parameter uncertainty is covered (e.g., α_E), but model-form uncertainty is weakly represented.

**Recommended pragmatic approach**
- introduce a model discrepancy term:
  - additive acceleration bias in drag direction (stochastic process), or
  - multiplicative correction with heavier tails than OU
- calibrate discrepancy using:
  - SLR residual structure and/or accelerometer comparisons
  - DSMC/truth-model comparisons for selected attitudes

**Deliverable**
- `docs/model_discrepancy.md` + implementation toggle (e.g., `DISCREPANCY_ON`).

---

## 4) Numerical reliability (parallel): Integrator + stochastic discretization

### 4.1 Integrator documentation + tolerances
Show numerical error is below estimation/UQ noise.

**Deliverables**
- `docs/numerics.md` documenting:
  - integrator type (DOP853/RK4/etc.), tolerances, step rejection
  - OU discretization method (exact discretization preferred)
  - STM propagation method and consistency checks

**Sanity checks**
- tolerance sweep: convergence of pos/vel/attitude
- OU stats: simulated mean/variance/autocorrelation vs theory

---

## 5) Scenario progression plan (what to run, in what order)

### Phase A — single-sat verification (ATT + OD basics)
1. ATT_A1 detumble (DET + MC; check torque authority realism)
2. ATT_A2 nadir pointing (UT vs MC comparison)
3. OD_C1 / OD_C2 GNSS-only (baseline OD consistency)
4. OD_C3 GNSS+accel (drag/parameter decoupling)

### Phase B — robustness and stress tests
1. ATT_A7 eclipse transition
2. OD_C9 measurement error stress tests
3. OD_C8 weather-limited SLR availability
4. FORM_B6 storm transient formation keeping

### Phase C — formation mission questions
1. FORM_B2 along-track keeping (2-sat)
2. FORM_B3 LVLH constraint (3-sat)
3. FORM_B5 close-approach risk probability
4. FORM_B7 cluster sparse contact

**Deliverable**
- `docs/experiment_matrix.md` (scenario × env toggles × estimator × UQ method × seeds)

---

## 6) Publication and dissertation alignment

### 6.1 Candidate paper splits
- Paper 1 (methods): coupled stochastic environment + GSI parameter estimation using GNSS/SLR
- Paper 2 (controls): aerodynamic panel attitude control robustness under thermospheric uncertainty
- Paper 3 (formation): differential-drag formation feasibility and close-approach risk under uncertainty

### 6.2 Thesis chapter mapping
- Ch2: rarefied aero + GSI + environment uncertainty
- Ch3: coupled dynamics + measurement models + estimation
- Ch4: UQ methods (MC/UT/STM + variance decomposition)
- Ch5: attitude control scenarios
- Ch6: formation + ground segment + operational constraints
- Ch7: conclusions + limitations + future work

---

## 7) Concrete “next commits” checklist

**Docs**
- [ ] docs/architecture.md
- [ ] docs/state_and_coupling.md
- [ ] docs/measurement_models.md
- [ ] docs/numerics.md
- [ ] docs/experiment_matrix.md

**Code**
- [ ] regression tests for propagator + aero adapter
- [ ] OU exact discretization validation test
- [ ] cycle-slip/outlier injection module for GNSS (C9)
- [ ] station weather availability model for SLR (C8)
- [ ] Sobol/sensitivity postprocessing

**CI/HPC**
- [ ] add “small deterministic run” in CI (if feasible)
- [ ] standardize output schema (one JSON summary per run)
- [ ] enforce metadata: git hash, scenario_id, seed, toggles

---

## 8) Open questions to resolve (explicit decisions)

1. Exact 18-state definition (ordering + meaning + units + frames)?
2. Baseline GSI model in code: Maxwell-only, CLL, or Tuttas generalized model?
3. Primary OD estimator baseline: EKF vs EnKF vs batch LS?
4. Definition of “truth” for validation:
   - SLR/GNSS residuals only?
   - accelerometer-assisted truth?
   - DSMC subset truth?

---
End of document.
