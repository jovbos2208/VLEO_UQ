# Experiment Matrix

This matrix captures the paper-prep campaign order from `Plan_progression.md` with executable settings.

## 1) Phase A (single-sat verification)

| Scenario | Objective | Estimator | UQ Methods | Seed Set | Key Toggles |
|---|---|---|---|---|---|
| `ATT_A1_DETUMBLE_AERO_ONLY` | detumble authority realism | attitude filters + maneuver diagnostics | DET + MC (+UT optional) | `42, 43, 44` | `VLEO_MC_SAFE_MODE=1` |
| `ATT_A2_NADIR_POINTING` | nadir tracking robustness | attitude filters | DET + MC + UT | `42, 43, 44` | env sources on/off comparison |
| `OD_C1_GNSS_CONTINUOUS_BESTCASE` | best-case OD baseline | batch + EnKF | DET + MC + UT + POD | `42, 43, 44` | `pod_skip_slr=1` |
| `OD_C2_GNSS_REALISTIC_OUTAGES` | GNSS outage robustness | batch + EnKF | DET + MC + UT + POD | `42, 43, 44` | FOV/dropout/cycle-slip |
| `OD_C3_GNSS_PLUS_ACCEL` | GNSS+accel decoupling study | batch + EnKF | DET + MC + UT + POD | `42, 43, 44` | `SENS_ACCEL` scenario |

## 2) Phase B (stress and robustness)

| Scenario | Objective | Estimator | UQ Methods | Seed Set | Key Toggles |
|---|---|---|---|---|---|
| `ATT_A7_ECLIPSE_THERMAL_TRANSITION` | eclipse transition robustness | attitude filters | DET + MC + UT | `52, 53, 54` | thermal day/night wall model |
| `OD_C9_MEAS_ERROR_STRESS_TESTS` | outlier/slip stress | batch + EnKF | DET + MC + UT + POD | `52, 53, 54` | outlier + cycle-slip enabled |
| `OD_C8_STATION_WEATHER_AVAILABILITY` | weather-limited SLR | batch + EnKF | DET + MC + UT + POD | `52, 53, 54` | weather Markov SLR mask |
| `FORM_B6_ROBUST_FORMATION_STORM` | storm transient robustness | batch OD for analysis | DET + MC + UT | `52, 53, 54` | storm pulse + OU wind |

## 3) Phase C (formation mission studies)

| Scenario | Objective | Estimator | UQ Methods | Seed Set | Key Toggles |
|---|---|---|---|---|---|
| `FORM_B2_ALONGTRACK_2SAT_DD_KEEPING` | along-track differential-drag keeping | GNSS-driven OD | DET + MC + UT | `62, 63, 64` | OU density + OU wind |
| `FORM_B3_LVLH_BOX_CONSTRAINT_3SAT` | LVLH box maintenance | GNSS-driven OD | DET + MC + UT | `62, 63, 64` | OU density + OU wind |
| `FORM_B5_CLOSE_APPROACH_RISK` | close-approach probability | GNSS-driven OD | DET + MC + UT | `62, 63, 64` | sparse contact and uncertainty stress |
| `FORM_B7_CLUSTER_SPARSE_GROUND_CONTACT` | sparse contact operations | GNSS + SLR | DET + MC + UT + POD | `62, 63, 64` | SLR + sparse contact |

## 4) Required Outputs Per Scenario

- `results/<run_id>/<scenario>/summary.json`
- `results/<run_id>/<scenario>/mc_ut_stm.npz`
- `results/<run_id>/<scenario>/pod/summary.json` (when POD enabled)
- gate reports:
  - `maneuver_gate_summary.json` (where applicable)
  - `measurement_realism_gate.json` (POD scenarios)
  - `pod_harness_validation.json` (POD scenarios)

## 5) Minimal Execution Template

1. Build config from catalog (`scripts/catalog_to_case_config.py`) or use a fixed config JSON.
2. Run phase-1 DET/STM (`scripts/run_det_stm_phase.py`).
3. Run MC/UT shards (`scripts/run_mc_job_array.py`, `scripts/run_ut_job_array.py`) and merge.
4. Run postprocess + POD (`scripts/postprocess_uq_pod.py`).
5. Run validation gates and archive reports with run metadata.
