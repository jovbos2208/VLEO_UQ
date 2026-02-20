````markdown
# Scenario Catalog v2 — VLEO Aerodynamics, Attitude/Mission/Orbit Control & UQ Framework

This catalog is designed to produce **meaningful Uncertainty Quantification (UQ)** results across four research areas:

- **ATT** — Attitude control (aerodynamic panels, torque authority, robustness)
- **MIS** — Mission control (formation flying, operational constraints, contact gaps)
- **ORB** — Orbit control (drag compensation, station keeping, maneuver execution under uncertainty)
- **AERO** — Aerodynamics (GSI physics, geometry/model-form uncertainty, calibration/validation)

Each scenario is a **base configuration** intended to be converted into runnable JSON by your existing catalog→config tooling (e.g., `catalog_to_case_config.py`). Numeric values are defaults (sweepable).

---

## Global conventions

### Frames and time
- Inertial: **GCRS/ECI**
- Earth-fixed: **ITRS/ECEF**
- Body: **B**
- Relative motion: **LVLH/R** (Hill frame)
- Time input: **UTC**, converted using your existing leap-second + EOP machinery.

### Default orbit family (unless overridden)
- Near-circular VLEO
- `alt_km: 350`
- `e: 0.001`
- `i_deg: 51.6`

### Default spacecraft family (unless overridden)
- Mesh/panel geometry with shadowing
- Aerodynamic hinge rotations enabled if requested
- Mass/inertia fixed within a batch (unless explicitly swept)

### UQ philosophy
Every scenario specifies:
1) **QoIs** (quantities of interest) + decision metrics  
2) **Uncertainty sources** to include (environment, GSI, geometry, sensors, numerics)  
3) **Propagation method** (MC/UT/STM/MFMC) with recommended settings

---

## Composable toggles (orthogonal modifiers)

Implement these toggles as additive blocks so scenarios combine without duplication.

### ENV toggles (thermosphere + winds)
- **ENV_Q1_QUIET**: quiet indices; no stochastic density factor
- **ENV_Q2_SCALE_UP**: fixed density factor (e.g., +50%)
- **ENV_Q3_SCALE_DOWN**: fixed density factor (e.g., −30%)
- **ENV_S1_OU_LOGRHO**: OU process in log-density factor `delta_logrho`
  - params: `tau_rho_s`, `sigma_logrho`
- **ENV_S2_STORM_PULSE**: transient density pulse in log-density
  - params: `t0_s`, `duration_s`, `amplitude_logrho`, `shape: {step|gaussian}`
- **ENV_W1_WIND_DET**: deterministic winds only
- **ENV_W2_WIND_OU_BIAS**: OU wind-bias state `delta_w` (vector OU)
- **ENV_COMP_O_N2**: composition perturbation/sensitivity mode (if available)

### AERO/GSI toggles
- **GSI_M1_MAXWELL**
- **GSI_C1_CLL** (if implemented)
- **GSI_S1_SENTMAN** (if implemented)
- **GSI_TUTTAS** (if implemented)
- **GSI_P1_EST_ALPHAE** / **GSI_P2_FIXED_ALPHAE**
- **GSI_T1_FIXED_TW** / **GSI_T2_DAYNIGHT_TW**
- **AERO_DISC_ON**: model discrepancy term enabled (bias process on drag accel/torque)
- **GEO_SHADOW_ON** / **GEO_SHADOW_OFF**
- **GEO_MISALIGN_ON**: panel normal/hinge misalignment uncertainty
- **GEO_ROUGHNESS_ON**: aging/roughness surrogate (time-varying accommodation)

### CONTROL toggles
- **CTL_NONE**
- **CTL_ATT_AERO_RATE**: detumble / rate damping
- **CTL_ATT_AERO_POINT**: pointing (nadir/target)
- **CTL_ATT_DRAG_MIN**: pointing with drag-min objective (multi-objective)
- **CTL_DD_ALONGTRACK**: differential-drag formation keeping
- **CTL_DD_REPHASE**: differential-drag rephasing maneuver
- **CTL_ORB_DRAG_COMP**: orbit control by drag-area scheduling (no propellant)
- **CTL_ORB_THRUST_SMA**: semi-major-axis hold by thrust pulses (if propulsion modeled)

### SENSOR/EST toggles (for nav/OD-in-the-loop)
- Sensors:
  - **SENS_GNSS_RAW_DUAL**, **SENS_GNSS_RAW_SINGLE** (raw code+carrier; Doppler optional)
  - **SENS_SLR**
  - **SENS_ACCEL**
  - **SENS_ATT_ST**, **SENS_ATT_GYRO**
- Estimation:
  - **EST_EKF**, **EST_ENKF**, **EST_BLS**, **EST_RTS**
- Ground networks:
  - **GRD_GLB_20**, **GRD_REG_EU**, **GRD_ONE_MID**
  - **GRD_WEATHER_ON**, **GRD_WEATHER_OFF**

### UQ/Propagation toggles
- **PROP_DET**: deterministic run
- **PROP_STM**: sensitivity/STM propagation
- **PROP_UT**: unscented transform
- **PROP_MC**: Monte Carlo
- **PROP_MFMC**: multi-fidelity MC (if available)

Suggested defaults:
- UT: `alpha=1e-3, beta=2, kappa=0`
- MC: `n_mc=50..300` depending on tail-risk QoIs
- STM: ensure finite-diff steps and tolerances are documented and consistent

---

## Scenario schema

```yaml
scenario_id: <string>
area: <ATT|MIS|ORB|AERO>
purpose: <string>
duration_s: <float>
dt_control_s: <float|null>
dt_meas_s: <float|null>

orbit:
  type: <kepler|eci_state>
  alt_km: <float|null>
  e: <float>
  i_deg: <float>
  raan_deg: <float|null>
  argp_deg: <float|null>
  M_deg: <float|null>

spacecraft:
  mass_kg: <float>
  inertia_kgm2: [Ixx,Iyy,Izz]
  geometry_ref: <string>

aero:
  model: <string>
  toggles: [ ... ]

environment:
  model: <string>
  toggles: [ ... ]

sensors:
  toggles: [ ... ]

ground:
  toggles: [ ... ]

estimation:
  toggles: [ ... ]

control:
  toggles: [ ... ]

uq:
  propagator: <PROP_DET|PROP_MC|PROP_UT|PROP_STM|PROP_MFMC>
  n_mc: <int|null>
  ut_alpha: <float|null>
  ut_beta: <float|null>
  ut_kappa: <float|null>

outputs:
  qois: [ ... ]
  metrics: [ ... ]
````

---

# ATTITUDE CONTROL SCENARIOS (ATT)

## ATT_01_DETUMBLE_AERO_ONLY

```yaml
scenario_id: ATT_01_DETUMBLE_AERO_ONLY
area: ATT
purpose: "Detumble from high initial body rates using aerodynamic panels only; quantify time-to-detumble uncertainty."
duration_s: 21600
dt_control_s: 1.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GSI_T1_FIXED_TW, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]}

sensors: {toggles: [SENS_ATT_GYRO, SENS_ATT_ST]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ATT_AERO_RATE]}
uq: {propagator: PROP_MC, n_mc: 200}

outputs:
  qois: [time_to_detumble_s, max_body_rate, settling_time_s]
  metrics: [rms_body_rate, p95_body_rate, aero_torque_p95, saturation_time_fraction]
```

## ATT_02_NADIR_POINTING_ROBUST

```yaml
scenario_id: ATT_02_NADIR_POINTING_ROBUST
area: ATT
purpose: "Maintain nadir pointing under density/wind uncertainty; quantify pointing error distribution and drag penalty."
duration_s: 43200
dt_control_s: 2.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_SHADOW_ON]}
environment: {model: DTM2013, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_ATT_ST, SENS_ATT_GYRO]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ATT_AERO_POINT]}
uq: {propagator: PROP_UT, ut_alpha: 1e-3, ut_beta: 2.0, ut_kappa: 0.0}

outputs:
  qois: [pointing_error_deg, drag_accel_mps2, torque_cmd_Nm]
  metrics: [rms_pointing_error, p95_pointing_error, mean_drag_accel, drag_penalty_vs_uncontrolled]
```

## ATT_03_YAW_STEERING_WIND_BIAS

```yaml
scenario_id: ATT_03_YAW_STEERING_WIND_BIAS
area: ATT
purpose: "Yaw steering/beta-angle regulation under wind-bias uncertainty; quantify robustness and coupling to wind states."
duration_s: 21600
dt_control_s: 2.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_W2_WIND_OU_BIAS, ENV_S1_OU_LOGRHO]}

sensors: {toggles: [SENS_ATT_ST, SENS_ATT_GYRO]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ATT_AERO_POINT]}
uq: {propagator: PROP_MC, n_mc: 120}

outputs:
  qois: [yaw_error_deg, wind_bias_state, drag_accel_mps2]
  metrics: [rms_yaw_error, p95_yaw_error, wind_bias_rmse, stability_margin_proxy]
```

## ATT_04_AUTHORITY_FEASIBILITY_MAP

```yaml
scenario_id: ATT_04_AUTHORITY_FEASIBILITY_MAP
area: ATT
purpose: "Torque authority feasibility under low density, high inertia, and geometry misalignment; quantify failure probability."
duration_s: 14400
dt_control_s: 1.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 450, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.25, 0.22, 0.40], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON, GEO_MISALIGN_ON]}
environment: {model: DTM2013, toggles: [ENV_Q3_SCALE_DOWN, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_ATT_GYRO, SENS_ATT_ST]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ATT_AERO_POINT]}
uq: {propagator: PROP_MC, n_mc: 150}

outputs:
  qois: [feasible_flag, torque_margin, pointing_error_deg]
  metrics: [feasible_fraction, torque_margin_p05, saturation_time_fraction, failure_mode_counts]
```

## ATT_05_ECLIPSE_THERMAL_MODE_SWITCH

```yaml
scenario_id: ATT_05_ECLIPSE_THERMAL_MODE_SWITCH
area: ATT
purpose: "Robustness through eclipse/terminator transitions (time-varying Tw); quantify transient pointing excursions."
duration_s: 32400
dt_control_s: 2.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GSI_T2_DAYNIGHT_TW, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_ATT_ST, SENS_ATT_GYRO]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ATT_AERO_POINT]}
uq: {propagator: PROP_UT, ut_alpha: 1e-3, ut_beta: 2.0, ut_kappa: 0.0}

outputs:
  qois: [max_pointing_excursion_deg, torque_spike_Nm, drag_variation]
  metrics: [max_pointing_excursion_deg, settling_time_after_eclipse_s, torque_spike_p99, drag_variation_p95]
```

## ATT_06_DRAG_MIN_POINTING_TRADE

```yaml
scenario_id: ATT_06_DRAG_MIN_POINTING_TRADE
area: ATT
purpose: "Multi-objective control: pointing with drag minimization; quantify the trade-space under uncertain density and GSI."
duration_s: 43200
dt_control_s: 5.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE, GEO_SHADOW_ON]}
environment: {model: DTM2013, toggles: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]}

sensors: {toggles: [SENS_ATT_ST, SENS_ATT_GYRO]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ATT_DRAG_MIN]}
uq: {propagator: PROP_MC, n_mc: 120}

outputs:
  qois: [pointing_error_deg, mean_drag_accel, alphaE_posterior]
  metrics: [p95_pointing_error, mean_drag_accel, drag_reduction_vs_CTL_ATT_AERO_POINT, alphaE_posterior_std]
```

---

# MISSION CONTROL SCENARIOS (MIS)

## MIS_01_SINGLE_DRAG_MANAGEMENT

```yaml
scenario_id: MIS_01_SINGLE_DRAG_MANAGEMENT
area: MIS
purpose: "Single-satellite drag management (mission ops): regulate decay rate via attitude/area scheduling; quantify lifetime uncertainty."
duration_s: 172800
dt_control_s: 30.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: DTM2013, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ORB_DRAG_COMP]}
uq: {propagator: PROP_MC, n_mc: 100}

outputs:
  qois: [delta_altitude_m, da_dt, lifetime_proxy]
  metrics: [da_dt_p95, lifetime_proxy_p05, control_duty_cycle, estimator_consistency]
```

## MIS_02_DIFF_DRAG_KEEPING_2SAT

```yaml
scenario_id: MIS_02_DIFF_DRAG_KEEPING_2SAT
area: MIS
purpose: "Two-satellite along-track formation keeping via differential drag; quantify separation error and constraint risk."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_ENKF]}

control: {toggles: [CTL_DD_ALONGTRACK]}
uq: {propagator: PROP_MC, n_mc: 80}

outputs:
  qois: [alongtrack_sep_m, control_schedule]
  metrics: [sep_error_rms, sep_error_p95, duty_cycle, min_sep_margin]
```

## MIS_03_REPHASE_2SAT

```yaml
scenario_id: MIS_03_REPHASE_2SAT
area: MIS
purpose: "Rephase along-track separation by prescribed delta_s using differential drag; quantify time-to-target and overshoot probability."
duration_s: 432000
dt_control_s: 120.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: DTM2013, toggles: [ENV_S2_STORM_PULSE, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_BLS]}

control: {toggles: [CTL_DD_REPHASE]}
uq: {propagator: PROP_MC, n_mc: 60}

outputs:
  qois: [time_to_target_s, overshoot_flag]
  metrics: [time_to_target_p50, overshoot_prob, post_rephase_stability]
```

## MIS_04_LVLH_BOX_3SAT

```yaml
scenario_id: MIS_04_LVLH_BOX_3SAT
area: MIS
purpose: "Three-satellite LVLH box constraint keeping; quantify box-violation probability and minimum separation under uncertainty."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_ENKF]}

control: {toggles: [CTL_DD_ALONGTRACK]}
uq: {propagator: PROP_MC, n_mc: 80}

outputs:
  qois: [box_violation_flag, min_inter_sat_range_m]
  metrics: [box_violation_prob, lvhl_error_rms, min_inter_sat_range_p05]
```

## MIS_05_CLOSE_APPROACH_RISK

```yaml
scenario_id: MIS_05_CLOSE_APPROACH_RISK
area: MIS
purpose: "Close-approach risk in formation; quantify P(min_range < threshold) including tail sensitivity (requires higher MC)."
duration_s: 172800
dt_control_s: 60.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: JB2008, toggles: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_ENKF]}

control: {toggles: [CTL_DD_ALONGTRACK]}
uq: {propagator: PROP_MC, n_mc: 250}

outputs:
  qois: [min_inter_sat_range_m, close_approach_flag]
  metrics: [p_min_range_below_threshold, min_range_p01, estimator_divergence_prob, recovery_time_p90]
```

## MIS_06_CLUSTER_SPARSE_CONTACT

```yaml
scenario_id: MIS_06_CLUSTER_SPARSE_CONTACT
area: MIS
purpose: "4–8 satellite cluster with sparse contact; quantify navigation/control degradation vs contact gaps and network choice."
duration_s: 604800
dt_control_s: 120.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: DTM2013, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL, SENS_SLR]}
ground: {toggles: [GRD_REG_EU, GRD_WEATHER_ON]}   # swap GRD_GLB_20 for comparison
estimation: {toggles: [EST_EKF, EST_RTS]}

control: {toggles: [CTL_DD_ALONGTRACK]}
uq: {propagator: PROP_MC, n_mc: 100}

outputs:
  qois: [contact_gap_stats, sep_error_m, pos_rms_m]
  metrics: [sep_error_p95, pos_rms_p95, reconvergence_time, network_resilience_index]
```

## MIS_07_STORM_RESPONSE_FORMATION

```yaml
scenario_id: MIS_07_STORM_RESPONSE_FORMATION
area: MIS
purpose: "Formation robustness to storm-like density transient; quantify constraint violation and recovery probability."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: JB2008, toggles: [ENV_S2_STORM_PULSE, ENV_W2_WIND_OU_BIAS]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_ENKF]}

control: {toggles: [CTL_DD_ALONGTRACK]}
uq: {propagator: PROP_MC, n_mc: 120}

outputs:
  qois: [constraint_violation_flag, recovery_time_s]
  metrics: [constraint_violation_prob, recovery_prob, recovery_time_p90, sep_error_p95]
```

---

# ORBIT CONTROL SCENARIOS (ORB)

## ORB_01_SMA_HOLD_BY_DRAG

```yaml
scenario_id: ORB_01_SMA_HOLD_BY_DRAG
area: ORB
purpose: "Altitude-band holding using drag-area scheduling (no propellant); quantify station-keeping performance under environment/GSI uncertainty."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 300, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ORB_DRAG_COMP]}
uq: {propagator: PROP_MC, n_mc: 150}

outputs:
  qois: [altitude_error_m, da_dt, control_duty_cycle]
  metrics: [altitude_band_violation_prob, da_dt_p95, control_duty_cycle, expected_lifetime_gain]
```

## ORB_02_SMA_HOLD_BY_THRUST

```yaml
scenario_id: ORB_02_SMA_HOLD_BY_THRUST
area: ORB
purpose: "Altitude/SMA holding via thrust pulses; quantify performance under thrust magnitude/direction/timing and mass uncertainty."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 300, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ORB_THRUST_SMA]}
uq: {propagator: PROP_MC, n_mc: 150}

outputs:
  qois: [altitude_error_m, delta_v_used_mps, targeting_error]
  metrics: [altitude_band_violation_prob, delta_v_p95, targeting_error_stats, robustness_score]
```

## ORB_03_RECOVERY_AFTER_STORM

```yaml
scenario_id: ORB_03_RECOVERY_AFTER_STORM
area: ORB
purpose: "Storm-induced decay + recovery control; quantify control cost and recovery probability (nav + environment jointly uncertain)."
duration_s: 432000
dt_control_s: 120.0
dt_meas_s: 10.0

orbit: {type: kepler, alt_km: 320, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: JB2008, toggles: [ENV_S2_STORM_PULSE, ENV_W2_WIND_OU_BIAS]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL, SENS_SLR]}
ground: {toggles: [GRD_GLB_20, GRD_WEATHER_ON]}
estimation: {toggles: [EST_EKF, EST_RTS]}

control: {toggles: [CTL_ORB_DRAG_COMP]}    # also run thrust variant if available
uq: {propagator: PROP_MC, n_mc: 200}

outputs:
  qois: [recovery_time_s, altitude_recovered_flag, control_cost]
  metrics: [recovery_prob, recovery_time_p90, control_cost_p95, nav_gap_sensitivity]
```

## ORB_04_MANEUVER_EXECUTION_UNC

```yaml
scenario_id: ORB_04_MANEUVER_EXECUTION_UNC
area: ORB
purpose: "Single maneuver execution uncertainty (timing/pointing); quantify post-maneuver orbit state uncertainty."
duration_s: 21600
dt_control_s: 1.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: DTM2013, toggles: [ENV_Q1_QUIET, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ORB_THRUST_SMA]}
uq: {propagator: PROP_UT, ut_alpha: 1e-3, ut_beta: 2.0, ut_kappa: 0.0}

outputs:
  qois: [delta_a_km, delta_e, delta_i_deg]
  metrics: [post_maneuver_pos_sigma, post_maneuver_vel_sigma, targeting_error_stats]
```

## ORB_05_NAV_OUTAGE_STATIONKEEPING

```yaml
scenario_id: ORB_05_NAV_OUTAGE_STATIONKEEPING
area: ORB
purpose: "Station keeping under GNSS outages + sparse SLR; quantify control degradation driven by navigation uncertainty."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL, SENS_SLR]}
ground: {toggles: [GRD_REG_EU, GRD_WEATHER_ON]}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ORB_DRAG_COMP]}
uq: {propagator: PROP_MC, n_mc: 120}

outputs:
  qois: [altitude_error_m, nav_error, control_cost]
  metrics: [altitude_band_violation_prob, nav_outage_sensitivity, recovery_time_p90, control_cost_p95]
```

---

# AERODYNAMICS SCENARIOS (AERO)

## AERO_01_GSI_PARAMETER_CALIBRATION

```yaml
scenario_id: AERO_01_GSI_PARAMETER_CALIBRATION
area: AERO
purpose: "Calibrate GSI parameters (alpha_E, Tw) using GNSS+SLR(+accelerometer); quantify posterior uncertainty and identifiability."
duration_s: 172800
dt_control_s: null
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE, GSI_T2_DAYNIGHT_TW, GEO_SHADOW_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL, SENS_SLR, SENS_ACCEL]}
ground: {toggles: [GRD_GLB_20, GRD_WEATHER_ON]}
estimation: {toggles: [EST_EKF, EST_RTS]}

control: {toggles: [CTL_NONE]}
uq: {propagator: PROP_UT, ut_alpha: 1e-3, ut_beta: 2.0, ut_kappa: 0.0}

outputs:
  qois: [alphaE_posterior, Tw_posterior, logrho_posterior]
  metrics: [alphaE_posterior_std, Tw_posterior_std, parameter_correlation_matrix, residual_whiteness_tests]
```

## AERO_02_MODEL_FORM_ENSEMBLE

```yaml
scenario_id: AERO_02_MODEL_FORM_ENSEMBLE
area: AERO
purpose: "Model-form UQ using ensemble runs (Maxwell vs CLL/Sentman/Tuttas); quantify spread in Cd/Ct and orbit decay."
duration_s: 86400
dt_control_s: null
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GEO_SHADOW_ON]}   # implementer: run multiple variants with different GSI_* toggles
environment: {model: DTM2013, toggles: [ENV_Q1_QUIET, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_NONE]}
uq: {propagator: PROP_MC, n_mc: 80}

outputs:
  qois: [Cd_eff, Ct_eff, da_dt]
  metrics: [model_spread_Cd, model_spread_da_dt, model_rank_by_residual_if_available]
```

## AERO_03_GEOMETRY_AND_ALIGNMENT_UQ

```yaml
scenario_id: AERO_03_GEOMETRY_AND_ALIGNMENT_UQ
area: AERO
purpose: "Geometry/shadowing/misalignment UQ; quantify impact on force/torque and derived CdA."
duration_s: 43200
dt_control_s: 5.0
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_SHADOW_ON, GEO_MISALIGN_ON]}
environment: {model: NRLMSISE00, toggles: [ENV_Q1_QUIET, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_ATT_ST, SENS_ATT_GYRO]}
ground: {toggles: []}
estimation: {toggles: [EST_EKF]}

control: {toggles: [CTL_ATT_AERO_POINT]}
uq: {propagator: PROP_MC, n_mc: 150}

outputs:
  qois: [force_body, torque_body, CdA_eff]
  metrics: [CdA_std, torque_std, shadowing_sensitivity_index, misalignment_effect_size]
```

## AERO_04_DISCREPANCY_TERM_TRANSITIONAL

```yaml
scenario_id: AERO_04_DISCREPANCY_TERM_TRANSITIONAL
area: AERO
purpose: "Transitional-regime/model discrepancy channel: enable discrepancy term and test whether it absorbs structured residuals."
duration_s: 172800
dt_control_s: null
dt_meas_s: 1.0

orbit: {type: kepler, alt_km: 250, e: 0.001, i_deg: 51.6}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, AERO_DISC_ON, GEO_SHADOW_ON]}
environment: {model: JB2008, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_GNSS_RAW_DUAL, SENS_SLR]}
ground: {toggles: [GRD_GLB_20, GRD_WEATHER_ON]}
estimation: {toggles: [EST_EKF, EST_RTS]}

control: {toggles: [CTL_NONE]}
uq: {propagator: PROP_MC, n_mc: 120}

outputs:
  qois: [discrepancy_state, residual_structure]
  metrics: [residual_whiteness_tests, discrepancy_state_std, improvement_vs_no_discrepancy]
```

## AERO_05_SURFACE_AGING_AO_EROSION

```yaml
scenario_id: AERO_05_SURFACE_AGING_AO_EROSION
area: AERO
purpose: "Time-varying accommodation/roughness surrogate (AO erosion): quantify long-arc drift in CdA and OD impact."
duration_s: 604800
dt_control_s: null
dt_meas_s: 5.0

orbit: {type: kepler, alt_km: 300, e: 0.001, i_deg: 97.0}
spacecraft: {mass_kg: 12.0, inertia_kgm2: [0.15, 0.12, 0.20], geometry_ref: baseline_mesh_v1}

aero: {model: mesh_panel_aero, toggles: [GSI_M1_MAXWELL, GEO_ROUGHNESS_ON, GEO_SHADOW_ON]}
environment: {model: DTM2013, toggles: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]}

sensors: {toggles: [SENS_SLR]}
ground: {toggles: [GRD_REG_EU, GRD_WEATHER_ON]}
estimation: {toggles: [EST_BLS]}

control: {toggles: [CTL_NONE]}
uq: {propagator: PROP_MC, n_mc: 80}

outputs:
  qois: [CdA_drift_rate, drift_parameter_proxy]
  metrics: [CdA_drift_rate_p95, long_arc_stability, station_bias_confounding_score]
```

---

## Recommended bundles (thesis-ready)

### Minimal publishable set (balanced coverage)

* ATT: **ATT_01, ATT_02, ATT_05**
* MIS: **MIS_02, MIS_05, MIS_06**
* ORB: **ORB_01, ORB_03, ORB_05**
* AERO: **AERO_01, AERO_03, AERO_04**

### Stress-test bundle (tail risk + ops fragility)

* MIS_05 (higher MC), MIS_07
* ORB_03, ORB_05
* AERO_04 + ENV_S2_STORM_PULSE

### Calibration/validation bundle (data-informed)

* AERO_01 (GNSS+SLR+accel), AERO_04 (discrepancy)
* Always run residual diagnostics + whiteness tests

---

## Scenario index (stable IDs)

### ATT

* ATT_01_DETUMBLE_AERO_ONLY
* ATT_02_NADIR_POINTING_ROBUST
* ATT_03_YAW_STEERING_WIND_BIAS
* ATT_04_AUTHORITY_FEASIBILITY_MAP
* ATT_05_ECLIPSE_THERMAL_MODE_SWITCH
* ATT_06_DRAG_MIN_POINTING_TRADE

### MIS

* MIS_01_SINGLE_DRAG_MANAGEMENT
* MIS_02_DIFF_DRAG_KEEPING_2SAT
* MIS_03_REPHASE_2SAT
* MIS_04_LVLH_BOX_3SAT
* MIS_05_CLOSE_APPROACH_RISK
* MIS_06_CLUSTER_SPARSE_CONTACT
* MIS_07_STORM_RESPONSE_FORMATION

### ORB

* ORB_01_SMA_HOLD_BY_DRAG
* ORB_02_SMA_HOLD_BY_THRUST
* ORB_03_RECOVERY_AFTER_STORM
* ORB_04_MANEUVER_EXECUTION_UNC
* ORB_05_NAV_OUTAGE_STATIONKEEPING

### AERO

* AERO_01_GSI_PARAMETER_CALIBRATION
* AERO_02_MODEL_FORM_ENSEMBLE
* AERO_03_GEOMETRY_AND_ALIGNMENT_UQ
* AERO_04_DISCREPANCY_TERM_TRANSITIONAL
* AERO_05_SURFACE_AGING_AO_EROSION


