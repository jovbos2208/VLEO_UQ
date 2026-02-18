# Scenario Catalog — VLEO UQ Satellite Aerodynamics Project

This file defines a **scenario catalog** for Codex to implement within the existing project.  
Each scenario is described in a structured, implementation-oriented way:

- **ID**: stable identifier (use as config `scenario_id`)
- **Area**: research area (ATT = attitude/aero control, FORM = formation/mission control, OD = orbit determination/ground)
- **Purpose**: what the run is meant to test
- **Core outputs**: metrics to log
- **Config schema**: a consistent block of fields Codex can map to your project’s config format

> Notes for implementers (Codex):
> - Treat each scenario as a *base config* that can be combined with environment / sensor / GSI toggles (see “Composable toggles”).
> - Where numeric values are given, they are **defaults**; implementers may sweep them.
> - “Optional” blocks can be omitted if unsupported; do not break the scenario ID list.

---

## Global conventions

### Time and frames
- Inertial frame: **GCRS/ECI**
- Earth-fixed frame: **ITRS/ECEF**
- Body frame: **B**
- LVLH for relative motion: **R** (Hill frame)
- Time scales: **UTC input**, convert internally using leap seconds + EOP as implemented.

### Default orbit family (unless overridden)
- Type: near-circular VLEO
- Altitude: **350 km**
- Inclination: **51.6°**
- e: **0.001**
- RAAN/argp/M: as needed (randomized per seed if unspecified)

### Default vehicle family (unless overridden)
- Mass: 5–20 kg (pick your platform; keep fixed per batch)
- Geometry: panelized surface model with shadowing option
- Actuation: aerodynamic panels/canted surfaces if scenario is ATT/FORM control

### Default environment model (unless overridden)
- Neutral density model: one of {NRLMSISE-00, DTM2013, JB2008}
- Winds: HWM14 or model used in code
- Uncertainty: multiplicative density factor (log-domain), optionally OU

---

## Composable toggles (apply to many scenarios)

Implement these as **orthogonal modifiers** so scenarios can be combined without duplicating definitions.

### ENV toggles (environment regimes)
- **ENV_Q1_QUIET**: quiet indices, no stochastic density factor
- **ENV_Q2_SCALE_UP**: fixed density scale factor (e.g., +50% in rho)
- **ENV_Q3_SCALE_DOWN**: fixed density scale factor (e.g., −30% in rho)
- **ENV_S1_OU_LOGRHO**: OU process in log-density factor  
  - params: `tau_rho` [s], `sigma_logrho` [-]
- **ENV_S2_STORM_PULSE**: transient density pulse (step or Gaussian bump) in log-density factor  
  - params: `t0`, `duration`, `amplitude_logrho`
- **ENV_W1_WIND_DET**: deterministic wind model only
- **ENV_W2_WIND_OU_BIAS**: OU wind-bias state added (vector OU)

### GSI toggles (gas-surface interaction / aero model)
- **GSI_M1_MAXWELL**: Maxwell model with accommodation parameters
- **GSI_C1_CLL**: CLL model (if implemented)
- **GSI_P1_EST_ALPHAE**: estimate energy accommodation `alpha_E`
- **GSI_P2_FIXED_ALPHAE**: fixed `alpha_E`
- **GSI_T1_FIXED_TW**: fixed wall temperature `T_w`
- **GSI_T2_DAYNIGHT_TW**: day/night (or sinusoidal) `T_w(t)` model
- **GEO_S1_SHADOWING_ON**: self-shadowing / visibility logic enabled
- **GEO_S2_SHADOWING_OFF**: no self-shadowing

### SENSOR toggles (measurement channels)
- **SENS_GNSS_RAW_DUAL**: raw GNSS code+carrier(+Doppler), dual-frequency
- **SENS_GNSS_RAW_SINGLE**: raw GNSS code+carrier, single-frequency
- **SENS_SLR**: SLR two-way range
- **SENS_ACCEL**: onboard accelerometer channel (non-grav specific force)
- **SENS_ATT_ST**: star tracker attitude measurement (quat)
- **SENS_ATT_GYRO**: gyro angular-rate measurement

### GROUND toggles (ground network)
- **GRD_GLB_20**: global ~20-station ILRS-like network
- **GRD_REG_EU**: Europe-only regional network
- **GRD_ONE_MID**: single mid-lat station
- **GRD_WEATHER_ON**: stochastic weather availability per station
- **GRD_WEATHER_OFF**: always-available (visibility mask only)

### EST toggles (estimation method)
- **EST_EKF**: EKF sequential
- **EST_ENKF**: ensemble Kalman filter
- **EST_BLS**: batch least squares / MAP per arc
- **EST_RTS**: RTS smoother post-processing (if available)

### CONTROL toggles
- **CTL_NONE**: no control (pure propagation)
- **CTL_ATT_AERO_RATE**: aero detumble / rate damping
- **CTL_ATT_AERO_POINT**: aero pointing controller (nadir/target)
- **CTL_DD_ALONGTRACK**: differential drag along-track formation keeping
- **CTL_DD_REPHASE**: differential drag rephasing maneuver

---

## Scenario definitions

Each scenario below uses this schema (fields can be mapped to your actual config keys):

```yaml
scenario_id: <string>
area: <ATT|FORM|OD>
purpose: <string>
duration_s: <float>
dt_control_s: <float|null>
dt_meas_s: <float|null>
orbit:
  type: <kepler|eci_state>
  a_km: <float|null>
  alt_km: <float|null>
  e: <float>
  i_deg: <float>
  raan_deg: <float|null>
  argp_deg: <float|null>
  M_deg: <float|null>
spacecraft:
  mass_kg: <float>
  inertia_kgm2: [Ixx,Iyy,Izz]        # or full tensor if supported
  panels: <vehicle geometry reference>
aero:
  model: <project model name>
  gsi_toggle: [ ... ]
environment:
  model: <NRLMSISE00|DTM2013|JB2008|...>
  env_toggle: [ ... ]
sensors:
  sensor_toggle: [ ... ]
ground:
  ground_toggle: [ ... ]
estimation:
  est_toggle: [ ... ]
control:
  control_toggle: <...>
uq:
  method: <MC|UT|MC+UT|...>
  n_mc: <int|null>
  ut_alpha: <float|null>
  ut_beta: <float|null>
  ut_kappa: <float|null>
outputs:
  metrics: [ ... ]
```

---

# ATTITUDE CONTROL SCENARIOS (ATT)

## ATT_A1_DETUMBLE_AERO_ONLY
- **Purpose:** detumble from high initial angular rates using aero torque authority only.
- **Recommended toggles:** `ENV_S1_OU_LOGRHO`, `ENV_W2_WIND_OU_BIAS`, `GEO_S1_SHADOWING_ON`, `GSI_P1_EST_ALPHAE`
- **Core outputs:** time-to-detumble, control effort, attitude error, torque margin.

```yaml
scenario_id: ATT_A1_DETUMBLE_AERO_ONLY
area: ATT
purpose: "Detumble from high initial rates using aerodynamic panels only."
duration_s: 21600         # 6 hours
dt_control_s: 1.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GSI_T1_FIXED_TW, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]
sensors:
  sensor_toggle: [SENS_ATT_GYRO, SENS_ATT_ST]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_RATE
uq:
  method: "MC"
  n_mc: 200
outputs:
  metrics: [time_to_detumble, max_body_rate, rms_body_rate, rms_pointing_error, aero_torque_rms, aero_torque_p95]
```

## ATT_A2_NADIR_POINTING
- **Purpose:** maintain nadir pointing within a bounded error using aero control.
- **Core outputs:** pointing error CDF, control effort distribution, drag penalty.

```yaml
scenario_id: ATT_A2_NADIR_POINTING
area: ATT
purpose: "Maintain nadir-pointing via aerodynamic panel control."
duration_s: 43200         # 12 hours
dt_control_s: 2.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_ST, SENS_ATT_GYRO]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
uq:
  method: "UT"
  ut_alpha: 1e-3
  ut_beta: 2.0
  ut_kappa: 0.0
outputs:
  metrics: [rms_pointing_error, p95_pointing_error, mean_drag, drag_penalty_vs_uncontrolled, control_effort]
```

## ATT_A3_YAW_STEERING_WIND_UNC
- **Purpose:** yaw steering / β-angle regulation under wind mismatch.
- **Core outputs:** yaw error statistics, robustness vs wind OU parameters.

```yaml
scenario_id: ATT_A3_YAW_STEERING_WIND_UNC
area: ATT
purpose: "Yaw steering (beta-angle) regulation under wind uncertainty."
duration_s: 21600
dt_control_s: 2.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_W2_WIND_OU_BIAS, ENV_S1_OU_LOGRHO]
sensors:
  sensor_toggle: [SENS_ATT_ST, SENS_ATT_GYRO]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
uq:
  method: "MC"
  n_mc: 100
outputs:
  metrics: [rms_yaw_error, p95_yaw_error, wind_bias_est_error, stability_margin_proxy]
```

## ATT_A4_SPIN_STAB_VS_3AXIS
- **Purpose:** compare passive spin-stabilized mode vs 3-axis aero control (stability vs drag penalty vs OD impact).
- **Core outputs:** stability metrics, drag increase, pointing performance.

```yaml
scenario_id: ATT_A4_SPIN_STAB_VS_3AXIS
area: ATT
purpose: "Compare passive spin stabilization vs 3-axis aerodynamic control."
duration_s: 43200
dt_control_s: 2.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_GYRO, SENS_ATT_ST]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
uq:
  method: "MC"
  n_mc: 100
outputs:
  metrics: [rms_pointing_error, max_pointing_error, mean_drag, drag_penalty_vs_uncontrolled, stability_mode_indicator]
```

## ATT_A5_CONTROL_AUTHORITY_LIMIT_MAP
- **Implementation note for Codex:** run a paired variant with `control_toggle: CTL_NONE` and an initial spin rate (e.g. 2–10 deg/s about max inertia axis) to represent “spin-stab”.
- **Purpose:** identify regimes where aero torque becomes insufficient (authority / feasibility map).
- **Core outputs:** feasible fraction, saturation time, failure modes.

```yaml
scenario_id: ATT_A5_CONTROL_AUTHORITY_LIMIT_MAP
area: ATT
purpose: "Map aerodynamic control authority limits across low-density/high-inertia cases."
duration_s: 14400
dt_control_s: 1.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 450, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.25, 0.22, 0.40]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_Q3_SCALE_DOWN, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_GYRO, SENS_ATT_ST]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [feasible_fraction, max_aero_torque, control_saturation_time_frac, p95_pointing_error, failure_mode_counts]
```


## ATT_A6_GSI_SENSITIVITY_CONTROL
- **Purpose:** quantify control authority sensitivity to accommodation / wall temperature / model form.
- **Implementation:** sweep `alpha_E`, `T_w`, and GSI model toggle.

```yaml
scenario_id: ATT_A6_GSI_SENSITIVITY_CONTROL
area: ATT
purpose: "Sensitivity of aero control authority to GSI parameters (alpha_E, T_w) and model selection."
duration_s: 14400
dt_control_s: 2.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE, GSI_T2_DAYNIGHT_TW, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_Q1_QUIET, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_ST, SENS_ATT_GYRO]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
uq:
  method: "MC"
  n_mc: 150
outputs:
  metrics: [torque_authority_mean, torque_authority_p05, pointing_error_cdf, alphaE_posterior_std, Tw_effect_size]
```

## ATT_A7_ECLIPSE_THERMAL_TRANSITION
- **Purpose:** robustness through eclipse/terminator thermal transitions (affecting reflected temperature).
- **Core outputs:** transient pointing excursions, controller stability, torque spikes.

```yaml
scenario_id: ATT_A7_ECLIPSE_THERMAL_TRANSITION
area: ATT
purpose: "Controller robustness through eclipse/terminator transitions via time-varying wall temperature."
duration_s: 32400
dt_control_s: 2.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 97.0}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GSI_T2_DAYNIGHT_TW, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_ST, SENS_ATT_GYRO]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
uq:
  method: "UT"
  ut_alpha: 1e-3
  ut_beta: 2.0
  ut_kappa: 0.0
outputs:
  metrics: [max_pointing_excursion, settling_time_after_eclipse, torque_spike_p99, drag_variation]
```

## ATT_A8_ROLL_ROTATION_BASELINE
- **Purpose:** baseline aerodynamic roll maneuver from zero attitude to half/full turn.
- **Core outputs:** time to 180°/360° roll, overshoot, peak roll rate, cross-axis coupling.

```yaml
scenario_id: ATT_A8_ROLL_ROTATION_BASELINE
area: ATT
purpose: "Baseline roll maneuver: start at 0 deg attitude and execute aerodynamic roll half/full turn."
duration_s: 5545          # ~1 orbit at 350 km
dt_control_s: 1.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_Q1_QUIET, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_ST, SENS_ATT_GYRO]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
maneuver:
  initial_attitude_euler_deg_zyx: [0.0, 0.0, 0.0]
  axis: roll
  half_turn_deg: 180.0
  full_turn_deg: 360.0
  wing_initial_deg: 8.0
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [time_to_roll_180_deg, time_to_roll_360_deg, overshoot_roll_360_deg, peak_roll_rate_deg_s, rms_pitch_excursion_deg, rms_yaw_excursion_deg]
```

## ATT_A9_PITCH_ROTATION_BASELINE
- **Purpose:** baseline aerodynamic pitch maneuver from zero attitude to half/full turn.
- **Core outputs:** time to 180°/360° pitch, overshoot, peak pitch rate, cross-axis coupling.

```yaml
scenario_id: ATT_A9_PITCH_ROTATION_BASELINE
area: ATT
purpose: "Baseline pitch maneuver: start at 0 deg attitude and execute aerodynamic pitch half/full turn."
duration_s: 5545          # ~1 orbit at 350 km
dt_control_s: 1.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_Q1_QUIET, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_ST, SENS_ATT_GYRO]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
maneuver:
  initial_attitude_euler_deg_zyx: [0.0, 0.0, 0.0]
  axis: pitch
  half_turn_deg: 180.0
  full_turn_deg: 360.0
  wing_initial_deg: 8.0
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [time_to_pitch_180_deg, time_to_pitch_360_deg, overshoot_pitch_360_deg, peak_pitch_rate_deg_s, rms_roll_excursion_deg, rms_yaw_excursion_deg]
```

## ATT_A10_YAW_ROTATION_BASELINE
- **Purpose:** baseline aerodynamic yaw maneuver from zero attitude to half/full turn.
- **Core outputs:** time to 180°/360° yaw, overshoot, peak yaw rate, cross-axis coupling.

```yaml
scenario_id: ATT_A10_YAW_ROTATION_BASELINE
area: ATT
purpose: "Baseline yaw maneuver: start at 0 deg attitude and execute aerodynamic yaw half/full turn."
duration_s: 5545          # ~1 orbit at 350 km
dt_control_s: 1.0
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_Q1_QUIET, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_ATT_ST, SENS_ATT_GYRO]
ground: {ground_toggle: []}
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
maneuver:
  initial_attitude_euler_deg_zyx: [0.0, 0.0, 0.0]
  axis: yaw
  half_turn_deg: 180.0
  full_turn_deg: 360.0
  wing_initial_deg: 8.0
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [time_to_yaw_180_deg, time_to_yaw_360_deg, overshoot_yaw_360_deg, peak_yaw_rate_deg_s, rms_roll_excursion_deg, rms_pitch_excursion_deg]
```

---

# FORMATION / MISSION CONTROL SCENARIOS (FORM)

## FORM_B1_SINGLE_DRAG_MANAGEMENT
- **Purpose:** regulate mean decay rate / altitude loss by modulating effective drag.
- **Core outputs:** altitude dispersion, mean semi-major axis decay, control schedule.

```yaml
scenario_id: FORM_B1_SINGLE_DRAG_MANAGEMENT
area: FORM
purpose: "Single-satellite drag management via scheduling effective CdA."
duration_s: 172800       # 48 hours
dt_control_s: 30.0
dt_meas_s: 10.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []      # GNSS space-based, no ground network needed
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_ATT_AERO_POINT
uq:
  method: "MC"
  n_mc: 100
outputs:
  metrics: [delta_altitude, mean_da_dt, cdA_schedule, estimation_consistency, rho_scale_est_error]
```

## FORM_B2_ALONGTRACK_2SAT_DD_KEEPING
- **Purpose:** 2-satellite along-track separation holding via differential drag.
- **Core outputs:** separation error, control duty cycle, robustness vs density OU.

```yaml
scenario_id: FORM_B2_ALONGTRACK_2SAT_DD_KEEPING
area: FORM
purpose: "Two-satellite along-track formation keeping using differential drag."
duration_s: 259200       # 72 hours
dt_control_s: 60.0
dt_meas_s: 10.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_ENKF]
control:
  control_toggle: CTL_DD_ALONGTRACK
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [sep_error_rms, sep_error_p95, reconfig_time, duty_cycle, collision_margin_min]
```




## FORM_B3_LVLH_BOX_CONSTRAINT_3SAT
- **Implementation note for Codex:** treat this as a **sweep scenario** over density scale (`ENV_Q3_SCALE_DOWN` magnitude), inertia scaling, and/or panel lever-arm scale.
- **Purpose:** constrain relative motion in LVLH (bounded “box”) for 3 satellites using differential drag + passive dynamics.
- **Core outputs:** box violation probability, min separation, reconfiguration events.

```yaml
scenario_id: FORM_B3_LVLH_BOX_CONSTRAINT_3SAT
area: FORM
purpose: "Maintain 3-satellite LVLH box constraints using differential drag."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 10.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_ENKF]
control:
  control_toggle: CTL_DD_ALONGTRACK
uq:
  method: "MC"
  n_mc: 60
outputs:
  metrics: [box_violation_prob, lvhl_rms_error, min_inter_sat_range, reconfig_event_count, duty_cycle]
```


## FORM_B4_REPHASE_DIFFERENTIAL_DRAG
- **Purpose:** re-phase along-track by a prescribed Δs without propellant.
- **Core outputs:** time-to-target, overshoot probability, density bias sensitivity.

```yaml
scenario_id: FORM_B4_REPHASE_DIFFERENTIAL_DRAG
area: FORM
purpose: "Re-phase along-track separation using differential drag (no propellant)."
duration_s: 432000       # 5 days
dt_control_s: 120.0
dt_meas_s: 10.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S2_STORM_PULSE, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_BLS]
control:
  control_toggle: CTL_DD_REPHASE
uq:
  method: "MC"
  n_mc: 60
outputs:
  metrics: [time_to_target_sep, overshoot_prob, fuel_equivalent_metric, post_rephase_stability]
```

## FORM_B5_CLOSE_APPROACH_RISK
- **Purpose:** quantify probability of violating minimum separation under worst-case environment uncertainty + measurement gaps.
- **Core outputs:** P(min range < threshold), divergence probability, recovery time.

```yaml
scenario_id: FORM_B5_CLOSE_APPROACH_RISK
area: FORM
purpose: "Close-approach risk under density/wind uncertainty and measurement gaps."
duration_s: 172800
dt_control_s: 60.0
dt_meas_s: 10.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "JB2008"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_ENKF]
control:
  control_toggle: CTL_DD_ALONGTRACK
uq:
  method: "MC"
  n_mc: 200
outputs:
  metrics: [p_min_range_below_threshold, min_range_distribution, close_approach_count, estimator_divergence_prob, recovery_time]
```

## FORM_B6_ROBUST_FORMATION_STORM
- **Purpose:** test robustness during storm-like density transients.
- **Core outputs:** probability of constraint violation, estimator divergence rate.

```yaml
scenario_id: FORM_B6_ROBUST_FORMATION_STORM
area: FORM
purpose: "Robust formation keeping under storm-like density transient and measurement gaps."
duration_s: 259200
dt_control_s: 60.0
dt_meas_s: 10.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "JB2008"
  env_toggle: [ENV_S2_STORM_PULSE, ENV_W2_WIND_OU_BIAS]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_ENKF]
control:
  control_toggle: CTL_DD_ALONGTRACK
uq:
  method: "MC"
  n_mc: 100
outputs:
  metrics: [constraint_violation_prob, estimator_divergence_prob, sep_error_p95, recovery_time]
```

## FORM_B7_CLUSTER_SPARSE_GROUND_CONTACT
- **Purpose:** cluster (4–8 sats) with sparse ground contact; determine minimal contact plan for stable OD + formation control.
- **Core outputs:** estimation error growth vs contact gaps, formation violation probability.

```yaml
scenario_id: FORM_B7_CLUSTER_SPARSE_GROUND_CONTACT
area: FORM
purpose: "4–8 sat cluster with sparse ground contact: minimal ground plan for stable OD + formation control."
duration_s: 604800
dt_control_s: 120.0
dt_meas_s: 10.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL, SENS_SLR]
ground:
  ground_toggle: [GRD_REG_EU, GRD_WEATHER_ON]
estimation:
  est_toggle: [EST_EKF, EST_RTS]
control:
  control_toggle: CTL_DD_ALONGTRACK
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [contact_gap_stats, estimation_error_growth, reconvergence_time, formation_violation_prob, station_plan_efficiency]
```


---

# ORBIT DETERMINATION & GROUND SEGMENT SCENARIOS (OD)

## OD_C1_GNSS_CONTINUOUS_BESTCASE
- **Purpose:** best-case OD with continuous raw GNSS observables.
- **Core outputs:** orbit RMS, drag parameter posteriors, residual stats.

```yaml
scenario_id: OD_C1_GNSS_CONTINUOUS_BESTCASE
area: OD
purpose: "Best-case OD with continuous raw GNSS code/carrier(+Doppler)."
duration_s: 86400
dt_control_s: null
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_NONE
uq:
  method: "UT"
  ut_alpha: 1e-3
  ut_beta: 2.0
  ut_kappa: 0.0
outputs:
  metrics: [pos_rms_eci, vel_rms_eci, cdA_posterior_std, alphaE_posterior_std, gnss_residual_stats]
```

## OD_C2_GNSS_REALISTIC_OUTAGES
- **Purpose:** OD performance under GNSS outages and antenna constraints.
- **Core outputs:** error growth vs outages, re-convergence behavior.

```yaml
scenario_id: OD_C2_GNSS_REALISTIC_OUTAGES
area: OD
purpose: "OD with raw GNSS under realistic outages (sky mask, eclipse, antenna constraints)."
duration_s: 86400
dt_control_s: null
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 97.0}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_NONE
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [pos_error_growth_rate, reconvergence_time, outage_sensitivity, filter_consistency]
```

## OD_C3_GNSS_PLUS_ACCEL
- **Purpose:** OD with raw GNSS augmented by onboard accelerometer to reduce drag/parameter ambiguity.
- **Core outputs:** posterior variance of density scale and GSI params, accel bias stability.

```yaml
scenario_id: OD_C3_GNSS_PLUS_ACCEL
area: OD
purpose: "OD + parameter estimation with raw GNSS + accelerometer."
duration_s: 86400
dt_control_s: null
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL, SENS_ACCEL]
ground:
  ground_toggle: []
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_NONE
uq:
  method: "UT"
  ut_alpha: 1e-3
  ut_beta: 2.0
  ut_kappa: 0.0
outputs:
  metrics: [pos_rms_eci, logrho_state_std, alphaE_posterior_std, accel_bias_est, drag_decomposition_quality]
```

## OD_C4_SLR_ONLY_SPARSE
- **Purpose:** OD with SLR-only sparse passes (precision but gaps).
- **Core outputs:** between-pass divergence, bias sensitivity, station geometry dependence.

```yaml
scenario_id: OD_C4_SLR_ONLY_SPARSE
area: OD
purpose: "SLR-only OD with sparse passes; quantify divergence between passes and bias sensitivity."
duration_s: 604800        # 7 days
dt_control_s: null
dt_meas_s: 5.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_SLR]
ground:
  ground_toggle: [GRD_GLB_20, GRD_WEATHER_ON]
estimation:
  est_toggle: [EST_BLS]
control:
  control_toggle: CTL_NONE
uq:
  method: "MC"
  n_mc: 50
outputs:
  metrics: [pos_rms_eci, between_pass_divergence, station_bias_est, slr_residual_stats]
```

## OD_C5_GNSS_SLR_HYBRID_ANCHOR
- **Purpose:** hybrid OD with raw GNSS + SLR anchor to control drifts/biases.
- **Core outputs:** improvement vs GNSS-only in parameter posteriors and long-arc stability.

```yaml
scenario_id: OD_C5_GNSS_SLR_HYBRID_ANCHOR
area: OD
purpose: "Hybrid OD: raw GNSS + SLR anchor for drift/bias control and robustness."
duration_s: 172800
dt_control_s: null
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W2_WIND_OU_BIAS]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL, SENS_SLR]
ground:
  ground_toggle: [GRD_GLB_20, GRD_WEATHER_ON]
estimation:
  est_toggle: [EST_EKF, EST_RTS]
control:
  control_toggle: CTL_NONE
uq:
  method: "UT"
  ut_alpha: 1e-3
  ut_beta: 2.0
  ut_kappa: 0.0
outputs:
  metrics: [pos_rms_eci, alphaE_posterior_std, logrho_state_std, long_arc_stability, residual_whiteness_tests]
```

## OD_C6_SLR_SINGLE_STATION_WORSTCASE
- **Purpose:** what OD is possible with a single station + elevation mask.
- **Core outputs:** observability degradation, along-track uncertainty growth.

```yaml
scenario_id: OD_C6_SLR_SINGLE_STATION_WORSTCASE
area: OD
purpose: "Worst-case: SLR-only with a single mid-lat station; quantify observability limits."
duration_s: 604800
dt_control_s: null
dt_meas_s: 5.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO]
sensors:
  sensor_toggle: [SENS_SLR]
ground:
  ground_toggle: [GRD_ONE_MID, GRD_WEATHER_OFF]
estimation:
  est_toggle: [EST_BLS]
control:
  control_toggle: CTL_NONE
uq:
  method: "MC"
  n_mc: 40
outputs:
  metrics: [along_track_sigma, cross_track_sigma, radial_sigma, between_pass_divergence, station_geometry_sensitivity]
```

## OD_C7_NETWORK_COMPARISON_REGIONAL_VS_GLOBAL
- **Purpose:** compare OD performance across different ground station networks.
- **Implementation:** run the same scenario with `GRD_REG_EU` vs `GRD_GLB_20`.
- **Core outputs:** OD accuracy and parameter posteriors vs network.

```yaml
scenario_id: OD_C7_NETWORK_COMPARISON_REGIONAL_VS_GLOBAL
area: OD
purpose: "Compare OD/parameter estimation with regional vs global ground station networks."
duration_s: 259200
dt_control_s: null
dt_meas_s: 5.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_SLR]
ground:
  ground_toggle: [GRD_REG_EU, GRD_WEATHER_ON]   # implementer: swap GRD_GLB_20 for comparison
estimation:
  est_toggle: [EST_BLS]
control:
  control_toggle: CTL_NONE
uq:
  method: "MC"
  n_mc: 40
outputs:
  metrics: [pos_rms_eci, alphaE_posterior_std, network_effect_size, pass_count, residual_stats]
```

## OD_C8_STATION_WEATHER_AVAILABILITY
- **Purpose:** OD robustness under weather-limited station availability; quantify expected error under stochastic visibility losses.
- **Core outputs:** expected orbit error, pass-loss statistics, resilience index.

```yaml
scenario_id: OD_C8_STATION_WEATHER_AVAILABILITY
area: OD
purpose: "OD robustness under stochastic station weather availability (SLR)."
duration_s: 604800
dt_control_s: null
dt_meas_s: 5.0
orbit: {type: kepler, alt_km: 400, e: 0.001, i_deg: 51.6}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P2_FIXED_ALPHAE]
environment:
  model: "DTM2013"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_SLR]
ground:
  ground_toggle: [GRD_GLB_20, GRD_WEATHER_ON]
estimation:
  est_toggle: [EST_BLS]
control:
  control_toggle: CTL_NONE
uq:
  method: "MC"
  n_mc: 80
outputs:
  metrics: [expected_pos_rms, percentile_pos_rms, pass_loss_stats, between_pass_divergence, network_resilience_index]
```

## OD_C9_MEAS_ERROR_STRESS_TESTS
- **Purpose:** stress-test measurement models (GNSS multipath/outliers, cycle slips, SLR timing biases).
- **Core outputs:** robustness metrics, false-event rates, residual tail behavior.

```yaml
scenario_id: OD_C9_MEAS_ERROR_STRESS_TESTS
area: OD
purpose: "Stress-test GNSS/SLR measurement errors: outliers, cycle slips, timing biases."
duration_s: 172800
dt_control_s: null
dt_meas_s: 1.0
orbit: {type: kepler, alt_km: 350, e: 0.001, i_deg: 97.0}
spacecraft:
  mass_kg: 12.0
  inertia_kgm2: [0.15, 0.12, 0.20]
  panels: "baseline_panels_v1"
aero:
  model: "panel_gsi"
  gsi_toggle: [GSI_M1_MAXWELL, GSI_P1_EST_ALPHAE, GEO_S1_SHADOWING_ON]
environment:
  model: "NRLMSISE00"
  env_toggle: [ENV_S1_OU_LOGRHO, ENV_W1_WIND_DET]
sensors:
  sensor_toggle: [SENS_GNSS_RAW_DUAL, SENS_SLR]
ground:
  ground_toggle: [GRD_REG_EU, GRD_WEATHER_ON]
estimation:
  est_toggle: [EST_EKF]
control:
  control_toggle: CTL_NONE
uq:
  method: "MC"
  n_mc: 120
outputs:
  metrics: [outlier_rejection_rate, cycle_slip_detection_rate, bias_est_drift, false_density_event_rate, residual_heavy_tail_score]
```

---

## Scenario index (stable IDs)

### ATT
- ATT_A1_DETUMBLE_AERO_ONLY
- ATT_A2_NADIR_POINTING
- ATT_A3_YAW_STEERING_WIND_UNC
- ATT_A4_SPIN_STAB_VS_3AXIS
- ATT_A5_CONTROL_AUTHORITY_LIMIT_MAP
- ATT_A6_GSI_SENSITIVITY_CONTROL
- ATT_A7_ECLIPSE_THERMAL_TRANSITION
- ATT_A8_ROLL_ROTATION_BASELINE
- ATT_A9_PITCH_ROTATION_BASELINE
- ATT_A10_YAW_ROTATION_BASELINE

### FORM
- FORM_B1_SINGLE_DRAG_MANAGEMENT
- FORM_B2_ALONGTRACK_2SAT_DD_KEEPING
- FORM_B3_LVLH_BOX_CONSTRAINT_3SAT
- FORM_B4_REPHASE_DIFFERENTIAL_DRAG
- FORM_B5_CLOSE_APPROACH_RISK
- FORM_B6_ROBUST_FORMATION_STORM
- FORM_B7_CLUSTER_SPARSE_GROUND_CONTACT

### OD
- OD_C1_GNSS_CONTINUOUS_BESTCASE
- OD_C2_GNSS_REALISTIC_OUTAGES
- OD_C3_GNSS_PLUS_ACCEL
- OD_C4_SLR_ONLY_SPARSE
- OD_C5_GNSS_SLR_HYBRID_ANCHOR
- OD_C6_SLR_SINGLE_STATION_WORSTCASE
- OD_C7_NETWORK_COMPARISON_REGIONAL_VS_GLOBAL
- OD_C8_STATION_WEATHER_AVAILABILITY
- OD_C9_MEAS_ERROR_STRESS_TESTS

---

## Minimal “selected set” for thesis runs (recommended)

If you need a compact set with strong coverage:

- **ATT:** ATT_A1, ATT_A2, ATT_A6, ATT_A7  
- **FORM:** FORM_B2, FORM_B4, FORM_B6  
- **OD:** OD_C2, OD_C4, OD_C5, OD_C7  

Pair each with environment toggles: `ENV_Q1_QUIET`, `ENV_S1_OU_LOGRHO`, `ENV_S2_STORM_PULSE`.
