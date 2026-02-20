# Measurement Models

This note documents the implemented measurement models used by the current OD/POD and attitude workflows.

## 1) GNSS Raw Measurements

Implementation reference: `python/vleo_uq/pod_uq.py` (`simulate_gnss_measurements`, `batch_od_measurements`, `run_pod_uq_measurements*`).

Per visible satellite and epoch:

- Code:
  - `y_code = rho + b_rx + b_tropo + b_sat + eps_code`
- Carrier:
  - `y_carrier = rho + b_rx + b_tropo + b_sat + lambda*N + eps_carrier`

Where:

- `rho`: geometric range between spacecraft and GNSS satellite.
- `b_rx`: receiver clock bias state (random walk in truth generation and estimated nuisance state in OD).
- `b_tropo`: optional scalar nuisance term (kept optional and default-off in current runs).
- `b_sat`: optional satellite clock bias from SP3 clocks when enabled.
- `N`: integer ambiguity per carrier link/frequency.
- `eps_*`: zero-mean Gaussian measurement noise, configured by `code_sigma_m` and `carrier_sigma_m`.

### Implemented realism options

- LOS and Earth-occultation screening.
- FOV half-angle mask.
- Random dropout probability.
- Cycle slips (`cycle_slip_gap_s`, `cycle_slip_prob_per_min`).
- Code/carrier outlier injection (`*_outlier_prob`, `*_outlier_sigma_scale`).

### OD treatment

- Batch OD linearizes around propagated states and solves weighted least squares with nuisance states.
- EnKF path jointly updates dynamic and nuisance states.

## 2) SLR Measurements

Implementation reference: `python/vleo_uq/pod_uq.py` (`simulate_slr_measurements`, `slr_availability_mask`, `batch_od_measurements`).

Per station and epoch:

- `y_slr = rho_station + b_station + eps_slr`

Where:

- `rho_station`: geometric range spacecraft to station.
- `b_station`: station range bias (sampled truth bias and optionally estimated in OD).
- `eps_slr`: Gaussian white noise with configured `sigma_m`.

### Implemented availability realism

- Elevation and LOS gating.
- Optional weather Markov masking (`weather_enabled`, clear/block transition probabilities).

## 3) Attitude Sensor Measurements

Implementation reference: `python/vleo_uq/attitude_filter.py`.

- Gyro:
  - `y_gyro = omega_true + b_g + eps_g`
  - `b_g` modeled as RW/GM-style bias in truth generation and filter process model.
- Star tracker:
  - quaternion attitude samples with configurable cadence and angular noise.

Filter implementations include MEKF, UKF, and EnKF variants.

## 4) Accelerometer Status

`SENS_ACCEL` scenario toggles exist in catalog/config, but a full accelerometer measurement equation is not yet wired into the POD OD state update path in `python/vleo_uq/pod_uq.py`.

Current implication:

- GNSS/SLR realism is production-ready for OD stress scenarios.
- Accelerometer-assisted OD remains a documented gap for the next implementation block.

## 5) Explicitly Not Modeled (Current Baseline)

- Full ionosphere dual-frequency correction modeling in OD state (measurement generation uses noise/bias abstractions).
- Full light-time/relativistic corrections.
- Antenna phase-center offsets and detailed lever-arm corrections.

These omissions are intentional for current campaign throughput; they should be revisited before final publication-grade absolute orbit accuracy claims.
