# State and Coupling Definition

This document is the explicit state-vector and coupling reference for the current C++/Python runtime.

## 1) Conventions

- Translational frame: ECI (inertial) for `r_I`, `v_I`, force balance.
- Attitude quaternion: `q_BI = [w, x, y, z]`, mapping inertial vectors to body frame.
- Angular rate: `w_BI_B` in body frame.
- State size: 18.

Implementation references:

- `cpp/include/vleo/propagator.hpp`
- `cpp/src/propagator.cpp`
- `cpp/bindings/bindings.cpp`

## 2) State vector (18D)

- `x[0:3]`: `r_I` [m]
- `x[3:6]`: `v_I` [m/s]
- `x[6:10]`: `q_BI` [unit quaternion, wxyz]
- `x[10:13]`: `w_BI_B` [rad/s]
- `x[13]`: `log_rho_fast` [dimensionless]
- `x[14]`: `log_rho_bias` [dimensionless]
- `x[15:18]`: `wind_bias_I` [m/s]

Latent indexing in code:

- `kIdxLogRhoFast = 13`
- `kIdxLogRhoBias = 14`
- `kIdxWind = 15` (3 components)

## 3) Dynamics coupling

### 3.1 Translational dynamics

- `r_dot = v_I`
- `v_dot = a_grav + a_aero`
- `a_grav = -mu * r_I / ||r_I||^3`
- `a_aero = R_BI^T * F_B / m`

Where aerodynamic force/torque come from `AeroAdapter::computeFT(...)` using environment plus latent states.

### 3.2 Attitude dynamics

- `q_dot = quat_derivative_BI(q_BI, w_BI_B)`
- `w_dot = I^{-1} (tau_B - w x (I w))`

If `freeze_attitude=true`, angular propagation is disabled and evaluation uses `w=0` for force/torque calls.

### 3.3 Latent atmosphere/wind coupling

- Effective density:
  - `rho_eff = env.density * exp(log_rho_fast + log_rho_bias)`
- Effective inertial wind:
  - `wind_eff_I = env.wind_I + wind_bias_I`

These effective values are passed to aero force/torque evaluation.

## 4) Latent process model

Latent states are advanced by exact OU-style discrete updates per integration step:

- For `tau > 0`:
  - `x <- exp(-dt/tau) * x + sigma * sqrt(1 - exp(-2dt/tau)) * N(0,1)`
- For `tau <= 0` and `sigma > 0`:
  - random-walk fallback `x <- x + sigma * sqrt(dt) * N(0,1)`

Applies to:

- `log_rho_fast` via `rho_fast_tau_s`, `rho_fast_sigma`
- `log_rho_bias` via `rho_bias_tau_s`, `rho_bias_sigma`
- `wind_bias_I` components via `wind_tau_s`, `wind_sigma`

## 5) Error-state mapping used by UT/STM paths

The current code uses a reduced error-state dimension (`17`) with quaternion scalar exclusion for finite-difference/Jacobian mappings.

- Full-state index mapping array in code: `kErrToFull`.
- Scaling factors in code: `kErrScale`.

This mapping must remain synchronized with any future state-layout changes.

## 6) Measurement coupling summary

Current POD measurement models consume propagated trajectory states as follows:

- GNSS/SLR geometry depends primarily on spacecraft position (`x[0:3]`) and timing.
- RTN error/sigma metrics derive from propagated `r_I`, `v_I` basis.
- Additional POD nuisance states (clock, tropo, ambiguities, station bias) are estimator-side variables, not part of the 18D propagated dynamical state.

Key implementation:

- `python/vleo_uq/pod_uq.py`

## 7) Event coupling summary

Discrete events are applied through event utilities and can perturb:

- velocity (`delta_v`),
- attitude/rates (`attitude_reset`),
- and covariance handling around event boundaries.

Key implementation:

- `python/vleo_uq/events.py`
- `scripts/run_case_studies.py` (`propagate_all_with_events`, event-to-covariance logic)
