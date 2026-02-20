# Model Discrepancy Channel

This document defines the runtime model-discrepancy channels implemented for atmospheric and aerodynamic mismatch studies.

## 1) Purpose

Parameter uncertainty (mass/inertia/density scale/OU params) is already represented.
These channels add explicit model-form uncertainty as stochastic corrections applied directly to environment inputs used by propagation:

- density discrepancy (`rho_model -> rho_eff`),
- composition discrepancy (`m_particle_model -> m_particle_eff`),
- storm jump regime on top of density,
- AO erosion surrogate on aerodynamic accommodation angles (`eta1`, `eta2`),
- plasma/charging surrogate on magnetic-field scaling,
- structural-flex surrogate on SRP scale and accommodation angles.
- rarefied-regime switching surrogate on `temperature_ratio_method`.

## 2) Mathematical Form

### 2.1 Density discrepancy (continuous OU/RW)

- `rho_eff(t) = rho_model(t) * s(t)`
- `s(t) = exp(x(t))`

`x(t)` is a zero-mean first-order process:

- OU form when `tau_s > 0`.
- Random-walk limit when `tau_s = 0`.

Innovation distribution:

- Student-t with `df` degrees of freedom (scaled to unit variance when `df > 2`).
- Gaussian fallback when `df <= 2`.

The resulting scale can be clipped to:

- `[1 - clip_rel, 1 + clip_rel]`

to prevent pathological excursions during large campaign runs.

### 2.2 Composition discrepancy (continuous OU/RW surrogate)

We optionally apply a parallel stochastic channel to effective particle mass used in aero evaluation:

- `m_eff(t) = m_model(t) * exp(x_c(t))`

with OU/RW dynamics and the same heavy-tail option as density.
This provides an explicit composition-uncertainty surrogate without adding a new propagated state dimension.

### 2.3 Storm jump regime (piecewise jumps)

We optionally superimpose a jump process on density log-scale:

- `rho_eff(t) = rho_eff,base(t) * exp(j(t))`

where `j(t)` is piecewise constant during storm windows:

- storm starts are sampled from a Poisson-rate model (`rate_per_day`),
- each event holds for `duration_s`,
- jump amplitude is sampled from a heavy-tailed distribution in relative-log space.

This is a practical disturbance-regime model for storm-like thermosphere mismatch.

### 2.4 AO erosion surrogate

We optionally apply a stochastic scale to aerodynamic accommodation angles:

- `eta_eff(t) = clip(eta_model(t) * s_ao(t), [0, \pi/2])`

where `s_ao(t)` combines an OU/RW stochastic component and an optional linear drift term
(`ao_erosion_rate_rel_per_day`) to emulate cumulative AO-driven surface property change over time.

### 2.5 Plasma/charging surrogate

We optionally apply a stochastic scale to magnetic field inputs used for residual dipole torque:

- `B_eff(t) = B_model(t) * s_plasma(t)`

with OU/RW + heavy-tail options analogous to density/composition channels.

### 2.6 Structural flex surrogate

We optionally apply a stochastic scale to `srp_scale` and a square-root scaled perturbation on
`eta1/eta2` to represent flexible-structure induced cross-sectional and aero-coupling variability.

### 2.7 Rarefied-regime switching surrogate

We optionally apply a stochastic switching process for aerodynamic model regime by toggling
`temperature_ratio_method` during active windows. This provides an explicit regime-switch
uncertainty representation without introducing a new propagated Knudsen state.

## 3) Runtime Controls

Implementation reference: `scripts/model_discrepancy.py`.

### 3.1 Density discrepancy keys

- `model_discrepancy_on` (bool)
- `model_discrepancy_sigma_rel` (float, relative sigma)
- `model_discrepancy_tau_s` (float)
- `model_discrepancy_df` (float)
- `model_discrepancy_clip_rel` (float)
- `model_discrepancy_seed_offset` (int)

Environment overrides:

- `VLEO_MODEL_DISCREPANCY_ON`
- `VLEO_MODEL_DISCREPANCY_SIGMA_REL`
- `VLEO_MODEL_DISCREPANCY_TAU_S`
- `VLEO_MODEL_DISCREPANCY_DF`
- `VLEO_MODEL_DISCREPANCY_CLIP_REL`
- `VLEO_MODEL_DISCREPANCY_SEED_OFFSET`

### 3.2 Composition discrepancy keys

- `composition_discrepancy_on` (bool)
- `composition_sigma_rel` (float)
- `composition_tau_s` (float)
- `composition_df` (float)
- `composition_clip_rel` (float)
- `composition_seed_offset` (int)

Environment overrides:

- `VLEO_COMPOSITION_DISCREPANCY_ON`
- `VLEO_COMPOSITION_DISCREPANCY_SIGMA_REL`
- `VLEO_COMPOSITION_DISCREPANCY_TAU_S`
- `VLEO_COMPOSITION_DISCREPANCY_DF`
- `VLEO_COMPOSITION_DISCREPANCY_CLIP_REL`
- `VLEO_COMPOSITION_DISCREPANCY_SEED_OFFSET`

### 3.3 Storm jump keys

- `storm_jump_on` (bool)
- `storm_jump_rate_per_day` (float)
- `storm_jump_duration_s` (float)
- `storm_jump_sigma_rel` (float)
- `storm_jump_mean_rel` (float)
- `storm_jump_df` (float)
- `storm_jump_clip_rel` (float)
- `storm_jump_seed_offset` (int)

Environment overrides:

- `VLEO_STORM_JUMP_ON`
- `VLEO_STORM_JUMP_RATE_PER_DAY`
- `VLEO_STORM_JUMP_DURATION_S`
- `VLEO_STORM_JUMP_SIGMA_REL`
- `VLEO_STORM_JUMP_MEAN_REL`
- `VLEO_STORM_JUMP_DF`
- `VLEO_STORM_JUMP_CLIP_REL`
- `VLEO_STORM_JUMP_SEED_OFFSET`

### 3.4 AO erosion surrogate keys

- `ao_erosion_on` (bool)
- `ao_erosion_sigma_rel` (float)
- `ao_erosion_tau_s` (float)
- `ao_erosion_df` (float)
- `ao_erosion_clip_rel` (float)
- `ao_erosion_rate_rel_per_day` (float)
- `ao_erosion_seed_offset` (int)

Environment overrides:

- `VLEO_AO_EROSION_ON`
- `VLEO_AO_EROSION_SIGMA_REL`
- `VLEO_AO_EROSION_TAU_S`
- `VLEO_AO_EROSION_DF`
- `VLEO_AO_EROSION_CLIP_REL`
- `VLEO_AO_EROSION_RATE_REL_PER_DAY`
- `VLEO_AO_EROSION_SEED_OFFSET`

### 3.5 Plasma/charging surrogate keys

- `plasma_charging_on` (bool)
- `plasma_charging_sigma_rel` (float)
- `plasma_charging_tau_s` (float)
- `plasma_charging_df` (float)
- `plasma_charging_clip_rel` (float)
- `plasma_charging_seed_offset` (int)

Environment overrides:

- `VLEO_PLASMA_CHARGING_ON`
- `VLEO_PLASMA_CHARGING_SIGMA_REL`
- `VLEO_PLASMA_CHARGING_TAU_S`
- `VLEO_PLASMA_CHARGING_DF`
- `VLEO_PLASMA_CHARGING_CLIP_REL`
- `VLEO_PLASMA_CHARGING_SEED_OFFSET`

### 3.6 Structural-flex surrogate keys

- `structural_flex_on` (bool)
- `structural_flex_sigma_rel` (float)
- `structural_flex_tau_s` (float)
- `structural_flex_df` (float)
- `structural_flex_clip_rel` (float)
- `structural_flex_seed_offset` (int)

Environment overrides:

- `VLEO_STRUCTURAL_FLEX_ON`
- `VLEO_STRUCTURAL_FLEX_SIGMA_REL`
- `VLEO_STRUCTURAL_FLEX_TAU_S`
- `VLEO_STRUCTURAL_FLEX_DF`
- `VLEO_STRUCTURAL_FLEX_CLIP_REL`
- `VLEO_STRUCTURAL_FLEX_SEED_OFFSET`

### 3.7 Rarefied-regime switch keys

- `rarefied_regime_switch_on` (bool)
- `rarefied_switch_rate_per_hour` (float)
- `rarefied_switch_mean_duration_s` (float)
- `rarefied_switch_alt_temperature_ratio_method` (int)
- `rarefied_switch_seed_offset` (int)

Environment overrides:

- `VLEO_RAREFIED_REGIME_SWITCH_ON`
- `VLEO_RAREFIED_SWITCH_RATE_PER_HOUR`
- `VLEO_RAREFIED_SWITCH_MEAN_DURATION_S`
- `VLEO_RAREFIED_SWITCH_ALT_TEMPERATURE_RATIO_METHOD`
- `VLEO_RAREFIED_SWITCH_SEED_OFFSET`

## 4) Seed and Reproducibility

The effective RNG seed is deterministically mixed from:

- run seed
- stable hash of scenario name/id
- per-channel seed offsets

so the same scenario + seed reproduces identical discrepancy series across runs.

## 5) Integration Points

The discrepancy channels are applied to environment entries in:

- `scripts/run_case_studies.py`
- `scripts/mc_scenarios.py`
- `scripts/ut_scenarios.py`
- `scripts/run_det_stm_phase.py`
- `scripts/postprocess_uq_pod.py`

Each case summary records discrepancy metadata including per-channel settings and applied scale statistics.

## 6) Current Limits

- Channels are currently shared per scenario run (not independent per MC particle).
- Calibration against external truth sources (e.g., accelerometer/DSMC residual structure) is still pending and should be handled in the next calibration sprint.
