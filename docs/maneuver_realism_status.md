# Maneuver Realism Status (P1..P8)

Snapshot date: 2026-02-18

This document records current implementation status for maneuver realism items from `Plan_progression.md`.

## Status summary

- `P1`: complete
- `P2`: complete
- `P3`: complete
- `P4`: complete
- `P5`: complete
- `P6`: complete
- `P7`: complete
- `P8`: complete

## Evidence by point

### P1 Actuator realism limits

Status: complete

Evidence:

- Runtime extraction and defaults for limits in `scripts/run_mc_job_array.py`:
  - `eta_max_deg` (`110` default)
  - `eta_rate_max_deg_s` (`5` default)
  - `eta_accel_max_deg_s2` (`1` default)
- Config/env override behavior implemented in the same function.

### P2 Finite-rate actuator evolution

Status: complete

Evidence:

- Rate and acceleration-limited command tracking loop in `scripts/run_mc_job_array.py`:
  - target phase logic (`phase=0/1/2`),
  - bounded desired rate,
  - bounded rate-change (`eta_accel_max_deg_s2`),
  - bounded command (`eta_max_deg`).
- The same function stores commanded profile and rates in diagnostics.

### P3 Control semantics beyond waveform fallback

Status: complete

Evidence:

- Catalog conversion currently maps control toggle mainly to waveform parameter defaults (`wing_amp_deg`, `wing_period_s`) in `scripts/catalog_to_case_config.py`.
- `scripts/run_mc_job_array.py` now resolves control mode from catalog toggles (`CTL_ATT_AERO_RATE`, `CTL_ATT_AERO_POINT`) and applies explicit closed-loop maneuver command synthesis in `_propagate_member_with_axis_runtime(...)`:
  - `point` mode: attitude-pointing style command from angle/rate feedback,
  - `rate` mode: axis-rate tracking command from rate error.
- The runtime keeps finite actuator-rate/acceleration limits and writes control-mode diagnostics for auditability.

### P4 Maneuver robustness diagnostics

Status: complete

Evidence:

- `scripts/postprocess_mc_only.py` `build_maneuver_actuator_stats(...)` computes:
  - `control_switch_count`
  - command jump and jump-rate maxima
  - commanded-rate maxima
  - saturation fraction
  - half/full-turn reachability and time-to-full-turn
  - reached-fraction parsing support

### P5 Wing-sweep physical-domain guidance

Status: complete

Evidence:

- `scripts/summarize_att_wing_sweep.py` emits `in_recommended_range = (|wing_deg| <= 60)`.
- `scripts/build_att_wing_sweep_config.py` now enforces recommended range by default
  (`--recommended_abs_max_deg`, default `60`) and requires explicit opt-in
  (`--allow_out_of_recommended_range`) for out-of-range sweeps.
- `ATT_baseline.slurm` defaults sweep angles to `0..60 deg` and forwards the explicit
  override toggle (`ATT_ALLOW_OUT_OF_RECOMMENDED_RANGE=1`) only when intentionally requested.

### P6 Static authority map gate before time-domain sweeps

Status: complete

Evidence:

- Static authority-map pre-gate implemented in `scripts/validate_maneuver_realism_gate.py` (`--stage pre`):
  - computes axis torque map vs commanded wing angle at scenario initial state,
  - requires bidirectional authority and minimum torque floor.
- Integrated into ATT launch flow before time-domain sweeps in `run_all_A_scenarios.slurm`
  (`VG-4 pre-gate authority check` block).

### P7 Tighten maneuver-study initial attitude-rate spread

Status: complete

Evidence:

- Scenario-aware initial-rate covariance tightening in:
  - `scripts/mc_scenarios.py`
  - `scripts/ut_scenarios.py`
- Both use `initial_rate_sigma_deg_s` with default `0.5` for maneuver-axis scenarios.

### P8 Formal validation gates (pass/fail hard stops)

Status: complete

Evidence:

- Post-run pass/fail enforcement implemented in `scripts/validate_maneuver_realism_gate.py` (`--stage post`):
  - checks saturation fraction, command-rate limit adherence, and turn reachability from `summary.json`.
  - enforces that pre-gate authority result exists/passes (unless explicitly disabled).
- Integrated into `run_all_A_scenarios.slurm` in both branches:
  - MC-only flow (after `postprocess_mc_only.py`),
  - full flow (after phase-2 postprocess).

## Action items to close open/partial points

- Tune VG-4 thresholds using accumulated ATT_A8/A9/A10 campaign statistics.
