# Numerics

## Integrator

- Propagation uses adaptive high-order integration in `cpp/src/propagator.cpp` (`step_dop853` + step-doubling error control).
- Acceptance/rejection uses relative/absolute tolerances (`rtol`, `atol`) from `PropagatorConfig`.
- Quaternion states are renormalized after accepted steps.

## Stochastic Discretization

- Latent OU/random-walk channels are discretized in `apply_ou_update(...)`:
  - OU exact update with `phi = exp(-dt/tau)` and stationary variance term.
  - Random-walk fallback when `tau <= 0`.
- UT process-noise augmentation and STM process-noise `Q` use the same OU/RW variance model (`build_noise_spec(...)`).

## Validation Evidence

### OU statistics (VG-2)

- `tests/test_ou_discretization.py` validates:
  - mean and variance behavior vs theory,
  - lag-1 autocorrelation consistency,
  - random-walk fallback variance scaling,
  - seeded reproducibility.

### Integrator tolerance sweep (D-5 / VG-1)

- `scripts/validate_tolerance_sweep.py` runs a deterministic tolerance sweep for mission and attitude cases.
- `tests/test_tolerance_sweep.py` enforces that tighter tolerances reduce final state error versus a strict reference.

Run manually:

```bash
source /home/jovan/venv/bin/activate
python3 scripts/validate_tolerance_sweep.py --out /tmp/vleo_tolerance_sweep.json
```

Pass criterion:

- For mission and attitude cases, finest tolerance must have lower final position and velocity error than coarsest tolerance (and lower final attitude error for attitude case).

### UT/STM coverage realism vs MC (VG-3)

- `scripts/validate_coverage_short_arc.py` evaluates 1-sigma/2-sigma/3-sigma marginal coverage using MC final-state samples as reference.
- `tests/test_coverage_validation.py` enforces automated pass/fail for UT and STM coverage consistency.

Run manually:

```bash
source /home/jovan/venv/bin/activate
python3 scripts/validate_coverage_short_arc.py --out /tmp/vleo_coverage_short_arc.json
```
