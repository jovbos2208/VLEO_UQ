Below is a concrete Python + C++ (pybind11) “Codex plan” tailored to your setup:

* You already have **`fore_moment_module` (C++)** that computes forces and torques from an object model.
* You want **highest precision** and **comparisons between propagation/UQ methods**.
* You want **POD UQ**, **attitude UQ**, **payload-resolution impact**, and **source-by-source sensitivity**.

I’ll structure this as (A) architecture, (B) C++ interfaces to wrap your force/torque module cleanly, (C) binding/build system, (D) propagation methods implementation plan, (E) UQ sources representation, (F) POD + attitude + payload impact workflow, (G) sensitivity analysis workflow, (H) validation/performance checklist, plus skeleton code.

---

## A) Architecture: what goes in C++ vs Python (high precision + high throughput)

### Put in C++ (performance-critical)

1. **Numerical integration** (orbit + attitude + latent stochastic states).
2. **Batch propagation** for ensembles / sigma points (vectorized loops, OpenMP).
3. Optional: **variational equations / STM** and numerical Jacobians (expensive but contained).
4. Random process stepping (OU / random walk) with reproducible RNG.

### Keep in Python (orchestration + analysis)

1. Scenario configuration (YAML/JSON), experiment grids.
2. UQ source registry (distributions, grouping, toggles).
3. Sensitivity (Morris/Sobol) bookkeeping, plotting, reporting.
4. Measurement simulation + OD pipeline control (but heavy computations can be delegated to C++ if needed).

The key design point: **C++ should accept “environment drivers” as precomputed time series** (F10.7, Ap, density scale factors, wind perturbations, eclipse flags), so you avoid Python callbacks in the inner propagation loop.

---

## B) C++ interfaces: wrap `fore_moment_module` into a clean “model” ABI

### B1) State layout (single-satellite)

Use a strict, contiguous layout; keep it consistent for all propagators:

* Orbit: `r(3), v(3)` in ECI
* Attitude: quaternion `q(4)` (scalar-last or scalar-first; pick one and enforce), body rates `w(3)`
* Latent processes (as states):

  * `delta_log_rho_fast` (OU)
  * `delta_log_rho_bias` (random walk or constant parameter)
  * `delta_wind(3)` (OU)
  * `beta_drag` (optional combined scale)
  * sensor biases if doing attitude/POD filters

Define:

```cpp
constexpr int NX = 6 + 7 + N_LATENT; // r,v + q,w + latent
```

### B2) Force/torque evaluation wrapper

Create a small adapter around `fore_moment_module` so the propagator calls exactly one function:

```cpp
struct ForceMoment {
  Eigen::Vector3d accel_eci;   // m/s^2
  Eigen::Vector3d torque_body; // N*m (body frame)
};

struct EnvInputs {
  double f107;
  double ap;
  Eigen::Vector3d wind_eci;     // or wind in ECEF/ECI; be explicit
  double srp_scale;             // optional
  bool in_eclipse;
  // Optional: density model outputs if you precompute them externally
  double rho_model;             // kg/m^3 at sat location/time (if provided)
};

struct VehicleParams {
  double mass;
  Eigen::Matrix3d inertia;
  // Any surface/optical/aero parameters that fore_moment_module needs
};

class IForceMomentModel {
public:
  virtual ~IForceMomentModel() = default;
  virtual ForceMoment evaluate(double t, const Eigen::VectorXd& x,
                               const EnvInputs& env,
                               const VehicleParams& p) const = 0;
};
```

Then implement:

```cpp
class ForeMomentAdapter final : public IForceMomentModel {
public:
  explicit ForeMomentAdapter(std::shared_ptr<fore_moment_module::Model> model);
  ForceMoment evaluate(double t, const Eigen::VectorXd& x,
                       const EnvInputs& env,
                       const VehicleParams& p) const override;
private:
  std::shared_ptr<fore_moment_module::Model> model_;
};
```

**Thread safety requirement:** ensure `fore_moment_module::Model` is either:

* immutable during `evaluate()`, or
* cloneable per thread (recommended for OpenMP).

If it uses caches, make them thread-local or lock-free per-thread copies.

---

## C) Python bindings + build system (pybind11 + scikit-build-core)

### C1) Binding tool choice

Use **pybind11**. It is the standard for high-performance Python/C++ bindings and works cleanly with Eigen.

### C2) Packaging/build

Use:

* **CMake**
* **scikit-build-core** (PEP517) for pip builds
* Build `fore_moment_module` as a CMake target and link it in.

Recommended repository layout:

```
vleo_uq/
  cpp/
    CMakeLists.txt
    include/vleo/...
    src/...
    bindings/bindings.cpp
  python/
    vleo_uq/__init__.py
    vleo_uq/scenarios/...
  pyproject.toml
```

### C3) Minimal `pyproject.toml`

```toml
[build-system]
requires = ["scikit-build-core>=0.9", "pybind11>=2.11", "numpy"]
build-backend = "scikit_build_core.build"

[project]
name = "vleo-uq"
version = "0.1.0"
requires-python = ">=3.10"
dependencies = ["numpy", "scipy", "pydantic", "xarray", "h5py", "SALib"]
```

### C4) CMake sketch

```cmake
cmake_minimum_required(VERSION 3.21)
project(vleo_uq LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

# fore_moment_module as subproject or external
add_subdirectory(fore_moment_module)

add_library(vleo_core
  src/propagator.cpp
  src/integrators.cpp
  src/models.cpp
  src/processes.cpp
  src/variational.cpp
)
target_include_directories(vleo_core PUBLIC include)
target_link_libraries(vleo_core PUBLIC Eigen3::Eigen fore_moment_module_target)

pybind11_add_module(_vleo_uq bindings/bindings.cpp)
target_link_libraries(_vleo_uq PRIVATE vleo_core)
```

---

## D) Propagation methods to implement (highest precision + comparison)

You said you want “highest precision” and have compute. Therefore:

### D1) Truth baseline

* **Special perturbations (SP)** propagation of augmented state.
* Integrator: **DOP853** or high-order RK8(7) with tight tolerances.
* Discrete events: maneuvers, mode switches, eclipse transitions.
* Latent stochastic processes integrated as states (OU) or via exact discrete updates.

### D2) Methods to compare

1. **MC ensemble** (benchmark)
2. **UT / sigma points** (augmented state)
3. **STM/covariance** (variational equations + process noise)
4. Optional: **PCE / sparse collocation**
5. Optional: **DA** (if you add a DA library later)

In code, make this a clean interface:

```cpp
class IPropagator {
public:
  virtual ~IPropagator() = default;
  virtual void propagate(
    const Eigen::MatrixXd& X0,      // (N, NX) ensemble; N=1 for deterministic
    const Eigen::VectorXd& t_grid,  // (Nt,)
    const std::vector<EnvInputs>& env_series,
    Eigen::MatrixXd* X_out          // (N*Nt, NX) or (N,Nt,NX)
  ) const = 0;
};
```

Then implement:

* `DeterministicPropagator` (N=1)
* `EnsemblePropagatorMC`
* `SigmaPointPropagatorUT`
* `STMCovPropagator` (outputs `P(t)` not samples)

---

## E) UQ sources: representation and how they enter propagation (code-friendly)

Make each uncertainty source one of:

### E1) Uncertain initial state

* `x0 ~ N(mu0, P0)` or mixture
* Implement sampling in C++ for speed (NumPy -> C++ also fine).

### E2) Uncertain parameters (constant per arc)

Examples:

* mass, inertia, center-of-pressure offset
* optical coefficients for SRP
* GSI/accommodation parameters
* gravity parameter perturbations (if included)

Represent as a “param vector” `theta` sampled per run.

### E3) Slowly varying parameters

Examples:

* Cd drift, accommodation drift, mass depletion bias
  Model as random walk / long-time GM process.

### E4) Fast stochastic processes (dominant in VLEO)

Examples:

* log-density fluctuation OU
* wind OU
* residual acceleration OU (RTN)
* residual torque OU (body frame)

These should be **latent states** so propagation methods treat them consistently.

### E5) Discrete random events

Examples:

* ΔV execution errors (scale/direction/time)
* attitude mode transitions time/state
* eclipse boundary uncertainty (optional)

Events are a list of `Event{t0, type, params}` with uncertainties.

---

## F) POD UQ + attitude UQ + payload resolution: workflow design

### F1) POD UQ (requires OD loop)

You want “relatable POD UQ numbers.” That means:

1. Generate truth trajectory (from one MC draw).
2. Simulate measurements (GNSS code/carrier, SLR, etc.).
3. Run OD (batch least squares / smoothing) to recover estimated states + covariance.
4. Repeat over many runs; report:

   * posterior RTN sigmas
   * RMS error vs covariance (coverage)

Implementation plan:

* Keep measurement simulation in Python initially, but make the measurement model C++ callable for speed later.
* Implement batch OD in Python (SciPy least squares) or in C++ (Ceres) depending on complexity.

Given your compute and precision goals, I recommend:

* **Batch smoother** (multi-arc) for POD metrics
* Use your high-fidelity propagation as the dynamics in the estimator.

### F2) Attitude UQ

Two levels:

* Fast: commanded attitude + stochastic error-state OU + optional jitter PSD model.
* Full: propagate attitude dynamics + sensor models + MEKF/UKF.

Since you want “highest precision,” implement the full chain eventually, but start with the fast one for sensitivity sweeps.

### F3) Payload impact (resolution degradation)

Implement a module that maps orbit+attitude uncertainty to measurement error:

Examples you likely want:

* ground geolocation error from attitude + orbit (nadir imager)
* smear over exposure time from pointing jitter
* for altimetry: mispointing-to-range bias
* for SAR: simplified pointing-to-Doppler centroid shift (proxy)

These should be pure Python at first (fast enough), fed by MC outputs.

---

## G) Sensitivity study: many sources, large compute

Do it hierarchically:

1. Group-level global sensitivity (Atmosphere vs Aero vs Attitude vs Radiation vs Maneuvers vs Gravity/Frames vs Measurement)
2. Expand the dominant groups into sub-parameters.

Implement:

* Morris screening to reduce dimension
* Sobol indices (Saltelli) for final quantification

To make this tractable with time series outputs:

* Define scalar outputs at key horizons (e.g., RTN σ at 6h/24h/3d; payload RMS at 10s/1s/0.1s integration windows)
* Also compute time-dependent indices for a small set of outputs (optional)

---

## H) Validation + performance checklist (so results are publishable)

1. **Integrator convergence**: tighten tolerances until statistics stop changing.
2. **MC convergence**: show variance stabilization vs N.
3. **Thread safety**: no shared mutable state in force/torque model.
4. **No-double-counting policies** enforced via configuration:

   * either density-centric or drag-scale-centric experiments
5. **Covariance realism** metrics for STM/UT vs MC:

   * coverage fractions at 1σ/2σ/3σ
6. **Ablation runs**: enable one source group at a time; confirm dominance logic.

---

## Concrete binding skeleton (pybind11)

### `bindings/bindings.cpp` (sketch)

```cpp
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "vleo/propagators.hpp"
#include "vleo/models.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_vleo_uq, m) {
  m.doc() = "VLEO UQ propagation core";

  py::class_<vleo::EnvInputs>(m, "EnvInputs")
    .def(py::init<>())
    .def_readwrite("f107", &vleo::EnvInputs::f107)
    .def_readwrite("ap", &vleo::EnvInputs::ap)
    .def_readwrite("wind_eci", &vleo::EnvInputs::wind_eci)
    .def_readwrite("srp_scale", &vleo::EnvInputs::srp_scale)
    .def_readwrite("in_eclipse", &vleo::EnvInputs::in_eclipse)
    .def_readwrite("rho_model", &vleo::EnvInputs::rho_model);

  py::class_<vleo::VehicleParams>(m, "VehicleParams")
    .def(py::init<>())
    .def_readwrite("mass", &vleo::VehicleParams::mass)
    .def_readwrite("inertia", &vleo::VehicleParams::inertia);

  py::class_<vleo::ForeMomentAdapter, std::shared_ptr<vleo::ForeMomentAdapter>>(m, "ForeMomentAdapter")
    .def(py::init<std::shared_ptr<fore_moment_module::Model>>());

  py::class_<vleo::DeterministicPropagator>(m, "DeterministicPropagator")
    .def(py::init<std::shared_ptr<vleo::IForceMomentModel>, vleo::VehicleParams>())
    .def("propagate",
      [](const vleo::DeterministicPropagator& self,
         const Eigen::VectorXd& x0,
         const Eigen::VectorXd& t_grid,
         const std::vector<vleo::EnvInputs>& env) {
           Eigen::MatrixXd X;
           self.propagate(x0, t_grid, env, &X);
           return X;
      });

  py::class_<vleo::EnsemblePropagatorMC>(m, "EnsemblePropagatorMC")
    .def(py::init<std::shared_ptr<vleo::IForceMomentModel>, vleo::VehicleParams>())
    .def("propagate",
      [](const vleo::EnsemblePropagatorMC& self,
         const Eigen::MatrixXd& X0,
         const Eigen::VectorXd& t_grid,
         const std::vector<vleo::EnvInputs>& env) {
           Eigen::MatrixXd X;
           self.propagate(X0, t_grid, env, &X);
           return X;
      })
    .def("set_seed", &vleo::EnsemblePropagatorMC::set_seed)
    .def("set_num_threads", &vleo::EnsemblePropagatorMC::set_num_threads);
}
```

### Python façade API (`python/vleo_uq/__init__.py`)

Wrap `_vleo_uq` into a stable Python interface, add scenario/config parsing, and analysis tools.

---

## Implementation milestones (Codex tickets)

### Ticket 1: Build & link `fore_moment_module` into a Python extension

* Create CMake superbuild
* Confirm `pip install -e .` builds `_vleo_uq` and imports.

### Ticket 2: Deterministic propagator calling `ForeMomentAdapter`

* Implement high-order integrator (DOP853)
* Validate numerics on:

  * two-body (no forces)
  * drag-only sanity case

### Ticket 3: Batch ensemble propagation + OpenMP

* Memory layout: `X0` as `(N, NX)` contiguous
* Output as `(N, Nt, NX)` (expose as NumPy array)

### Ticket 4: Latent OU processes as state components

* Exact discrete OU update or SDE integration
* Seeded RNG per particle for reproducibility

### Ticket 5: Sigma-point UT propagator

* Sigma point generation in augmented error-state
* Compare mean/cov vs MC on controlled cases

### Ticket 6: STM/covariance (optional early; easier after MC/UT)

* Numerical Jacobians initially (finite differences)
* Later: add AD if you want (CppAD, autodiff)

### Ticket 7: POD UQ harness

* Truth → measurement sim → batch OD → posterior stats
* Output RTN sigmas, radial RMS, coverage

### Ticket 8: Attitude UQ + payload impact

* MEKF/UKF optional
* Payload mapping functions and smear metrics

### Ticket 9: Sensitivity framework

* Morris screening → Sobol
* Grouped then expanded

---

## Key engineering choices for your case

1. **Thread safety of `fore_moment_module`** is non-negotiable if you want MC at scale. If it’s not thread-safe, implement `clone()` and use one model instance per thread.
2. **Avoid Python callbacks** in `evaluate()`; pass all environment inputs as arrays/time series.
3. **Use deterministic seeding per particle** so results are reproducible across thread counts.
4. **Start with MC + UT first**. STM requires Jacobians; you can add it once you have a validated reference.

