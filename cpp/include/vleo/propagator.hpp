#ifndef VLEO_PROPAGATOR_HPP
#define VLEO_PROPAGATOR_HPP

#include "aero_adapter.h"
#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <vector>

namespace vleo {

struct EnvInputs {
    double density = 0.0;
    double temperature_K = 0.0;
    double particles_mass_kg = 0.0;
    Eigen::Vector3d wind_I = Eigen::Vector3d::Zero();
    double eta1_rad = 0.0;
    double eta2_rad = 0.0;
    int temperature_ratio_method = 1;
};

struct VehicleParams {
    double mass_kg = 1.0;
    Eigen::Matrix3d inertia_B = Eigen::Matrix3d::Identity();
};

struct PropagatorConfig {
    double mu_earth_m3_s2 = 3.986004418e14; // WGS-84 GM
    double rtol = 1e-4;
    double atol = 1e-6;
    double min_step_s = 1e-6;
    double max_step_s = 510.0;
    double step_safety = 0.9;
    double min_step_factor = 0.2;
    double max_step_factor = 5.0;
    int max_steps = 100000;
    // Latent OU processes (set sigma to 0.0 to disable noise).
    double rho_fast_tau_s = 0.0;
    double rho_fast_sigma = 0.0;
    double rho_bias_tau_s = 0.0;
    double rho_bias_sigma = 0.0;
    double wind_tau_s = 0.0;
    double wind_sigma = 0.0;
    std::uint64_t rng_seed = 0;
    bool verbose = false;
    int progress_stride = 50;
    int particle_stride = 10;
    bool debug_state = false;
    int debug_stride = 50;
    bool freeze_attitude = false;
};

class DeterministicPropagator {
public:
    DeterministicPropagator(std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
                            const VehicleParams& vehicle,
                            const PropagatorConfig& config = PropagatorConfig());

    void propagate(const Eigen::VectorXd& x0,
                   const Eigen::VectorXd& t_grid,
                   const std::vector<EnvInputs>& env,
                   Eigen::MatrixXd* X_out) const;

    static constexpr int kBaseStateSize = 13;
    static constexpr int kLatentSize = 5;
    static constexpr int kStateSize = kBaseStateSize + kLatentSize;

private:
    std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero_;
    VehicleParams vehicle_;
    PropagatorConfig config_;
};

class EnsemblePropagatorMC {
public:
    EnsemblePropagatorMC(std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
                         const VehicleParams& vehicle,
                         const PropagatorConfig& config = PropagatorConfig());

    void propagate(const Eigen::MatrixXd& X0,
                   const Eigen::VectorXd& t_grid,
                   const std::vector<EnvInputs>& env,
                   double* out) const;

    static constexpr int kStateSize = DeterministicPropagator::kStateSize;

private:
    std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero_;
    VehicleParams vehicle_;
    PropagatorConfig config_;
};

class SigmaPointPropagatorUT {
public:
    SigmaPointPropagatorUT(std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
                           const VehicleParams& vehicle,
                           const PropagatorConfig& config = PropagatorConfig(),
                           double alpha = 1.0,
                           double beta = 2.0,
                           double kappa = 0.0);

    void propagate(const Eigen::VectorXd& x0,
                   const Eigen::MatrixXd& P0,
                   const Eigen::VectorXd& t_grid,
                   const std::vector<EnvInputs>& env,
                   double* mean_out,
                   double* cov_out,
                   bool augment_process_noise = false) const;

    static constexpr int kStateSize = DeterministicPropagator::kStateSize;

private:
    std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero_;
    VehicleParams vehicle_;
    PropagatorConfig config_;
    double alpha_;
    double beta_;
    double kappa_;
};

class StmPropagator {
public:
    StmPropagator(std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
                  const VehicleParams& vehicle,
                  const PropagatorConfig& config = PropagatorConfig(),
                  double fd_eps_rel = 1e-6,
                  double fd_eps_abs = 1e-8);

    void propagate(const Eigen::VectorXd& x0,
                   const Eigen::MatrixXd& P0,
                   const Eigen::VectorXd& t_grid,
                   const std::vector<EnvInputs>& env,
                   double* mean_out,
                   double* cov_out,
                   double* stm_out = nullptr) const;

    static constexpr int kStateSize = DeterministicPropagator::kStateSize;

private:
    std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero_;
    VehicleParams vehicle_;
    PropagatorConfig config_;
    double fd_eps_rel_;
    double fd_eps_abs_;
};

} // namespace vleo

#endif // VLEO_PROPAGATOR_HPP
