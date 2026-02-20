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
    Eigen::Vector3d sun_position_I_m = Eigen::Vector3d::Zero();
    Eigen::Vector3d moon_position_I_m = Eigen::Vector3d::Zero();
    Eigen::Vector3d magnetic_field_I_T = Eigen::Vector3d::Zero();
    double srp_scale = 1.0;
    double albedo_ir_scale = 1.0;
    Eigen::Vector3d tide_loading_accel_I_m_s2 = Eigen::Vector3d::Zero();
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
    bool use_j2_perturbation = false;
    double j2_earth = 1.08262668e-3; // WGS-84 J2
    bool use_j3_perturbation = false;
    double j3_earth = -2.53215306e-6; // WGS-84 J3
    bool use_j4_perturbation = false;
    double j4_earth = -1.61098761e-6; // WGS-84 J4
    double gravity_fd_step_m = 10.0;
    double earth_equatorial_radius_m = 6378137.0;
    bool use_sun_third_body = false;
    bool use_moon_third_body = false;
    double mu_sun_m3_s2 = 1.32712440018e20;
    double mu_moon_m3_s2 = 4.9048695e12;
    double sun_ephemeris_scale = 1.0;
    double moon_ephemeris_scale = 1.0;
    bool use_srp_acceleration = false;
    double srp_cr = 1.2;
    double srp_area_m2 = 0.0;
    double solar_pressure_1au_n_m2 = 4.56e-6;
    double astronomical_unit_m = 149597870700.0;
    bool use_albedo_ir_acceleration = false;
    double albedo_ir_cr = 1.0;
    double albedo_ir_area_m2 = 0.0;
    double albedo_pressure_n_m2 = 1.2e-6;
    double earth_ir_pressure_n_m2 = 1.0e-6;
    bool use_tide_loading_acceleration = false;
    double tide_loading_scale = 1.0;
    bool use_magnetic_torque = false;
    Eigen::Vector3d residual_dipole_B_A_m2 = Eigen::Vector3d::Zero();
    double residual_dipole_scale = 1.0;
    double magnetic_field_scale = 1.0;
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
