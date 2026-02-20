#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "aero_adapter.h"
#include "env_config.h"
#include "vleo/propagator.hpp"

#include <stdexcept>

namespace py = pybind11;
using vleo_aerodynamics_core::AeroAdapter;
using vleo_aerodynamics_core::AeroConfig;
using vleo_aerodynamics_core::GeometryConfig;
using vleo_aerodynamics_core::MassProperties;
using vleo_aerodynamics_core::ShuttlecockDefaults;
using vleo_aerodynamics_core::SpaceWeatherConfig;
using vleo::DeterministicPropagator;
using vleo::EnsemblePropagatorMC;
using vleo::SigmaPointPropagatorUT;
using vleo::StmPropagator;
using vleo::EnvInputs;
using vleo::PropagatorConfig;
using vleo::VehicleParams;

PYBIND11_MODULE(_vleo_uq, m) {
    m.doc() = "VLEO UQ bindings for force/torque evaluation";

    py::class_<MassProperties>(m, "MassProperties")
        .def(py::init<>())
        .def_readwrite("mass_kg", &MassProperties::mass_kg)
        .def_readwrite("inertia_B", &MassProperties::inertia_B);

    py::class_<GeometryConfig>(m, "GeometryConfig")
        .def(py::init<>())
        .def_readwrite("object_files", &GeometryConfig::object_files)
        .def_readwrite("hinge_points_CAD", &GeometryConfig::hinge_points_CAD)
        .def_readwrite("hinge_axes_CAD", &GeometryConfig::hinge_axes_CAD)
        .def_readwrite("temperatures_K", &GeometryConfig::temperatures_K)
        .def_readwrite("eac", &GeometryConfig::eac)
        .def_readwrite("DCM_B_from_CAD", &GeometryConfig::DCM_B_from_CAD)
        .def_readwrite("CoM_CAD", &GeometryConfig::CoM_CAD);

    py::class_<AeroConfig>(m, "AeroConfig")
        .def(py::init<>())
        .def_readwrite("temperature_ratio_method", &AeroConfig::temperature_ratio_method)
        .def_readwrite("use_nrlmsis_particle_mass", &AeroConfig::use_nrlmsis_particle_mass)
        .def_readwrite("fixed_particle_mass_kg", &AeroConfig::fixed_particle_mass_kg);

    py::class_<SpaceWeatherConfig>(m, "SpaceWeatherConfig")
        .def(py::init<>())
        .def_readwrite("f107a", &SpaceWeatherConfig::f107a)
        .def_readwrite("f107", &SpaceWeatherConfig::f107)
        .def_readwrite("ap", &SpaceWeatherConfig::ap);

    m.def("default_geometry", &ShuttlecockDefaults::geometry,
          "Return the default shuttlecock geometry configuration.");

    py::class_<AeroAdapter, std::shared_ptr<AeroAdapter>>(m, "AeroAdapter")
        .def(py::init<>())
        .def("init", &AeroAdapter::init, py::arg("geometry"))
        .def_property("verbose", &AeroAdapter::verbose, &AeroAdapter::set_verbose)
        .def(
            "compute_ft",
            [](AeroAdapter& self,
               const Eigen::Vector4d& q_wxyz,
               const Eigen::Vector3d& omega_BI_B,
               const Eigen::Vector3d& v_I,
               const Eigen::Vector3d& wind_I,
               double density,
               double temperature_K,
               double particles_mass_kg,
               double eta1_rad,
               double eta2_rad,
               int temperature_ratio_method) {
                if (q_wxyz.norm() == 0.0) {
                    throw std::invalid_argument("q_wxyz must be non-zero");
                }
                Eigen::Quaterniond q(q_wxyz[0], q_wxyz[1], q_wxyz[2], q_wxyz[3]);
                q.normalize();
                auto out = self.computeFT(q, omega_BI_B, v_I, wind_I,
                                          density, temperature_K, particles_mass_kg,
                                          eta1_rad, eta2_rad, temperature_ratio_method);
                return py::make_tuple(out.first, out.second);
            },
            py::arg("q_wxyz"),
            py::arg("omega_BI_B"),
            py::arg("v_I"),
            py::arg("wind_I"),
            py::arg("density"),
            py::arg("temperature_K"),
            py::arg("particles_mass_kg"),
            py::arg("eta1_rad"),
            py::arg("eta2_rad"),
            py::arg("temperature_ratio_method") = 1,
            "Compute body-frame force/torque. Quaternion is [w, x, y, z].");

    py::class_<EnvInputs>(m, "EnvInputs")
        .def(py::init<>())
        .def_readwrite("density", &EnvInputs::density)
        .def_readwrite("temperature_K", &EnvInputs::temperature_K)
        .def_readwrite("particles_mass_kg", &EnvInputs::particles_mass_kg)
        .def_readwrite("wind_I", &EnvInputs::wind_I)
        .def_readwrite("sun_position_I_m", &EnvInputs::sun_position_I_m)
        .def_readwrite("moon_position_I_m", &EnvInputs::moon_position_I_m)
        .def_readwrite("magnetic_field_I_T", &EnvInputs::magnetic_field_I_T)
        .def_readwrite("srp_scale", &EnvInputs::srp_scale)
        .def_readwrite("albedo_ir_scale", &EnvInputs::albedo_ir_scale)
        .def_readwrite("tide_loading_accel_I_m_s2", &EnvInputs::tide_loading_accel_I_m_s2)
        .def_readwrite("eta1_rad", &EnvInputs::eta1_rad)
        .def_readwrite("eta2_rad", &EnvInputs::eta2_rad)
        .def_readwrite("temperature_ratio_method", &EnvInputs::temperature_ratio_method);

    py::class_<VehicleParams>(m, "VehicleParams")
        .def(py::init<>())
        .def_readwrite("mass_kg", &VehicleParams::mass_kg)
        .def_readwrite("inertia_B", &VehicleParams::inertia_B);

    py::class_<PropagatorConfig>(m, "PropagatorConfig")
        .def(py::init<>())
        .def_readwrite("mu_earth_m3_s2", &PropagatorConfig::mu_earth_m3_s2)
        .def_readwrite("use_j2_perturbation", &PropagatorConfig::use_j2_perturbation)
        .def_readwrite("j2_earth", &PropagatorConfig::j2_earth)
        .def_readwrite("use_j3_perturbation", &PropagatorConfig::use_j3_perturbation)
        .def_readwrite("j3_earth", &PropagatorConfig::j3_earth)
        .def_readwrite("use_j4_perturbation", &PropagatorConfig::use_j4_perturbation)
        .def_readwrite("j4_earth", &PropagatorConfig::j4_earth)
        .def_readwrite("gravity_fd_step_m", &PropagatorConfig::gravity_fd_step_m)
        .def_readwrite("earth_equatorial_radius_m", &PropagatorConfig::earth_equatorial_radius_m)
        .def_readwrite("use_sun_third_body", &PropagatorConfig::use_sun_third_body)
        .def_readwrite("use_moon_third_body", &PropagatorConfig::use_moon_third_body)
        .def_readwrite("mu_sun_m3_s2", &PropagatorConfig::mu_sun_m3_s2)
        .def_readwrite("mu_moon_m3_s2", &PropagatorConfig::mu_moon_m3_s2)
        .def_readwrite("sun_ephemeris_scale", &PropagatorConfig::sun_ephemeris_scale)
        .def_readwrite("moon_ephemeris_scale", &PropagatorConfig::moon_ephemeris_scale)
        .def_readwrite("use_srp_acceleration", &PropagatorConfig::use_srp_acceleration)
        .def_readwrite("srp_cr", &PropagatorConfig::srp_cr)
        .def_readwrite("srp_area_m2", &PropagatorConfig::srp_area_m2)
        .def_readwrite("solar_pressure_1au_n_m2", &PropagatorConfig::solar_pressure_1au_n_m2)
        .def_readwrite("astronomical_unit_m", &PropagatorConfig::astronomical_unit_m)
        .def_readwrite("use_albedo_ir_acceleration", &PropagatorConfig::use_albedo_ir_acceleration)
        .def_readwrite("albedo_ir_cr", &PropagatorConfig::albedo_ir_cr)
        .def_readwrite("albedo_ir_area_m2", &PropagatorConfig::albedo_ir_area_m2)
        .def_readwrite("albedo_pressure_n_m2", &PropagatorConfig::albedo_pressure_n_m2)
        .def_readwrite("earth_ir_pressure_n_m2", &PropagatorConfig::earth_ir_pressure_n_m2)
        .def_readwrite("use_tide_loading_acceleration", &PropagatorConfig::use_tide_loading_acceleration)
        .def_readwrite("tide_loading_scale", &PropagatorConfig::tide_loading_scale)
        .def_readwrite("use_magnetic_torque", &PropagatorConfig::use_magnetic_torque)
        .def_readwrite("residual_dipole_B_A_m2", &PropagatorConfig::residual_dipole_B_A_m2)
        .def_readwrite("residual_dipole_scale", &PropagatorConfig::residual_dipole_scale)
        .def_readwrite("magnetic_field_scale", &PropagatorConfig::magnetic_field_scale)
        .def_readwrite("rtol", &PropagatorConfig::rtol)
        .def_readwrite("atol", &PropagatorConfig::atol)
        .def_readwrite("min_step_s", &PropagatorConfig::min_step_s)
        .def_readwrite("max_step_s", &PropagatorConfig::max_step_s)
        .def_readwrite("step_safety", &PropagatorConfig::step_safety)
        .def_readwrite("min_step_factor", &PropagatorConfig::min_step_factor)
        .def_readwrite("max_step_factor", &PropagatorConfig::max_step_factor)
        .def_readwrite("max_steps", &PropagatorConfig::max_steps)
        .def_readwrite("rho_fast_tau_s", &PropagatorConfig::rho_fast_tau_s)
        .def_readwrite("rho_fast_sigma", &PropagatorConfig::rho_fast_sigma)
        .def_readwrite("rho_bias_tau_s", &PropagatorConfig::rho_bias_tau_s)
        .def_readwrite("rho_bias_sigma", &PropagatorConfig::rho_bias_sigma)
        .def_readwrite("wind_tau_s", &PropagatorConfig::wind_tau_s)
        .def_readwrite("wind_sigma", &PropagatorConfig::wind_sigma)
        .def_readwrite("rng_seed", &PropagatorConfig::rng_seed)
        .def_readwrite("verbose", &PropagatorConfig::verbose)
        .def_readwrite("progress_stride", &PropagatorConfig::progress_stride)
        .def_readwrite("particle_stride", &PropagatorConfig::particle_stride)
        .def_readwrite("debug_state", &PropagatorConfig::debug_state)
        .def_readwrite("debug_stride", &PropagatorConfig::debug_stride)
        .def_readwrite("freeze_attitude", &PropagatorConfig::freeze_attitude);

    py::class_<DeterministicPropagator>(m, "DeterministicPropagator")
        .def(py::init<std::shared_ptr<AeroAdapter>, const VehicleParams&, const PropagatorConfig&>(),
             py::arg("aero"),
             py::arg("vehicle"),
             py::arg("config") = PropagatorConfig())
        .def_property_readonly_static("state_size", [](py::object) {
            return DeterministicPropagator::kStateSize;
        })
        .def(
            "propagate",
            [](const DeterministicPropagator& self,
               const Eigen::VectorXd& x0,
               const Eigen::VectorXd& t_grid,
               const std::vector<EnvInputs>& env) {
                Eigen::MatrixXd X;
                self.propagate(x0, t_grid, env, &X);
                return X;
            },
            py::arg("x0"),
            py::arg("t_grid"),
            py::arg("env"),
            "Propagate state. Layout: r(3), v(3), q_wxyz(4), w_BI_B(3), log_rho_fast, log_rho_bias, wind_I(3).");

    py::class_<EnsemblePropagatorMC>(m, "EnsemblePropagatorMC")
        .def(py::init<std::shared_ptr<AeroAdapter>, const VehicleParams&, const PropagatorConfig&>(),
             py::arg("aero"),
             py::arg("vehicle"),
             py::arg("config") = PropagatorConfig())
        .def_property_readonly_static("state_size", [](py::object) {
            return EnsemblePropagatorMC::kStateSize;
        })
        .def(
            "propagate",
            [](const EnsemblePropagatorMC& self,
               const Eigen::MatrixXd& X0,
               const Eigen::VectorXd& t_grid,
               const std::vector<EnvInputs>& env) {
                const int nt = static_cast<int>(t_grid.size());
                const int n = static_cast<int>(X0.rows());
                py::array_t<double> out({nt, n, EnsemblePropagatorMC::kStateSize});
                self.propagate(X0, t_grid, env, out.mutable_data());
                return out;
            },
            py::arg("X0"),
            py::arg("t_grid"),
            py::arg("env"),
            "Propagate ensemble. Output layout: (Nt, N, state_size).");

    py::class_<SigmaPointPropagatorUT>(m, "SigmaPointPropagatorUT")
        .def(py::init<std::shared_ptr<AeroAdapter>, const VehicleParams&, const PropagatorConfig&,
                      double, double, double>(),
             py::arg("aero"),
             py::arg("vehicle"),
             py::arg("config") = PropagatorConfig(),
             py::arg("alpha") = 1.0,
             py::arg("beta") = 2.0,
             py::arg("kappa") = 0.0)
        .def_property_readonly_static("state_size", [](py::object) {
            return SigmaPointPropagatorUT::kStateSize;
        })
        .def(
            "propagate",
            [](const SigmaPointPropagatorUT& self,
               const Eigen::VectorXd& x0,
               const Eigen::MatrixXd& P0,
               const Eigen::VectorXd& t_grid,
               const std::vector<EnvInputs>& env,
               bool augment_process_noise) {
                const int nt = static_cast<int>(t_grid.size());
                py::array_t<double> mean({nt, SigmaPointPropagatorUT::kStateSize});
                py::array_t<double> cov({nt, SigmaPointPropagatorUT::kStateSize, SigmaPointPropagatorUT::kStateSize});
                self.propagate(x0, P0, t_grid, env, mean.mutable_data(), cov.mutable_data(),
                               augment_process_noise);
                return py::make_tuple(mean, cov);
            },
            py::arg("x0"),
            py::arg("P0"),
            py::arg("t_grid"),
            py::arg("env"),
            py::arg("augment_process_noise") = false,
            "Propagate sigma points and return (mean, covariance).");

    py::class_<StmPropagator>(m, "StmPropagator")
        .def(py::init<std::shared_ptr<AeroAdapter>, const VehicleParams&, const PropagatorConfig&,
                      double, double>(),
             py::arg("aero"),
             py::arg("vehicle"),
             py::arg("config") = PropagatorConfig(),
             py::arg("fd_eps_rel") = 1e-6,
             py::arg("fd_eps_abs") = 1e-8)
        .def_property_readonly_static("state_size", [](py::object) {
            return StmPropagator::kStateSize;
        })
        .def(
            "propagate",
            [](const StmPropagator& self,
               const Eigen::VectorXd& x0,
               const Eigen::MatrixXd& P0,
               const Eigen::VectorXd& t_grid,
               const std::vector<EnvInputs>& env,
               bool return_stm) -> py::object {
                const int nt = static_cast<int>(t_grid.size());
                py::array_t<double> mean({nt, StmPropagator::kStateSize});
                py::array_t<double> cov({nt, StmPropagator::kStateSize, StmPropagator::kStateSize});
                if (return_stm) {
                    py::array_t<double> stm({nt, StmPropagator::kStateSize, StmPropagator::kStateSize});
                    self.propagate(x0, P0, t_grid, env, mean.mutable_data(), cov.mutable_data(),
                                   stm.mutable_data());
                    return py::make_tuple(mean, cov, stm);
                }
                self.propagate(x0, P0, t_grid, env, mean.mutable_data(), cov.mutable_data(), nullptr);
                return py::make_tuple(mean, cov);
            },
            py::arg("x0"),
            py::arg("P0"),
            py::arg("t_grid"),
            py::arg("env"),
            py::arg("return_stm") = false,
            "Propagate mean/cov via finite-difference STM. If return_stm, returns (mean, cov, stm).");
}
