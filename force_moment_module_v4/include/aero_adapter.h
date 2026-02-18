#ifndef AERO_ADAPTER_H
#define AERO_ADAPTER_H

#include "env_config.h"
#include "vleo_aerodynamics.h"
#include <Eigen/Dense>
#include <utility>

namespace vleo_aerodynamics_core {

// Thin wrapper that owns imported bodies and invokes vleoAerodynamics
class AeroAdapter {
public:
    AeroAdapter() = default;

    // Initialize geometry from config (imports bodies once)
    void init(const GeometryConfig& geom_cfg);

    // Compute aerodynamic force/torque in body frame using the black-box function
    // angles are eta1, eta2 for shuttlecock; mapping applied internally
    std::pair<Eigen::Vector3d, Eigen::Vector3d>
    computeFT(const Eigen::Quaterniond& q_BI,
              const Eigen::Vector3d& omega_BI_B,
              const Eigen::Vector3d& V_I,
              const Eigen::Vector3d& wind_I,
              double density,
              double temperature_K,
              double particles_mass_kg,
              double eta1_rad,
              double eta2_rad,
              int temperature_ratio_method = 1) const;

    void set_verbose(bool verbose) { verbose_ = verbose; }
    bool verbose() const { return verbose_; }

private:
    std::vector<Body> bodies_;
    bool verbose_ = false;
};

} // namespace vleo_aerodynamics_core

#endif // AERO_ADAPTER_H
