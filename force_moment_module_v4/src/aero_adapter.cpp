#include "aero_adapter.h"
#include "body_importer.h"
#include <stdexcept>
#include <iostream>
#include <cstdlib>

namespace vleo_aerodynamics_core {

void AeroAdapter::init(const GeometryConfig& geom_cfg) {
    bodies_.clear();
    bodies_ = importMultipleBodies(geom_cfg.object_files,
                                   geom_cfg.hinge_points_CAD,
                                   geom_cfg.hinge_axes_CAD,
                                   geom_cfg.temperatures_K,
                                   geom_cfg.eac,
                                   geom_cfg.DCM_B_from_CAD,
                                   geom_cfg.CoM_CAD);
    if (bodies_.size() != geom_cfg.object_files.size()) {
        throw std::runtime_error("AeroAdapter: imported bodies count mismatch");
    }
}

std::pair<Eigen::Vector3d, Eigen::Vector3d>
AeroAdapter::computeFT(const Eigen::Quaterniond& q_BI,
                       const Eigen::Vector3d& omega_BI_B,
                       const Eigen::Vector3d& V_I,
                       const Eigen::Vector3d& wind_I,
                       double density,
                       double temperature_K,
                       double particles_mass_kg,
                       double eta1_rad,
                       double eta2_rad,
                       int temperature_ratio_method) const {
    // Angles vector must match bodies_ order: [MainBody, WingRight, WingTop, WingLeft, WingBottom]
    std::vector<double> angles(bodies_.size(), 0.0);
    if (angles.size() >= 5) {
        angles[1] = -eta1_rad; // WingRight
        angles[3] = +eta1_rad; // WingLeft
        angles[2] = +eta2_rad; // WingTop
        angles[4] = -eta2_rad; // WingBottom
        if (verbose_) {
            std::cout << "WingRight: " << angles[1] << '\n';
            std::cout << "WingLeft: " << angles[3] << '\n';
            std::cout << "WingTop: " << angles[2] << '\n';
            std::cout << "WingBottom: " << angles[4] << '\n';
        }
    }
    std::string debug_vtk_path;
    const std::string* vtk_output_path = nullptr;
    if (const char* env = std::getenv("FM_DEBUG_VTK"); env && env[0] != '\0') {
        debug_vtk_path = env;
        vtk_output_path = &debug_vtk_path;
    }

    AeroForceAndTorque out = vleoAerodynamics(q_BI, omega_BI_B, V_I, wind_I,
                                              density, temperature_K, particles_mass_kg,
                                              bodies_, angles, temperature_ratio_method,
                                              vtk_output_path);
    return {out.force, out.torque};
}

} // namespace vleo_aerodynamics_core
