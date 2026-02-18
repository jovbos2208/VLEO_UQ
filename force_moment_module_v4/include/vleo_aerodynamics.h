#ifndef VLEO_AERODYNAMICS_H
#define VLEO_AERODYNAMICS_H

#include "body_importer.h"
#include "aerodynamics.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace vleo_aerodynamics_core {

AeroForceAndTorque vleoAerodynamics(const Eigen::Quaterniond& attitude_quaternion_BI,
                                    const Eigen::Vector3d& rotational_velocity_BI_B,
                                    const Eigen::Vector3d& velocity_I_I,
                                    const Eigen::Vector3d& wind_velocity_I_I,
                                    double density,
                                    double temperature,
                                    double particles_mass,
                                    const std::vector<Body>& bodies,
                                    const std::vector<double>& bodies_rotation_angles_rad,
                                    int temperature_ratio_method,
                                    const std::string* vtk_output_path = nullptr);

}

#endif // VLEO_AERODYNAMICS_H
