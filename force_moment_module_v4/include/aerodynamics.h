#ifndef AERODYNAMICS_H
#define AERODYNAMICS_H

#include <Eigen/Dense>
#include <vector>

namespace vleo_aerodynamics_core {

struct AeroForceAndTorque {
    Eigen::Vector3d force;
    Eigen::Vector3d torque;
};

AeroForceAndTorque calcAeroForceAndTorque(const Eigen::Ref<const Eigen::VectorXd>& areas,
                                          const Eigen::Ref<const Eigen::Matrix3Xd>& normals,
                                          const Eigen::Ref<const Eigen::Matrix3Xd>& centroids,
                                          const Eigen::Ref<const Eigen::Matrix3Xd>& v_rels,
                                          double density,
                                          double gas_temperature,
                                          const Eigen::Ref<const Eigen::VectorXd>& surface_temperatures,
                                          const Eigen::Ref<const Eigen::VectorXd>& energy_accommodation_coefficients,
                                          double particles_mass,
                                          int temperature_ratio_method,
                                          std::vector<Eigen::Vector3d>* out_pressures = nullptr);

}

#endif // AERODYNAMICS_H
