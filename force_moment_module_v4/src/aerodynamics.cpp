#include "aerodynamics.h"
#include <iostream>
#include <cmath>
#include <vector>

namespace vleo_aerodynamics_core {

const double KB = 1.38064852e-23; // Boltzmann constant

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
                                          std::vector<Eigen::Vector3d>* out_pressures) {

    const int count = static_cast<int>(areas.size());
    if (count == 0) {
        return {Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()};
    }

    const double cm = std::sqrt(2.0 * KB * gas_temperature / particles_mass);
    const double prefactor = 0.5 * density * cm * cm;
    const double inv_sqrt_pi = 1.0 / std::sqrt(M_PI);
    const double sqrt_pi_over_two = std::sqrt(M_PI) / 2.0;

    Eigen::ArrayXd V = v_rels.colwise().norm().array();
    Eigen::ArrayXd invV = Eigen::ArrayXd::Zero(count);
    Eigen::Array<bool, Eigen::Dynamic, 1> active = (V > 1e-12);
    invV = active.select(V.cwiseInverse(), 0.0);
    Eigen::ArrayXd s = V / cm;
    Eigen::ArrayXd s2 = s.square();

    Eigen::Matrix3Xd vhat = v_rels;
    Eigen::RowVectorXd invV_row = invV.matrix().transpose();
    for (int row = 0; row < vhat.rows(); ++row) {
        vhat.row(row).array() *= invV_row.array();
    }
    Eigen::ArrayXd cosdelta = -(vhat.array() * normals.array()).colwise().sum();
    Eigen::ArrayXd scosdelta = s * cosdelta;
    Eigen::ArrayXd e = (-scosdelta.square()).exp();
    Eigen::ArrayXd erfcterm = (-scosdelta).unaryExpr([](double x) { return std::erfc(x); });
    Eigen::ArrayXd G1 = scosdelta * inv_sqrt_pi * e + (0.5 + scosdelta.square()) * erfcterm;
    Eigen::ArrayXd G2 = inv_sqrt_pi * e + scosdelta * erfcterm;

    Eigen::ArrayXd T_rat(count);
    Eigen::ArrayXd eac = energy_accommodation_coefficients.array();
    Eigen::ArrayXd invV2 = invV.square();
    Eigen::ArrayXd temps = surface_temperatures.array();
    switch (temperature_ratio_method) {
        case 1: {
            Eigen::ArrayXd enum_term = scosdelta * erfcterm;
            Eigen::ArrayXd denom = inv_sqrt_pi * e + enum_term;
            Eigen::ArrayXd term1 = eac * (2.0 * KB * temps * invV2 * s2);
            Eigen::ArrayXd term2 = (1.0 - eac) * (1.0 + s2 / 2.0 + 0.25 * enum_term / denom);
            T_rat = term1 + term2;
            break;
        }
        case 2: {
            Eigen::ArrayXd inner = (4.0 * KB * temps * invV2) - 1.0;
            T_rat = 0.5 * s2 * (1.0 + eac * inner) + 1.25 * (1.0 - eac);
            break;
        }
        case 3: {
            Eigen::ArrayXd inner = (4.0 * KB * temps * invV2) - 1.0;
            T_rat = 0.5 * s2 * (1.0 + eac * inner);
            break;
        }
        default:
            throw std::runtime_error("Invalid temperature ratio method");
    }

    Eigen::ArrayXd sqrtTr = T_rat.sqrt();
    Eigen::ArrayXd coeff_normals = -(G1 + sqrt_pi_over_two * sqrtTr * G2);
    Eigen::Matrix3Xd pressures = normals;
    Eigen::RowVectorXd coeff_row = coeff_normals.matrix().transpose();
    for (int row = 0; row < pressures.rows(); ++row) {
        pressures.row(row).array() *= coeff_row.array();
    }
    Eigen::Matrix3Xd scaled_normals = normals;
    Eigen::RowVectorXd cos_row = cosdelta.matrix().transpose();
    for (int row = 0; row < scaled_normals.rows(); ++row) {
        scaled_normals.row(row).array() *= cos_row.array();
    }
    Eigen::Matrix3Xd term = vhat + scaled_normals;
    Eigen::ArrayXd sG2 = s * G2;
    Eigen::RowVectorXd sG2_row = sG2.matrix().transpose();
    for (int row = 0; row < pressures.rows(); ++row) {
        pressures.row(row).array() += term.row(row).array() * sG2_row.array();
    }
    pressures *= prefactor;
    Eigen::RowVectorXd mask_row = active.cast<double>().matrix().transpose();
    for (int row = 0; row < pressures.rows(); ++row) {
        pressures.row(row).array() *= mask_row.array();
    }

    if (out_pressures) {
        out_pressures->resize(count);
        for (int i = 0; i < count; ++i) {
            (*out_pressures)[i] = pressures.col(i);
        }
    }

    Eigen::Matrix3Xd forces = pressures;
    Eigen::RowVectorXd areas_row = areas.transpose();
    for (int row = 0; row < forces.rows(); ++row) {
        forces.row(row).array() *= areas_row.array();
    }
    Eigen::Vector3d total_force = forces.rowwise().sum();

    Eigen::Vector3d total_torque = Eigen::Vector3d::Zero();
    for (int i = 0; i < count; ++i) {
        total_torque += centroids.col(i).cross(forces.col(i));
    }

    return {total_force, total_torque};
}

}
