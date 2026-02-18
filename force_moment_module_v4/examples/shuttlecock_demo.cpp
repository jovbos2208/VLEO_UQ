#include "shuttlecock_aero.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    using namespace vleo_aerodynamics_core;

    const double pi = std::acos(-1.0);
    const double deg2rad = pi / 180.0;
    // const double eta1_rad = 10.0 * deg2rad;
    // const double eta2_rad = -10.0 * deg2rad;
    const double rho = 1e-11;
    // const double temperature_K = 1000.0;
    // const double speed_ratio = 5.0;

    if(argc != 10){
        std::cout << "Input Error" << '\n';
        std::cout << "Input should be ./shuttlecock_demo Vx Vy Vz eta1_deg eta2_deg eta3_deg eta4_deg s T_K" << '\n';
        return 1;
    }
    
    auto [force, torque] = shuttlecock_aero(
        std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]),
        (std::stof(argv[4])* deg2rad), (std::stof(argv[5])* deg2rad),(std::stof(argv[6])* deg2rad), (std::stof(argv[7])* deg2rad),
        rho,
        std::stof(argv[9]),
        std::stof(argv[8]));

    // std::cout << "Force (N): " << force.transpose() << '\n';
    // std::cout << "Torque (Nm): " << torque.transpose() << '\n';
    std::cout << force.transpose() << " " << torque.transpose() << '\n';

    return 0;
}
