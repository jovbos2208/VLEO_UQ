#include "shuttlecock_aero.h"

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
using json = nlohmann::json;


int main() {
    // getting input from python from json and saving to vector
    std::string input;
    std::getline(std::cin, input);

    json parameters = json::parse(input);
    std::vector<double> vec = parameters["data"];

    //mapping input vector so that it is a 9xn Matrix with n being the number of runs
    Eigen::Map<const Eigen::Matrix<double, 9, Eigen::Dynamic>> inputMatrix(vec.data(), 9, vec.size() / 9);

    using namespace vleo_aerodynamics_core;
    const double deg2rad = std::acos(-1.0) / 180.0;

    Eigen::Vector3d accumulated_force = Eigen::Vector3d::Zero();
    Eigen::Vector3d accumulated_torque = Eigen::Vector3d::Zero();

    const double rho = 1e-11;

    for (int c = 0; c < inputMatrix.cols(); ++c) {
        // Access column c as a vector
        auto col = inputMatrix.col(c);
        auto [force, torque] = shuttlecock_aero(
        col[0], col[1], col[2],
        (col[3]* deg2rad), (col[4]* deg2rad),(col[5]* deg2rad), (col[6]* deg2rad),
        rho,
        col[8],
        col[7]);

        std::cout << col[0] << " " << col[1] << " " << col[2] << " " << col[3] << " " << col[4] << " " << col[5] << " " << col[6] << " " << col[7] << " " << col[8] << " " << force.transpose() << " " << torque.transpose() << '\n';
    }
}