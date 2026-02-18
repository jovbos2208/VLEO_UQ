#ifndef SHUTTLECOCK_AERO_H
#define SHUTTLECOCK_AERO_H

#include <Eigen/Dense>
#include <string>

namespace vleo_aerodynamics_core {

struct ShuttlecockDebugOptions {
    bool enable_debug_vtk = false;
    std::string vtk_filename = "shuttlecock_debug.vtk";
};

// Computes aerodynamic force and torque for a shuttlecock-like satellite
// with a box body and four wings. Wings along +X/-X are deflected by eta1,
// wings along +Y/-Y by eta2, rotating outward about the box edges (hinge lines)
// parallel to +Z at x=±box_half and y=±box_half.
// Inputs:
//  Vx,Vy,Vz: components of inertial-relative velocity [m/s]
//  eta1_rad: deflection angle for ±X wings [rad]
//  eta2_rad: deflection angle for ±Y wings [rad]
//  rho: density [kg/m^3]
//  T_K: gas temperature [K]
//  s: speed ratio = |V| / cm, cm = sqrt(2 kB T / m), used to derive particle mass
// Returns force and torque in body frame.
std::pair<Eigen::Vector3d, Eigen::Vector3d> shuttlecock_aero(double Vx, double Vy, double Vz,
                                                             double eta1_rad, double eta2_rad,
                                                             double eta3_rad, double eta4_rad,
                                                             double rho, double T_K, double s,
                                                             const ShuttlecockDebugOptions& debug_opts = {});

// Pre-loads shuttlecock geometry (OBJ meshes, hinge data, etc.) so subsequent
// calls to shuttlecock_aero() only adjust the pre-existing bodies.
void preload_shuttlecock_geometry();

}

#endif // SHUTTLECOCK_AERO_H
