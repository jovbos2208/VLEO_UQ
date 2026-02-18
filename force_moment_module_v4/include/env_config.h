#ifndef ENV_CONFIG_H
#define ENV_CONFIG_H

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <array>

namespace vleo_aerodynamics_core {

struct MassProperties {
    double mass_kg = 1.0;                 // vehicle mass
    Eigen::Matrix3d inertia_B = Eigen::Matrix3d::Identity(); // body inertia (kg m^2)
};

struct GeometryConfig {
    // Five-part shuttlecock default: [MainBody, WingRight, WingTop, WingLeft, WingBottom]
    std::vector<std::string> object_files;
    Eigen::Matrix3Xd hinge_points_CAD;    // 3xN, meters in CAD frame
    Eigen::Matrix3Xd hinge_axes_CAD;      // 3xN, unit vectors in CAD frame
    std::vector<std::vector<double>> temperatures_K; // per-body temperature lists
    std::vector<std::vector<double>> eac;            // per-body energy accommodation
    Eigen::Matrix3d DCM_B_from_CAD = Eigen::Matrix3d::Identity();
    Eigen::Vector3d CoM_CAD = Eigen::Vector3d::Zero();
};

struct AeroConfig {
    int temperature_ratio_method = 1;   // keep 1 by default
    // Particle mass selection:
    // If use_nrlmsis_particle_mass = true, compute mean molecular mass from NRLMSIS composition.
    // Otherwise, use fixed_particle_mass_kg.
    bool use_nrlmsis_particle_mass = true;
    double fixed_particle_mass_kg = 16.0 * 1.6605390689252e-27; // 16u by default
};

struct ShuttlecockDefaults {
    // Helper to populate a default five-part shuttlecock configuration
    static GeometryConfig geometry();
};

struct SpaceWeatherConfig {
    float f107a = 150.0f;
    float f107 = 150.0f;
    std::array<float,7> ap = {15.0f,12.0f,12.0f,12.0f,12.0f,12.0f,12.0f};
};

struct EnvConfig {
    GeometryConfig geometry = ShuttlecockDefaults::geometry();
    AeroConfig aero;
    MassProperties mass_props;
    SpaceWeatherConfig space_weather;
    // Initial conditions
    double init_altitude_m = 300e3;  // if range_min == range_max, use this exact value
    double init_altitude_min_m = 300e3;
    double init_altitude_max_m = 300e3;
    double init_incl_min_rad = 0.0;
    double init_incl_max_rad = 0.0;
    bool randomize_raan = true;
    bool randomize_argp = true;
    bool randomize_true_anomaly = true;
    double init_raan_rad = 0.0;
    double init_argp_rad = 0.0;
    double init_true_anomaly_rad = 0.0;
    // Targets
    double target_altitude_m = 300e3;
    // Control limits
    double eta_limit_rad = 100.0 * M_PI/180.0;
    double thrust_min_N = 0.0;
    double thrust_max_N = 0.0;
    double total_impulse_budget_Ns = 0.0; // total thrust integral budget per episode
    // Episode control
    int N_orbits_per_episode = 10;
    double orbit_period_min_s = 3000.0;  // clamps
    double orbit_period_max_s = 20000.0;
    // Integration defaults
    double dt_initial_s = 1.0;
    double dp54_atol = 1e-6;
    double dp54_rtol = 1e-6;
    double dp54_min_dt = 1e-4;
    double dp54_max_dt = 10.0;
    double dp54_safety = 0.9;
    // Atmosphere sampling cache (disabled by default)
    // If period > 0, reuse NRLMSIS results for up to this duration and small geo changes
    double atmo_cache_period_s = 0.0;      // max age of cache in seconds (0 => disabled)
    double atmo_cache_alt_tol_m = 0.0;     // altitude tolerance for cache reuse
    double atmo_cache_latlon_tol_deg = 0.0;// lat/lon tolerance for cache reuse
    // Rewards
    double w_alt = 1e-6; // per meter (tracking)
    double w_att_ao = 1.0; // AoA/AoS tracking weight
    double w_att_roll = 1.0; // roll tracking weight
    double w_effort = 1e-2; // control effort penalty
    double w_dwell_alt = 1.0; // bonus per-orbit for staying in altitude band
    double w_dwell_att = 1.0; // bonus per-orbit for staying in attitude band
    double w_spin = 0.1; // penalty weight for angular rates (rad/s) across all three axes

    // Bands
    double alt_band_halfwidth_m = 1000.0; // +/- around target_altitude_m
    double aoa_band_rad = 5.0 * M_PI/180.0; // +/-
    double aos_band_rad = 5.0 * M_PI/180.0; // +/-
    double roll_band_rad = 5.0 * M_PI/180.0; // +/-
    double omega_band_rad_s = 5.0 * M_PI/180.0; // angular rate magnitude band

    // Initialization
    bool align_attitude_to_velocity = false; // if true, set body +X along initial Vrel
};

} // namespace vleo_aerodynamics_core

#endif // ENV_CONFIG_H
