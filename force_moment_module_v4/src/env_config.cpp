#include "env_config.h"
#include <filesystem>

namespace vleo_aerodynamics_core {

static std::string resolve_example_path_impl(const std::string& name) {
    namespace fs = std::filesystem;
    fs::path here = fs::path(__FILE__).parent_path();
    std::vector<fs::path> candidates = {
        here / ".." / "data" / name,
        fs::current_path() / "data" / name,
        here / name
    };
    for (const auto& cand : candidates) {
        fs::path abs = fs::absolute(cand);
        if (fs::exists(abs)) {
            return abs.string();
        }
    }
    return fs::absolute(here / ".." / "data" / name).string();
}



static std::string choose_first_existing(const std::initializer_list<const char*>& names) {
    for (auto* n : names) {
        std::string p = resolve_example_path_impl(n);
        if (std::filesystem::exists(p)) return p;
    }
    // Return first resolved even if not existing (vleo-aero tools may create later)
    auto it = names.begin();
    return resolve_example_path_impl(*it);
}

GeometryConfig ShuttlecockDefaults::geometry() {
    GeometryConfig cfg;
    // Files order must match hinge/axis columns: [MainBody, Wing1, Wing2, Wing3, Wing4]
    // Keep this consistent with shuttlecock_aero.cpp so Python and C++ paths use the same assembly.
    cfg.object_files = {
        choose_first_existing({"mainBody.obj", "MainBody.obj", "shuttle_box.obj", "box.obj", "cube.obj"}),
        choose_first_existing({"wing_1.obj", "wing_right.obj", "wing_pos_x.obj"}),
        choose_first_existing({"wing_2.obj", "wing_pos_y.obj"}),
        choose_first_existing({"wing_3.obj", "wing_left.obj", "wing_neg_x.obj"}),
        choose_first_existing({"wing_4.obj", "wing_neg_y.obj"}),
    };

    cfg.hinge_points_CAD.resize(3, 5);
    cfg.hinge_axes_CAD.resize(3, 5);
    // Same hinge layout used in shuttlecock_aero.cpp.
    const double box_half = 0.05;
    cfg.hinge_points_CAD.col(0) << 0.0, 0.0, 0.0;
    cfg.hinge_points_CAD.col(1) << -0.34965, 0.0, -box_half;
    cfg.hinge_points_CAD.col(2) << -0.34965, -box_half, 0.0;
    cfg.hinge_points_CAD.col(3) << -0.34965, 0.0, +box_half;
    cfg.hinge_points_CAD.col(4) << -0.34965, +box_half, 0.0;

    cfg.hinge_axes_CAD.col(0) = Eigen::Vector3d::UnitZ();
    cfg.hinge_axes_CAD.col(1) = Eigen::Vector3d::UnitY();
    cfg.hinge_axes_CAD.col(2) = Eigen::Vector3d::UnitZ();
    cfg.hinge_axes_CAD.col(3) = Eigen::Vector3d::UnitY();
    cfg.hinge_axes_CAD.col(4) = Eigen::Vector3d::UnitZ();

    cfg.temperatures_K = std::vector<std::vector<double>>(5, std::vector<double>{300.0});
    cfg.eac            = std::vector<std::vector<double>>(5, std::vector<double>{0.9});

    // Keep CAD frame identical to body frame for this shuttlecock assembly.
    cfg.DCM_B_from_CAD = Eigen::Matrix3d::Identity();
    cfg.CoM_CAD = Eigen::Vector3d::Zero();

    return cfg;
}

} // namespace vleo_aerodynamics_core
