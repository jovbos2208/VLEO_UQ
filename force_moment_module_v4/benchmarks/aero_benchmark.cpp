#include "shuttlecock_aero.h"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

struct BenchmarkConfig {
    int iterations = 100;
    double rho = 1e-11;
    double temperature_K = 1000.0;
    double speed_ratio = 5.0;
};

BenchmarkConfig parseArgs(int argc, char** argv) {
    BenchmarkConfig cfg;
    if (argc >= 2) {
        cfg.iterations = std::max(1, std::stoi(argv[1]));
    }
    if (argc >= 3) {
        cfg.rho = std::stod(argv[2]);
    }
    if (argc >= 4) {
        cfg.temperature_K = std::stod(argv[3]);
    }
    if (argc >= 5) {
        cfg.speed_ratio = std::stod(argv[4]);
    }
    return cfg;
}

class ScopedOstreamSilencer {
public:
    explicit ScopedOstreamSilencer(std::ostream& stream)
        : stream_(stream), old_buf_(stream.rdbuf(sink_.rdbuf())) {}
    ScopedOstreamSilencer(const ScopedOstreamSilencer&) = delete;
    ScopedOstreamSilencer& operator=(const ScopedOstreamSilencer&) = delete;
    ~ScopedOstreamSilencer() {
        stream_.rdbuf(old_buf_);
    }
private:
    std::ostream& stream_;
    std::streambuf* old_buf_;
    std::ostringstream sink_;
};

} // namespace

int main(int argc, char** argv) {
    using namespace vleo_aerodynamics_core;
    const BenchmarkConfig cfg = parseArgs(argc, argv);
    const double deg2rad = std::acos(-1.0) / 180.0;

    Eigen::Vector3d accumulated_force = Eigen::Vector3d::Zero();
    Eigen::Vector3d accumulated_torque = Eigen::Vector3d::Zero();

    auto runIteration = [&](int i) {
        const double base_velocity = 7600.0;
        const double vx = base_velocity;
        const double vy = 10.0 * std::sin(0.1 * static_cast<double>(i));
        const double vz = 5.0 * std::cos(0.07 * static_cast<double>(i));

        const double sweep = 8.0 * std::sin(0.05 * static_cast<double>(i));
        const double eta1 = (sweep + 2.0) * deg2rad;
        const double eta2 = (-sweep + 1.0) * deg2rad;
        const double eta3 = (0.5 * sweep) * deg2rad;
        const double eta4 = (-0.25 * sweep) * deg2rad;

        auto [force, torque] = shuttlecock_aero(vx, vy, vz,
                                                eta1, eta2, eta3, eta4,
                                                cfg.rho, cfg.temperature_K, cfg.speed_ratio);
        accumulated_force += force;
        accumulated_torque += torque;
    };

    std::chrono::steady_clock::time_point t_start;
    std::chrono::steady_clock::time_point t_end;
    {
        ScopedOstreamSilencer silencer(std::cout);
        // Warm-up iteration to amortize one-time allocations before timing
        runIteration(0);

        t_start = std::chrono::steady_clock::now();
        for (int i = 0; i < cfg.iterations; ++i) {
            runIteration(i);
        }
        t_end = std::chrono::steady_clock::now();
    }

    const std::chrono::duration<double> elapsed = t_end - t_start;
    const double avg_seconds = elapsed.count() / static_cast<double>(cfg.iterations);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Iterations       : " << cfg.iterations << '\n';
    std::cout << "Total time [s]   : " << elapsed.count() << '\n';
    std::cout << "Avg / call  [s]  : " << avg_seconds << '\n';
    std::cout << "Accum. force [N] : " << accumulated_force.transpose() << '\n';
    std::cout << "Accum. torque[Nm]: " << accumulated_torque.transpose() << '\n';
    return 0;
}
