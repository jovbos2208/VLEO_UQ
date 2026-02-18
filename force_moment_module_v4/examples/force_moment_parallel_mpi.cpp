#include "shuttlecock_aero.h"

#include <Eigen/Dense>
#include <mpi.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>

namespace {

struct RunConfig {
    bool coarse_sweep = true;   // matches TESTLAUF_100K behaviour from the Python script
    std::string output_file = "results.csv";
    std::size_t chunk_rows = 256;  // number of rows sent per MPI message / write flush
    double orbital_speed = 7500.0; // m/s
    std::uint64_t max_cases = 0;   // 0 => no limit
};

struct SweepValues {
    std::vector<int> speed_ratios;
    std::vector<int> temperatures;
    std::vector<int> angles;
    std::vector<int> wing_deflections;
};

struct SweepSteps {
    int speed_ratio = 1;
    int temperature = 50;
    int angle = 1;
    int wing = 1;
};

struct DecodedIndices {
    std::size_t speed_ratio = 0;
    std::size_t temperature = 0;
    std::size_t aoa = 0;
    std::size_t aos = 0;
    std::size_t eta2 = 0;
    std::size_t eta1 = 0;
};

struct WorkRange {
    std::uint64_t start = 0;
    std::uint64_t count = 0;
};

struct SweepCache {
    std::size_t angle_span = 0;
    std::vector<Eigen::Vector3d> velocities;
};

constexpr int kRowWidth = 15;
constexpr double kEta3Deg = 0.0;
constexpr double kEta4Deg = 0.0;
constexpr double kRho = 1e-11;
constexpr int kTagRows = 1;
constexpr int kTagData = 2;

std::vector<int> make_range(int start, int stop_exclusive, int step) {
    if (step <= 0) {
        throw std::invalid_argument("step must be positive");
    }
    std::vector<int> values;
    for (int value = start; value < stop_exclusive; value += step) {
        values.push_back(value);
    }
    return values;
}

SweepSteps steps_from_config(const RunConfig& config) {
    if (config.coarse_sweep) {
        return {10, 100, 11, 10};
    }
    return {1, 50, 1, 1};
}

SweepValues build_sweep_values(const RunConfig& config) {
    const auto steps = steps_from_config(config);
    SweepValues values;
    values.speed_ratios = make_range(7, 11, steps.speed_ratio);
    values.temperatures = make_range(800, 1001, steps.temperature);
    values.angles = make_range(-90, 91, steps.angle);
    values.wing_deflections = make_range(0, 101, steps.wing);
    return values;
}

std::uint64_t total_cases(const SweepValues& values) {
    const auto multiply = [](std::uint64_t a, std::size_t b) {
        return a * static_cast<std::uint64_t>(b);
    };
    return multiply(
        multiply(
            multiply(
                multiply(
                    multiply(values.speed_ratios.size(), values.temperatures.size()),
                    values.angles.size()),
                values.angles.size()),
            values.wing_deflections.size()),
        values.wing_deflections.size());
}

DecodedIndices decode_indices(std::uint64_t index, const SweepValues& values) {
    DecodedIndices indices{};
    auto remaining = index;

    const std::size_t eta_span = values.wing_deflections.size();
    const std::size_t angle_span = values.angles.size();
    const std::size_t temperature_span = values.temperatures.size();

    indices.eta1 = remaining % eta_span;
    remaining /= eta_span;

    indices.eta2 = remaining % eta_span;
    remaining /= eta_span;

    indices.aos = remaining % angle_span;
    remaining /= angle_span;

    indices.aoa = remaining % angle_span;
    remaining /= angle_span;

    indices.temperature = remaining % temperature_span;
    remaining /= temperature_span;

    indices.speed_ratio = remaining % values.speed_ratios.size();
    return indices;
}

Eigen::Vector3d rotate_speed_x(double speed_x, double aoa_deg, double aos_deg) {
    constexpr double deg_to_rad = M_PI / 180.0;
    const double theta_y = aoa_deg * deg_to_rad;
    const double theta_z = aos_deg * deg_to_rad;

    Eigen::Matrix3d Ry;
    Ry << std::cos(theta_y), 0.0, std::sin(theta_y),
          0.0, 1.0, 0.0,
          -std::sin(theta_y), 0.0, std::cos(theta_y);

    Eigen::Matrix3d Rz;
    Rz << std::cos(theta_z), -std::sin(theta_z), 0.0,
          std::sin(theta_z), std::cos(theta_z), 0.0,
          0.0, 0.0, 1.0;

    return Rz * (Ry * Eigen::Vector3d(speed_x, 0.0, 0.0));
}

SweepCache build_sweep_cache(const SweepValues& values, double orbital_speed) {
    SweepCache cache;
    cache.angle_span = values.angles.size();
    cache.velocities.resize(cache.angle_span * cache.angle_span);

    for (std::size_t aoa_idx = 0; aoa_idx < values.angles.size(); ++aoa_idx) {
        for (std::size_t aos_idx = 0; aos_idx < values.angles.size(); ++aos_idx) {
            const std::size_t linear = aoa_idx * cache.angle_span + aos_idx;
            cache.velocities[linear] = rotate_speed_x(orbital_speed,
                                                      static_cast<double>(values.angles[aoa_idx]),
                                                      static_cast<double>(values.angles[aos_idx]));
        }
    }
    return cache;
}

std::array<double, kRowWidth> compute_row(std::uint64_t index, const SweepValues& values,
                                           const SweepCache& cache) {
    using vleo_aerodynamics_core::shuttlecock_aero;
    constexpr double deg_to_rad = M_PI / 180.0;

    const auto indices = decode_indices(index, values);
    const double speed_ratio = static_cast<double>(values.speed_ratios[indices.speed_ratio]);
    const double temperature = static_cast<double>(values.temperatures[indices.temperature]);
    const double eta2_deg = static_cast<double>(values.wing_deflections[indices.eta2]);
    const double eta1_deg = static_cast<double>(values.wing_deflections[indices.eta1]);

    const std::size_t vel_idx = indices.aoa * cache.angle_span + indices.aos;
    const Eigen::Vector3d& velocity = cache.velocities[vel_idx];
    const auto [force, torque] = shuttlecock_aero(
        velocity.x(), velocity.y(), velocity.z(),
        eta1_deg * deg_to_rad,
        eta2_deg * deg_to_rad,
        kEta3Deg * deg_to_rad,
        kEta4Deg * deg_to_rad,
        kRho,
        temperature,
        speed_ratio);

    return {velocity.x(), velocity.y(), velocity.z(),
            eta1_deg, eta2_deg, kEta3Deg, kEta4Deg,
            speed_ratio, temperature,
            force.x(), force.y(), force.z(),
            torque.x(), torque.y(), torque.z()};
}

void write_header(std::ofstream& csv) {
    csv << "Vx,Vy,Vz,eta1_deg,eta2_deg,eta3_deg,eta4_deg,s,T_K,Fx,Fy,Fz,Mx,My,Mz\n";
}

void write_rows(std::ofstream& csv, const std::vector<double>& buffer) {
    csv << std::setprecision(12) << std::scientific;
    for (std::size_t offset = 0; offset < buffer.size(); offset += kRowWidth) {
        for (int c = 0; c < kRowWidth; ++c) {
            csv << buffer[offset + c];
            if (c + 1 < kRowWidth) {
                csv << ',';
            }
        }
        csv << '\n';
    }
}

void process_range(std::uint64_t start, std::uint64_t count,
                   const RunConfig& config, const SweepValues& values,
                   const SweepCache& cache,
                   std::ofstream& csv, bool write_header_flag) {
    if (write_header_flag) {
        write_header(csv);
    }

    if (count == 0) {
        return;
    }

    std::vector<double> buffer;
    buffer.reserve(config.chunk_rows * kRowWidth);

    const std::uint64_t end = start + count;
    for (std::uint64_t index = start; index < end; ++index) {
        const auto row = compute_row(index, values, cache);
        buffer.insert(buffer.end(), row.begin(), row.end());

        if (buffer.size() / kRowWidth >= config.chunk_rows) {
            write_rows(csv, buffer);
            buffer.clear();
        }
    }

    if (!buffer.empty()) {
        write_rows(csv, buffer);
    }
}

WorkRange compute_work_range(std::uint64_t total, int world_size, int rank) {
    WorkRange range{};
    const std::uint64_t base = total / static_cast<std::uint64_t>(world_size);
    const std::uint64_t remainder = total % static_cast<std::uint64_t>(world_size);
    range.count = base + (static_cast<std::uint64_t>(rank) < remainder ? 1 : 0);
    range.start = base * static_cast<std::uint64_t>(rank) +
                  std::min(static_cast<std::uint64_t>(rank), remainder);
    return range;
}

std::string temp_filename(const RunConfig& config, int rank) {
    return config.output_file + ".rank" + std::to_string(rank) + ".tmp";
}

void merge_temp_files(const RunConfig& config, int world_size) {
    std::ofstream final_csv(config.output_file, std::ios::trunc);
    if (!final_csv) {
        throw std::runtime_error("failed to open output file '" + config.output_file + "'");
    }
    write_header(final_csv);

    for (int rank = 0; rank < world_size; ++rank) {
        const std::string tmp = temp_filename(config, rank);
        if (!std::filesystem::exists(tmp)) {
            continue;
        }
        std::ifstream input(tmp);
        if (!input) {
            throw std::runtime_error("failed to read temporary file '" + tmp + "'");
        }
        final_csv << input.rdbuf();
        input.close();
        std::filesystem::remove(tmp);
    }
}

RunConfig parse_arguments(int argc, char** argv, int rank) {
    RunConfig config;
    bool show_help = false;
    std::string error_message;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--full-sweep") {
            config.coarse_sweep = false;
        } else if (arg == "--test-sweep") {
            config.coarse_sweep = true;
        } else if (arg.rfind("--output=", 0) == 0) {
            config.output_file = arg.substr(9);
        } else if (arg.rfind("--chunk-size=", 0) == 0) {
            try {
                config.chunk_rows = std::stoull(arg.substr(13));
                if (config.chunk_rows == 0) {
                    error_message = "chunk size must be greater than zero";
                }
            } catch (const std::exception&) {
                error_message = "invalid chunk size";
            }
        } else if (arg.rfind("--max-cases=", 0) == 0) {
            try {
                config.max_cases = std::stoull(arg.substr(12));
            } catch (const std::exception&) {
                error_message = "invalid max-cases value";
            }
        } else if (arg == "-h" || arg == "--help") {
            show_help = true;
        } else {
            error_message = "unknown option: " + arg;
        }
    }

    if (show_help) {
        if (rank == 0) {
            std::cout << "Usage: mpirun -n <ranks> force_moment_parallel_mpi [options]\n"
                      << "Options:\n"
                      << "  --full-sweep       Run the full parameter sweep (default is test sweep).\n"
                      << "  --test-sweep       Run the reduced 100k sweep.\n"
                      << "  --output=FILE      Write CSV output to FILE (default: results.csv).\n"
                      << "  --chunk-size=N     Rows per batch flush to disk (default: 256).\n"
                      << "  --max-cases=N      Limit the total number of parameter combinations.\n"
                      << "  -h, --help         Show this help message.\n";
        }
        MPI_Finalize();
        std::exit(0);
    }

    if (!error_message.empty()) {
        if (rank == 0) {
            std::cerr << error_message << '\n';
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return config;
}

} // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank = 0;
    int world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    vleo_aerodynamics_core::preload_shuttlecock_geometry();

    const RunConfig config = parse_arguments(argc, argv, world_rank);
    const SweepValues values = build_sweep_values(config);

    if (values.speed_ratios.empty() || values.temperatures.empty() ||
        values.angles.empty() || values.wing_deflections.empty()) {
        if (world_rank == 0) {
            std::cerr << "Invalid sweep definition (one of the ranges is empty)." << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const SweepCache cache = build_sweep_cache(values, config.orbital_speed);

    std::uint64_t total = total_cases(values);
    if (config.max_cases > 0 && config.max_cases < total) {
        total = config.max_cases;
    }
    if (total == 0) {
        if (world_rank == 0) {
            std::cerr << "No parameter combinations to evaluate." << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (world_rank == 0) {
        std::cout << "Evaluating " << total << " combinations using " << world_size
                  << " MPI rank(s)." << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    const auto start_time = std::chrono::steady_clock::now();

    const WorkRange my_range = compute_work_range(total, world_size, world_rank);

    if (world_size == 1) {
        std::ofstream csv(config.output_file, std::ios::trunc);
        if (!csv) {
            std::cerr << "Failed to open output file '" << config.output_file << "'." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        process_range(my_range.start, my_range.count, config, values, cache, csv, true);
    } else {
        if (my_range.count > 0) {
            const std::string tmp_name = temp_filename(config, world_rank);
            std::ofstream tmp(tmp_name, std::ios::trunc);
            if (!tmp) {
                std::cerr << "Failed to open temporary file '" << tmp_name << "'." << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            process_range(my_range.start, my_range.count, config, values, cache, tmp, false);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == 0) {
            try {
                merge_temp_files(config, world_size);
            } catch (const std::exception& ex) {
                std::cerr << ex.what() << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {
        const auto elapsed = std::chrono::steady_clock::now() - start_time;
        const double seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed).count();
        std::cout << "Finished in " << seconds << " s." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
