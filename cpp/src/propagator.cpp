#include "vleo/propagator.hpp"
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <random>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace vleo {

namespace {

constexpr int kStages = 12;
constexpr int kStateSize = DeterministicPropagator::kStateSize;
constexpr int kErrStateSize = kStateSize - 1;
constexpr int kOrder = 8;
constexpr int kBaseStateSize = DeterministicPropagator::kBaseStateSize;
constexpr int kIdxLatentStart = DeterministicPropagator::kBaseStateSize;
constexpr int kIdxLogRhoFast = kIdxLatentStart;
constexpr int kIdxLogRhoBias = kIdxLatentStart + 1;
constexpr int kIdxWind = kIdxLatentStart + 2;
using State = Eigen::Matrix<double, kStateSize, 1>;

const std::array<int, kErrStateSize> kErrToFull = {
    0, 1, 2,
    3, 4, 5,
    7, 8, 9,
    10, 11, 12,
    13, 14, 15, 16, 17
};

const std::array<double, kErrStateSize> kErrScale = {
    1.0, 1.0, 1.0,
    1.0, 1.0, 1.0,
    2.0, 2.0, 2.0,
    1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0
};

const double kC[kStages] = {
    0.0,
    0.05260015195876773,
    0.0789002279381516,
    0.1183503419072274,
    0.2816496580927726,
    0.3333333333333333,
    0.25,
    0.3076923076923077,
    0.6512820512820513,
    0.6,
    0.8571428571428571,
    1.0};

const double kA[kStages][kStages] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.05260015195876773, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0197250569845379, 0.0591751709536137, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.02958758547680685, 0.0, 0.08876275643042054, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.2413651341592667, 0.0, -0.8845494793282861, 0.924834003261792, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.037037037037037035, 0.0, 0.0, 0.17082860872947386, 0.12546768756682242, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.037109375, 0.0, 0.0, 0.17025221101954405, 0.06021653898045596, -0.017578125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.03709200011850479, 0.0, 0.0, 0.17038392571223998, 0.10726203044637328, -0.015319437748624402, 0.008273789163814023, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.6241109587160757, 0.0, 0.0, -3.3608926294469414, -0.868219346841726, 27.59209969944671, 20.154067550477894, -43.48988418106996, 0.0, 0.0, 0.0, 0.0},
    {0.47766253643826434, 0.0, 0.0, -2.4881146199716677, -0.590290826836843, 21.230051448181193, 15.279233632882423, -33.28821096898486, -0.020331201708508627, 0.0, 0.0, 0.0},
    {-0.9371424300859873, 0.0, 0.0, 5.186372428844064, 1.0914373489967295, -8.149787010746927, -18.52006565999696, 22.739487099350505, 2.4936055526796523, -3.0467644718982196, 0.0, 0.0},
    {2.273310147516538, 0.0, 0.0, -10.53449546673725, -2.0008720582248625, -17.9589318631188, 27.94888452941996, -2.8589982771350235, -8.87285693353063, 12.360567175794303, 0.6433927460157636, 0.0}};

const double kB[kStages] = {
    0.054293734116568765,
    0.0,
    0.0,
    0.0,
    0.0,
    4.450312892752409,
    1.8915178993145003,
    -5.801203960010585,
    0.3111643669578199,
    -0.1521609496625161,
    0.20136540080403034,
    0.04471061572777259};

Eigen::Vector4d quat_derivative_BI(const Eigen::Quaterniond& q_BI,
                                   const Eigen::Vector3d& w_BI_B) {
    // q_BI maps inertial vectors into body; sign follows Rdot = -[w]x R convention.
    Eigen::Quaterniond omega(0.0, w_BI_B.x(), w_BI_B.y(), w_BI_B.z());
    Eigen::Quaterniond qdot = q_BI * omega;
    return Eigen::Vector4d(-0.5 * qdot.w(), -0.5 * qdot.x(), -0.5 * qdot.y(), -0.5 * qdot.z());
}

State state_derivative(const State& x,
                       const EnvInputs& env,
                       const VehicleParams& vehicle,
                       const PropagatorConfig& config,
                       const vleo_aerodynamics_core::AeroAdapter& aero) {
    Eigen::Vector3d r_I = x.segment<3>(0);
    Eigen::Vector3d v_I = x.segment<3>(3);
    Eigen::Quaterniond q_BI(x[6], x[7], x[8], x[9]);
    Eigen::Vector3d w_BI_B = x.segment<3>(10);
    Eigen::Vector3d w_eval = w_BI_B;
    if (config.freeze_attitude) {
        w_eval.setZero();
    }

    if (q_BI.norm() == 0.0) {
        throw std::invalid_argument("quaternion must be non-zero");
    }
    q_BI.normalize();

    const double log_rho = x[kIdxLogRhoFast] + x[kIdxLogRhoBias];
    const double density = env.density * std::exp(log_rho);
    const Eigen::Vector3d wind_I = env.wind_I + x.segment<3>(kIdxWind);

    Eigen::Vector3d force_B = Eigen::Vector3d::Zero();
    Eigen::Vector3d torque_B = Eigen::Vector3d::Zero();
    if (density > 0.0) {
        auto ft = aero.computeFT(q_BI, w_eval, v_I, wind_I,
                                 density, env.temperature_K,
                                 env.particles_mass_kg, env.eta1_rad,
                                 env.eta2_rad, env.temperature_ratio_method);
        force_B = ft.first;
        torque_B = ft.second;
    }

    const Eigen::Matrix3d R_BI = q_BI.toRotationMatrix();
    const Eigen::Vector3d force_I = R_BI.transpose() * force_B;

    const double r_norm = r_I.norm();
    if (r_norm == 0.0) {
        throw std::invalid_argument("position magnitude must be non-zero");
    }

    const Eigen::Vector3d accel_grav = -config.mu_earth_m3_s2 * r_I / (r_norm * r_norm * r_norm);
    const Eigen::Vector3d accel_I = accel_grav + force_I / vehicle.mass_kg;

    Eigen::Vector3d wdot_B = Eigen::Vector3d::Zero();
    if (!config.freeze_attitude) {
        const Eigen::Matrix3d inertia_B = vehicle.inertia_B;
        wdot_B = inertia_B.inverse() *
            (torque_B - w_BI_B.cross(inertia_B * w_BI_B));
    }

    State dx = State::Zero();
    dx.segment<3>(0) = v_I;
    dx.segment<3>(3) = accel_I;
    if (!config.freeze_attitude) {
        dx.segment<4>(6) = quat_derivative_BI(q_BI, w_BI_B);
        dx.segment<3>(10) = wdot_B;
    }

    return dx;
}

State step_dop853(const State& x,
                  double dt,
                  const EnvInputs& env,
                  const VehicleParams& vehicle,
                  const PropagatorConfig& config,
                  const vleo_aerodynamics_core::AeroAdapter& aero) {
    std::array<State, kStages> k;
    k[0] = state_derivative(x, env, vehicle, config, aero);

    for (int i = 1; i < kStages; ++i) {
        State x_stage = x;
        for (int j = 0; j < i; ++j) {
            if (kA[i][j] != 0.0) {
                x_stage += dt * kA[i][j] * k[j];
            }
        }
        k[i] = state_derivative(x_stage, env, vehicle, config, aero);
    }

    State x_next = x;
    for (int j = 0; j < kStages; ++j) {
        if (kB[j] != 0.0) {
            x_next += dt * kB[j] * k[j];
        }
    }

    Eigen::Vector4d q_next = x_next.segment<4>(6);
    if (q_next.norm() == 0.0) {
        throw std::invalid_argument("quaternion must be non-zero");
    }
    q_next.normalize();
    x_next.segment<4>(6) = q_next;

    return x_next;
}

double scaled_error(const State& err,
                    const State& ref,
                    const PropagatorConfig& config) {
    double max_norm = 0.0;
    for (int i = 0; i < kStateSize; ++i) {
        const double scale = config.atol + config.rtol * std::abs(ref[i]);
        const double denom = scale > 0.0 ? scale : 1.0;
        const double ratio = std::abs(err[i]) / denom;
        if (ratio > max_norm) {
            max_norm = ratio;
        }
    }
    return max_norm;
}

void apply_ou_update(State& x,
                     double dt,
                     const PropagatorConfig& config,
                     std::mt19937_64* rng) {
    if (dt <= 0.0) {
        return;
    }

    const bool use_noise = rng != nullptr;
    std::normal_distribution<double> dist(0.0, 1.0);

    auto ou_step = [&](double& val, double tau, double sigma) {
        if (tau > 0.0) {
            const double phi = std::exp(-dt / tau);
            val *= phi;
            if (use_noise && sigma > 0.0) {
                const double var = 1.0 - phi * phi;
                if (var > 0.0) {
                    val += sigma * std::sqrt(var) * dist(*rng);
                }
            }
        } else if (use_noise && sigma > 0.0) {
            val += sigma * std::sqrt(dt) * dist(*rng);
        }
    };

    ou_step(x[kIdxLogRhoFast], config.rho_fast_tau_s, config.rho_fast_sigma);
    ou_step(x[kIdxLogRhoBias], config.rho_bias_tau_s, config.rho_bias_sigma);

    for (int i = 0; i < 3; ++i) {
        ou_step(x[kIdxWind + i], config.wind_tau_s, config.wind_sigma);
    }
}

std::string state_to_string(const State& x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << "[";
    for (int i = 0; i < kStateSize; ++i) {
        oss << x[i];
        if (i + 1 < kStateSize) {
            oss << ", ";
        }
    }
    oss << "]";
    return oss.str();
}

State integrate_interval(const State& x0,
                         double t0,
                         double t1,
                         const EnvInputs& env,
                         const VehicleParams& vehicle,
                         const PropagatorConfig& config,
                         const vleo_aerodynamics_core::AeroAdapter& aero,
                         std::mt19937_64* rng,
                         int* steps_out) {
    if (t1 <= t0) {
        if (steps_out) {
            *steps_out = 0;
        }
        return x0;
    }

    State x = x0;
    double t = t0;
    double h = std::min(config.max_step_s, t1 - t0);
    if (h <= 0.0) {
        if (steps_out) {
            *steps_out = 0;
        }
        return x;
    }

    int steps = 0;
    while (t < t1) {
        if (++steps > config.max_steps) {
            throw std::runtime_error("adaptive integrator exceeded max_steps");
        }
        if (t + h > t1) {
            h = t1 - t;
        }

        const State x_full = step_dop853(x, h, env, vehicle, config, aero);
        State x_half = step_dop853(x, 0.5 * h, env, vehicle, config, aero);
        x_half = step_dop853(x_half, 0.5 * h, env, vehicle, config, aero);

        const State err = (x_half - x_full) / (std::pow(2.0, kOrder) - 1.0);
        const double err_norm = scaled_error(err, x_half, config);

        if (err_norm <= 1.0 || h <= config.min_step_s) {
            t += h;
            x = x_half;
            apply_ou_update(x, h, config, rng);

            double factor = config.max_step_factor;
            if (err_norm > 0.0) {
                factor = config.step_safety * std::pow(1.0 / err_norm, 1.0 / (kOrder + 1));
                factor = std::clamp(factor, config.min_step_factor, config.max_step_factor);
            }
            h = std::min(config.max_step_s, h * factor);
        } else {
            double factor = config.step_safety * std::pow(1.0 / err_norm, 1.0 / (kOrder + 1));
            factor = std::clamp(factor, config.min_step_factor, config.max_step_factor);
            h = std::max(config.min_step_s, h * factor);
        }
    }

    if (steps_out) {
        *steps_out = steps;
    }
    return x;
}

void write_state(double* out, int nt, int n, int t_idx, int p_idx, const State& x) {
    const std::size_t offset = (static_cast<std::size_t>(t_idx) * n + p_idx) * kStateSize;
    for (int k = 0; k < kStateSize; ++k) {
        out[offset + static_cast<std::size_t>(k)] = x[k];
    }
}

void normalize_quaternion_state(State& x) {
    Eigen::Vector4d q = x.segment<4>(6);
    if (q.norm() == 0.0) {
        throw std::invalid_argument("quaternion must be non-zero");
    }
    q.normalize();
    x.segment<4>(6) = q;
}

void align_quaternion_sign(State& x, const State& ref) {
    const double dot = x.segment<4>(6).dot(ref.segment<4>(6));
    if (dot < 0.0) {
        x.segment<4>(6) *= -1.0;
    }
}

Eigen::Quaterniond quat_from_small_angle(const Eigen::Vector3d& dtheta) {
    const double angle = dtheta.norm();
    if (angle < 1e-12) {
        Eigen::Quaterniond dq(1.0,
                              0.5 * dtheta.x(),
                              0.5 * dtheta.y(),
                              0.5 * dtheta.z());
        dq.normalize();
        return dq;
    }
    const Eigen::Vector3d axis = dtheta / angle;
    const double half = 0.5 * angle;
    const double s = std::sin(half);
    return Eigen::Quaterniond(std::cos(half),
                              axis.x() * s,
                              axis.y() * s,
                              axis.z() * s);
}

Eigen::Vector3d quat_error_vector(const Eigen::Quaterniond& q,
                                  const Eigen::Quaterniond& q_ref) {
    Eigen::Quaterniond dq = q * q_ref.conjugate();
    if (dq.w() < 0.0) {
        dq.coeffs() *= -1.0;
    }
    const double w = std::clamp(dq.w(), -1.0, 1.0);
    const double angle = 2.0 * std::acos(w);
    const double s = std::sqrt(std::max(0.0, 1.0 - w * w));
    if (s < 1e-8) {
        return 2.0 * dq.vec();
    }
    return angle * dq.vec() / s;
}

Eigen::MatrixXd full_to_error_cov(const Eigen::MatrixXd& P_full) {
    Eigen::MatrixXd P_err(kErrStateSize, kErrStateSize);
    for (int i = 0; i < kErrStateSize; ++i) {
        const int fi = kErrToFull[i];
        const double si = kErrScale[i];
        for (int j = 0; j < kErrStateSize; ++j) {
            const int fj = kErrToFull[j];
            const double sj = kErrScale[j];
            P_err(i, j) = si * sj * P_full(fi, fj);
        }
    }
    return P_err;
}

Eigen::Matrix<double, kStateSize, kStateSize> error_to_full_cov(const Eigen::MatrixXd& P_err) {
    Eigen::Matrix<double, kStateSize, kStateSize> P_full =
        Eigen::Matrix<double, kStateSize, kStateSize>::Zero();
    for (int i = 0; i < kErrStateSize; ++i) {
        const int fi = kErrToFull[i];
        const double si = kErrScale[i];
        for (int j = 0; j < kErrStateSize; ++j) {
            const int fj = kErrToFull[j];
            const double sj = kErrScale[j];
            const double scale = (si * sj);
            if (scale > 0.0) {
                P_full(fi, fj) = P_err(i, j) / scale;
            }
        }
    }
    return P_full;
}

Eigen::Vector4d average_quaternion_markley(const std::vector<State>& sigma_states,
                                           double w0m,
                                           double wi,
                                           const Eigen::Vector4d& ref) {
    Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
    for (std::size_t i = 0; i < sigma_states.size(); ++i) {
        Eigen::Vector4d q = sigma_states[i].segment<4>(6);
        if (q.norm() == 0.0) {
            continue;
        }
        q.normalize();
        const double w = (i == 0) ? w0m : wi;
        M += w * (q * q.transpose());
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(M);
    if (solver.info() != Eigen::Success) {
        return ref;
    }
    Eigen::Vector4d q_mean = solver.eigenvectors().col(3);
    if (q_mean.dot(ref) < 0.0) {
        q_mean *= -1.0;
    }
    if (q_mean.norm() == 0.0) {
        return ref;
    }
    q_mean.normalize();
    return q_mean;
}

Eigen::MatrixXd expand_covariance(const Eigen::MatrixXd& P,
                                  int target_size) {
    if (P.rows() == target_size && P.cols() == target_size) {
        return P;
    }
    if (P.rows() != P.cols()) {
        throw std::invalid_argument("P0 must be square");
    }
    if (P.rows() > target_size) {
        throw std::invalid_argument("P0 larger than target size");
    }
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(target_size, target_size);
    out.topLeftCorner(P.rows(), P.cols()) = P;
    return out;
}

Eigen::MatrixXd robust_cholesky(const Eigen::MatrixXd& P, double scale) {
    Eigen::MatrixXd sym = 0.5 * (P + P.transpose());
    Eigen::MatrixXd A = scale * sym;
    double jitter = 0.0;
    for (int attempt = 0; attempt < 6; ++attempt) {
        Eigen::LLT<Eigen::MatrixXd> llt(A);
        if (llt.info() == Eigen::Success) {
            return llt.matrixL();
        }
        const double diag_max = A.diagonal().cwiseAbs().maxCoeff();
        jitter = (jitter == 0.0 ? 1e-10 : jitter * 10.0) * std::max(1.0, diag_max);
        A = scale * sym + jitter * Eigen::MatrixXd::Identity(P.rows(), P.cols());
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(sym);
    if (solver.info() != Eigen::Success) {
        throw std::invalid_argument("P0 is not positive definite");
    }
    Eigen::VectorXd evals = solver.eigenvalues();
    const double max_eval = evals.maxCoeff();
    const double eps = 1e-12 * std::max(1.0, max_eval);
    evals = evals.cwiseMax(eps);
    Eigen::MatrixXd sqrt_d = evals.cwiseSqrt().asDiagonal();
    return std::sqrt(scale) * solver.eigenvectors() * sqrt_d;
}

Eigen::MatrixXd stabilize_covariance(const Eigen::MatrixXd& P) {
    Eigen::MatrixXd sym = 0.5 * (P + P.transpose());
    if (!sym.allFinite()) {
        throw std::runtime_error("covariance contains non-finite values");
    }
    double jitter = 0.0;
    const double diag_max = sym.diagonal().cwiseAbs().maxCoeff();
    for (int attempt = 0; attempt < 6; ++attempt) {
        Eigen::MatrixXd A = sym + jitter * Eigen::MatrixXd::Identity(P.rows(), P.cols());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A);
        if (solver.info() == Eigen::Success) {
            Eigen::VectorXd evals = solver.eigenvalues();
            const double max_eval = evals.maxCoeff();
            const double eps = 1e-12 * std::max(1.0, max_eval);
            evals = evals.cwiseMax(eps);
            return solver.eigenvectors() * evals.asDiagonal() * solver.eigenvectors().transpose();
        }
        jitter = (jitter == 0.0 ? 1e-12 : jitter * 10.0) * std::max(1.0, diag_max);
    }
    Eigen::MatrixXd L = robust_cholesky(sym, 1.0);
    return L * L.transpose();
}

struct NoiseSpec {
    std::vector<int> indices;
    Eigen::VectorXd variances;
};

NoiseSpec build_noise_spec(double dt, const PropagatorConfig& config) {
    NoiseSpec spec;
    std::vector<double> vars;

    auto push_var = [&](double tau, double sigma, int idx) {
        if (sigma <= 0.0) {
            return;
        }
        double var = 0.0;
        if (tau > 0.0) {
            const double phi = std::exp(-dt / tau);
            var = sigma * sigma * (1.0 - phi * phi);
        } else {
            var = sigma * sigma * dt;
        }
        if (var > 0.0) {
            spec.indices.push_back(idx);
            vars.push_back(var);
        }
    };

    push_var(config.rho_fast_tau_s, config.rho_fast_sigma, kIdxLogRhoFast);
    push_var(config.rho_bias_tau_s, config.rho_bias_sigma, kIdxLogRhoBias);
    push_var(config.wind_tau_s, config.wind_sigma, kIdxWind + 0);
    push_var(config.wind_tau_s, config.wind_sigma, kIdxWind + 1);
    push_var(config.wind_tau_s, config.wind_sigma, kIdxWind + 2);

    spec.variances = Eigen::VectorXd::Zero(static_cast<int>(vars.size()));
    for (int i = 0; i < static_cast<int>(vars.size()); ++i) {
        spec.variances[i] = vars[i];
    }
    return spec;
}

void apply_ou_noise_only(State& x,
                         double dt,
                         const PropagatorConfig& config,
                         const Eigen::VectorXd& noise,
                         const NoiseSpec& spec) {
    if (spec.indices.empty()) {
        return;
    }
    if (noise.size() != static_cast<int>(spec.indices.size())) {
        throw std::invalid_argument("noise vector dimension mismatch");
    }

    auto add_noise = [&](int idx, double tau, double sigma, double w) {
        if (sigma <= 0.0) {
            return;
        }
        double coeff = 0.0;
        if (tau > 0.0) {
            const double phi = std::exp(-dt / tau);
            coeff = sigma * std::sqrt(1.0 - phi * phi);
        } else {
            coeff = sigma * std::sqrt(dt);
        }
        x[idx] += coeff * w;
    };

    int cursor = 0;
    if (config.rho_fast_sigma > 0.0) {
        add_noise(kIdxLogRhoFast, config.rho_fast_tau_s, config.rho_fast_sigma, noise[cursor++]);
    }
    if (config.rho_bias_sigma > 0.0) {
        add_noise(kIdxLogRhoBias, config.rho_bias_tau_s, config.rho_bias_sigma, noise[cursor++]);
    }
    if (config.wind_sigma > 0.0) {
        add_noise(kIdxWind + 0, config.wind_tau_s, config.wind_sigma, noise[cursor++]);
        add_noise(kIdxWind + 1, config.wind_tau_s, config.wind_sigma, noise[cursor++]);
        add_noise(kIdxWind + 2, config.wind_tau_s, config.wind_sigma, noise[cursor++]);
    }
}

} // namespace

DeterministicPropagator::DeterministicPropagator(
    std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
    const VehicleParams& vehicle,
    const PropagatorConfig& config)
    : aero_(std::move(aero)), vehicle_(vehicle), config_(config) {
    if (!aero_) {
        throw std::invalid_argument("aero adapter must be non-null");
    }
}

void DeterministicPropagator::propagate(const Eigen::VectorXd& x0,
                                        const Eigen::VectorXd& t_grid,
                                        const std::vector<EnvInputs>& env,
                                        Eigen::MatrixXd* X_out) const {
    if (!X_out) {
        throw std::invalid_argument("X_out must be non-null");
    }
    if (x0.size() != kBaseStateSize && x0.size() != kStateSize) {
        throw std::invalid_argument("x0 size must be 13 or 18");
    }
    if (t_grid.size() < 2) {
        throw std::invalid_argument("t_grid must have at least 2 entries");
    }
    if (static_cast<int>(env.size()) != t_grid.size()) {
        throw std::invalid_argument("env size must match t_grid size");
    }

    const int nt = static_cast<int>(t_grid.size());
    X_out->resize(nt, kStateSize);

    State x = State::Zero();
    if (x0.size() == kBaseStateSize) {
        x.segment(0, kBaseStateSize) = x0;
    } else {
        x = x0;
    }
    X_out->row(0) = x.transpose();

    const bool use_noise = (config_.rho_fast_sigma > 0.0 ||
                            config_.rho_bias_sigma > 0.0 ||
                            config_.wind_sigma > 0.0);
    std::mt19937_64 rng(config_.rng_seed);
    std::mt19937_64* rng_ptr = use_noise ? &rng : nullptr;

    if (config_.verbose) {
        std::cout << "[vleo] det start nt=" << nt << std::endl;
    }

    for (int i = 0; i < nt - 1; ++i) {
        const double dt = t_grid[i + 1] - t_grid[i];
        if (dt <= 0.0) {
            throw std::invalid_argument("t_grid must be strictly increasing");
        }

        const EnvInputs& env_i = env[i];
        int substeps = 0;
        int* steps_ptr = config_.debug_state ? &substeps : nullptr;
        x = integrate_interval(x, t_grid[i], t_grid[i + 1], env_i, vehicle_, config_, *aero_, rng_ptr, steps_ptr);

        X_out->row(i + 1) = x.transpose();

        if (config_.verbose) {
            const bool last = (i == nt - 2);
            if (config_.progress_stride > 0 && (i % config_.progress_stride == 0 || last)) {
                std::cout << std::fixed << std::setprecision(2)
                          << "[vleo] det step " << (i + 1) << "/" << (nt - 1)
                          << " t=" << t_grid[i + 1] << " dt=" << dt << std::endl;
            }
        }
        if (config_.debug_state) {
            const bool last = (i == nt - 2);
            if (config_.debug_stride > 0 && (i % config_.debug_stride == 0 || last)) {
                std::cout << "[vleo] det debug step=" << (i + 1)
                          << " substeps=" << substeps
                          << " x=" << state_to_string(x) << std::endl;
            }
        }
    }
}

EnsemblePropagatorMC::EnsemblePropagatorMC(
    std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
    const VehicleParams& vehicle,
    const PropagatorConfig& config)
    : aero_(std::move(aero)), vehicle_(vehicle), config_(config) {
    if (!aero_) {
        throw std::invalid_argument("aero adapter must be non-null");
    }
}

void EnsemblePropagatorMC::propagate(const Eigen::MatrixXd& X0,
                                     const Eigen::VectorXd& t_grid,
                                     const std::vector<EnvInputs>& env,
                                     double* out) const {
    if (!out) {
        throw std::invalid_argument("out must be non-null");
    }
    if (X0.cols() != kBaseStateSize && X0.cols() != kStateSize) {
        throw std::invalid_argument("X0 must have 13 or 18 columns");
    }
    if (t_grid.size() < 2) {
        throw std::invalid_argument("t_grid must have at least 2 entries");
    }
    if (static_cast<int>(env.size()) != t_grid.size()) {
        throw std::invalid_argument("env size must match t_grid size");
    }

    const int nt = static_cast<int>(t_grid.size());
    const int n = static_cast<int>(X0.rows());

    for (int i = 0; i < nt - 1; ++i) {
        const double dt = t_grid[i + 1] - t_grid[i];
        if (dt <= 0.0) {
            throw std::invalid_argument("t_grid must be strictly increasing");
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int p = 0; p < n; ++p) {
        const bool do_particle_log =
            config_.verbose && (config_.particle_stride > 0) && (p % config_.particle_stride == 0);
#ifdef _OPENMP
#pragma omp critical
#endif
        if (do_particle_log) {
            std::cout << "[vleo] mc particle " << (p + 1) << "/" << n << " start" << std::endl;
        }

        const bool do_step_log = config_.verbose && (p == 0);
        const bool do_debug = config_.debug_state &&
                              (config_.particle_stride > 0) &&
                              (p % config_.particle_stride == 0);
        State x = State::Zero();
        if (X0.cols() == kBaseStateSize) {
            x.segment(0, kBaseStateSize) = X0.row(p).transpose();
        } else {
            x = X0.row(p).transpose();
        }
        write_state(out, nt, n, 0, p, x);

        const bool use_noise = (config_.rho_fast_sigma > 0.0 ||
                                config_.rho_bias_sigma > 0.0 ||
                                config_.wind_sigma > 0.0);
        const std::uint64_t seed = config_.rng_seed + static_cast<std::uint64_t>(p) * 0x9e3779b97f4a7c15ULL;
        std::mt19937_64 rng(seed);
        std::mt19937_64* rng_ptr = use_noise ? &rng : nullptr;

        for (int i = 0; i < nt - 1; ++i) {
            const EnvInputs& env_i = env[i];
            int substeps = 0;
            int* steps_ptr = (config_.debug_state && do_debug) ? &substeps : nullptr;
            x = integrate_interval(x, t_grid[i], t_grid[i + 1], env_i, vehicle_, config_, *aero_, rng_ptr, steps_ptr);
            write_state(out, nt, n, i + 1, p, x);

            if (do_step_log) {
                const bool last = (i == nt - 2);
                if (config_.progress_stride > 0 && (i % config_.progress_stride == 0 || last)) {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        const double dt = t_grid[i + 1] - t_grid[i];
                        std::cout << std::fixed << std::setprecision(2)
                                  << "[vleo] mc step " << (i + 1) << "/" << (nt - 1)
                                  << " t=" << t_grid[i + 1] << " dt=" << dt << std::endl;
                    }
                }
            }

            if (do_debug) {
                const bool last = (i == nt - 2);
                if (config_.debug_stride > 0 && (i % config_.debug_stride == 0 || last)) {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        std::cout << "[vleo] mc debug p=" << (p + 1)
                                  << " step=" << (i + 1)
                                  << " substeps=" << substeps
                                  << " x=" << state_to_string(x) << std::endl;
                    }
                }
            }
        }

#ifdef _OPENMP
#pragma omp critical
#endif
        if (do_particle_log) {
            std::cout << "[vleo] mc particle " << (p + 1) << "/" << n << " done" << std::endl;
        }
    }
}

SigmaPointPropagatorUT::SigmaPointPropagatorUT(
    std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
    const VehicleParams& vehicle,
    const PropagatorConfig& config,
    double alpha,
    double beta,
    double kappa)
    : aero_(std::move(aero)),
      vehicle_(vehicle),
      config_(config),
      alpha_(alpha),
      beta_(beta),
      kappa_(kappa) {
    if (!aero_) {
        throw std::invalid_argument("aero adapter must be non-null");
    }
    if (alpha_ <= 0.0) {
        throw std::invalid_argument("alpha must be positive");
    }
}

void SigmaPointPropagatorUT::propagate(const Eigen::VectorXd& x0,
                                       const Eigen::MatrixXd& P0,
                                       const Eigen::VectorXd& t_grid,
                                       const std::vector<EnvInputs>& env,
                                       double* mean_out,
                                       double* cov_out,
                                       bool augment_process_noise) const {
    if (!mean_out || !cov_out) {
        throw std::invalid_argument("mean_out and cov_out must be non-null");
    }
    if (x0.size() != kBaseStateSize && x0.size() != kStateSize) {
        throw std::invalid_argument("x0 size must be 13 or 18");
    }
    if (t_grid.size() < 2) {
        throw std::invalid_argument("t_grid must have at least 2 entries");
    }
    if (static_cast<int>(env.size()) != t_grid.size()) {
        throw std::invalid_argument("env size must match t_grid size");
    }
    if (P0.rows() != P0.cols()) {
        throw std::invalid_argument("P0 must be square");
    }
    if (P0.rows() != x0.size()) {
        throw std::invalid_argument("P0 size must match x0 size");
    }

    const int n_state = kStateSize;
    const int nt = static_cast<int>(t_grid.size());

    State mean = State::Zero();
    if (x0.size() == kBaseStateSize) {
        mean.segment(0, kBaseStateSize) = x0;
    } else {
        mean = x0;
    }
    normalize_quaternion_state(mean);

    Eigen::MatrixXd P_full = expand_covariance(P0, n_state);
    P_full = 0.5 * (P_full + P_full.transpose());
    if (config_.freeze_attitude) {
        P_full.block(6, 0, 7, n_state).setZero();
        P_full.block(0, 6, n_state, 7).setZero();
    }

    auto write_outputs = [&](int t_idx, const State& x_mean, const Eigen::MatrixXd& cov) {
        const std::size_t mean_offset = static_cast<std::size_t>(t_idx) * n_state;
        for (int i = 0; i < n_state; ++i) {
            mean_out[mean_offset + static_cast<std::size_t>(i)] = x_mean[i];
        }
        const std::size_t cov_offset = static_cast<std::size_t>(t_idx) * n_state * n_state;
        for (int i = 0; i < n_state; ++i) {
            for (int j = 0; j < n_state; ++j) {
                cov_out[cov_offset + static_cast<std::size_t>(i) * n_state + j] = cov(i, j);
            }
        }
    };

    write_outputs(0, mean, P_full);

    if (config_.verbose) {
        std::cout << "[vleo] ut start nt=" << nt << " augment=" << (augment_process_noise ? "true" : "false")
                  << std::endl;
    }

    for (int step = 0; step < nt - 1; ++step) {
        const double dt = t_grid[step + 1] - t_grid[step];
        if (dt <= 0.0) {
            throw std::invalid_argument("t_grid must be strictly increasing");
        }

        Eigen::MatrixXd P_err = full_to_error_cov(P_full);
        if (config_.freeze_attitude) {
            P_err.block(6, 0, 6, kErrStateSize).setZero();
            P_err.block(0, 6, kErrStateSize, 6).setZero();
        }
        P_err = stabilize_covariance(P_err);
        const Eigen::Quaterniond q_mean_ref(mean[6], mean[7], mean[8], mean[9]);

        NoiseSpec noise_spec;
        if (augment_process_noise) {
            noise_spec = build_noise_spec(dt, config_);
        }
        const int n_noise = static_cast<int>(noise_spec.indices.size());
        const int n_aug = kErrStateSize + n_noise;

        Eigen::VectorXd mean_aug = Eigen::VectorXd::Zero(n_aug);

        Eigen::MatrixXd P_aug = Eigen::MatrixXd::Zero(n_aug, n_aug);
        P_aug.topLeftCorner(kErrStateSize, kErrStateSize) = P_err;
        if (n_noise > 0) {
            P_aug.bottomRightCorner(n_noise, n_noise) = noise_spec.variances.asDiagonal();
        }

        const double lambda = alpha_ * alpha_ * (n_aug + kappa_) - n_aug;
        const double scale = n_aug + lambda;
        if (scale <= 0.0) {
            throw std::invalid_argument("alpha/kappa produce non-positive scaling");
        }

        const Eigen::MatrixXd L = robust_cholesky(P_aug, scale);
        const int n_sigma = 2 * n_aug + 1;

        const double w0m = lambda / scale;
        const double w0c = w0m + (1.0 - alpha_ * alpha_ + beta_);
        const double wi = 1.0 / (2.0 * scale);

        std::vector<State> sigma_states(n_sigma);
        std::vector<int> sigma_steps;
        if (config_.debug_state) {
            sigma_steps.assign(n_sigma, 0);
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int s = 0; s < n_sigma; ++s) {
            Eigen::VectorXd sigma = mean_aug;
            if (s > 0) {
                const int idx = (s - 1) % n_aug;
                const double sign = (s <= n_aug) ? 1.0 : -1.0;
                sigma += sign * L.col(idx);
            }

            const Eigen::VectorXd err = sigma.head(kErrStateSize);
            State x_state = mean;
            x_state.segment<3>(0) += err.segment<3>(0);
            x_state.segment<3>(3) += err.segment<3>(3);
            if (!config_.freeze_attitude) {
                const Eigen::Vector3d dtheta = err.segment<3>(6);
                const Eigen::Quaterniond dq = quat_from_small_angle(dtheta);
                const Eigen::Quaterniond q_sigma = q_mean_ref * dq;
                x_state[6] = q_sigma.w();
                x_state[7] = q_sigma.x();
                x_state[8] = q_sigma.y();
                x_state[9] = q_sigma.z();
            }
            x_state.segment<3>(10) += err.segment<3>(9);
            x_state.segment<5>(kIdxLatentStart) += err.segment<5>(12);
            normalize_quaternion_state(x_state);

            int* steps_ptr = config_.debug_state ? &sigma_steps[s] : nullptr;
            x_state = integrate_interval(
                x_state, t_grid[step], t_grid[step + 1], env[step], vehicle_, config_, *aero_, nullptr, steps_ptr);

            if (n_noise > 0) {
                const Eigen::VectorXd noise = sigma.tail(n_noise);
                apply_ou_noise_only(x_state, dt, config_, noise, noise_spec);
            }
            sigma_states[s] = x_state;
        }

        State mean_next = w0m * sigma_states[0];
        for (int s = 1; s < n_sigma; ++s) {
            mean_next += wi * sigma_states[s];
        }
        if (config_.freeze_attitude) {
            mean_next.segment<4>(6) = mean.segment<4>(6);
            mean_next.segment<3>(10) = mean.segment<3>(10);
        } else {
            const Eigen::Vector4d q_mean = average_quaternion_markley(
                sigma_states, w0m, wi, mean.segment<4>(6));
            mean_next.segment<4>(6) = q_mean;
        }
        normalize_quaternion_state(mean_next);

        auto sigma_diff = [&](int s) -> Eigen::VectorXd {
            Eigen::VectorXd diff = Eigen::VectorXd::Zero(kErrStateSize);
            diff.segment<3>(0) = sigma_states[s].segment<3>(0) - mean_next.segment<3>(0);
            diff.segment<3>(3) = sigma_states[s].segment<3>(3) - mean_next.segment<3>(3);
            diff.segment<3>(9) = sigma_states[s].segment<3>(10) - mean_next.segment<3>(10);
            diff.segment<5>(12) = sigma_states[s].segment<5>(kIdxLatentStart) -
                                  mean_next.segment<5>(kIdxLatentStart);
            if (!config_.freeze_attitude) {
                const Eigen::Quaterniond q_sigma(sigma_states[s][6],
                                                 sigma_states[s][7],
                                                 sigma_states[s][8],
                                                 sigma_states[s][9]);
                const Eigen::Quaterniond q_ref(mean_next[6],
                                               mean_next[7],
                                               mean_next[8],
                                               mean_next[9]);
                diff.segment<3>(6) = quat_error_vector(q_sigma, q_ref);
            }
            return diff;
        };

        Eigen::MatrixXd P_err_next =
            w0c * sigma_diff(0) * sigma_diff(0).transpose();
        for (int s = 1; s < n_sigma; ++s) {
            const Eigen::VectorXd diff = sigma_diff(s);
            P_err_next += wi * diff * diff.transpose();
        }

        mean = mean_next;
        P_err_next = stabilize_covariance(P_err_next);
        P_full = error_to_full_cov(P_err_next);
        write_outputs(step + 1, mean, P_full);

        if (config_.verbose) {
            const bool last = (step == nt - 2);
            if (config_.progress_stride > 0 && (step % config_.progress_stride == 0 || last)) {
                std::cout << std::fixed << std::setprecision(2)
                          << "[vleo] ut step " << (step + 1) << "/" << (nt - 1)
                          << " t=" << t_grid[step + 1] << " dt=" << dt << std::endl;
            }
        }
        if (config_.debug_state) {
            const bool last = (step == nt - 2);
            if (config_.debug_stride > 0 && (step % config_.debug_stride == 0 || last)) {
                int min_steps = 0;
                int max_steps = 0;
                double mean_steps = 0.0;
                if (!sigma_steps.empty()) {
                    min_steps = sigma_steps[0];
                    max_steps = sigma_steps[0];
                    int sum = 0;
                    for (int val : sigma_steps) {
                        min_steps = std::min(min_steps, val);
                        max_steps = std::max(max_steps, val);
                        sum += val;
                    }
                    mean_steps = static_cast<double>(sum) / static_cast<double>(sigma_steps.size());
                }
                std::cout << "[vleo] ut debug step=" << (step + 1)
                          << " substeps_min=" << min_steps
                          << " substeps_mean=" << std::fixed << std::setprecision(1) << mean_steps
                          << " substeps_max=" << max_steps
                          << " x=" << state_to_string(mean) << std::endl;
            }
        }
    }
}

StmPropagator::StmPropagator(std::shared_ptr<vleo_aerodynamics_core::AeroAdapter> aero,
                             const VehicleParams& vehicle,
                             const PropagatorConfig& config,
                             double fd_eps_rel,
                             double fd_eps_abs)
    : aero_(std::move(aero)),
      vehicle_(vehicle),
      config_(config),
      fd_eps_rel_(fd_eps_rel),
      fd_eps_abs_(fd_eps_abs) {
    if (!aero_) {
        throw std::invalid_argument("aero adapter must be non-null");
    }
    if (fd_eps_rel_ <= 0.0 || fd_eps_abs_ <= 0.0) {
        throw std::invalid_argument("finite-difference eps must be positive");
    }
}

void StmPropagator::propagate(const Eigen::VectorXd& x0,
                              const Eigen::MatrixXd& P0,
                              const Eigen::VectorXd& t_grid,
                              const std::vector<EnvInputs>& env,
                              double* mean_out,
                              double* cov_out,
                              double* stm_out) const {
    if (!mean_out || !cov_out) {
        throw std::invalid_argument("mean_out and cov_out must be non-null");
    }
    if (x0.size() != kBaseStateSize && x0.size() != kStateSize) {
        throw std::invalid_argument("x0 size must be 13 or 18");
    }
    if (t_grid.size() < 2) {
        throw std::invalid_argument("t_grid must have at least 2 entries");
    }
    if (static_cast<int>(env.size()) != t_grid.size()) {
        throw std::invalid_argument("env size must match t_grid size");
    }
    if (P0.rows() != P0.cols()) {
        throw std::invalid_argument("P0 must be square");
    }
    if (P0.rows() != x0.size()) {
        throw std::invalid_argument("P0 size must match x0 size");
    }

    const int n_state = kStateSize;
    const int nt = static_cast<int>(t_grid.size());

    State mean = State::Zero();
    if (x0.size() == kBaseStateSize) {
        mean.segment(0, kBaseStateSize) = x0;
    } else {
        mean = x0;
    }
    normalize_quaternion_state(mean);

    Eigen::MatrixXd P_init = expand_covariance(P0, n_state);
    P_init = 0.5 * (P_init + P_init.transpose());
    Eigen::Matrix<double, kStateSize, kStateSize> P =
        P_init.cast<double>();

    Eigen::Matrix<double, kStateSize, kStateSize> Phi =
        Eigen::Matrix<double, kStateSize, kStateSize>::Identity();

    auto write_outputs = [&](int t_idx,
                             const State& x_mean,
                             const Eigen::Matrix<double, kStateSize, kStateSize>& cov,
                             const Eigen::Matrix<double, kStateSize, kStateSize>* stm) {
        const std::size_t mean_offset = static_cast<std::size_t>(t_idx) * n_state;
        for (int i = 0; i < n_state; ++i) {
            mean_out[mean_offset + static_cast<std::size_t>(i)] = x_mean[i];
        }
        const std::size_t cov_offset = static_cast<std::size_t>(t_idx) * n_state * n_state;
        for (int i = 0; i < n_state; ++i) {
            for (int j = 0; j < n_state; ++j) {
                cov_out[cov_offset + static_cast<std::size_t>(i) * n_state + j] = cov(i, j);
            }
        }
        if (stm) {
            const std::size_t stm_offset = static_cast<std::size_t>(t_idx) * n_state * n_state;
            for (int i = 0; i < n_state; ++i) {
                for (int j = 0; j < n_state; ++j) {
                    stm_out[stm_offset + static_cast<std::size_t>(i) * n_state + j] = (*stm)(i, j);
                }
            }
        }
    };

    write_outputs(0, mean, P, stm_out ? &Phi : nullptr);

    if (config_.verbose) {
        std::cout << "[vleo] stm start nt=" << nt << std::endl;
    }

    for (int step = 0; step < nt - 1; ++step) {
        const double dt = t_grid[step + 1] - t_grid[step];
        if (dt <= 0.0) {
            throw std::invalid_argument("t_grid must be strictly increasing");
        }

        const EnvInputs& env_i = env[step];
        int substeps = 0;
        int* steps_ptr = config_.debug_state ? &substeps : nullptr;
        const State x_next = integrate_interval(mean, t_grid[step], t_grid[step + 1],
                                                env_i, vehicle_, config_, *aero_, nullptr, steps_ptr);

        Eigen::Matrix<double, kStateSize, kStateSize> F;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j < n_state; ++j) {
            const double scale = std::max(1.0, std::abs(mean[j]));
            const double eps = fd_eps_abs_ + fd_eps_rel_ * scale;

            State x_plus = mean;
            x_plus[j] += eps;
            normalize_quaternion_state(x_plus);

            State x_minus = mean;
            x_minus[j] -= eps;
            normalize_quaternion_state(x_minus);

            State x_plus_next = integrate_interval(x_plus, t_grid[step], t_grid[step + 1],
                                                   env_i, vehicle_, config_, *aero_, nullptr, nullptr);
            State x_minus_next = integrate_interval(x_minus, t_grid[step], t_grid[step + 1],
                                                    env_i, vehicle_, config_, *aero_, nullptr, nullptr);

            align_quaternion_sign(x_plus_next, x_next);
            align_quaternion_sign(x_minus_next, x_next);

            F.col(j) = (x_plus_next - x_minus_next) / (2.0 * eps);
        }

        Eigen::Matrix<double, kStateSize, kStateSize> Q =
            Eigen::Matrix<double, kStateSize, kStateSize>::Zero();
        const NoiseSpec noise_spec = build_noise_spec(dt, config_);
        for (int i = 0; i < static_cast<int>(noise_spec.indices.size()); ++i) {
            const int idx = noise_spec.indices[i];
            Q(idx, idx) = noise_spec.variances[i];
        }

        Eigen::Matrix<double, kStateSize, kStateSize> P_next = F * P * F.transpose() + Q;
        P_next = 0.5 * (P_next + P_next.transpose());

        if (stm_out) {
            Phi = F * Phi;
        }

        mean = x_next;
        P = P_next;

        write_outputs(step + 1, mean, P, stm_out ? &Phi : nullptr);

        if (config_.verbose) {
            const bool last = (step == nt - 2);
            if (config_.progress_stride > 0 && (step % config_.progress_stride == 0 || last)) {
                std::cout << std::fixed << std::setprecision(2)
                          << "[vleo] stm step " << (step + 1) << "/" << (nt - 1)
                          << " t=" << t_grid[step + 1] << " dt=" << dt << std::endl;
            }
        }
        if (config_.debug_state) {
            const bool last = (step == nt - 2);
            if (config_.debug_stride > 0 && (step % config_.debug_stride == 0 || last)) {
                std::cout << "[vleo] stm debug step=" << (step + 1)
                          << " substeps=" << substeps
                          << " x=" << state_to_string(mean) << std::endl;
            }
        }
    }
}

} // namespace vleo
