#include "PoseLib/solvers/p1p2ll.h"
#include "PoseLib/solvers/p2p1ll.h"
#include "problem_generator.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>

// Reference sources are extracted by compare_mixed_pose.py, never linked into PoseLib.
namespace poselib {
int reference_p1p2ll(const std::vector<Eigen::Vector3d> &, const std::vector<Eigen::Vector3d> &,
                     const std::vector<Eigen::Vector3d> &, const std::vector<Eigen::Vector3d> &,
                     const std::vector<Eigen::Vector3d> &, std::vector<CameraPose> *);
int reference_p2p1ll(const std::vector<Eigen::Vector3d> &, const std::vector<Eigen::Vector3d> &,
                     const std::vector<Eigen::Vector3d> &, const std::vector<Eigen::Vector3d> &,
                     const std::vector<Eigen::Vector3d> &, std::vector<CameraPose> *);
} // namespace poselib

namespace {
using namespace poselib;
using Solver = decltype(&p1p2ll);

// Identical optional Newton polishing for both methods to measure the runtime /
// accuracy frontier. These six geometric equations are independent of either solver.
void polish(const AbsolutePoseProblemInstance &x, CameraPose &pose) {
    Eigen::Matrix<double, 6, 6> J;
    Eigen::Matrix<double, 6, 1> residual;
    int row = 0;
    for (size_t i = 0; i < x.X_point_.size(); ++i) {
        const Eigen::Vector3d P = pose.rotate(x.X_point_[i]);
        const Eigen::Vector3d a = x.x_point_[i].unitOrthogonal();
        const Eigen::Vector3d b = x.x_point_[i].cross(a);
        for (const Eigen::Vector3d &n : {a, b}) {
            J.block<1, 3>(row, 0) = P.cross(n).transpose();
            J.block<1, 3>(row, 3) = n.transpose();
            residual(row++) = n.dot(P + pose.t);
        }
    }
    for (size_t i = 0; i < x.X_line_line_.size(); ++i) {
        const Eigen::Vector3d n = x.l_line_line_[i].normalized();
        const Eigen::Vector3d P = pose.rotate(x.X_line_line_[i]);
        const Eigen::Vector3d v = pose.rotate(x.V_line_line_[i]);
        J.block<1, 3>(row, 0) = P.cross(n).transpose();
        J.block<1, 3>(row, 3) = n.transpose();
        residual(row++) = n.dot(P + pose.t);
        J.block<1, 3>(row, 0) = v.cross(n).transpose();
        J.block<1, 3>(row, 3).setZero();
        residual(row++) = n.dot(v);
    }
    const Eigen::Matrix<double, 6, 1> step = J.partialPivLu().solve(-residual);
    if (step.allFinite()) {
        pose.q = quat_multiply(quat_exp(step.head<3>()), pose.q).normalized();
        pose.t += step.tail<3>();
    }
}

void solve(Solver solver, const AbsolutePoseProblemInstance &x, int iterations, CameraPoseVector &solutions) {
    solver(x.x_point_, x.X_point_, x.l_line_line_, x.X_line_line_, x.V_line_line_, &solutions);
    // Compare the same physical solution set. Include reference filtering in
    // timing; production mixed solvers select positive lifts before reconstruction.
    const auto nonpositive = [&](const CameraPose &pose) {
        for (size_t i = 0; i < x.X_point_.size(); ++i)
            if (!(x.x_point_[i].dot(pose.apply(x.X_point_[i])) > 0.0))
                return true;
        return false;
    };
    if (solver == reference_p1p2ll || solver == reference_p2p1ll)
        solutions.erase(std::remove_if(solutions.begin(), solutions.end(), nonpositive), solutions.end());
    for (int iter = 0; iter < iterations; ++iter)
        for (auto &pose : solutions)
            polish(x, pose);
}

// A finite-aware angular residual supplements the existing PoseLib validator,
// whose point residual is quadratic in the small angular error.
double residual(const AbsolutePoseProblemInstance &x, const CameraPose &pose) {
    if (!pose.q.allFinite() || !pose.t.allFinite())
        return std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < x.X_point_.size(); ++i)
        if (!(x.x_point_[i].dot(pose.apply(x.X_point_[i])) > 0.0))
            return std::numeric_limits<double>::infinity();
    double error = std::abs(pose.q.squaredNorm() - 1.0);
    for (size_t i = 0; i < x.X_point_.size(); ++i)
        error = std::max(error, x.x_point_[i].cross(pose.apply(x.X_point_[i]).normalized()).norm());
    for (size_t i = 0; i < x.X_line_line_.size(); ++i) {
        const Eigen::Vector3d n = x.l_line_line_[i].normalized();
        error = std::max(error, std::abs(n.dot(pose.apply(x.X_line_line_[i]).normalized())));
        error = std::max(error, std::abs(n.dot(pose.rotate(x.V_line_line_[i]).normalized())));
    }
    return error;
}

void transform(std::vector<AbsolutePoseProblemInstance> &data, const std::string &scene, unsigned int seed) {
    std::mt19937 random(seed);
    std::normal_distribution<double> normal;
    const auto normal3 = [&]() {
        const double x = normal(random), y = normal(random), z = normal(random);
        return Eigen::Vector3d(x, y, z);
    };
    for (auto &x : data) {
        if (scene == "paper") {
            const Eigen::Vector3d axis = normal3().normalized();
            const double angle = normal(random);
            const Eigen::Matrix3d R = Eigen::AngleAxisd(angle, axis).toRotationMatrix();
            const Eigen::Vector3d center = normal3().normalized();
            x.pose_gt = CameraPose(R, -R * center);
            const Eigen::Vector3d mean(0, 0, 5);
            for (auto &P : x.X_point_)
                P = mean + normal3();
            for (size_t i = 0; i < x.X_line_line_.size(); ++i) {
                x.X_line_line_[i] = mean + normal3();
                const Eigen::Vector3d endpoint = mean + normal3();
                x.V_line_line_[i] = (endpoint - x.X_line_line_[i]).normalized();
            }
        } else if (scene == "coplanar" || scene == "near_coplanar") {
            const double thickness = scene == "coplanar" ? 0.0 : 1e-6;
            for (auto &P : x.X_point_)
                P.z() *= thickness;
            for (auto &P : x.X_line_line_)
                P.z() *= thickness;
            for (auto &v : x.V_line_line_) {
                v.z() *= thickness;
                v.normalize();
            }
        } else if (scene == "incident") {
            Eigen::Vector3d P = x.pose_gt.apply(x.X_point_[0]);
            P -= x.l_line_line_[0].dot(P) * x.l_line_line_[0];
            x.X_point_[0] = x.pose_gt.apply_inverse(P);
        } else if (scene == "large" || scene == "small") {
            const double scale = scene == "large" ? 1e3 : 1e-3;
            const Eigen::Vector3d shift(3, -2, 1);
            for (auto &P : x.X_point_)
                P = scale * (P + shift);
            for (auto &P : x.X_line_line_)
                P = scale * (P + shift);
            x.pose_gt.t = scale * (x.pose_gt.t - x.pose_gt.rotate(shift));
        }
        if (scene != "generic") {
            for (size_t k = 0; k < x.X_point_.size(); ++k)
                x.x_point_[k] = x.pose_gt.apply(x.X_point_[k]).normalized();
            for (size_t k = 0; k < x.X_line_line_.size(); ++k)
                x.l_line_line_[k] =
                    x.pose_gt.apply(x.X_line_line_[k]).cross(x.pose_gt.rotate(x.V_line_line_[k])).normalized();
        }
    }
}

struct Result {
    std::vector<double> errors, residuals, times;
    size_t solutions = 0, valid = 0, found = 0, nonfinite = 0;
};

void measure_accuracy(Solver solver, const std::vector<AbsolutePoseProblemInstance> &data, int iterations, Result &r) {
    CameraPoseVector poses;
    for (const auto &x : data) {
        solve(solver, x, iterations, poses);
        r.solutions += poses.size();
        double error = std::numeric_limits<double>::infinity();
        for (const auto &pose : poses) {
            if (!pose.q.allFinite() || !pose.t.allFinite()) {
                ++r.nonfinite;
            } else {
                r.valid += CalibPoseValidator::is_valid(x, pose, 1.0, 1e-6);
                error = std::min(error, CalibPoseValidator::compute_pose_error(x, pose, 1.0));
            }
            r.residuals.push_back(residual(x, pose));
        }
        r.found += error < 1e-6;
        r.errors.push_back(error);
    }
    std::sort(r.errors.begin(), r.errors.end());
    std::sort(r.residuals.begin(), r.residuals.end());
}

double quantile(const std::vector<double> &v, double q) {
    return v.empty() ? std::numeric_limits<double>::infinity() : v[std::min(v.size() - 1, size_t(q * v.size()))];
}

// Stage timings are deliberately separate from the end-to-end measurement.
// Restore the same unpolished candidates outside the polishing-only timer on
// every batch. No per-pose clocks, allocations, validation, or copies are timed.
volatile double profiling_checksum = 0.0;
void profile_polishing(Solver solver, const std::vector<AbsolutePoseProblemInstance> &data, const std::string &scene,
                       int points, int method, int iterations, unsigned int seed) {
    std::vector<CameraPoseVector> initial(data.size()), working;
    std::srand(seed + 1);
    size_t candidates = 0;
    for (size_t i = 0; i < data.size(); ++i) {
        solve(solver, data[i], 0, initial[i]);
        candidates += initial[i].size();
    }
    std::vector<double> timings[3]; // solver, polishing only, complete pipeline
    CameraPoseVector scratch;
    for (int rep = -1; rep < 10; ++rep) {
        for (int index = 0; index < 3; ++index) {
            const int stage = (index + rep + 3) % 3;
            if (stage == 1)
                working = initial;
            std::srand(seed + 1);
            const auto start = std::chrono::steady_clock::now();
            if (stage == 1) {
                for (size_t i = 0; i < data.size(); ++i)
                    for (int iter = 0; iter < iterations; ++iter)
                        for (auto &pose : working[i])
                            polish(data[i], pose);
            } else {
                for (const auto &x : data)
                    solve(solver, x, stage == 0 ? 0 : iterations, scratch);
            }
            const double ns =
                std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count() /
                data.size();
            if (rep >= 0)
                timings[stage].push_back(ns);
            // Consume all independently polished poses after stopping the clock.
            // The observable checksum prevents dead-store elimination.
            double checksum = 0.0;
            if (stage == 1) {
                for (const auto &poses : working)
                    for (const auto &pose : poses)
                        checksum += pose.q.squaredNorm() + 1e-16 * pose.t.squaredNorm();
            } else {
                for (const auto &pose : scratch)
                    checksum += pose.q.squaredNorm() + 1e-16 * pose.t.squaredNorm();
            }
            profiling_checksum = checksum;
        }
    }
    std::cout << scene << ',' << points << ',' << (method == 0 ? "reference" : "new") << ',' << iterations << ','
              << data.size() << ',' << candidates;
    for (auto &stage : timings) {
        std::sort(stage.begin(), stage.end());
        std::cout << ',' << quantile(stage, 0.5) << ',' << stage.front() << ',' << stage.back();
    }
    std::cout << '\n';
}

void print(const std::string &scene, int points, int method, int iterations, Result &r) {
    std::sort(r.times.begin(), r.times.end());
    const double count = r.errors.size();
    const double mean = std::accumulate(r.errors.begin(), r.errors.end(), 0.0) / count;
    std::cout << scene << ',' << points << ',' << (method == 0 ? "reference" : "new") << ',' << iterations << ','
              << r.errors.size() << ',' << quantile(r.times, 0.5) << ',' << r.times.front() << ',' << r.times.back()
              << ',' << r.solutions / count << ',' << r.valid << ',' << r.solutions << ',' << r.found << ','
              << r.nonfinite << ',' << mean << ',' << quantile(r.errors, 0.5) << ',' << quantile(r.errors, 0.99) << ','
              << r.errors.back() << ',' << quantile(r.residuals, 0.99) << ',' << quantile(r.residuals, 1.0) << '\n';
}
} // namespace

int main(int argc, char **argv) {
    const int count = argc > 1 ? std::stoi(argv[1]) : 100000;
    const unsigned int seed = argc > 2 ? std::stoul(argv[2]) : 1;
    const std::string scene = argc > 3 ? argv[3] : "generic";
    const int max_iterations = argc > 4 ? std::stoi(argv[4]) : 2;
    const bool profiling = argc > 5 && std::string(argv[5]) == "--profile-polishing";
    if (count < 1 || max_iterations < 0 || max_iterations > 10 || (profiling && max_iterations == 0))
        return 1;
    if (scene != "generic" && scene != "paper" && scene != "coplanar" && scene != "near_coplanar" &&
        scene != "incident" && scene != "large" && scene != "small") {
        std::cerr << "Unknown scene: " << scene << '\n';
        return 1;
    }
    std::cout << std::setprecision(12);
    if (profiling)
        std::cout << "scene,points,method,polish,instances,candidates,solver_ns,solver_min_ns,solver_max_ns,"
                     "polishing_ns,polishing_min_ns,polishing_max_ns,pipeline_ns,pipeline_min_ns,pipeline_max_ns\n";
    else
        std::cout
            << "scene,points,method,polish,instances,ns,ns_min,ns_max,solutions_per_instance,valid,solutions,found,"
               "nonfinite,mean_error,median_error,p99_error,max_error,p99_residual,max_residual\n";
    for (int points : {2, 1}) {
        ProblemOptions options;
        options.n_point_point_ = points;
        options.n_line_line_ = 3 - points;
        std::vector<AbsolutePoseProblemInstance> data;
        std::srand(seed);
        generate_abspose_problems(count, &data, options);
        transform(data, scene, seed);
        const Solver solvers[2] = {points == 2 ? reference_p2p1ll : reference_p1p2ll, points == 2 ? p2p1ll : p1p2ll};
        if (profiling) {
            for (int iterations = 1; iterations <= max_iterations; ++iterations)
                for (int method = 0; method < 2; ++method)
                    profile_polishing(solvers[method], data, scene, points, method, iterations, seed);
            std::cout.flush();
            continue;
        }
        for (int iterations = 0; iterations <= max_iterations; ++iterations) {
            Result results[2];
            for (int method = 0; method < 2; ++method) {
                std::srand(seed + 1);
                measure_accuracy(solvers[method], data, iterations, results[method]);
            }
            // Alternate execution order; exclude validation, generation, and I/O.
            // Match the stock benchmark's ten batches and reused output storage.
            CameraPoseVector poses[2];
            for (int rep = -1; rep < 10; ++rep) {
                for (int index = 0; index < 2; ++index) {
                    const int method = (index + (rep & 1)) % 2;
                    const auto start = std::chrono::steady_clock::now();
                    for (const auto &x : data)
                        solve(solvers[method], x, iterations, poses[method]);
                    const double ns =
                        std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count() /
                        count;
                    if (rep >= 0)
                        results[method].times.push_back(ns);
                }
            }
            for (int method = 0; method < 2; ++method)
                print(scene, points, method, iterations, results[method]);
            std::cout.flush();
        }
    }
}
