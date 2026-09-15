#include "PoseLib/misc/univariate.h"
#include "PoseLib/solvers/p3p_common.h"
#include "test.h"

using namespace poselib;

namespace {
bool test_projective_quadratic() {
    // Zero roots, repeated roots, linear polynomials, and roots at infinity.
    for (double scale : {1e-200, 1.0, 1e200}) {
        for (const Eigen::Vector3d &coeff :
             {Eigen::Vector3d(1, -3, 2), Eigen::Vector3d(1, 0, -1), Eigen::Vector3d(1, -1, 0),
              Eigen::Vector3d(0, 1, -2), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(1, 2, 1), Eigen::Vector3d(1, 0, 0),
              Eigen::Vector3d(0, 0, 1)}) {
            Eigen::Vector2d roots[2];
            const int count =
                univariate::solve_quadratic_real(scale * coeff(0), scale * coeff(1), scale * coeff(2), roots);
            const double disc = coeff(1) * coeff(1) - 4 * coeff(0) * coeff(2);
            REQUIRE_EQ(count, disc == 0.0 ? 1 : 2);
            for (int i = 0; i < count; ++i) {
                REQUIRE(roots[i].allFinite());
                REQUIRE_SMALL(roots[i].norm() - 1.0, 1e-14);
                const double u = roots[i](0), v = roots[i](1);
                REQUIRE_SMALL(coeff(0) * u * u + coeff(1) * u * v + coeff(2) * v * v, 1e-14);
            }
            if (count == 2)
                REQUIRE(std::abs(roots[0](0) * roots[1](1) - roots[0](1) * roots[1](0)) > 0.1);
        }
    }
    Eigen::Vector2d roots[2];
    REQUIRE_EQ(univariate::solve_quadratic_real(0, 0, 0, roots), 0);
    REQUIRE_EQ(univariate::solve_quadratic_real(1, 0, 1, roots), 0);
    REQUIRE_EQ(univariate::solve_quadratic_real(1, std::numeric_limits<double>::quiet_NaN(), 1, roots), 0);
    double affine[2];
    REQUIRE_EQ(univariate::solve_quadratic_real(1, -2, 1, affine), 2);
    REQUIRE_EQ(affine[0], 1.0);
    REQUIRE_EQ(affine[1], 1.0);
    return true;
}

bool test_pivoted_conic_factorization() {
    // The first coordinate vanishes for the coplanar mixed-pose pencil.
    // Exercise every permutation, including zero first row/column.
    for (int k = 0; k < 3; ++k) {
        Eigen::Vector3d p = Eigen::Vector3d::Zero(), q = Eigen::Vector3d::Zero();
        p((k + 1) % 3) = 1;
        p((k + 2) % 3) = 2;
        q((k + 1) % 3) = 3;
        q((k + 2) % 3) = -1;
        const Eigen::Matrix3d C = (p * q.transpose() + q * p.transpose()) * 0.5;
        const auto factors = compute_pq(C, ConicFactorStrategy::MAX_ABS);
        REQUIRE(factors[0].allFinite() && factors[1].allFinite());
        Eigen::Matrix3d reconstructed =
            (factors[0] * factors[1].transpose() + factors[1] * factors[0].transpose()) * 0.5;
        const double scale = reconstructed.cwiseProduct(C).sum() / C.squaredNorm();
        REQUIRE(std::abs(scale) > 0);
        REQUIRE_SMALL((reconstructed / scale - C).norm(), 1e-13);
    }
    const Eigen::Matrix3d H = (Eigen::Matrix3d() << 1, 2, 3, 0, 4, 5, 1, 0, 6).finished();
    REQUIRE_SMALL((H * adjugate(H) - H.determinant() * Eigen::Matrix3d::Identity()).norm(), 1e-13);
    return true;
}
} // namespace

std::vector<Test> register_solver_primitives_test() {
    return {TEST(test_projective_quadratic), TEST(test_pivoted_conic_factorization)};
}
