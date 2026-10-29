// Unit tests for the optimiser's coordinate scaling and finite-difference
// gradient in optimization.hpp.
//
// L-BFGS-B works on log(rate), log(multiplier) and log(branch length)
// (to_internal / from_internal), and takes its gradient by central
// differences kept inside the bounds (fd_gradient). These tests pin the
// helpers down on functions whose gradient is known exactly.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "optimization.hpp"

using Catch::Approx;

TEST_CASE("to_internal and from_internal are inverses and keep the bounds in order", "[optimiser]"){
    CHECK(from_internal(to_internal(0.001, true), true) == Approx(0.001));
    CHECK(from_internal(to_internal(37.5, true), true) == Approx(37.5));
    CHECK(to_internal(0.37, false) == 0.37);
    CHECK(from_internal(0.37, false) == 0.37);
    CHECK(to_internal(MIN_MRATE, true) < to_internal(MAX_MRATE, true));
    CHECK(to_internal(M_MIN, true) < to_internal(1.0, true));
    CHECK(to_internal(1.0, true) == 0.0);
}

TEST_CASE("fd_gradient takes central differences inside the bounds and one-sided at them", "[optimiser]"){
    // f(z) = (z0 - 1)^2 + exp(z1) + 3 z2, gradient (2 (z0 - 1), exp(z1), 3).
    auto f = [](const double* z){
        return (z[0] - 1.0) * (z[0] - 1.0) + std::exp(z[1]) + 3.0 * z[2];
    };
    double z[3] = {0.5, 0.0, 2.0};
    double l[3] = {0.0, -1.0, 2.0};   // z2 sits on its lower bound
    double u[3] = {1.0, 1.0, 5.0};
    double g[3] = {0.0, 0.0, 0.0};

    double fz = fd_gradient(f, 3, z, l, u, g);

    CHECK(fz == Approx(0.25 + 1.0 + 6.0));
    CHECK(g[0] == Approx(-1.0).margin(1e-7));
    CHECK(g[1] == Approx(1.0).margin(1e-7));
    CHECK(g[2] == Approx(3.0).margin(1e-7));
    // the point is left as it was
    CHECK(z[0] == 0.5);
    CHECK(z[1] == 0.0);
    CHECK(z[2] == 2.0);
}

TEST_CASE("fd_gradient is more accurate than a forward difference on a curved function", "[optimiser]"){
    // exp(z) at z = 5: a forward step of 5e-4 is off by about half the step
    // times the curvature (0.04); central differences are off by the step
    // squared (about 6e-6).
    auto f = [](const double* z){ return std::exp(z[0]); };
    double z[1] = {5.0};
    double l[1] = {0.0};
    double u[1] = {10.0};
    double g[1] = {0.0};

    fd_gradient(f, 1, z, l, u, g);

    double exact = std::exp(5.0);
    CHECK(std::fabs(g[0] - exact) < 1e-4 * exact);
    double h = FD_REL_STEP * 5.0;
    double forward = (std::exp(5.0 + h) - exact) / h;
    CHECK(std::fabs(forward - exact) > 10 * std::fabs(g[0] - exact));
}

TEST_CASE("fd_gradient never evaluates outside the bounds", "[optimiser]"){
    double l[2] = {2.0, -3.0};
    double u[2] = {2.0 + 1e-6, 3.0};   // the first box is narrower than the step
    double z[2] = {2.0, 3.0};         // the second sits on its upper bound
    int outside = 0;
    auto f = [&](const double* zz){
        for(int i = 0; i < 2; i++){
            if(zz[i] < l[i] || zz[i] > u[i]) outside++;
        }
        return 5.0 * zz[0] - 2.0 * zz[1];
    };
    double g[2] = {0.0, 0.0};

    fd_gradient(f, 2, z, l, u, g);

    CHECK(outside == 0);
    CHECK(g[0] == Approx(5.0));
    CHECK(g[1] == Approx(-2.0));
}

TEST_CASE("fd_gradient reports a zero gradient for a fixed variable", "[optimiser]"){
    auto f = [](const double* z){ return z[0] * z[0]; };
    double z[1] = {1.5};
    double l[1] = {1.5};
    double u[1] = {1.5};
    double g[1] = {7.0};

    fd_gradient(f, 1, z, l, u, g);

    CHECK(g[0] == 0.0);
}
