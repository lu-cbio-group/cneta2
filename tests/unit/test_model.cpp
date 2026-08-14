// Unit tests for the copy-number evolution model in model.cpp.
//
// These pin down the mathematical invariants that every rate/transition
// matrix must satisfy, independent of the particular rates used. They are
// the cheapest guard we have against a refactor quietly breaking the
// likelihood.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <vector>

#include "model.hpp"

namespace {

// The matrices are stored column-major: m[i + j * n] is row i, column j.
double row_sum(const double* m, int n, int row){
    double sum = 0.0;
    for(int j = 0; j < n; ++j){
        sum += m[row + j * n];
    }
    return sum;
}

// A generator (Q) matrix must have every row summing to zero.
//
// model.hpp offers check_matrix_row_sum() for this, but it compares against
// exactly 0.0 and prints a line per row, so it is unusable as an assertion.
// We do the same check with a tolerance instead.
void require_generator(const double* m, int n){
    for(int i = 0; i < n; ++i){
        REQUIRE_THAT(row_sum(m, n, i),
                     Catch::Matchers::WithinAbs(0.0, 1e-12));
    }
}

}  // namespace

TEST_CASE("compute_haplotype_change_dim follows the parity rule", "[model]"){
    // Odd n gives (2n+1)^2, even n gives (2n)^2.
    REQUIRE(compute_haplotype_change_dim(1) == 9);     // (2*1+1)^2
    REQUIRE(compute_haplotype_change_dim(3) == 49);    // (2*3+1)^2
    REQUIRE(compute_haplotype_change_dim(2) == 16);    // (2*2)^2
    REQUIRE(compute_haplotype_change_dim(4) == 64);    // (2*4)^2
}

TEST_CASE("get_rate_matrix_bounded is a valid generator", "[model]"){
    const int cn_max = 4;
    const int n = cn_max + 1;
    std::vector<double> q(n * n, 0.0);

    get_rate_matrix_bounded(q.data(), 0.01, 0.02, cn_max);

    require_generator(q.data(), n);

    SECTION("off-diagonal entries are non-negative, diagonal non-positive"){
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                if(i == j){
                    REQUIRE(q[i + j * n] <= 0.0);
                }else{
                    REQUIRE(q[i + j * n] >= 0.0);
                }
            }
        }
    }

    SECTION("copy number zero is absorbing"){
        // Nothing can be duplicated or deleted once the copy number is 0.
        for(int j = 0; j < n; ++j){
            REQUIRE(q[0 + j * n] == 0.0);
        }
    }

    SECTION("the maximum copy number can only be deleted"){
        // At the ceiling there is no duplication left to do, so the only
        // non-zero off-diagonal is the step down.
        REQUIRE(q[cn_max + (cn_max - 1) * n] > 0.0);
        REQUIRE(q[cn_max + cn_max * n] < 0.0);
    }

    SECTION("rates scale linearly with copy number"){
        // A cell with copy number i has 2i copies available to duplicate.
        const double dup_rate = 0.01;
        for(int i = 1; i < cn_max; ++i){
            REQUIRE_THAT(q[i + (i + 1) * n],
                         Catch::Matchers::WithinRel(2.0 * i * dup_rate, 1e-12));
        }
    }
}

TEST_CASE("get_rate_matrix_haplotype_specific is a valid generator", "[model]"){
    const int cn_max = 4;
    const int n = (cn_max + 1) * (cn_max + 2) / 2;   // 15 states
    std::vector<double> q(n * n, 0.0);

    get_rate_matrix_haplotype_specific(q.data(), 0.01, 0.02, cn_max);

    REQUIRE(n == 15);
    require_generator(q.data(), n);
}

TEST_CASE("get_rate_matrix_wgd is a valid generator", "[model]"){
    const int max_wgd = 2;
    const int n = max_wgd + 1;
    std::vector<double> q(n * n, 0.0);

    get_rate_matrix_wgd(q.data(), 0.03, max_wgd);

    require_generator(q.data(), n);

    SECTION("WGD only ever moves the count upwards"){
        REQUIRE_THAT(q[0 + 1 * n], Catch::Matchers::WithinRel(0.03, 1e-12));
        REQUIRE_THAT(q[1 + 2 * n], Catch::Matchers::WithinRel(0.03, 1e-12));
        REQUIRE(q[1 + 0 * n] == 0.0);
        REQUIRE(q[2 + 1 * n] == 0.0);
    }

    SECTION("the maximum number of WGD events is absorbing"){
        for(int j = 0; j < n; ++j){
            REQUIRE(q[max_wgd + j * n] == 0.0);
        }
    }
}

TEST_CASE("get_transition_matrix_bounded produces a stochastic matrix", "[model]"){
    const int cn_max = 4;
    const int n = cn_max + 1;
    std::vector<double> q(n * n, 0.0);
    std::vector<double> p(n * n, 0.0);

    get_rate_matrix_bounded(q.data(), 0.01, 0.02, cn_max);

    SECTION("rows sum to one and entries are probabilities"){
        for(double t : {0.5, 5.0, 50.0}){
            get_transition_matrix_bounded(q.data(), p.data(), t, n);

            for(int i = 0; i < n; ++i){
                REQUIRE_THAT(row_sum(p.data(), n, i),
                             Catch::Matchers::WithinAbs(1.0, 1e-9));
                for(int j = 0; j < n; ++j){
                    REQUIRE(p[i + j * n] >= 0.0);
                    REQUIRE(p[i + j * n] <= 1.0 + 1e-12);
                }
            }
        }
    }

    SECTION("zero elapsed time gives the identity"){
        get_transition_matrix_bounded(q.data(), p.data(), 0.0, n);

        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                const double expected = (i == j) ? 1.0 : 0.0;
                REQUIRE_THAT(p[i + j * n],
                             Catch::Matchers::WithinAbs(expected, 1e-12));
            }
        }
    }

    SECTION("over long times everything drains into copy number zero"){
        // State 0 is absorbing and reachable from every other state, so the
        // chain must eventually concentrate there.
        get_transition_matrix_bounded(q.data(), p.data(), 1.0e5, n);

        for(int i = 0; i < n; ++i){
            REQUIRE_THAT(p[i + 0 * n], Catch::Matchers::WithinAbs(1.0, 1e-6));
        }
    }
}

TEST_CASE("get_transition_prob implements the Mk model", "[model]"){
    const double mu = 0.1;

    SECTION("probabilities over the five states sum to one"){
        for(double blen : {0.0, 1.0, 10.0, 100.0}){
            for(int from = 0; from < 5; ++from){
                double total = 0.0;
                for(int to = 0; to < 5; ++to){
                    total += get_transition_prob(mu, blen, from, to);
                }
                REQUIRE_THAT(total, Catch::Matchers::WithinAbs(1.0, 1e-12));
            }
        }
    }

    SECTION("a zero-length branch cannot change state"){
        REQUIRE_THAT(get_transition_prob(mu, 0.0, 2, 2),
                     Catch::Matchers::WithinAbs(1.0, 1e-12));
        REQUIRE_THAT(get_transition_prob(mu, 0.0, 2, 3),
                     Catch::Matchers::WithinAbs(0.0, 1e-12));
    }

    SECTION("over long branches every state becomes equally likely"){
        for(int to = 0; to < 5; ++to){
            REQUIRE_THAT(get_transition_prob(mu, 1.0e6, 0, to),
                         Catch::Matchers::WithinAbs(0.2, 1e-9));
        }
    }

    SECTION("staying is always at least as likely as any particular change"){
        REQUIRE(get_transition_prob(mu, 1.0, 1, 1)
                > get_transition_prob(mu, 1.0, 1, 4));
    }
}
