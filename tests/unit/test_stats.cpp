// Unit tests for the small numeric helpers in stats.hpp.

#include <catch2/catch_test_macros.hpp>

#include "stats.hpp"

TEST_CASE("pairInteger is symmetric", "[stats]"){
    // The pairing function is used to key cached site patterns, so it must
    // give the same answer whichever way round the arguments arrive.
    for(int a = 0; a <= 12; ++a){
        for(int b = 0; b <= 12; ++b){
            REQUIRE(pairInteger(a, b) == pairInteger(b, a));
        }
    }
}

TEST_CASE("pairInteger is injective over unordered pairs", "[stats]"){
    // Distinct unordered pairs must not collide, or two different site
    // patterns would share a cache slot.
    std::vector<int> seen;
    for(int a = 0; a <= 20; ++a){
        for(int b = a; b <= 20; ++b){
            seen.push_back(pairInteger(a, b));
        }
    }

    std::sort(seen.begin(), seen.end());
    REQUIRE(std::adjacent_find(seen.begin(), seen.end()) == seen.end());
}

TEST_CASE("pairInteger matches the Cantor pairing definition", "[stats]"){
    // Spot values from T(a+b) + max(a,b), where T(n) = n(n+1)/2.
    // Note it is the *larger* argument that is added, not the smaller: both
    // branches of the implementation add whichever of the two is bigger.
    REQUIRE(pairInteger(0, 0) == 0);    // T(0) + 0
    REQUIRE(pairInteger(0, 1) == 2);    // T(1) + 1
    REQUIRE(pairInteger(1, 1) == 4);    // T(2) + 1
    REQUIRE(pairInteger(2, 3) == 18);   // T(5) + 3
}

TEST_CASE("fact computes small factorials", "[stats]"){
    REQUIRE(fact(0) == 1);
    REQUIRE(fact(1) == 1);
    REQUIRE(fact(2) == 2);
    REQUIRE(fact(5) == 120);
    // 12! is the largest factorial that fits in a 32-bit int; anything above
    // this overflows, which is why callers must keep n small.
    REQUIRE(fact(12) == 479001600);
}

TEST_CASE("rchoose picks the only outcome with non-zero weight", "[stats]"){
    // Weighted sampling is inherently random, but a degenerate weight vector
    // makes the outcome deterministic, so this is safe to assert exactly.
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 12345);

    SECTION("all weight on the first element"){
        std::vector<double> rates{1.0, 0.0, 0.0, 0.0};
        for(int i = 0; i < 50; ++i){
            REQUIRE(rchoose(r, rates) == 0);
        }
    }

    SECTION("all weight on the last element"){
        std::vector<double> rates{0.0, 0.0, 0.0, 1.0};
        for(int i = 0; i < 50; ++i){
            REQUIRE(rchoose(r, rates) == 3);
        }
    }

    SECTION("a spread of weights stays in range"){
        std::vector<double> rates{0.1, 0.2, 0.3, 0.4};
        for(int i = 0; i < 200; ++i){
            int k = rchoose(r, rates);
            REQUIRE(k >= 0);
            REQUIRE(k < static_cast<int>(rates.size()));
        }
    }

    gsl_rng_free(r);
}

TEST_CASE("runiform stays within its bounds", "[stats]"){
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 12345);

    for(int i = 0; i < 500; ++i){
        double x = runiform(r, 2.0, 5.0);
        REQUIRE(x >= 2.0);
        REQUIRE(x <= 5.0);
    }

    gsl_rng_free(r);
}
