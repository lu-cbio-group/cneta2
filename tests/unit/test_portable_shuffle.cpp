// Unit tests for portable_shuffle (portable_shuffle.hpp), the replacement for
// std::random_shuffle.
//
// The pinned sequences are what libstdc++'s std::random_shuffle produces on
// the same generator (checked with g++ 16). They are the point of the
// helper: every seeded result in cneta depends on this exact order, so a
// change here is a change in every simulation and every tree search.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <numeric>
#include <vector>

#include "portable_shuffle.hpp"

namespace {

// The generator the reference sequences were produced with: a 64-bit LCG
// whose high bits are reduced into [0, n), with the fp_myrng signature.
struct Lcg {
    unsigned long long state;
    int calls = 0;
    explicit Lcg(unsigned long long seed) : state(seed) {}
    unsigned long operator()(unsigned long n){
        ++calls;
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return (unsigned long)((state >> 33) % n);
    }
};

std::vector<int> iota_vector(int n){
    std::vector<int> v(n);
    std::iota(v.begin(), v.end(), 0);
    return v;
}

std::vector<int> shuffled(int n, unsigned long long seed = 12345){
    std::vector<int> v = iota_vector(n);
    Lcg rng(seed);
    portable_shuffle(v.begin(), v.end(), rng);
    return v;
}

// A generator with the exact type cneta passes: long unsigned (*)(long unsigned).
unsigned long draw_last(unsigned long n){ return n - 1; }

}  // namespace

TEST_CASE("portable_shuffle reproduces libstdc++'s random_shuffle order",
          "[portable_shuffle]"){
    // std::random_shuffle(v.begin(), v.end(), Lcg(12345)) under g++ 16.
    CHECK(shuffled(0) == std::vector<int>{});
    CHECK(shuffled(1) == std::vector<int>{0});
    CHECK(shuffled(2) == std::vector<int>{1, 0});
    CHECK(shuffled(3) == std::vector<int>{1, 2, 0});
    CHECK(shuffled(5) == std::vector<int>{1, 4, 3, 0, 2});
    CHECK(shuffled(10) == std::vector<int>{1, 8, 6, 7, 5, 2, 3, 0, 9, 4});
    CHECK(shuffled(23) == std::vector<int>{16, 8, 12, 10, 21, 18, 3, 0, 9, 17,
                                           7, 13, 6, 20, 2, 14, 19, 4, 15, 1,
                                           11, 22, 5});
}

TEST_CASE("portable_shuffle permutes the range", "[portable_shuffle]"){
    for(int n : {2, 7, 64, 1000}){
        std::vector<int> v = shuffled(n, 99);
        std::sort(v.begin(), v.end());
        CHECK(v == iota_vector(n));
    }
}

TEST_CASE("portable_shuffle draws n - 1 numbers, one per position after the first",
          "[portable_shuffle]"){
    // The number and order of draws ties every later random number to the
    // seed: a helper that drew differently would still permute, but would
    // shift everything that follows it.
    for(int n : {0, 1, 2, 10, 64}){
        std::vector<int> v = iota_vector(n);
        Lcg rng(1);
        portable_shuffle(v.begin(), v.end(), rng);
        CHECK(rng.calls == std::max(0, n - 1));
    }
}

TEST_CASE("portable_shuffle is deterministic for a given generator state",
          "[portable_shuffle]"){
    CHECK(shuffled(50, 7) == shuffled(50, 7));
    CHECK(shuffled(50, 7) != shuffled(50, 8));
}

TEST_CASE("portable_shuffle accepts a plain function pointer like fp_myrng",
          "[portable_shuffle]"){
    unsigned long (*fp)(unsigned long) = &draw_last;
    std::vector<int> v = iota_vector(5);
    portable_shuffle(v.begin(), v.end(), fp);
    // Always drawing the last index means j == i at every step, so nothing
    // moves.
    CHECK(v == iota_vector(5));
}
