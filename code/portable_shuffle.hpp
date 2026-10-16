#ifndef PORTABLE_SHUFFLE_HPP
#define PORTABLE_SHUFFLE_HPP

#include <algorithm>   // std::iter_swap

/**
 * @brief Shuffle [first, last) in place, drawing from rng(n), which must
 * return a value in [0, n).
 *
 * This replaces std::random_shuffle, which C++17 removed and which, worse,
 * libstdc++ (GCC, the HPC) and libc++ (Clang, macOS) implement differently:
 * the same seed gave different trees and simulated data on the two
 * platforms. With this helper the sequence of draws is the same everywhere.
 *
 * The loop is libstdc++'s, so results produced with GCC before this helper
 * existed can still be reproduced from their seeds; macOS output changed
 * once when it was introduced. Do not "improve" the loop: any change to the
 * sequence of rng() calls changes every seeded result.
 *
 * @param first, last: random-access iterators
 * @param rng: callable, rng(n) in [0, n); typically an Rng (stats.hpp) over a GSL generator
 */
template <typename RandomIt, typename Rng>
void portable_shuffle(RandomIt first, RandomIt last, Rng&& rng){
    if(first == last) return;
    for(RandomIt i = first + 1; i != last; ++i){
        RandomIt j = first + rng((i - first) + 1);
        if(i != j) std::iter_swap(i, j);
    }
}

#endif
