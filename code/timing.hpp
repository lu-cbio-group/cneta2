// timing.hpp
#ifndef TIMING_HPP
#define TIMING_HPP

#include <chrono>

namespace timing {

// Use a unified Clock type (steady_clock)
using Clock     = std::chrono::steady_clock;
using TimePoint = Clock::time_point;

// Current time point
inline TimePoint now() {
    return Clock::now();
}

// Elapsed seconds between two time points
double elapsed_seconds(TimePoint start, TimePoint end);

} // namespace timing

#endif // TIMING_HPP
