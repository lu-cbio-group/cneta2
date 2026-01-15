// timing.cpp
#include "timing.hpp"

namespace timing {

double elapsed_seconds(TimePoint start, TimePoint end) {
    return std::chrono::duration<double>(end - start).count();
}

} // namespace timing
