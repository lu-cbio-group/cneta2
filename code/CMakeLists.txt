cmake_minimum_required(VERSION 3.10)
project(cneta VERSION 2.0 LANGUAGES CXX)

# Tell CMake to use the modern Boost search behavior and silence CMP0167
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# 1. Set Global C++ Standards
# C++14 is the minimum for Catch2 v3 (see tests/unit/). The sources build
# cleanly at this level; only std::random_shuffle warns as deprecated.
# Moving to C++17 additionally requires replacing the two random_shuffle
# calls in tree_op.cpp and genome.cpp, which changes RNG draw order and so
# shifts fixed-seed simulation output -- deliberately left for later.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Forces standard C++ over compiler-specific extensions

# All three user-facing tools land in <repo>/bin, regardless of where the
# build tree lives. The run-*.sh scripts and the test suite look there.
#
# Guarded so that -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=... still wins: that is how
# a throwaway build (CI, the build-system tests) writes its executables
# somewhere else instead of overwriting the ones in the working tree.
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)
endif()

# Apply the Boost flags from the legacy makefile globally
add_compile_definitions(BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_LIB)

# 2. Find External Dependencies
# CMake will look for these automatically
find_package(Boost REQUIRED COMPONENTS program_options filesystem)
find_package(GSL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP) # Optional, used in cnetml

# 3. Handle Local Sub-Libraries
# lbfgsb already uses CMake, so we can just include it natively
add_subdirectory(lbfgsb)

# gzstream uses an old Makefile. Instead of calling make, it's easier to just compile
# its single source file directly into a CMake library.
add_library(gzstream gzstream/gzstream.C)
target_include_directories(gzstream PUBLIC gzstream)
target_link_libraries(gzstream PRIVATE ZLIB::ZLIB)

# 4. Create the core scientific library
# We compile all the shared evolutionary math into one library.
add_library(cneta_core
        matexp/matrix_exponential.cpp
        matexp/r8lib.cpp
        evo_tree.cpp
        tree_op.cpp
        genome.cpp
        model.cpp
        parse_cn.cpp
        stats.cpp
        likelihood.cpp
        nni.cpp
        optimization.cpp
        state.cpp
        timing.cpp
)

# Tell the core library where to find headers and link external libraries
target_include_directories(cneta_core PUBLIC ./ lbfgsb)
target_link_libraries(cneta_core
        PUBLIC
        Boost::program_options
        Boost::filesystem
        GSL::gsl
        GSL::gslcblas
        lbfgsb
        gzstream
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(cneta_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# 5. Build the Executables
# Now, the main programs can be built:
add_executable(cnets cnets.cpp)
target_link_libraries(cnets cneta_core)

add_executable(cnetml cnetml.cpp)
target_link_libraries(cnetml cneta_core)

add_executable(cnetmcmc cnetmcmc.cpp)
target_link_libraries(cnetmcmc cneta_core)

# 6. Tests (opt-in)
# OFF by default, so `./build.sh local` behaves exactly as it always has and
# needs no network access. Turn on with:
#     cmake -S code -B code/build -DCNETA_BUILD_TESTING=ON
option(CNETA_BUILD_TESTING "Build the Catch2 unit tests" OFF)
if(CNETA_BUILD_TESTING)
    enable_testing()
    # The tests live outside this source tree, so the binary directory has
    # to be named explicitly.
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../tests/unit ${CMAKE_BINARY_DIR}/tests-unit)
endif()