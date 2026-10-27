cmake_minimum_required(VERSION 3.16)

# The version lives in one place, the VERSION file at the repository root
# (MAJOR.MINOR.PATCH, see docs/source/developer-guide/releases.md). The
# programs print it with --version and write it into their outputs.
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION" CNETA_VERSION LIMIT_COUNT 1)
string(STRIP "${CNETA_VERSION}" CNETA_VERSION)
if(NOT CNETA_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    message(FATAL_ERROR "The VERSION file must hold MAJOR.MINOR.PATCH, e.g. 2.1.0; found '${CNETA_VERSION}'")
endif()
project(cneta VERSION ${CNETA_VERSION} LANGUAGES CXX)

# Which commit the build came from, so results produced between releases can
# be traced too: "v2.1.0" on a tagged release, "v2.1.0-5-gabc1234" five
# commits later, "-dirty" with uncommitted changes. Taken when CMake
# configures (build.sh configures every time); empty outside a git checkout.
set(CNETA_GIT_DESCRIBE "")
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../.git")
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --tags --always --dirty
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.."
        OUTPUT_VARIABLE CNETA_GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
endif()
message(STATUS "cneta version ${CNETA_VERSION} (git: ${CNETA_GIT_DESCRIBE})")

# Tell CMake to use the modern Boost search behavior and silence CMP0167
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# 1. Set Global C++ Standards
# -----------------------------
# C++20 (GCC 10, Clang 10 or Apple clang 12 and newer; the HPC's GCC 13 and
# the conda compilers qualify). std::random_shuffle, which blocked anything
# past C++14, was replaced by portable_shuffle.hpp. Catch2 v3 (tests/unit/)
# needs at least C++14; the legacy code/makefile uses the same standard.

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Forces standard C++ over compiler-specific extensions

# Floating-point contraction (fused multiply-add) is on by default in Clang
# and off in GCC under a strict -std=c++XX, as used here. Fused and unfused
# results differ in the last bit, which can be enough to send an optimiser
# down a different path, so turn it off everywhere: the same seed should give
# the same tree on a Mac and on the HPC.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-ffp-contract=off)
endif()

# Build optimised by default. With no build type CMake adds no -O flag at all,
# so build.sh produced unoptimised programs while the legacy makefile uses
# -O3. Match the makefile (-O3, assertions kept) so the two builds compute the
# same thing; pass -DCMAKE_BUILD_TYPE=Debug for an unoptimised build.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Where the finished programs (cnets, cnetml, cnetmcmc) are written.
#
# Default: <repo>/bin, wherever the build directory is, next to the
# bin/run-*.sh scripts. The scripts and the tests look there.
#
# To send them somewhere else, pass -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=<dir> when
# configuring. The build-system tests do this so a throwaway build does not
# overwrite the executables in your working tree.
#
# The NOT DEFINED check is what makes that work: -D options are set before this
# file runs, and an unconditional set() here would silently replace them.
# The setting applies to every executable defined below, which is why
# tests/unit/CMakeLists.txt redirects its own test binary.
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)
endif()

# Apply the Boost flags from the legacy makefile globally
add_compile_definitions(BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_LIB)

# 2. Find External Dependencies
# ------------------------------
# CMake will look for these automatically
find_package(Boost REQUIRED COMPONENTS program_options filesystem)
find_package(GSL REQUIRED)
find_package(ZLIB REQUIRED)
# OpenMP is optional: it parallelises cnetml's tree search, and without it the
# pragmas are ignored and cnetml simply runs single-threaded.
find_package(OpenMP QUIET)
if(NOT OpenMP_CXX_FOUND AND APPLE)
    # Apple's clang ships without OpenMP. Homebrew's libomp provides it, but is
    # keg-only, so FindOpenMP cannot see it without being told where it is
    # (OpenMP_ROOT is not enough for AppleClang). /opt/homebrew is Apple
    # Silicon, /usr/local is Intel.
    foreach(_libomp_prefix /opt/homebrew/opt/libomp /usr/local/opt/libomp)
        if(EXISTS "${_libomp_prefix}/include/omp.h" AND EXISTS "${_libomp_prefix}/lib/libomp.dylib")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${_libomp_prefix}/include")
            set(OpenMP_CXX_LIB_NAMES omp)
            set(OpenMP_omp_LIBRARY "${_libomp_prefix}/lib/libomp.dylib")
            find_package(OpenMP QUIET)
            break()
        endif()
    endforeach()
endif()
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP: enabled (${OpenMP_CXX_FLAGS})")
else()
    message(STATUS "OpenMP: not found, cnetml will run single-threaded"
                   " (macOS: brew install libomp)")
endif()

# 3. Handle Local Sub-Libraries
# ------------------------------
# lbfgsb already uses CMake, so we can just include it natively
add_subdirectory(lbfgsb)

# gzstream uses an old Makefile. Instead of calling make, it's easier to just compile
# its single source file directly into a CMake library.
add_library(gzstream gzstream/gzstream.C)
target_include_directories(gzstream PUBLIC gzstream)
target_link_libraries(gzstream PRIVATE ZLIB::ZLIB)

# 4. Create the core scientific library
# --------------------------------------
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

# The version and commit, for common.hpp (PUBLIC: the programs and the unit
# tests see the same values).
target_compile_definitions(cneta_core PUBLIC
        CNETA_VERSION="${CNETA_VERSION}"
        CNETA_GIT_DESCRIBE="${CNETA_GIT_DESCRIBE}")

# 5. Build the Executables
# -------------------------
# Now, the main programs can be built:
add_executable(cnets cnets.cpp)
target_link_libraries(cnets cneta_core)

add_executable(cnetml cnetml.cpp)
target_link_libraries(cnetml cneta_core)

add_executable(cnetmcmc cnetmcmc.cpp)
target_link_libraries(cnetmcmc cneta_core)

# 6. Tests (opt-in for unit testing)
# -----------------------------------
# OFF by default, so `./build.sh local` needs no network access. To use them:
#     cmake -S code -B code/build -DCNETA_BUILD_TESTING=ON
#     cmake --build code/build -j 4
#     ctest --test-dir code/build --output-on-failure
# (Configuring alone builds nothing.) Catch2 comes from the system if
# installed, otherwise it is downloaded once into <build>/_deps.
option(CNETA_BUILD_TESTING "Build the Catch2 unit tests" OFF)
if(CNETA_BUILD_TESTING)
    enable_testing()
    # The tests live outside this source tree, so the binary directory has
    # to be named explicitly.
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../tests/unit ${CMAKE_BINARY_DIR}/tests-unit)
endif()