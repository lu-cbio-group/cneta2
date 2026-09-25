cmake_minimum_required(VERSION 3.10)
project(cneta VERSION 2.0 LANGUAGES CXX)

# Tell CMake to use the modern Boost search behavior and silence CMP0167
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# 1. Set Global C++ Standards
# -----------------------------
# C++14 is the minimum for Catch2 v3 (see tests/unit/). The sources seem to
# build cleanly at this level; only std::random_shuffle warns as deprecated.
# TODO: Moving to C++17
# C++17 additionally requires replacing the two random_shuffle
# calls in tree_op.cpp and genome.cpp, which changes RNG draw order and so
# shifts fixed-seed simulation output.

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Forces standard C++ over compiler-specific extensions

# Where the finished programs (cnets, cnetml, cnetmcmc) are written.
#
# Default: <repo>/bin, wherever the build directory is. The run-*.sh scripts
# and the tests look there.
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