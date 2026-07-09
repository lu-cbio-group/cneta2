cmake_minimum_required(VERSION 3.10)
project(cneta VERSION 2.0 LANGUAGES CXX)

# Tell CMake to use the modern Boost search behavior and silence CMP0167
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# 1. Set Global C++ Standards
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Forces standard C++ over compiler-specific extensions

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