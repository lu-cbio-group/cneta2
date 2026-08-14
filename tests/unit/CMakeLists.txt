# Catch2 unit tests for the shared cneta library.
#
# Enabled from code/CMakeLists.txt with -DCNETA_BUILD_TESTING=ON.
#
# These link against `cneta_core`, the static library that already holds all
# the shared science. When that target becomes `libcneta`, only the
# target_link_libraries line below needs to change.

if(CMAKE_VERSION VERSION_LESS 3.16)
    message(FATAL_ERROR
        "CNETA_BUILD_TESTING requires CMake >= 3.16 (found ${CMAKE_VERSION}). "
        "The main build still works with 3.10; only the tests need a newer CMake.")
endif()

# Prefer a Catch2 that is already installed (Homebrew, apt, an HPC module).
# Fall back to fetching a pinned release. On an offline machine you can point
# at a local checkout instead:
#     cmake ... -DFETCHCONTENT_SOURCE_DIR_CATCH2=/path/to/Catch2
find_package(Catch2 3 QUIET)

if(NOT Catch2_FOUND)
    message(STATUS "System Catch2 not found - fetching v3.7.1")
    include(FetchContent)
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.7.1
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
else()
    message(STATUS "Using system Catch2 ${Catch2_VERSION}")
endif()

add_executable(cneta_unit_tests
    test_stats.cpp
    test_model.cpp
    test_parse_cn.cpp
    test_tree_op.cpp
)

target_link_libraries(cneta_unit_tests PRIVATE cneta_core Catch2::Catch2WithMain)

# Fixture files are read at run time, so the tests need to know where they are
# regardless of the working directory ctest happens to use.
target_compile_definitions(cneta_unit_tests PRIVATE
    CNETA_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data")

# bin/ is reserved for the three user-facing tools; test binaries stay in the
# build tree so `ls bin/` keeps meaning "the programs a user runs".
set_target_properties(cneta_unit_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

include(Catch)
catch_discover_tests(cneta_unit_tests)
