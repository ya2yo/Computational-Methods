cmake_minimum_required(VERSION 3.20)

project(computational_methods VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

option(COMPUTATIONAL_METHODS_BUILD_TESTS "Build automated tests" ON)
option(COMPUTATIONAL_METHODS_ENABLE_WARNINGS "Enable compiler warnings" ON)
option(COMPUTATIONAL_METHODS_ENABLE_COVERAGE "Enable gcov coverage instrumentation" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(computational_methods
    src/curve_fit.cpp
    src/direct_solvers.cpp
    src/ordinary_differential_equations.cpp
    src/root_finding.cpp
    src/interpolation.cpp
    src/iterative_solvers.cpp
    src/numerical_integration.cpp
)

target_include_directories(computational_methods
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/computational_methods>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/computational_methods>"
)

if(COMPUTATIONAL_METHODS_ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(computational_methods PRIVATE --coverage -O0 -g)
    target_link_options(computational_methods PUBLIC --coverage)
endif()

if(COMPUTATIONAL_METHODS_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(computational_methods PRIVATE /W4 /permissive-)
    else()
        target_compile_options(computational_methods PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

add_executable(computational_methods_demo examples/main.cpp)
target_link_libraries(computational_methods_demo PRIVATE computational_methods)

include(CTest)
if(BUILD_TESTING AND COMPUTATIONAL_METHODS_BUILD_TESTS)
    add_executable(computational_methods_tests tests/test_algorithms.cpp)
    target_link_libraries(computational_methods_tests PRIVATE computational_methods)
    add_test(NAME computational_methods_tests COMMAND computational_methods_tests)
endif()

set_target_properties(computational_methods computational_methods_demo
    PROPERTIES
        CXX_EXTENSIONS OFF
)
