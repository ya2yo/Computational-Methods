cmake_minimum_required(VERSION 3.20)

project(computational_methods VERSION 1.0.0 LANGUAGES CXX)

option(COMPUTATIONAL_METHODS_BUILD_TESTS "Build automated tests" ON)
option(COMPUTATIONAL_METHODS_ENABLE_WARNINGS "Enable compiler warnings" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(computational_methods
    "Computational Methods/CurveFit.cpp"
    "Computational Methods/Direct.cpp"
    "Computational Methods/Ordinary_Differential_Equations.cpp"
    "Computational Methods/Unary_Equation.cpp"
    "Computational Methods/function_insert.cpp"
    "Computational Methods/linear_equations_iterative.cpp"
    "Computational Methods/numerical_integration.cpp"
)

target_include_directories(computational_methods PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/Computational Methods"
)

if(COMPUTATIONAL_METHODS_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(computational_methods PRIVATE /W4 /permissive-)
    else()
        target_compile_options(computational_methods PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

add_executable(computational_methods_demo "Computational Methods/main.cpp")
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
