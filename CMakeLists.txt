cmake_minimum_required(VERSION 3.14)

project(binning-analysis VERSION 0.1.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(BA_BUILD_EXAMPLES "Build examples" OFF)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++20" COMPILER_SUPPORTS_CXX20)
if(NOT COMPILER_SUPPORTS_CXX20)
    message(FATAL_ERROR "binning-analysis requires a compiler with C++20 support")
endif()

add_library(binning_analysis_lib INTERFACE)
target_include_directories(binning_analysis_lib INTERFACE "${PROJECT_SOURCE_DIR}/include")


if (BA_BUILD_EXAMPLES)
    add_executable(1d_ising_mcmc "examples/1d_ising_mcmc.cpp")
    target_link_libraries(1d_ising_mcmc PRIVATE binning_analysis)
    set_property(TARGET 1d_ising_mcmc PROPERTY CXX_STANDARD 20)
endif()

export(TARGETS binning_analysis_lib NAMESPACE BinningAnalysis:: FILE BinningAnalysis.cmake)
