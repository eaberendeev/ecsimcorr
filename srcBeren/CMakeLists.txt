cmake_minimum_required(VERSION 3.18)

set(CMAKE_CXX_COMPILER "g++")
project(beren3d LANGUAGES CXX)

option(BUILD_TESTS "Build test executables" OFF)
option(USE_TIMERS "Build with flat and tree timers" OFF)
option(ENABLE_COVERAGE "Instrument project for GCC gcov coverage" OFF)
set(MEMORY_TIMERS_LEVEL 0 CACHE STRING "Level of memory profiling: 0(disable), 1 (*llocators), 2 (+ memcpy, memset)")
set(USE_MPI OFF)
set(USE_AMGCL ON)
set(BIN_PATH ${PROJECT_BINARY_DIR}/bin)

if(NOT USE_TIMERS AND MEMORY_TIMERS_LEVEL GREATER 0)
    message(FATAL_ERROR "Memory profiling requires USE_TIMERS=ON")
endif()

find_package(OpenMP REQUIRED)

# Common project compiler settings propagate only through explicitly linked targets.
add_library(beren_build_options INTERFACE)
target_compile_features(beren_build_options INTERFACE cxx_std_20)
target_compile_options(beren_build_options INTERFACE -march=native -Wall -Wextra -fmax-errors=5)
target_compile_definitions(beren_build_options INTERFACE
    $<$<BOOL:${USE_TIMERS}>:USE_TIMERS>
    MEMORY_TIMERS_LEVEL=${MEMORY_TIMERS_LEVEL})
if(ENABLE_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "ENABLE_COVERAGE currently requires GCC/gcov")
    endif()
    target_compile_options(beren_build_options INTERFACE --coverage -O0)
    target_link_options(beren_build_options INTERFACE --coverage)
endif()

add_library(beren_eigen INTERFACE)
target_include_directories(beren_eigen SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/3rd_party/eigen)
add_library(beren_json INTERFACE)
target_include_directories(beren_json SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/3rd_party)
add_library(beren_amgcl INTERFACE)
target_link_libraries(beren_amgcl INTERFACE beren_eigen)
target_include_directories(beren_amgcl SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/3rd_party/amgcl)
target_compile_definitions(beren_amgcl INTERFACE AMGCL_NO_BOOST)

add_subdirectory(core)
add_subdirectory(grid)
add_subdirectory(runtime)
add_subdirectory(utils)
add_subdirectory(config)
add_subdirectory(fields)
add_subdirectory(particles)
add_subdirectory(collisions)
add_subdirectory(coupling)
add_subdirectory(io)
add_subdirectory(simulation)
add_subdirectory(diagnostics)

add_executable(beren3d app/main.cpp)
# OBJECT dependencies supply usage requirements transitively, but their object
# files must be named directly by each executable. Preserve the existing archives.
target_link_libraries(beren3d PRIVATE simulation diagnostics fields particles coupling utils
    simulation_boundaries collisions beren_io beren_runtime configuration)
install(TARGETS beren3d DESTINATION ${BIN_PATH})

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

message(STATUS "Project ${PROJECT_NAME} configured successfully")
message(STATUS "Binary will be installed to: ${BIN_PATH}")
