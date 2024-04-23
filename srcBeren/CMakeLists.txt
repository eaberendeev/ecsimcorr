cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-Wall -fopenmp -Wextra -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_VERBOSE_MAKEFILE on)

project(beren3d LANGUAGES CXX)

set(LIBRARIES Diagnostics Fields Particles Utils)
set(INCLUDES constants diagnostics fields particles simulation utils .)
set(BIN_PATH ${PROJECT_BINARY_DIR}/bin)
#set(TESTS_PATH ${BIN_PATH}/tests)

include_directories(SYSTEM ${PATH_TO_EIGEN})
include_directories(${INCLUDES})

add_subdirectory(diagnostics)
add_subdirectory(fields)
add_subdirectory(particles)
add_subdirectory(utils)

add_subdirectory(${CURRENT_SIMULATION})

#if (NOT ${ENABLE_TESTS} STREQUAL "NONE")
#add_subdirectory(testing)
#endif()
