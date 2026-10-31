cmake_minimum_required(VERSION 3.18)

set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_STANDARD 20)
project(beren3d LANGUAGES CXX)

#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
option(BUILD_TESTS "Build test executables" OFF)
option(USE_TIMERS "Build with flat and tree timers" OFF)
option(ENABLE_COVERAGE "Instrument project for GCC gcov coverage" OFF)
set(MEMORY_TIMERS_LEVEL 0 CACHE STRING "Level of memory profiling: 0(disable), 1 (*llocators), 2 (+ memcpy, memset)")
set(USE_MPI OFF)
set(USE_AMGCL ON)

if(NOT ${USE_TIMERS} AND ${MEMORY_TIMERS_LEVEL} GREATER 0)
    message(FATAL_ERROR "Your are request profiling memory with disable profile")
endif()

add_compile_definitions($<$<BOOL:${USE_TIMERS}>:USE_TIMERS>)
add_compile_definitions("MEMORY_TIMERS_LEVEL=${MEMORY_TIMERS_LEVEL}")
if (${MEMORY_TIMERS_LEVEL} GREATER 0)
    message(STATUS "Common memory profiling")
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--wrap realloc -Wl,--wrap malloc -Wl,--wrap calloc -Wl,--wrap reallocarray -Wl,--wrap free ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--wrap _Znwm -Wl,--wrap _Znam -Wl,--wrap _ZdlPvm ${CMAKE_EXE_LINKER_FLAGS}")
endif()
if (${MEMORY_TIMERS_LEVEL} GREATER 1)
    message(STATUS "Extended memory profiling")
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--wrap memcpy -Wl,--wrap memset  ${CMAKE_EXE_LINKER_FLAGS}")
endif()

find_package(OpenMP REQUIRED)
link_libraries(OpenMP::OpenMP_CXX)
add_compile_options(-march=native -Wall -Wextra)
add_compile_options(-fmax-errors=5)
include_directories(include algorithms diagnostics fields particles utils simulation vahedi_surendra_advance .)

include_directories(SYSTEM 3rd_party)
include_directories(SYSTEM 3rd_party/eigen)

if(${USE_AMGCL})
    add_compile_definitions(AMGCL_NO_BOOST)
    include_directories(SYSTEM 3rd_party/amgcl)
endif()

if(ENABLE_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "ENABLE_COVERAGE currently requires GCC/gcov")
    endif()
    add_compile_options(--coverage -O0)
    add_link_options(--coverage)
endif()

# Pure value types, geometry and JSON loading do not need the numerical runtime.
add_subdirectory(core)
add_subdirectory(grid)
add_subdirectory(config)

set(LIBRARIES Diagnostics Fields Particles Utils Simulation advVS94)
set(BIN_PATH ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(diagnostics)
add_subdirectory(fields)
add_subdirectory(particles)
add_subdirectory(utils)
add_subdirectory(simulation)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
add_subdirectory(vahedi_surendra_advance)

foreach(grid_client IN ITEMS Diagnostics Fields Particles Utils Simulation)
    target_link_libraries(${grid_client} PUBLIC BerenGrid)
endforeach()
target_link_libraries(Simulation PRIVATE Configuration)
foreach(core_client IN ITEMS advVS94 collision_processing)
    target_link_libraries(${core_client} PUBLIC BerenCore)
endforeach()

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS ${PROJECT_NAME} DESTINATION ${BIN_PATH})

message(STATUS "Project ${PROJECT_NAME} configured successfully")
message(STATUS "Binary will be installed to: ${BIN_PATH}")
