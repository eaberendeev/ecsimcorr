cmake_minimum_required(VERSION 3.16)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(USE_MPI OFF)
set(USE_AMGCL ON)

if(${USE_MPI})
    set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx.mpich")
else()
    set(CMAKE_CXX_COMPILER "g++")
endif()

#set(CMAKE_CXX_FLAGS "-Wall -g -fsanitize=address -fno-omit-frame-pointer -O1  -march=native -fopenmp -Wextra -std=c++17")
set(CMAKE_CXX_FLAGS "-Wall -g -O3 -march=native -fopenmp -Wextra -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_VERBOSE_MAKEFILE on)

project(beren3d LANGUAGES CXX)

set(LIBRARIES Diagnostics Fields Particles Utils Simulation advVS94)
set(INCLUDES constants diagnostics fields particles utils simulation vahedi_surendra_advance .)
set(BIN_PATH ${PROJECT_BINARY_DIR}/bin)
#set(TESTS_PATH ${BIN_PATH}/tests)

include_directories(SYSTEM ${PATH_TO_EIGEN})

if(${USE_AMGCL})
include_directories(SYSTEM ${PATH_TO_AMGCL})
endif()

include_directories(${INCLUDES})

if(${USE_MPI})
    add_library(PETSC UNKNOWN IMPORTED)
    find_package(MPI REQUIRED)
    find_library(petsc libpetsc.so /home/master/soft/install/petsc/lib)
    message(STATUS "${petsc}")
    set_target_properties(PETSC PROPERTIES
    IMPORTED_LOCATION ${petsc}
    )
    target_include_directories(PETSC INTERFACE /home/master/soft/install/petsc/include)
    include_directories(SYSTEM /usr/lib/x86_64-linux-gnu/mpich/include)
endif()

add_subdirectory(diagnostics)
add_subdirectory(fields)
add_subdirectory(particles)
add_subdirectory(utils)
add_subdirectory(simulation)
add_subdirectory(tests)
add_subdirectory(vahedi_surendra_advance)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS ${PROJECT_NAME} DESTINATION ${BIN_PATH})

message(STATUS "Project ${PROJECT_NAME} configured successfully")
message(STATUS "Binary will be installed to: ${BIN_PATH}")
