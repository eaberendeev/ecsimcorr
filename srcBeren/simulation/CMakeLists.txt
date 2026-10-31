project(Simulation LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)

set(SOURCE_LIB  simulation.cpp
                recovery.cpp)


add_library(Simulation ${SOURCE_LIB})
add_subdirectory(ecsim)
add_subdirectory(ecsim_corr)
# CMake 3.16 reads source properties in the directory that creates Simulation;
# TARGET_DIRECTORY, which would allow setting them in ecsim/, needs CMake 3.18.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(ecsim/response/response_kernels.cpp PROPERTIES COMPILE_FLAGS "-Ofast -funroll-loops")
endif()
target_include_directories(Simulation
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ecsim
    ${CMAKE_CURRENT_SOURCE_DIR}/ecsim_corr
)
