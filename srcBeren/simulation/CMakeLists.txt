project(Simulation LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)

set(SOURCE_LIB  simulation.cpp
                recovery.cpp)


add_library(Simulation ${SOURCE_LIB})
add_subdirectory(ecsim)
add_subdirectory(ecsim_corr)
target_include_directories(Simulation
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ecsim
    ${CMAKE_CURRENT_SOURCE_DIR}/ecsim_corr
)
