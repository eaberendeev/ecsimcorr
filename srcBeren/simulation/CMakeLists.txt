project(simulation LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)

set(SOURCE_LIB  simulation.cpp
                recovery.cpp)


# Boundary laws are also instantiated by configuration. Link their objects
# directly so resolving the configuration archive does not depend on link order.
add_library(simulation_boundaries OBJECT
                boundary_conditions/boundary_condition_handler.cpp
                boundary_conditions/open_boundary_conditions.cpp
                boundary_conditions/periodic_boundary_conditions.cpp
                boundary_conditions/er0_condition.cpp
                boundary_conditions/electron_reflection_condition.cpp)
target_link_libraries(simulation_boundaries PUBLIC beren_grid)

add_library(simulation ${SOURCE_LIB})
add_subdirectory(ecsim)
add_subdirectory(ecsim_corr)
target_include_directories(simulation
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ecsim
    ${CMAKE_CURRENT_SOURCE_DIR}/ecsim_corr
)
