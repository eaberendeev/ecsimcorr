add_library(simulation_boundaries OBJECT
    boundary_conditions/boundary_condition_handler.cpp
    boundary_conditions/open_boundary_conditions.cpp
    boundary_conditions/periodic_boundary_conditions.cpp
    boundary_conditions/er0_condition.cpp
    boundary_conditions/electron_reflection_condition.cpp)
target_link_libraries(simulation_boundaries PUBLIC particles fields)

add_library(simulation STATIC simulation.cpp)
# Keep the archive limited to scheme implementation. Its OBJECT providers are
# linked explicitly by executables, rather than copied into this archive.
target_link_libraries(simulation PUBLIC beren_grid beren_json beren_eigen OpenMP::OpenMP_CXX
    PRIVATE configuration)
add_subdirectory(ecsim)
add_subdirectory(ecsim_corr)
