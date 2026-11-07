# Manual full-size benchmark; initialization and warm-up are excluded.
add_executable(step_timing main.cpp)
target_link_libraries(step_timing PRIVATE
    simulation
    diagnostics
    fields
    particles
    coupling
    utils
    simulation_boundaries
    collisions
    beren_io
    beren_runtime
    configuration
)
