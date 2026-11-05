set(UNIT_TARGETS
    collision_configuration
    field_configuration
    # domain проверяет использование только заголовков — намеренно без линковки библиотек.
    domain
    # Domain JSON must link without field solvers, particles or simulation.
    domain_config
    # Тесты Дирихле-маски условия er0 (требуют simulation: boundary_conditions).
    er0_dirichlet
    # Тесты политики обрезки траектории частицы (требуют simulation: boundary_conditions,
    # particles/particle_path).
    path_policy
    averaging
    diagnostic_outputs
    diagnostic_density
    reflected_current
    solver_invariants
    simulation_invariants
    simulation_scenarios
    wall_scenarios
    field_solve
    coulomb_invariants
    kernel_invariants
    particle_motion
    species_invariants
    environment_options
    profiler_window
    operator_invariants
    neutral_collision_invariants
    recovery_fields
    simulation_restart
)

foreach(test_name IN LISTS UNIT_TARGETS)
    add_executable(${test_name} ${test_name}.cpp)
    install(TARGETS ${test_name} DESTINATION ${BIN_PATH})
endforeach()

# Header-only geometry and its JSON decoder intentionally avoid numerical objects.
target_link_libraries(domain PRIVATE beren_grid)
target_link_libraries(domain_config PRIVATE configuration)

target_link_libraries(environment_options PRIVATE beren_core)
target_link_libraries(particle_motion PRIVATE beren_grid utils beren_runtime beren_json)
target_link_libraries(solver_invariants PRIVATE fields utils beren_runtime)
target_link_libraries(field_configuration PRIVATE fields utils beren_runtime configuration)
target_link_libraries(recovery_fields PRIVATE beren_io particles utils beren_runtime configuration)
target_link_libraries(species_invariants PRIVATE particles utils beren_runtime configuration)

foreach(test_name IN ITEMS collision_configuration coulomb_invariants neutral_collision_invariants)
    target_link_libraries(${test_name} PRIVATE collisions particles utils beren_runtime configuration)
endforeach()
foreach(test_name IN ITEMS er0_dirichlet path_policy)
    target_link_libraries(${test_name} PRIVATE simulation_boundaries fields particles utils beren_runtime configuration)
endforeach()
foreach(test_name IN ITEMS averaging diagnostic_outputs)
    target_link_libraries(${test_name} PRIVATE diagnostics coupling beren_io simulation_boundaries
        fields particles utils beren_runtime configuration)
endforeach()
target_sources(diagnostic_outputs PRIVATE diagnostic_output_validation.cpp)

# Field-system tests use scheme algebra without constructing the simulation loop.
target_link_libraries(operator_invariants PRIVATE simulation simulation_boundaries fields
    particles utils beren_runtime configuration)

# These tests execute a scheme or its fused mover. Name the full OBJECT closure:
# linking simulation's usage requirements alone does not link provider objects.
foreach(test_name IN ITEMS diagnostic_density reflected_current simulation_invariants simulation_scenarios
        wall_scenarios field_solve kernel_invariants profiler_window simulation_restart)
    target_link_libraries(${test_name} PRIVATE simulation diagnostics fields particles coupling utils
        simulation_boundaries collisions beren_io beren_runtime configuration)
endforeach()
