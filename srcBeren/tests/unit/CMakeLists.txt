set(UNIT_TARGETS
    collision_configuration
    field_configuration
    # domain проверяет использование только заголовков — намеренно без линковки библиотек.
    domain
    grid_nodes
    stencil_gather
    density_stencil
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
    eccapfim_invariants
    drift_kinetic_invariants
    halton_loading
    symplectic_invariants
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

# Numerical unit tests share one complete runtime link set. The two geometry
# clients deliberately prove that their API is usable without numerical objects.
foreach(test_name IN LISTS UNIT_TARGETS)
    add_executable(${test_name} ${test_name}.cpp)
    if(test_name STREQUAL "domain")
        target_link_libraries(${test_name} PRIVATE beren_grid)
    elseif(test_name STREQUAL "domain_config")
        target_link_libraries(${test_name} PRIVATE configuration)
    else()
        target_link_libraries(${test_name} PRIVATE
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
    endif()
    install(TARGETS ${test_name} DESTINATION ${BIN_PATH})
endforeach()

target_sources(diagnostic_outputs PRIVATE diagnostic_output_validation.cpp)
