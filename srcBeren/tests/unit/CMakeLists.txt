project(domain LANGUAGES CXX)

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
    operator_invariants
    neutral_collision_invariants
    recovery_fields
)

foreach(test_name IN LISTS UNIT_TARGETS)
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} beren_grid)
    if(test_name STREQUAL "domain_config")
        target_link_libraries(${test_name} configuration)
    elseif(NOT test_name STREQUAL "domain")
        target_link_libraries(${test_name} ${LIBRARIES})
    endif()
    install(TARGETS ${test_name} DESTINATION ${BIN_PATH})
endforeach()

target_link_libraries(simulation_invariants configuration)
target_link_libraries(collision_configuration configuration)
target_sources(diagnostic_outputs PRIVATE diagnostic_output_validation.cpp)

foreach(config_client IN ITEMS path_policy er0_dirichlet reflected_current field_configuration operator_invariants)
    target_link_libraries(${config_client} configuration)
endforeach()
