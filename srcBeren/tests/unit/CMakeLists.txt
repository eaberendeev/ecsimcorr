project(domain LANGUAGES CXX)

set(UNIT_TARGETS
    # domain проверяет использование только заголовков — намеренно без линковки библиотек.
    domain
    # Тесты Дирихле-маски условия er0 (требуют utils: boundary_conditions).
    er0_dirichlet
    # Тесты политики обрезки траектории частицы (требуют utils: boundary_conditions,
    # particle_path).
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
    if(NOT test_name STREQUAL "domain")
        target_link_libraries(${test_name} ${LIBRARIES})
    endif()
    install(TARGETS ${test_name} DESTINATION ${BIN_PATH})
endforeach()

target_link_libraries(simulation_invariants configuration)
target_sources(diagnostic_outputs PRIVATE diagnostic_output_validation.cpp)

foreach(config_client IN ITEMS path_policy er0_dirichlet reflected_current)
    target_link_libraries(${config_client} configuration)
endforeach()
