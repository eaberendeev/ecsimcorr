project(domain LANGUAGES CXX)

set(UNIT_TARGETS
    # domain проверяет использование только заголовков — намеренно без линковки библиотек.
    domain
    # Тесты Дирихле-маски условия er0 (требуют Utils: boundary_conditions).
    er0_dirichlet
    # Тесты политики обрезки траектории частицы (требуют Utils: boundary_conditions,
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
    environment_options
    operator_invariants
    neutral_collision_invariants
    recovery_fields
)

foreach(test_name IN LISTS UNIT_TARGETS)
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} BerenGrid)
    if(NOT test_name STREQUAL "domain")
        target_link_libraries(${test_name} ${LIBRARIES})
    endif()
    install(TARGETS ${test_name} DESTINATION ${BIN_PATH})
endforeach()

target_link_libraries(simulation_invariants Configuration)
target_sources(diagnostic_outputs PRIVATE diagnostic_output_validation.cpp)
