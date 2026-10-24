project(domain LANGUAGES CXX)

set(SOURCE_EXE domain.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})
install(TARGETS domain DESTINATION ${BIN_PATH})

# Тесты Дирихле-маски условия er0 (требуют Utils: boundary_conditions).
add_executable(er0_dirichlet er0_dirichlet.cpp)
target_link_libraries(er0_dirichlet ${LIBRARIES})
install(TARGETS er0_dirichlet DESTINATION ${BIN_PATH})

# Тесты политики обрезки траектории частицы (требуют Utils: boundary_conditions,
# particle_path).
add_executable(path_policy path_policy.cpp)
target_link_libraries(path_policy ${LIBRARIES})
install(TARGETS path_policy DESTINATION ${BIN_PATH})

# Тесты общего хелпера усреднения диагностических 2D-полей.
add_executable(averaging averaging.cpp)
target_link_libraries(averaging ${LIBRARIES})
install(TARGETS averaging DESTINATION ${BIN_PATH})

add_executable(diagnostic_outputs diagnostic_outputs.cpp)
target_link_libraries(diagnostic_outputs ${LIBRARIES})
install(TARGETS diagnostic_outputs DESTINATION ${BIN_PATH})

add_executable(reflected_current reflected_current.cpp)
target_link_libraries(reflected_current ${LIBRARIES})
install(TARGETS reflected_current DESTINATION ${BIN_PATH})

foreach(invariant_test solver_invariants simulation_invariants coulomb_invariants kernel_invariants operator_invariants neutral_collision_invariants recovery_fields)
    add_executable(${invariant_test} ${invariant_test}.cpp)
    target_link_libraries(${invariant_test} ${LIBRARIES})
    install(TARGETS ${invariant_test} DESTINATION ${BIN_PATH})
endforeach()
