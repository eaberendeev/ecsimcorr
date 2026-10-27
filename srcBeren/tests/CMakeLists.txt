add_subdirectory(collision)
add_subdirectory(coulomb)
add_subdirectory(unit)
add_subdirectory(particles_array)
add_subdirectory(sparse)
add_subdirectory(lmat2_assembly)

set(UNIT_TESTS
    domain
    er0_dirichlet
    path_policy
    averaging
    diagnostic_outputs
    reflected_current
    lmat2_assembly
    solver_invariants
    simulation_invariants
    coulomb_invariants
    kernel_invariants
    operator_invariants
    neutral_collision_invariants
    recovery_fields
)
set(VALIDATION_TESTS
    collision
    coulomb
    particles_array
    sparse
)

foreach(test_name IN LISTS UNIT_TESTS VALIDATION_TESTS)
    set(test_workdir "${CMAKE_BINARY_DIR}/test_results/${test_name}")
    file(MAKE_DIRECTORY "${test_workdir}")
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    set_tests_properties(${test_name} PROPERTIES
        WORKING_DIRECTORY "${test_workdir}"
        ENVIRONMENT "OMP_NUM_THREADS=2"
        PROCESSORS 2
    )
    if(test_name IN_LIST UNIT_TESTS)
        set_tests_properties(${test_name} PROPERTIES LABELS unit TIMEOUT 60)
    else()
        set_tests_properties(${test_name} PROPERTIES LABELS validation TIMEOUT 3600)
    endif()
    add_custom_target(run_${test_name}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "^${test_name}$"
        DEPENDS ${test_name}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        USES_TERMINAL
    )
endforeach()

set_tests_properties(collision coulomb PROPERTIES LABELS "validation;experiment")
set_tests_properties(sparse PROPERTIES LABELS "validation;benchmark")

configure_file(particles_array/config.json "${CMAKE_BINARY_DIR}/test_results/particles_array/config.json" COPYONLY)
configure_file(collision/test_config.json "${CMAKE_BINARY_DIR}/test_results/collision/test_config.json" COPYONLY)

add_test(NAME solver_double COMMAND $<TARGET_FILE:solver_invariants>)
set_tests_properties(solver_double PROPERTIES LABELS unit TIMEOUT 60
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test_results/solver_invariants"
    ENVIRONMENT "OMP_NUM_THREADS=2;USE_MIXED_PRECISION=0")

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME python_checks COMMAND ${Python3_EXECUTABLE} -m unittest discover
        -s "${CMAKE_SOURCE_DIR}/../tests/integration" -p "test_*.py")
    set_tests_properties(python_checks PROPERTIES LABELS unit TIMEOUT 30)
endif()

add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L unit
    DEPENDS ${UNIT_TESTS}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    USES_TERMINAL
)
add_custom_target(run_validation
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L validation
    DEPENDS ${VALIDATION_TESTS}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    USES_TERMINAL
)

