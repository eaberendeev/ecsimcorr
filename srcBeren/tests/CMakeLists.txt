find_package(Python3 COMPONENTS Interpreter)

add_subdirectory(coulomb)
add_subdirectory(unit)
add_subdirectory(particles_array)
add_subdirectory(sparse)
add_subdirectory(lmat2_assembly)

set(UNIT_TESTS
    collision_configuration
    field_configuration
    domain
    grid_nodes
    domain_config
    er0_dirichlet
    path_policy
    averaging
    diagnostic_outputs
    diagnostic_density
    reflected_current
    lmat2_assembly
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
set(VALIDATION_TESTS
    particles_array
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

# A separate process keeps cached environment options independent of the window test.
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/test_results/profiler_initialization")
add_test(NAME profiler_initialization COMMAND $<TARGET_FILE:profiler_window> --include-zero)
set_tests_properties(profiler_initialization PROPERTIES
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test_results/profiler_initialization"
    ENVIRONMENT "OMP_NUM_THREADS=2"
    PROCESSORS 2
    LABELS "unit;runtime"
    TIMEOUT 60
)

# Full experiments stay explicit. Coulomb data are checked numerically when
# Python is available; the bounded relaxation case is registered separately.
foreach(experiment IN ITEMS coulomb sparse)
    set(experiment_workdir "${CMAKE_BINARY_DIR}/test_results/${experiment}")
    file(MAKE_DIRECTORY "${experiment_workdir}")
    if(experiment STREQUAL "coulomb" AND Python3_Interpreter_FOUND)
        add_custom_target(run_coulomb
            COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=2 $<TARGET_FILE:coulomb>
            COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/coulomb/check_relaxation.py"
            DEPENDS coulomb
            WORKING_DIRECTORY "${experiment_workdir}"
            USES_TERMINAL
        )
    else()
        add_custom_target(run_${experiment}
            COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=2 $<TARGET_FILE:${experiment}>
            DEPENDS ${experiment}
            WORKING_DIRECTORY "${experiment_workdir}"
            USES_TERMINAL
        )
    endif()
endforeach()

set_tests_properties(particle_motion kernel_invariants reflected_current simulation_scenarios simulation_invariants
    coulomb_invariants neutral_collision_invariants PROPERTIES LABELS "unit;physics")
set_tests_properties(wall_scenarios PROPERTIES LABELS "unit;physics;geometry_boundary")
set_tests_properties(domain domain_config path_policy er0_dirichlet PROPERTIES LABELS "unit;geometry_boundary")
set_tests_properties(field_configuration operator_invariants solver_invariants lmat2_assembly field_solve PROPERTIES LABELS "unit;numerics")
set_tests_properties(averaging diagnostic_outputs recovery_fields simulation_restart PROPERTIES LABELS "unit;io")
set_tests_properties(diagnostic_density PROPERTIES LABELS "unit;physics;io")
set_tests_properties(species_invariants PROPERTIES LABELS "unit;physics;geometry_boundary;sources")
set_tests_properties(environment_options profiler_window PROPERTIES LABELS "unit;runtime")
set_tests_properties(collision_configuration coulomb_invariants neutral_collision_invariants
    PROPERTIES LABELS "unit;physics;collisions")
# Explicit physical names; the old targets remain for existing scripts.
add_custom_target(run_coulomb_relaxation_experiment DEPENDS run_coulomb)

configure_file(particles_array/config.json "${CMAKE_BINARY_DIR}/test_results/particles_array/config.json" COPYONLY)

add_test(NAME solver_double COMMAND $<TARGET_FILE:solver_invariants>)
set_tests_properties(solver_double PROPERTIES LABELS "unit;numerics" TIMEOUT 60
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test_results/solver_invariants"
    ENVIRONMENT "OMP_NUM_THREADS=2;USE_MIXED_PRECISION=0")

if(Python3_Interpreter_FOUND)
    add_test(NAME physics_scenarios
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/../tests/integration/validate_physics.py"
            --binary $<TARGET_FILE:beren3d> --profile smoke --threads 2)
    set_tests_properties(physics_scenarios PROPERTIES LABELS "unit;physics;integration" TIMEOUT 120 PROCESSORS 2)
    add_custom_target(run_physics_validation
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -V -R "^physics_scenarios$"
        DEPENDS beren3d WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" USES_TERMINAL)
    add_test(NAME architecture_dependencies
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/../tests/check_architecture.py"
            --source-dir "${CMAKE_SOURCE_DIR}")
    add_test(NAME architecture_checker
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/../tests/test_architecture.py")
    set_tests_properties(architecture_dependencies architecture_checker PROPERTIES
        LABELS "unit;architecture" TIMEOUT 30)
    add_custom_target(check_architecture
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/../tests/check_architecture.py"
            --source-dir "${CMAKE_SOURCE_DIR}"
        VERBATIM)

    set(relaxation_workdir "${CMAKE_BINARY_DIR}/test_results/coulomb_relaxation")
    file(MAKE_DIRECTORY "${relaxation_workdir}")
    add_test(NAME coulomb_relaxation
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/coulomb/check_relaxation.py"
            --binary $<TARGET_FILE:coulomb>)
    set_tests_properties(coulomb_relaxation PROPERTIES
        LABELS "validation;physics;collisions" TIMEOUT 180 PROCESSORS 2
        ENVIRONMENT "OMP_NUM_THREADS=2" WORKING_DIRECTORY "${relaxation_workdir}")
    add_custom_target(run_coulomb_relaxation
        COMMAND ${CMAKE_CTEST_COMMAND} --verbose -R "^coulomb_relaxation$"
        DEPENDS coulomb WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" USES_TERMINAL)
    add_custom_target(run_collisions
        COMMAND ${CMAKE_CTEST_COMMAND} --verbose -L collisions
        DEPENDS coulomb collision_configuration coulomb_invariants neutral_collision_invariants
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" USES_TERMINAL)

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

if(Python3_Interpreter_FOUND)
    add_dependencies(run_validation coulomb)
    add_dependencies(run_tests beren3d)
endif()
