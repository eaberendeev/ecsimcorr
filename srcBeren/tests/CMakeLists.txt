add_subdirectory(collision)
add_subdirectory(coulomb)
add_subdirectory(unit)
add_subdirectory(particles_array)
add_subdirectory(sparse)
add_subdirectory(solver)

add_custom_target(run_domain
    COMMAND ${BIN_PATH}/domain
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running domain unit tests"
)
add_custom_target(run_collision
    COMMAND ${BIN_PATH}/collision
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running collision tests"
)
add_custom_target(run_coulomb
    COMMAND ${BIN_PATH}/coulomb
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running Coulomb collision tests"
)
add_custom_target(run_particles_array
    COMMAND ${BIN_PATH}/particles_array
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running particles array tests"
)
add_custom_target(run_sparse
    COMMAND ${BIN_PATH}/sparse
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running test of sparse matrix utils"
)
add_custom_target(run_solver
    COMMAND ${BIN_PATH}/solver
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running solver stage-1 tests"
)
add_custom_target(run_tests
    DEPENDS run_domain run_collision run_solver
    COMMENT "Running all tests"
)

