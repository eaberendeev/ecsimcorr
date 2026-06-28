add_subdirectory(collision)
add_subdirectory(coulomb)
add_subdirectory(unit)
add_subdirectory(particles_array)

add_custom_target(run_domain
    COMMAND ${BIN_PATH}/domain.exe
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running domain unit tests"
)
add_custom_target(run_collision
    COMMAND ${BIN_PATH}/collision.exe
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running collision tests"
)
add_custom_target(run_coulomb
    COMMAND ${BIN_PATH}/coulomb.exe
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running Coulomb collision tests"
)
add_custom_target(run_particles_array
    COMMAND ${BIN_PATH}/particles_array.exe
    WORKING_DIRECTORY ${BIN_PATH}
    COMMENT "Running particles array tests"
)
add_custom_target(run_tests
    DEPENDS run_domain run_collision
    COMMENT "Running all tests"
)

