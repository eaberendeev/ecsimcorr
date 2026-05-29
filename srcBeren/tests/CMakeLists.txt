add_subdirectory(collision)
add_subdirectory(unit)

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
add_custom_target(run_tests
    DEPENDS run_domain run_collision
    COMMENT "Running all tests"
)

