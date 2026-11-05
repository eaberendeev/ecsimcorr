add_executable(coulomb main.cpp)
target_link_libraries(coulomb PRIVATE beren_build_options collisions particles utils beren_runtime configuration)
install(TARGETS coulomb DESTINATION ${BIN_PATH})
# Analysis scripts stay in the source folder; build.py runs the executable and
# scripts with test_results/coulomb as their working directory.
