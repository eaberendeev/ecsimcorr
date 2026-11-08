# Explicit experiment; never part of the fast CTest suite.
# The full-simulation runtime link set (mirrors the numerical unit tests)
# supplies the beren_core include root and the solver/diagnostic symbols.
set(WEIBEL_WHITNEY_LIBRARIES
    simulation
    diagnostics
    fields
    particles
    coupling
    utils
    simulation_boundaries
    collisions
    beren_io
    beren_runtime
    configuration
)
add_executable(weibel_whitney main.cpp)
target_link_libraries(weibel_whitney PRIVATE ${WEIBEL_WHITNEY_LIBRARIES})
# Separate binary so an experiment can finish with its original executable
# while a validated cell-assembly implementation is built alongside it.
add_executable(weibel_whitney_cell main.cpp)
target_link_libraries(weibel_whitney_cell PRIVATE ${WEIBEL_WHITNEY_LIBRARIES})
