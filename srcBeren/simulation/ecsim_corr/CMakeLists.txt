set(ECSIM_SOURCES
    simulation_ecsim_corr.cpp
    algorithms_ecsim_corr.cpp
)

target_sources(simulation PRIVATE ${ECSIM_SOURCES})

message(STATUS "ECSIM scheme sources added: ${ECSIM_SOURCES}")
