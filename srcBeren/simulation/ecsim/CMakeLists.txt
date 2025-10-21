set(ECSIM_SOURCES
    simulation_ecsim.cpp
)

target_sources(Simulation PRIVATE ${ECSIM_SOURCES})

target_include_directories(Simulation 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "ECSIM scheme sources added: ${ECSIM_SOURCES}")