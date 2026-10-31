set(ECSIM_SOURCES
    algorithms_ecsim.cpp
    simulation_ecsim.cpp
    field_system.cpp
    response/response.cpp
    response/response_assembly.cpp
    response/response_kernels.cpp
)

target_sources(Simulation PRIVATE ${ECSIM_SOURCES})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(response/response_kernels.cpp TARGET_DIRECTORY Simulation
        PROPERTIES COMPILE_FLAGS "-Ofast -funroll-loops")
endif()

target_include_directories(Simulation 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "ECSIM scheme sources added: ${ECSIM_SOURCES}")
