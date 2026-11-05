target_sources(simulation PRIVATE
    algorithms_ecsim.cpp
    simulation_ecsim.cpp
    simulation_ecsim_base.cpp
    field_system.cpp
    response/response.cpp
    response/response_assembly.cpp
    response/response_kernels.cpp
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set_source_files_properties(response/response_kernels.cpp TARGET_DIRECTORY simulation
        PROPERTIES COMPILE_FLAGS "-Ofast -funroll-loops")
endif()
