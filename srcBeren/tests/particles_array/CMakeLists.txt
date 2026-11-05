add_executable(particles_array main.cpp)
target_link_libraries(particles_array PRIVATE beren_build_options simulation coupling simulation_boundaries
    fields particles utils beren_runtime configuration)
install(TARGETS particles_array DESTINATION ${BIN_PATH})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/config.json DESTINATION ${BIN_PATH})
