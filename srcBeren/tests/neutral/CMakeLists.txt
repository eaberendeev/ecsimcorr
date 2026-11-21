add_executable(neutral_ionization main.cpp)
target_link_libraries(neutral_ionization PRIVATE collisions particles utils beren_runtime configuration)
install(TARGETS neutral_ionization DESTINATION ${BIN_PATH})
