add_executable(particles_array.exe main.cpp)

target_link_libraries(particles_array.exe ${LIBRARIES})
install(TARGETS particles_array.exe DESTINATION ${BIN_PATH})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/config.json DESTINATION ${BIN_PATH})

message(STATUS INSTALLED to ${BIN_PATH})
