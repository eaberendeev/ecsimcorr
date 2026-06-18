add_executable(sparse.exe main.cpp)

target_link_libraries(sparse.exe ${LIBRARIES})
install(TARGETS sparse.exe DESTINATION ${BIN_PATH})

message(STATUS INSTALLED to ${BIN_PATH})
