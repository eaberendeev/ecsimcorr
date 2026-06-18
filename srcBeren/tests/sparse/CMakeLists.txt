add_executable(sparse main.cpp)

target_link_libraries(sparse ${LIBRARIES})
install(TARGETS sparse DESTINATION ${BIN_PATH})

message(STATUS INSTALLED to ${BIN_PATH})
