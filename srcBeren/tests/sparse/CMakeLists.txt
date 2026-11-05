add_executable(sparse main.cpp)
target_link_libraries(sparse PRIVATE beren_build_options utils beren_runtime)
install(TARGETS sparse DESTINATION ${BIN_PATH})
