add_executable(lmat2_assembly main.cpp)
target_link_libraries(lmat2_assembly PRIVATE beren_build_options simulation particles utils beren_runtime configuration)
install(TARGETS lmat2_assembly DESTINATION ${BIN_PATH})
