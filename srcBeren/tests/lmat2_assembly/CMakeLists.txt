project(lmat2_assembly LANGUAGES CXX)

set(SOURCE_EXE main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS lmat2_assembly DESTINATION ${BIN_PATH})

message(STATUS INSTALLED to ${BIN_PATH})