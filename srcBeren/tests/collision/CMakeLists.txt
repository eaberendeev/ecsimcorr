project(collision LANGUAGES CXX)

set(SOURCE_EXE  main.cpp)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS collision DESTINATION ${BIN_PATH})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_config.json
     DESTINATION ${BIN_PATH})
