project(tracking LANGUAGES CXX)

set(SOURCE_EXE  main.cpp
                test_track_in_fieldB.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS tracking DESTINATION ${BIN_PATH})
