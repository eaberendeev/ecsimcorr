project(beren3d_implicit.exe LANGUAGES CXX)

set(SOURCE_EXE  main.cpp
                simulation_implicit.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS beren3d_implicit.exe DESTINATION ${BIN_PATH})
