project(beren3d_ecsim.exe LANGUAGES CXX)

set(SOURCE_EXE  main.cpp
                simulation_ecsim_corr.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS beren3d_ecsim.exe DESTINATION ${BIN_PATH})
