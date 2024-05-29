project(beren3d.exe LANGUAGES CXX)

set(SOURCE_EXE  main.cpp
                simulation.cpp
                recovery.cpp
                simulation_output.cpp
                simulation_init.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} PETSC mpich ${LIBRARIES})
install(TARGETS beren3d.exe DESTINATION ${BIN_PATH})
