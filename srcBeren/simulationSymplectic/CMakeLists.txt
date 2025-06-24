project(beren3d_symplectic.exe LANGUAGES CXX)

set(SOURCE_EXE  main.cpp
                particles_symplectic.cpp
                simulation_symplectic.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS beren3d_symplectic.exe DESTINATION ${BIN_PATH})
