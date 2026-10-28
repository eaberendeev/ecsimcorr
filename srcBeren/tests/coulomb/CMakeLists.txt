project(coulomb LANGUAGES CXX)

set(SOURCE_EXE  main.cpp)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS coulomb DESTINATION ${BIN_PATH})
# Analysis scripts stay in the source folder; build.py runs the executable and
# scripts with test_results/coulomb as their working directory.
