project(coulomb.exe LANGUAGES CXX)

set(SOURCE_EXE  main.cpp)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})

target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
install(TARGETS coulomb.exe DESTINATION ${BIN_PATH})
# compare_theory.py / plot_coulomb.py are run in-place from this source folder
# (build.py runs the test binary here too), so they are not copied to BIN_PATH.
