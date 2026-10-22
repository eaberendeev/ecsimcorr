project(domain LANGUAGES CXX)

set(SOURCE_EXE domain.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})
install(TARGETS domain DESTINATION ${BIN_PATH})

# Тесты Дирихле-маски условия er0 (требуют Utils: boundary_conditions).
add_executable(er0_dirichlet er0_dirichlet.cpp)
target_link_libraries(er0_dirichlet ${LIBRARIES})
install(TARGETS er0_dirichlet DESTINATION ${BIN_PATH})
