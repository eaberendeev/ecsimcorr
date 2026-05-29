project(domain.exe LANGUAGES CXX)

set(SOURCE_EXE domain.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_EXE})
install(TARGETS domain.exe DESTINATION ${BIN_PATH})
