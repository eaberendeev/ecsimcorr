project(utils)

get_filename_component(PARENT_DIR ../../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB utils.cpp)

add_library(utils STATIC ${SOURCE_LIB})