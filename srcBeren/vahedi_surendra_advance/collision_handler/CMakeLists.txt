project(collision_handler)

get_filename_component(PARENT_DIR ../../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB collision_handler.cpp)

add_library(collision_handler STATIC ${SOURCE_LIB})