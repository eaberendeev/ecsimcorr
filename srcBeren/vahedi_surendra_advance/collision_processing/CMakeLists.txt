project(collision_processing)

get_filename_component(PARENT_DIR ../../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB collision_processing.cpp)

add_library(collision_processing STATIC ${SOURCE_LIB})