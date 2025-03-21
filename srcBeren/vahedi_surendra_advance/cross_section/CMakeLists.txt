project(cross_section)

get_filename_component(PARENT_DIR ../../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB cross_section.cpp)

add_library(cross_section STATIC ${SOURCE_LIB})