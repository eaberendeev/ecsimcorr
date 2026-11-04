project(utils LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  sparse.cpp
                util.cpp
                voxel_traversal.cpp)

add_library(utils OBJECT ${SOURCE_LIB})
