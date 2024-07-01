project(Utils LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  parameters_map.cpp
                Read.cpp 
                service.cpp 
                output_util.cpp
                voxel_traversal.cpp 
                World.cpp)

add_library(Utils OBJECT ${SOURCE_LIB})
