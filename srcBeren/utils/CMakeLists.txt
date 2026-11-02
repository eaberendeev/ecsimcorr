project(utils LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  Read.cpp
                sparse.cpp
                output_util.cpp
                util.cpp
                voxel_traversal.cpp
                boundary_conditions/boundary_condition_handler.cpp
                $<$<BOOL:${USE_TIMERS}>:timer.cpp>)

add_library(utils OBJECT ${SOURCE_LIB})
