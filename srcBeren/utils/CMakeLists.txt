project(Utils LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  Read.cpp
                sparse/matrices.cpp
                sparse/sparse.cpp
                output_util.cpp
                util.cpp
                voxel_traversal.cpp
                boundary_conditions/boundary_condition_handler.cpp
                boundary_conditions/bphi_condition.cpp
                boundary_conditions/electron_reflection_condition.cpp
                boundary_conditions/er0_condition.cpp
                boundary_conditions/open_boundary_conditions.cpp
                boundary_conditions/periodic_boundary_conditions.cpp
                boundary_conditions/second_emission_condition.cpp
                $<$<BOOL:${USE_TIMERS}>:timer.cpp>)

add_library(Utils OBJECT ${SOURCE_LIB})
