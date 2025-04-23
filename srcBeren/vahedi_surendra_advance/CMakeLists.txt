

project(advVS94 LANGUAGES CXX)


get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

include_directories(collision_processing)
include_directories(utils)

include_directories(../utils)

add_subdirectory(./collision_processing)

set(SOURCE_LIB collisions_with_neutrals.cpp collision_handler.cpp cross_section.cpp collision_utils.cpp)

add_library(advVS94 STATIC ${SOURCE_LIB})

#add_library(advVS94 OBJECT ${SOURCE_LIB} ${LINK_LIBS} clsVec collision_handler collision_processing cross_section utils)
target_link_libraries(advVS94 
    PUBLIC 
    collision_processing
)