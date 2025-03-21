

project(advVS94 LANGUAGES CXX)


get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

include_directories(clsVec)
include_directories(collision_handler)
include_directories(collision_processing)
include_directories(cross_section)
include_directories(utils)

include_directories(../utils)

add_subdirectory(./clsVec)
add_subdirectory(./collision_handler)
add_subdirectory(./collision_processing)
add_subdirectory(./cross_section)
add_subdirectory(./utils)

set(SOURCE_LIB collisions_with_neutrals.cpp)

add_library(advVS94 STATIC ${SOURCE_LIB})

#add_library(advVS94 OBJECT ${SOURCE_LIB} ${LINK_LIBS} clsVec collision_handler collision_processing cross_section utils)
target_link_libraries(advVS94 
    PUBLIC 
    clsVec
    collision_handler
    collision_processing
    cross_section
    utils
)