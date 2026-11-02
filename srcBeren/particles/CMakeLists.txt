project(Particles LANGUAGES CXX)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
#include_directories(${PARENT_DIR}/service)

set(SOURCE_LIB  collision.cpp
                collision_manager.cpp
                particle_array.cpp
                ../sources/particle_sources.cpp
                particle_moments.cpp)

add_library(Particles OBJECT ${SOURCE_LIB})

target_link_libraries(Particles PRIVATE Configuration)
