project(particles LANGUAGES CXX)

set(SOURCE_LIB  collision.cpp
                collision_manager.cpp
                particle_array.cpp
                particle_boundaries.cpp
                particle_sources.cpp
                particle_moments.cpp)

add_library(particles OBJECT ${SOURCE_LIB})

target_link_libraries(particles PRIVATE configuration)
