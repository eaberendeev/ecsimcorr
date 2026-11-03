project(particles LANGUAGES CXX)

set(SOURCE_LIB  particle_array.cpp
                second_emission_model.cpp
                particle_sources.cpp
                particle_moments.cpp)

add_library(particles OBJECT ${SOURCE_LIB})

target_link_libraries(particles PRIVATE configuration)
