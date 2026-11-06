add_library(particles OBJECT
    particle_array.cpp
    second_emission_model.cpp
    particle_sources.cpp
    particle_moments.cpp
)
target_link_libraries(particles PUBLIC beren_grid utils beren_json PRIVATE configuration)
