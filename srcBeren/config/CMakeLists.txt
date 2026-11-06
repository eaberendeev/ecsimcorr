add_library(configuration STATIC
    domain_config.cpp
    boundary_config.cpp
    particle_config.cpp
    field_config.cpp
    collision_config.cpp
)
target_link_libraries(configuration PUBLIC beren_grid beren_json beren_eigen OpenMP::OpenMP_CXX)
# Particle construction still calls configuration, whose factories construct
# particle sources, boundary laws, fields and collisions. Do not embed those
# OBJECT providers into this archive: each executable lists its required objects.
# domain_config uses only its own archive member and remains a standalone client.
