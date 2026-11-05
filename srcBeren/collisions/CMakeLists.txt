add_subdirectory(neutral)
add_library(collisions OBJECT
    coulomb.cpp
    neutral_collision.cpp
)
target_link_libraries(collisions PUBLIC beren_build_options particles adv_vs94 beren_runtime)
