add_subdirectory(neutral)

add_library(collisions OBJECT legacy_collisions.cpp coulomb.cpp neutral_collision.cpp)
target_link_libraries(collisions PUBLIC beren_grid adv_vs94)
