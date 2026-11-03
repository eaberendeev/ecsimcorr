add_subdirectory(neutral)

add_library(collisions OBJECT collision.cpp coulomb.cpp neutral_collision.cpp collision_manager.cpp)
target_link_libraries(collisions PUBLIC beren_grid adv_vs94)
