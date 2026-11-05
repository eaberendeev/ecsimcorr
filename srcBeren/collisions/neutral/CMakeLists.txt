add_library(adv_vs94 STATIC collisions_with_neutrals.cpp cross_section.cpp collision_processing.cpp)
target_link_libraries(adv_vs94 PUBLIC beren_core PRIVATE OpenMP::OpenMP_CXX)
