add_library(utils OBJECT
    sparse.cpp
    util.cpp
    voxel_traversal.cpp
)
target_link_libraries(utils PUBLIC beren_build_options beren_core beren_eigen beren_runtime)
