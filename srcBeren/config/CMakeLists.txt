add_library(configuration STATIC domain_config.cpp boundary_config.cpp particle_config.cpp field_config.cpp)
target_link_libraries(configuration PUBLIC beren_grid)
target_include_directories(configuration SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/3rd_party)
