add_library(configuration STATIC domainconfig.cpp particle_config.cpp)
target_link_libraries(configuration PUBLIC beren_grid)
target_include_directories(configuration SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/3rd_party)
