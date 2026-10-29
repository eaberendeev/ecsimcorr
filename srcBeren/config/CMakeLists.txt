add_library(Configuration STATIC DomainConfig.cpp)
target_link_libraries(Configuration PUBLIC BerenGrid)
target_include_directories(Configuration SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/3rd_party)
