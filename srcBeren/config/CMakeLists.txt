add_library(Configuration STATIC domainconfig.cpp)
target_link_libraries(Configuration PUBLIC BerenGrid)
target_include_directories(Configuration SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/3rd_party)
