add_library(Coupling OBJECT density.cpp current.cpp)
target_link_libraries(Coupling PUBLIC BerenGrid)
