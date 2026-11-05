add_library(coupling OBJECT density.cpp current.cpp)
target_link_libraries(coupling PUBLIC particles utils)
