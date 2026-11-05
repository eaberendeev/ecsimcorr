add_library(beren_io OBJECT
    text_input.cpp
    field_io.cpp
    recovery.cpp
)
target_link_libraries(beren_io PUBLIC beren_build_options particles utils)
