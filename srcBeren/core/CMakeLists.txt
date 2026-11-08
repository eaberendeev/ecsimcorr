add_library(beren_core INTERFACE)
target_include_directories(beren_core INTERFACE ${PROJECT_SOURCE_DIR})

target_compile_definitions(beren_core INTERFACE BEREN_GHOST_CELLS=${BEREN_GHOST_CELLS})
