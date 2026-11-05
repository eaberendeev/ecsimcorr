add_library(beren_core INTERFACE)
target_include_directories(beren_core INTERFACE ${PROJECT_SOURCE_DIR})
target_compile_features(beren_core INTERFACE cxx_std_20)
