add_library(BerenCore INTERFACE)
target_include_directories(BerenCore INTERFACE ${PROJECT_SOURCE_DIR})
target_compile_features(BerenCore INTERFACE cxx_std_20)
