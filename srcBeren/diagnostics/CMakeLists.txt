project(Diagnostics LANGUAGES CXX)

set(SOURCE_LIB
    Diagnostic.cpp
    Timer.cpp
    OutputFactory.cpp
    AveragedField2D.cpp
    outputs/TableOutputs.cpp
    outputs/SliceOutputs.cpp
    outputs/SliceFiles.cpp
    outputs/StateOutputs.cpp
    outputs/ProbeOutput.cpp
    outputs/EnergySpectrumOutput.cpp
)

add_library(Diagnostics OBJECT ${SOURCE_LIB})
