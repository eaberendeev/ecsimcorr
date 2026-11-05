set(SOURCE_LIB
    Diagnostic.cpp
    DiagnosticMeasurements.cpp
    DensityMeasurements.cpp
    EnergyMeasurements.cpp
    FieldSlices.cpp
    OutputFactory.cpp
    AveragedField2D.cpp
    outputs/TableOutputs.cpp
    outputs/SliceOutputs.cpp
    outputs/SliceFiles.cpp
    outputs/StateOutputs.cpp
    outputs/ProbeOutput.cpp
    outputs/EnergySpectrumOutput.cpp
)

add_library(diagnostics OBJECT ${SOURCE_LIB})

target_link_libraries(diagnostics PUBLIC particles utils simulation_boundaries
    PRIVATE coupling beren_io configuration)
