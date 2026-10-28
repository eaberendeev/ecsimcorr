#pragma once

#include <vector>

#include "containers.h"

// Plane coordinates retain the existing coord/cell_size -> storage-index
// convention, including ghost nodes in the files.
struct DiagnosticPlanes {
    std::vector<double> x, y, z;
};

// Validated once by OutputFactory. Missing interval disables snapshots;
// average_steps == 0 disables averages. The two modes are independent.
struct SliceOutputSettings {
    DiagnosticPlanes planes;
    int interval = 0;
    int average_steps = 0;
};

// Read-only selection bound to fields owned by Simulation. E/B optionally add
// a prescribed field; the selection owns neither field data nor file state.
struct OutputField {
    const char* basename;
    const Field3d* field;
    const Field3d* extra = nullptr;
};
