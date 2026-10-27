#include "Diagnostic.h"

#include "util.h"

void Diagnostics::make_folders() const {
    ensure_directory(".//Fields");
    ensure_directory(".//Fields//Diag2D");
    ensure_directory(".//Fields//Diag3D");
    ensure_directory(".//Recovery");
    ensure_directory(".//Recovery//Fields");
    ensure_directory(".//Recovery//Particles");
    ensure_directory(".//Particles");
    for (const auto& name : particleNames) {
        ensure_directory(".//Particles//" + name);
        ensure_directory(".//Particles//" + name + "//Diag2D");
        ensure_directory(".//Recovery//Particles//" + name);
    }
}
