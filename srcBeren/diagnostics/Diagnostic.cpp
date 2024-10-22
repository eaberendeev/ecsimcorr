#include "Diagnostic.h"

#include "output_util.h"

void Diagnostics::make_folders() const {
    create_directory(".//Anime");
    create_directory(".//Fields");
    create_directory(".//Fields//Diag3D");
    create_directory(".//Fields//Diag2D");
    create_directory(".//Fields//Diag1D");
    create_directory(".//Recovery");
    create_directory(".//Recovery//Fields");
    create_directory(".//Recovery//Particles");
    create_directory(".//Particles");
    for (const auto& name : particleNames) {
        create_directory(".//Particles//" + name);
        create_directory(".//Particles//" + name + "//Diag3D");
        create_directory(".//Particles//" + name + "//Diag2D");
        create_directory(".//Particles//" + name + "//Diag1D");
        create_directory(".//Recovery//Particles//" + name);
    }
    std::cerr << "Folders for output has been created: SUCCESS\n";
}
