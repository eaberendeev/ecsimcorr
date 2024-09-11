// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <iostream>

#include "recovery.h"
#include "simulation.h"
#include "util.h"
#include "service.h"

void Simulation::output_all(const Mesh& mesh,
                            const std::vector<ParticlesArray>& species,
                            int timestep) {
    // Recovery 
    for (auto& sp : species) {
        write_particles_to_recovery(sp, timestep, parameters.get_int("RecoveryInterval"));
    }
    write_fields_to_recovery(mesh.fieldEn, mesh.fieldB, timestep,
                             parameters.get_int("RecoveryInterval"));
}

void Simulation::make_folders() const {
    create_directory(".//Anime");
    create_directory(".//Fields");
    create_directory(".//Fields//Diag3D");
    create_directory(".//Fields//Diag2D");
    create_directory(".//Fields//Diag1D");
    create_directory(".//Recovery");
    create_directory(".//Recovery//Fields");
    create_directory(".//Recovery//Particles");
    create_directory(".//Performance");
    create_directory(".//Particles");
    for (const auto& sp : species) {
        create_directory(".//Particles//" + sp.name());
        create_directory(".//Particles//" + sp.name() + "//Diag3D");
        create_directory(".//Particles//" + sp.name() + "//Diag2D");
        create_directory(".//Particles//" + sp.name() + "//Diag1D");
        create_directory(".//Recovery//Particles//" + sp.name());
    }
    std::cerr << "Folders for output has been created: SUCCESS\n";
}

    // void Simulation::output_fields2D(const int timestep) {
    //     const int timestepDelay2D =
    //         int_value(parameters.time_delay_diag2D() / parameters.dt());

    //     if (timestep % timestepDelay2D != 0)
    //         return;

    //     const int delay = std::max(0, timestep / timestepDelay2D);

    //     const std::string startFilename = ".//Data//Fields//Diag2D//";
    //     const std::string sNumber = to_string(delay, 3);

    //     for (auto coordX : parameters.slice_fields_planeX()) {
    //         const int indCoordX = int_value(coordX / parameters.dx());
    //         const std::string pos = to_string(indCoordX, 3);

    //         const int axis = 0;
    //         mesh.output_field(fieldE, indCoordX, axis, startFilename +
    //         "FieldE",
    //                           sNumber);
    //         mesh.output_field(fieldB, indCoordX, axis, startFilename +
    //         "FieldB",
    //                           sNumber);
    //     }
    //     for (auto coordY : parameters.slice_fields_planeY()) {
    //         const int indCoordY = int_value(coordY / parameters.dy());
    //         const std::string pos = to_string(indCoordY, 3);

    //         const int axis = 1;
    //         mesh.output_field(fieldE, indCoordY, axis, startFilename +
    //         "FieldE",
    //                           sNumber);
    //         mesh.output_field(fieldB, indCoordY, axis, startFilename +
    //         "FieldB",
    //                           sNumber);
    //     }
    //     for (auto coordZ : parameters.slice_fields_planeZ()) {
    //         const int indCoordZ = int_value(coordZ / parameters.dz());
    //         const std::string pos = to_string(indCoordZ, 3);

    //         const int axis = 2;
    //         mesh.output_field(fieldE, indCoordZ, axis, startFilename +
    //         "FieldE",
    //                           sNumber);
    //         mesh.output_field(fieldB, indCoordZ, axis, startFilename +
    //         "FieldB",
    //                           sNumber);
    //     }
    // }

    // void Simulation::output_particles2D(const int timestep) {
    //     const int timestepDelay2D =
    //         int_value(parameters.time_delay_diag2D() / parameters.dt());

    //     if (timestep % timestepDelay2D != 0)
    //         return;

    //     const int delay = std::max(0, timestep / timestepDelay2D);
    //     for (const auto& sp : species) {
    //         const std::string startFilename =
    //             ".//Data//Particles//" + sp.name() + "//Diag2D//";
    //         const std::string sNumber = to_string(delay, 3);

    //         for (auto coordX : parameters.slice_fields_planeX()) {
    //             const int indCoordX = int_value(coordX / parameters.dx());
    //             const std::string pos = to_string(indCoordX, 3);

    //             const int axis = 0;
    //             sp.output_density_on_grid(indCoordX, axis,
    //                                       startFilename + "Density",
    //                                       sNumber);
    //         }
    //         for (auto coordY : parameters.slice_fields_planeY()) {
    //             const int indCoordY = int_value(coordY / parameters.dy());
    //             const std::string pos = to_string(indCoordY, 3);

    //             const int axis = 1;
    //             sp.output_density_on_grid(indCoordY, axis,
    //                                       startFilename + "Density",
    //                                       sNumber);
    //         }
    //         for (auto coordZ : parameters.slice_fields_planeZ()) {
    //             const int indCoordZ = int_value(coordZ / parameters.dz());
    //             const std::string pos = to_string(indCoordZ, 3);

    //             const int axis = 2;
    //             sp.output_density_on_grid(indCoordZ, axis,
    //                                       startFilename + "Density",
    //                                       sNumber);
    //         }
    //     }
    // }
