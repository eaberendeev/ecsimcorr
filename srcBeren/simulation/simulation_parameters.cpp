// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_parameters.h"

#include <string.h>

#include <fstream>
#include <iostream>

#include "parameters_map.h"
#include "service.h"

SimulationParameters::SimulationParameters(const std::string& fParamsName) {
    const ParametersMap data = load_parameters(fParamsName);
    if (data.is_empty()) {
        std::cerr << "Reading parameters from file '" + fParamsName +
                         "': FAILED. \n";
        exit(0);
    }
    if (!set_parameters(data)) {
        std::cerr << "Loading parameters: FAILED. \n"
                  << "Parameters is not set!\n ";
        exit(0);
    }
}

bool SimulationParameters::set_parameters(const ParametersMap& data) {
    bool success = true;
    try {
        std::cerr << "Parameters: set NumCells...\n";
        mNumCells = Vector3i(int_value(data.get_double("CellsX", 0)),
                             int_value(data.get_double("CellsY", 0)),
                             int_value(data.get_double("CellsZ", 0)));
        
        std::cerr << "Parameters: set GhostCells...\n";
        mGhostCells = Vector3i(int_value(data.get_double("GhostCells", 0)),
                               int_value(data.get_double("GhostCells", 0)),
                               int_value(data.get_double("GhostCells", 0)));

        mTotalCells = mNumCells + 2 * mGhostCells;

        std::cerr << "Parameters: set Dx, Dy, Dz...\n";
        mDx = data.get_double("Dx", 0);
        mDy = data.get_double("Dy", 0);
        mDz = data.get_double("Dz", 0);
        mCellSize = Vector3d(mDx, mDy, mDz);

        std::cerr << "Parameters: set Dt...\n";
        mDt = data.get_double("Dt", 0);

        std::cerr << "Parameters: set MaxTimestep...\n";
        mMaxTimestep = int_value(
            data.get_double("MaxTime", 0) / data.get_double("Dt", 0) + 1);

        std::cerr << "Parameters: set TimeDelayDiag2D...\n";
        mTimeDelayDiag2D = data.get_double("TimeDelayDiag2D", 0);

        std::cerr << "Parameters: set UniformB0...\n";
        mUniformB0 = Vector3d(data.get_double("UniformB0", 0),
                              data.get_double("UniformB0", 1),
                              data.get_double("UniformB0", 2));

        std::cerr << "Parameters: set SliceFieldsPlaneX, Y, Z...\n";
        mSliceFieldsPlaneX =
            string_to_doubles(data.get_values("SliceFieldsPlaneX"));
        mSliceFieldsPlaneY =
            string_to_doubles(data.get_values("SliceFieldsPlaneX"));
        mSliceFieldsPlaneZ =
            string_to_doubles(data.get_values("SliceFieldsPlaneX"));

    } catch (const std::out_of_range& e) {
        success = false;
    }
    return success;
}
