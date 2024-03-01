// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include <map>
#include <string>
#include <vector>

#include "util.h"
#include "parameters_map.h"

class SimulationParameters {
   public:
    SimulationParameters(const std::string& fParamsName);
    bool set_parameters(const ParametersMap& data);
    const Vector3i& num_cells() const { return mNumCells; }
    const Vector3i& ghost_cells() const { return mGhostCells; }
    const Vector3i& total_cells() const { return mTotalCells; }
    double dx() const { return mDx; }
    double dy() const { return mDy; }
    double dz() const { return mDz; }
    const Vector3d& cell_size() const { return mCellSize; }

    double dt() const { return mDt; }
    int max_timestep() const { return mMaxTimestep; }
    double time_delay_diag2D() const { return mTimeDelayDiag2D; }
    const Vector3d& uniformB0() const { return mUniformB0; }
    const std::vector<double>& slice_fields_planeX() const {
        return mSliceFieldsPlaneX;
    }
    const std::vector<double>& slice_fields_planeY() const {
        return mSliceFieldsPlaneY;
    }
    const std::vector<double>& slice_fields_planeZ() const {
        return mSliceFieldsPlaneZ;
    }

   private:
    Vector3i mNumCells;
    Vector3i mGhostCells;
    Vector3i mTotalCells;
    double mDx, mDy, mDz;
    Vector3d mCellSize;
    double mDt;
    int mMaxTimestep;
    double mTimeDelayDiag2D;
    Vector3d mUniformB0;
    std::vector<double> mSliceFieldsPlaneX;
    std::vector<double> mSliceFieldsPlaneY;
    std::vector<double> mSliceFieldsPlaneZ;
};

#endif   // SIMULATION_PARAMETERS_H
