#ifndef DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#include <fstream>

#include "Mesh.h"
#include "ParticlesArray.h"
#include "Read.h"
#include "Timer.h"
#include "World.h"
#include "util.h"
#include "Tracker.h"
#include "service.h"
#include "output_util.h"

class Diagnostics {
   public:
    Diagnostics(const ParametersMap& outputParameters, const Domain& domain,
                const Species& species)
        : _domain(domain), outParams(outputParameters) {
        outParams.print();
        fEnergy.open("energy.txt");

        for (size_t i = 0; i < species.size(); ++i) {
            const auto& sp = species[i];
            particleNames.push_back(sp->name());
        }
        make_folders();
        sliceCoordsPlaneX = outParams.get_double_values("sliceFieldsPlaneX");
        sliceCoordsPlaneY = outParams.get_double_values("sliceFieldsPlaneY");
        sliceCoordsPlaneZ = outParams.get_double_values("sliceFieldsPlaneZ");
    }
    ~Diagnostics() {
        fEnergy.close();
    }
    void make_folders() const;

    void track_particles(const std::vector<ParticlesArray>& species,
                         int timestep);

    void write_energy(const ParametersMap& parameters, int timestep);
    void output_energy_spectrum(
        const EnergySpectrum& spectrum, int timestep,
        const std::string& output_dir) const;
    template <typename T>
    void output_fields2D(
        const int timestep,
        const std::vector<std::pair<const T&, std::string>>& fields);
    void output_array2D(
        const int timestep,
        const std::vector<std::pair<Array3D<double>&, std::string>>& fields);

    //   private:
    const Domain& _domain;

    ParametersMap outParams;
    std::map<std::string, double> energy;
    std::vector<std::string> energyOrder;
    std::ofstream fEnergy;
    std::vector<std::string> particleNames;
    std::vector<double> sliceCoordsPlaneX, sliceCoordsPlaneY, sliceCoordsPlaneZ;

    void addEnergy(const std::string& key, double value);
};

template <typename T>
void output_field(const T& field, const int indCoord, const int axis,
                  const std::string& filename, const std::string& sNumber) {
    int overlap = 2 * GHOST_CELLS + 1;
    int3 overlaps = int3(overlap, overlap, overlap);
    int3 sizes = field.sizes();
    int3 start = {0, 0, 0};
    int3 end = sizes; // - overlaps;
    output_field_plane(field, start, end, indCoord, axis, field.nd(), filename,
                       sNumber);
}

template <typename T>
void Diagnostics::output_fields2D(
    const int timestep,
    const std::vector<std::pair<const T&, std::string>>& fields) {
    // static_assert(
    //     std::is_same_v<T, Field3d> || std::is_same_v<T, Array3D<double>>,
    //     "T must be either Field3d or Array<double>");
    const int timestepDelay2D = outParams.get_int("TimeStepDelayDiag2D");

    if (timestep % timestepDelay2D != 0)
        return;

    const int delay = std::max(0, timestep / timestepDelay2D);

    const std::string sNumber = to_string(delay, 4);

    auto output_for_axis = [&](const auto& coords, int axis, const double dr) {
        for (auto coord : coords) {
            const int indCoord = int_value(coord / dr);
            for (const auto& [field, fieldName] : fields) {
                output_field(field, indCoord, axis,
                             fieldName,
                                 sNumber);
            }
        }
    };
    auto cellSize = _domain.cell_size();
    output_for_axis(sliceCoordsPlaneX, 0, cellSize.x());
    output_for_axis(sliceCoordsPlaneY, 1, cellSize.y());
    output_for_axis(sliceCoordsPlaneZ, 2, cellSize.z());
}

#endif
