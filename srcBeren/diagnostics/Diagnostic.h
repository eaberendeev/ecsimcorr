#ifndef DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#include <fstream>

#include "Mesh.h"
#include "ParticlesArray.h"
#include "Read.h"
#include "Timer.h"
#include "World.h"
#include "output_util.h"

class Diagnostics {
   public:
    Diagnostics(const nlohmann::json& diagnostic_config, const Domain& domain, const Species& species)
        : _domain(domain), diagnostic_config(diagnostic_config) {
        for (const auto& kv : species) {
            const auto& sp = kv.second;
            particleNames.push_back(sp->name());
        }
        make_folders();
    }

    void make_folders() const;

    // Energy tracking API
    void addEnergy(const std::string& key, double value);
    void addBoundary(const std::string& key, double value);

    // Data for outputs to read
    const Domain& _domain;
    nlohmann::json diagnostic_config;
    std::map<std::string, double> energy;
    std::vector<std::string> energyOrder;
    std::vector<std::string> boundaryOrder;
    std::vector<std::string> particleNames;
};

static inline double calc_energy_field(const Field3d& field, const IndexRange& range) {
    return 0.5 * dot_product_sum(field, field, range);
}

#endif
