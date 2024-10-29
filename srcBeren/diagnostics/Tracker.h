#ifndef TRACKER_H_
#define TRACKER_H_

#ifdef SET_PARTICLE_IDS

#include <fstream>

#include "ParticlesArray.h"
#include "Timer.h"
#include "World.h"
#include "util.h"

/// @brief  todo : may be use particlesName from Diagnostics
class ParticleTracker {
   public:
    ParticleTracker(const Species& species, int _count,
                    const std::string& dirname, const std::string& type = "") : count(_count) {
        create_directory(".//" + dirname);

        // open files for particle tracking
        particleTrackFiles.resize(species.size());

        for (size_t i = 0; i < species.size(); ++i) {
            const auto& sp = species[i];
            particleNames.push_back(sp->name());

            particleTrackFiles[i].reserve(count);
            for (int j = 0; j < count; ++j) {
                std::string filename = ".//" + dirname + "//" + sp->name() +
                                       "_track_" + type + "_"+ to_string(j, 3) +
                                       ".txt";
                particleTrackFiles[i].emplace_back(filename);

                if (!particleTrackFiles[i].back().is_open()) {
                    throw std::runtime_error("Failed to open file: " +
                                             filename);
                }
            }
        }
    }

    void track_particles(const Species& species,
                         int timestep) {
        for (size_t i = 0; i < species.size(); ++i) {
            for (auto k = 0; k < species[i]->size(); ++k) {
                for (auto& particle : species[i]->particlesData(k)) {
                    if (particle.id < count)
                        particleTrackFiles[i][particle.id] << timestep << " "
                                                           << particle << "\n";
                }
            }
        }
    }

    ~ParticleTracker() {
        for (auto& speciesFiles : particleTrackFiles) {
            for (auto& file : speciesFiles) {
                if (file.is_open()) {
                    file.close();
                }
            }
        }
    }

   private:
    std::vector<std::string> particleNames;
    std::vector<std::vector<std::ofstream>> particleTrackFiles;
    int count;
};

#endif   // SET_PARTICLE_IDS

#endif   // TRACKER_H_
