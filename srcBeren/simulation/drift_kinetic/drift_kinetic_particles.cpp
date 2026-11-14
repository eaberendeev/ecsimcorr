#include "drift_kinetic_particles.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <type_traits>

namespace driftKinetic {
namespace {
constexpr char kTag[8] = {'B', 'E', 'R', 'E', 'N', 'G', 'C', '1'};
constexpr std::int32_t kVersion = 1;
static_assert(std::is_trivially_copyable_v<GuidingCenter>, "guiding centres are stored as raw records");
}   // namespace

void write_guiding_centers(const DriftKineticParticles& particles, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary | std::ios::trunc);
    if (!file)
        throw std::runtime_error("drift_kinetic: cannot write recovery file " + filename);
    const std::int64_t count = particles.get_total_num_of_particles();
    const std::int32_t record = sizeof(GuidingCenter);
    file.write(kTag, sizeof(kTag));
    file.write(reinterpret_cast<const char*>(&kVersion), sizeof(kVersion));
    file.write(reinterpret_cast<const char*>(&record), sizeof(record));
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));
    for (int cell = 0; cell < particles.size(); ++cell) {
        const auto& markers = particles(cell);
        if (!markers.empty())
            file.write(reinterpret_cast<const char*>(markers.data()),
                       static_cast<std::streamsize>(markers.size() * sizeof(GuidingCenter)));
    }
    file.close();
    if (!file)
        throw std::runtime_error("drift_kinetic: failed to write recovery file " + filename);
}

void read_guiding_centers(DriftKineticParticles& particles, const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file)
        throw std::runtime_error("drift_kinetic: recovery file not found: " + filename);
    char tag[sizeof(kTag)];
    std::int32_t version = 0, record = 0;
    std::int64_t count = -1;
    file.read(tag, sizeof(tag));
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    file.read(reinterpret_cast<char*>(&record), sizeof(record));
    file.read(reinterpret_cast<char*>(&count), sizeof(count));
    if (!file || std::memcmp(tag, kTag, sizeof(kTag)) != 0 || version != kVersion ||
        record != static_cast<std::int32_t>(sizeof(GuidingCenter)) || count < 0)
        throw std::runtime_error("drift_kinetic: not a guiding-centre recovery file: " + filename);
    for (std::int64_t n = 0; n < count; ++n) {
        GuidingCenter p;
        file.read(reinterpret_cast<char*>(&p), sizeof(p));
        if (!file)
            throw std::runtime_error("drift_kinetic: truncated recovery file " + filename);
        particles.add(p);
    }
}

}   // namespace driftKinetic
