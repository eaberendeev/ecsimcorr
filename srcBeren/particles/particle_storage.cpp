#include "particle_storage.h"

#include <atomic>
#include <iterator>

#include "timer.h"

void ParticleStorage::add_particle(const Particle& particle) {
    Vector3I cell_id = domain_.get_cell_index(particle.coord);

    cells_(cell_id.x(), cell_id.y(), cell_id.z()).push_back(particle);
}

void ParticleStorage::add_particles(const std::vector<Particle>& particles) {
    for (const auto& particle : particles) {
        add_particle(particle);
    }
}

void ParticleStorage::update_cells(const Domain& domain) {
    RECORD_TIMER;

    constexpr int colorStep = 3;

    const int nx = cells_.size().x();
    const int ny = cells_.size().y();
    const int nz = cells_.size().z();

    constexpr int64_t locksCount = 1 << 16;
    std::vector<std::atomic<int64_t>> locks(locksCount);
#pragma omp parallel for
    for (int i = 0; i < std::ssize(locks); ++i) {
        locks[i].store(0);
    }

#pragma omp parallel
    {
        for (int colorX = 0; colorX < colorStep; ++colorX) {
            for (int colorY = 0; colorY < colorStep; ++colorY) {
                for (int colorZ = 0; colorZ < colorStep; ++colorZ) {
#pragma omp for schedule(static) collapse(3)
                    for (int ix = colorX; ix < nx; ix += colorStep) {
                        for (int iy = colorY; iy < ny; iy += colorStep) {
                            for (int iz = colorZ; iz < nz; iz += colorStep) {
                                std::vector<Particle>& cell_particles = cells_(ix, iy, iz);
                                int ip = 0;

                                while (ip < static_cast<int>(cell_particles.size())) {
                                    Particle particle = cell_particles[ip];
                                    const Vector3I cell_id = domain.get_cell_index(particle.coord);

                                    const int dx = std::abs(cell_id.x() - ix);
                                    const int dy = std::abs(cell_id.y() - iy);
                                    const int dz = std::abs(cell_id.z() - iz);
                                    if (dx == 0 && dy == 0 && dz == 0) {
                                        ++ip;
                                    } else if (dx > 1 || dy > 1 || dz > 1) [[unlikely]] {
                                        std::cerr << "particle move error " << particle << std::endl;
                                        exit(1);
                                    } else {
                                        std::swap(cell_particles[ip], cell_particles.back());
                                        cell_particles.pop_back();
                                        const auto [ix2, iy2, iz2] = cell_id.split();

                                        std::atomic<int64_t>& lock =
                                            locks.at(domain_.sind(ix2, iy2, iz2) & (locksCount - 1));
                                        int64_t expected = 0;
                                        while (!lock.compare_exchange_strong(expected, 1)) {
                                            expected = 0;
                                            asm volatile("pause");
                                        }
                                        cells_(ix2, iy2, iz2).push_back(particle);
                                        lock.store(0);
                                    }
                                }
                            }
                        }
                    }   // sync
                }
            }
        }
    }
}
