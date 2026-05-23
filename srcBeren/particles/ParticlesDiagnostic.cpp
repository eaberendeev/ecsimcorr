#include <functional>

#include "ParticlesArray.h"
#include "Shape.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "sgs.h"
#include "timer.h"

std::ostream& operator<<(std::ostream& out, const ParticleSimple& particle) {
    out << particle.coord << " " << particle.velocity;
    return out;
}
std::ostream& operator<<(std::ostream& out, const ParticleMass& particle) {
    out << particle.coord << " " << particle.velocity << " " << particle.mass;
    return out;
}

double PulseFromKev(double kev, double mass) {
    double gama = kev / SGS::MC2 + mass;
    return sqrt((gama * gama) - mass);
}

template <int size>
struct Buffer3d {
    static constexpr int dims = 1;

    static_assert(size > 0, "SMAX must be positive");
    alignas(64) std::array<double, dims * size * size * size> data_;

    static inline constexpr int idx(int n, int m, int k, int dim) noexcept {
        assert(0 <= n && n < size);
        assert(0 <= m && m < size);
        assert(0 <= k && k < size);
        assert(0 <= dim && dim < dims);
        return (n * (size * size) + m * size + k) * dims + dim;
    }

    double operator()(int n, int m, int k, int dim) const {
        return data_[idx(n, m, k, dim)];
    }
    double& operator()(int n, int m, int k, int dim) {
        return data_[idx(n, m, k, dim)];
    }

    Buffer3d& operator+=(const Buffer3d& other) noexcept {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }
    inline void zero() noexcept {
        data_.fill(0.0);
    }
};

// Template specializations need to be explicitly instantiated in the cpp file
template void ParticlesArray::density_on_grid_update_impl<Shape, 2>();
template void ParticlesArray::density_on_grid_update_impl<Shape2, 2>();
// Implementation of the template function
template <ParticlesArray::ShapeFunction ShapeFn, int ShapeSize>
void ParticlesArray::density_on_grid_update_impl() {
    RECORD_TIMER;

    constexpr auto SMAX = 2 * ShapeSize;
    densityOnGrid.setZero();
    // std::cout << "density_on_grid_update_impl " << SMAX << std::endl;

    const auto& lambdaRef = [&]() {
        Field3d res(densityOnGrid);

#pragma omp parallel
        {
            timer::flatTimer timer("ref OMP");
#pragma omp for schedule(dynamic, 64)
            for (auto j = 0; j < size(); ++j) {
                alignas(64) double sx[SMAX];
                alignas(64) double sy[SMAX];
                alignas(64) double sz[SMAX];

                for (const auto& particle : particlesData(j)) {
                    // Vectorizable coordinate calculations
                    const double ix = particle.coord.x() / domain_.cell_size().x();
                    const double iy = particle.coord.y() / domain_.cell_size().y();
                    const double iz = particle.coord.z() / domain_.cell_size().z();

                    const int xk = floor(ix);
                    const int yk = floor(iy);
                    const int zk = floor(iz);

// Vectorizable shape calculations
#pragma omp simd
                    for (int n = 0; n < SMAX; ++n) {
                        sx[n] = ShapeFn(-ix + double(xk - GHOST_CELLS + n));
                        sy[n] = ShapeFn(-iy + double(yk - GHOST_CELLS + n));
                        sz[n] = ShapeFn(-iz + double(zk - GHOST_CELLS + n));
                    }

                    const double weight = is_neutral() ? mpw_ : mpw_ * charge;

                    // Density accumulation with loop unrolling
                    // #pragma unroll
                    for (int n = 0; n < SMAX; ++n) {
                        const int indx = xk + n;
                        const double sxw = sx[n];

                        for (int m = 0; m < SMAX; ++m) {
                            const int indy = yk + m;
                            const double sxyw = sxw * sy[m];

                            for (int k = 0; k < SMAX; ++k) {
                                const int indz = zk + k;
#pragma omp atomic update
                                res(indx, indy, indz, 0) += weight * sxyw * sz[k];
                            }
                        }
                    }
                }
            }
        }

        return res;
    };

    const auto& lambdaTestV4 = [&]() {
        Field3d res(densityOnGrid);

#pragma omp parallel
        {
            timer::flatTimer timer("test V4 OMP");
#pragma omp for schedule(dynamic, 64)
            for (auto j = 0; j < size(); ++j) {
                alignas(64) double sx[SMAX];
                alignas(64) double sy[SMAX];
                alignas(64) double sz[SMAX];

                // std::cout<<"Start particles data block!"<<std::endl;

                Buffer3d<SMAX> localBuffer;
                localBuffer.zero();

                for (const auto& particle : particlesData(j)) {
                    // Vectorizable coordinate calculations
                    const double ix = particle.coord.x() / domain_.cell_size().x();
                    const double iy = particle.coord.y() / domain_.cell_size().y();
                    const double iz = particle.coord.z() / domain_.cell_size().z();

                    const int xk = floor(ix);
                    const int yk = floor(iy);
                    const int zk = floor(iz);

                    // std::cout<<"*x: "<<xk<<" "<<yk<<" "<<zk<<std::endl;

// Vectorizable shape calculations
#pragma omp simd
                    for (int n = 0; n < SMAX; ++n) {
                        sx[n] = ShapeFn(-ix + double(xk - GHOST_CELLS + n));
                        sy[n] = ShapeFn(-iy + double(yk - GHOST_CELLS + n));
                        sz[n] = ShapeFn(-iz + double(zk - GHOST_CELLS + n));
                    }

                    const double weight = is_neutral() ? mpw_ : mpw_ * charge;

                    // Density accumulation with loop unrolling
                    // #pragma unroll
                    for (int n = 0; n < SMAX; ++n) {
                        const double sxw = sx[n];
                        for (int m = 0; m < SMAX; ++m) {
                            const double sxyw = sxw * sy[m];
                            for (int k = 0; k < SMAX; ++k) {
                                // #pragma omp atomic update
                                localBuffer(n, m, k, 0) += weight * sxyw * sz[k];
                            }
                        }
                    }
                }

                if (particlesData(j).size() != 0) {
                    const Particle particle = particlesData(j)[0];

                    const double ix = particle.coord.x() / domain_.cell_size().x();
                    const double iy = particle.coord.y() / domain_.cell_size().y();
                    const double iz = particle.coord.z() / domain_.cell_size().z();

                    const int xk = floor(ix);
                    const int yk = floor(iy);
                    const int zk = floor(iz);

                    for (int n = 0; n < SMAX; ++n) {
                        const int indx = xk + n;
                        const double sxw = sx[n];

                        for (int m = 0; m < SMAX; ++m) {
                            const int indy = yk + m;
                            const double sxyw = sxw * sy[m];

                            for (int k = 0; k < SMAX; ++k) {
                                const int indz = zk + k;
#pragma omp atomic update
                                res(indx, indy, indz, 0) += localBuffer(n, m, k, 0);
                            }
                        }
                    }
                }

                // std::cout<<"End particles data block!"<<std::endl;
            }
        }

        return res;
    };

    timer::timer timerRef("reference");
    const Field3d ref = lambdaRef();
    timerRef.finish();

    timer::timer timerTest("test V4");
    const Field3d testV4 = lambdaTestV4();
    timerTest.finish();

    timer::timer timerOpt("test V2");
#pragma omp parallel
    {
        timer::flatTimer timer("test V2 OMP");
        Field3d densityOnGridLocal(densityOnGrid);
        densityOnGridLocal.setZero();
#pragma omp barrier
#pragma omp for schedule(dynamic, 64)
        for (auto j = 0; j < size(); ++j) {
            alignas(64) double sx[SMAX];
            alignas(64) double sy[SMAX];
            alignas(64) double sz[SMAX];

            for (const Particle& particle : particlesData(j)) {
                // Vectorizable coordinate calculations
                const double ix = particle.coord.x() / domain_.cell_size().x();
                const double iy = particle.coord.y() / domain_.cell_size().y();
                const double iz = particle.coord.z() / domain_.cell_size().z();

                const int xk = floor(ix);
                const int yk = floor(iy);
                const int zk = floor(iz);

// Vectorizable shape calculations
#pragma omp simd
                for (int n = 0; n < SMAX; ++n) {
                    sx[n] = ShapeFn(-ix + double(xk - GHOST_CELLS + n));
                    sy[n] = ShapeFn(-iy + double(yk - GHOST_CELLS + n));
                    sz[n] = ShapeFn(-iz + double(zk - GHOST_CELLS + n));
                }

                const double weight = is_neutral() ? mpw_ : mpw_ * charge;

                // Density accumulation with loop unrolling
                // #pragma unroll
                for (int n = 0; n < SMAX; ++n) {
                    const int indx = xk + n;
                    const double sxw = sx[n];

                    for (int m = 0; m < SMAX; ++m) {
                        const int indy = yk + m;
                        const double sxyw = sxw * sy[m];

                        for (int k = 0; k < SMAX; ++k) {
                            const int indz = zk + k;
                            densityOnGridLocal(indx, indy, indz, 0) += weight * sxyw * sz[k];
                        }
                    }
                }
            }
        }

#pragma omp critical
        densityOnGrid += densityOnGridLocal;
    }

    timerOpt.finish();

    const Field3d diff = ref - densityOnGrid;
    const Field3d diff2 = ref - testV4;

    std::cout << "Ref vs test V2, norm: " << diff.norm() << ", normalized: " << diff.norm() / ref.norm() << std::endl;
    std::cout << "Ref vs test V4, norm: " << diff2.norm() << ", normalized: " << diff2.norm() / ref.norm() << std::endl;
}

void ParticlesArray::density_on_grid_update_impl_ngp() {
    densityOnGrid.setZero();
#pragma omp parallel for schedule(dynamic, 64)
    for (auto j = 0; j < size(); ++j) {
        for (const auto& particle : particlesData(j)) {
            // Vectorizable coordinate calculations
            const double x = particle.coord.x() / domain_.cell_size().x();
            const double y = particle.coord.y() / domain_.cell_size().y();
            const double z = particle.coord.z() / domain_.cell_size().z();
            const int indx = ngp(x);
            const int indy = ngp(y);
            const int indz = ngp(z);

            const double weight = is_neutral() ? mpw_ : mpw_ * charge;

#pragma omp atomic update
            densityOnGrid(indx, indy, indz, 0) += weight;
        }
    }
}

// Public interface that selects appropriate implementation
void ParticlesArray::density_on_grid_update(ShapeType type) {
    switch (type) {
        case ShapeType::NGP:
            density_on_grid_update_impl_ngp();
            break;
        case ShapeType::Linear:
            density_on_grid_update_impl<Shape, 2>();
            break;
        case ShapeType::Quadratic:
            density_on_grid_update_impl<Shape2, 2>();
            break;
    }
}

double ParticlesArray::get_kinetic_energy() const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            energy += get_energy_particle(particle.velocity, mass_, mpw_);
        }
    }

    return energy;
}

double ParticlesArray::get_init_kinetic_energy() const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            energy += get_energy_particle(particle.initVelocity, mass_, mpw_);
        }
    }

    return energy;
}

double ParticlesArray::get_kinetic_energy(int dim) const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            Vector3R velocity(0, 0, 0);
            velocity[dim] = particle.velocity[dim];
            energy += get_energy_particle(velocity, mass_, mpw_);
        }
    }

    return energy;
}

Vector3R ParticlesArray::get_kinetic_energy_component() const {
    double enx = 0;
    double eny = 0;
    double enz = 0;
#pragma omp parallel for reduction(+ : enx) reduction(+ : eny) reduction(+ : enz)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            auto velocity = particle.velocity;
            enx += get_energy_particle(velocity.x(), mass_, mpw_);
            eny += get_energy_particle(velocity.y(), mass_, mpw_);
            enz += get_energy_particle(velocity.z(), mass_, mpw_);
        }
    }

    return Vector3R(enx, eny, enz);
}
double ParticlesArray::get_kinetic_energy(int dim1, int dim2) const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            Vector3R velocity(0, 0, 0);
            velocity[dim1] = particle.velocity[dim1];
            velocity[dim2] = particle.velocity[dim2];
            energy += get_energy_particle(velocity, mass_, mpw_);
        }
    }

    return energy;
}

EnergySpectrum ParticlesArray::calculate_energy_spectrum() const {
    const int num_bins = 1000;
    // Находим минимальную и максимальную энергии
    double min_energy = std::numeric_limits<double>::max();
    double max_energy = std::numeric_limits<double>::lowest();

#pragma omp parallel for reduction(min : min_energy) reduction(max : max_energy)
    for (int k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            double e = get_energy_particle(particle.velocity, mass_, mpw_);
            if (e < min_energy)
                min_energy = e;
            if (e > max_energy)
                max_energy = e;
        }
    }

    // Обработка случая с одинаковыми энергиями
    if (min_energy == max_energy) {
        min_energy -= 0.1;
        max_energy += 0.1;
    }
    double bin_width = (max_energy - min_energy) / num_bins;

    // Счётчики частиц в бинах
    std::vector<int> counts(num_bins, 0);

// Заполнение гистограммы
#pragma omp parallel
    {
        std::vector<int> local_counts(num_bins, 0);

#pragma omp for
        for (int k = 0; k < size(); ++k) {
            for (const auto& particle : particlesData(k)) {
                double e = get_energy_particle(particle.velocity, mass_, mpw_);
                int index = static_cast<int>((e - min_energy) / bin_width);
                index = std::clamp(index, 0, num_bins - 1);
                local_counts[index]++;
            }
        }

#pragma omp critical
        {
            for (int i = 0; i < num_bins; ++i) {
                counts[i] += local_counts[i];
            }
        }
    }
    return EnergySpectrum(counts, min_energy, max_energy);
}
