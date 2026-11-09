#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "grid/domain.h"
#include "particles/Particle.h"
#include "runtime/random.h"

// ---------------- position types ----------------
struct IPositionDistribution {
    virtual ~IPositionDistribution() = default;
    virtual Vector3R sample(LehmerEngine& rndEng) const = 0;
    virtual double get_volume() const = 0;
};

struct IVelocityDistribution {
    virtual ~IVelocityDistribution() = default;
    virtual Vector3R sample(LehmerEngine& /*rndEng*/) const = 0;

    virtual Vector3R sample(const Vector3R& /*position*/, LehmerEngine& /*rndEng*/) const {
        assert(false && "Every velocity distribution must override sample(position, rng)");
        return Vector3R(0.0);
    }
};

using PosDistPtr = std::shared_ptr<IPositionDistribution>;
using VelDistPtr = std::shared_ptr<IVelocityDistribution>;

struct RectangleDistribution : public IPositionDistribution {
    Vector3R center;
    Vector3R half_extent;

    RectangleDistribution(const Vector3R& c, const Vector3R& half_) : center(c), half_extent(half_) {
    }

    double get_volume() const override {
        return 8.0 * half_extent.x() * half_extent.y() * half_extent.z();
    }

    Vector3R sample(LehmerEngine& rndEng) const override {
        std::uniform_real_distribution distr(0.0, 1.0);
        const double ux = distr(rndEng);
        const double uy = distr(rndEng);
        const double uz = distr(rndEng);
        return Vector3R(center.x() + half_extent.x() * (1.0 - 2.0 * ux),
                        center.y() + half_extent.y() * (1.0 - 2.0 * uy),
                        center.z() + half_extent.z() * (1.0 - 2.0 * uz));
    }
};

struct CylinderAxisAlignedDistribution : public IPositionDistribution {
    Vector3R center;
    double radius;
    double half_length;
    Axis ax;

    CylinderAxisAlignedDistribution(const Vector3R& c, double r, double halfL, Axis ax)
        : center(c), radius(r), half_length(halfL), ax(ax) {
    }

    double get_volume() const override {
        double r = std::max(radius, 0.0);
        return M_PI * r * r * (2.0 * half_length);
    }

    Vector3R sample(LehmerEngine& rndEng) const override {
        std::uniform_real_distribution distr(0.0, 1.0);

        double c1, c2, l;
        do {
            c1 = (1.0 - 2.0 * distr(rndEng)) * radius;
            c2 = (1.0 - 2.0 * distr(rndEng)) * radius;
        } while (c1 * c1 + c2 * c2 > radius * radius);

        l = (1.0 - 2.0 * distr(rndEng)) * half_length;

        // Смещения по главной оси (длина цилиндра)
        double x_off = (ax == Axis::X ? l : 0.0);
        double y_off = (ax == Axis::Y ? l : 0.0);
        double z_off = (ax == Axis::Z ? l : 0.0);

        // Смещения по двум осям сечения
        double x_cross = 0, y_cross = 0, z_cross = 0;
        if (ax == Axis::X) {
            x_cross = 0.0;
            y_cross = c1;
            z_cross = c2;
        } else if (ax == Axis::Y) {
            x_cross = c2;
            y_cross = 0.0;
            z_cross = c1;
        } else {
            x_cross = c1;
            y_cross = c2;
            z_cross = 0.0;
        }

        return Vector3R(center.x() + x_off + x_cross, center.y() + y_off + y_cross, center.z() + z_off + z_cross);
    }
};

struct CylinderRingZDistribution : public IPositionDistribution {
    Vector3R center;
    double r1;
    double r2;
    double half_length_z;

    CylinderRingZDistribution(const Vector3R& c, double r1_, double r2_, double halfZ_)
        : center(c), r1(r1_), r2(r2_), half_length_z(halfZ_) {
    }

    double get_volume() const override {
        if (r2 > r1 && r2 > 0.0) {
            return M_PI * (r2 * r2 - r1 * r1) * (2.0 * half_length_z);
        } else {
            return 0.0;
        }
    }

    Vector3R sample(LehmerEngine& rndEng) const override {
        std::uniform_real_distribution distr(0.0, 1.0);

        double r = std::sqrt(distr(rndEng) * (r2 * r2 - r1 * r1) + r1 * r1);
        double phi = 2.0 * M_PI * distr(rndEng);
        double rx = r * std::cos(phi);
        double ry = r * std::sin(phi);
        double rz = (1.0 - 2.0 * distr(rndEng)) * half_length_z;
        return Vector3R(center.x() + rx, center.y() + ry, center.z() + rz);
    }
};

// ---------------- velocity ----------------
struct GaussianVelocity : public IVelocityDistribution {
    Vector3R mean;
    Vector3R sigma;

    GaussianVelocity() {
        mean = Vector3R(0, 0, 0);
        sigma = Vector3R(0, 0, 0);
    }

    void set(const Vector3R& m, const Vector3R& s) {
        mean = m;
        sigma = s;
    }

    GaussianVelocity(const Vector3R& m, const Vector3R& s) : mean(m), sigma(s) {
    }

    Vector3R sample(LehmerEngine& rndEng) const override {
        std::normal_distribution normalDistr;
        const Vector3R res = {
            normalDistr(rndEng, std::normal_distribution<double>::param_type(mean.x(), sigma.x())),
            normalDistr(rndEng, std::normal_distribution<double>::param_type(mean.y(), sigma.y())),
            normalDistr(rndEng, std::normal_distribution<double>::param_type(mean.z(), sigma.z())),
        };
        return res;
    }

    Vector3R sample(const Vector3R& /*position*/, LehmerEngine& rndEng) const override {
        return sample(rndEng);
    }
};

struct TangentialVelocityDistribution : public IVelocityDistribution {
    Vector3R center;          // точка на оси вращения (x0, y0, z0 – z не используется)
    double mean_speed;        // средняя линейная скорость по касательной
    double sigma_speed;       // если >0 – скорость выбирается из N(mean_speed,
                              // sigma_speed)
    Vector3R thermal_sigma;   // стандартные отклонения теплового шума (σx, σy, σz)

    TangentialVelocityDistribution(const Vector3R& rot_center, double v_mean, double v_sigma,
                                   const Vector3R& therm_sigma)
        : center(rot_center), mean_speed(v_mean), sigma_speed(v_sigma), thermal_sigma(therm_sigma) {
    }

    // Без позиции не имеет смысла
    Vector3R sample(LehmerEngine&) const override {
        assert(false && "Tangential velocity requires sample(position, rng)");
        return Vector3R(0.0);
    }

    Vector3R sample(const Vector3R& pos, LehmerEngine& rndEng) const override {
        double rx = pos.x() - center.x();
        double ry = pos.y() - center.y();
        double r = std::sqrt(rx * rx + ry * ry);

        // Единичный касательный вектор (для r=0 зададим произвольное
        // направление, чтобы не падать)
        Vector3R e_phi(0.0, 0.0, 0.0);
        if (r > 0.0) {
            e_phi = Vector3R(-ry / r, rx / r, 0.0);
        } else {
            // если частица попала точно на ось – направление не определено,
            // ставим ноль (или можно кинуть исключение)
        }

        std::normal_distribution normalDistr;

        // Определяем линейную скорость вдоль касательной
        double v_tang = mean_speed;
        if (sigma_speed > 0.0) {
            v_tang = normalDistr(rndEng, std::normal_distribution<double>::param_type(mean_speed, sigma_speed));
        }

        // Базовый вектор направленной скорости
        Vector3R v_base = e_phi * v_tang;

        // Добавляем тепловой шум
        v_base.x() += normalDistr(rndEng, std::normal_distribution<double>::param_type(0.0, thermal_sigma.x()));
        v_base.y() += normalDistr(rndEng, std::normal_distribution<double>::param_type(0.0, thermal_sigma.y()));
        v_base.z() += normalDistr(rndEng, std::normal_distribution<double>::param_type(0.0, thermal_sigma.z()));

        return v_base;
    }
};

struct RigidRotationVelocity : public IVelocityDistribution {
    Vector3R rotation_center;   // точка на оси вращения (используются только x,y)
    double omega;               // угловая скорость [рад/с]

    RigidRotationVelocity(const Vector3R& rot_center, double omega_val)
        : rotation_center(rot_center), omega(omega_val) {
    }

    Vector3R sample(LehmerEngine& /*rndEng*/) const override {
        assert(false && "Rigid rotation requires sample(position, rng)");
        return Vector3R(0.0);
    }

    // Основная версия, зависящая от позиции
    Vector3R sample(const Vector3R& position, LehmerEngine& /*rndEng*/) const override {
        double rx = position.x() - rotation_center.x();
        double ry = position.y() - rotation_center.y();
        // v = ω × r  => v_x = -ω * ry, v_y = ω * rx, v_z = 0
        return Vector3R(-omega * ry, omega * rx, 0.0);
    }
};

inline double van_der_corput(std::size_t index, std::size_t base) {
    double q = 0.0;
    double bk = 1.0 / static_cast<double>(base);
    while (index > 0) {
        q += static_cast<double>(index % base) * bk;
        index /= base;
        bk /= static_cast<double>(base);
    }
    return q;
}

class IDistribution {
   public:
    IDistribution(std::shared_ptr<const IPositionDistribution> position,
                  std::shared_ptr<const IVelocityDistribution> velocity, const std::string& type, double mass,
                  double mpw)
        : position_(position), velocity_(velocity), type_(type), mass_(mass), mpw_(mpw) {
    }

    virtual ~IDistribution() = default;

    virtual int get_count_to_inject() = 0;
    std::shared_ptr<const IPositionDistribution> position_;
    std::shared_ptr<const IVelocityDistribution> velocity_;
    std::string type_;
    double mass_;
    double mpw_;

    Vector3R sample_position(LehmerEngine& rndEng) {
        return position_->sample(rndEng);
    }

    Vector3R sample_velocity(const Vector3R& position, LehmerEngine& rndEng) {
        return velocity_->sample(position, rndEng);
    }

    Vector3R sample_velocity(LehmerEngine& rndEng) {
        return velocity_->sample(rndEng);
    }

    // Optional index-aware sampling for low-discrepancy sequences. The
    // default keeps the existing RNG path unchanged; HaltonDistribution
    // overrides it and uses the particle index directly.
    virtual void sample_particle(std::size_t index, LehmerEngine& space, LehmerEngine& momentum, Vector3R& position,
                                 Vector3R& velocity) {
        (void) index;
        position = sample_position(space);
        velocity = sample_velocity(position, momentum);
    }

    const std::string& get_type() const {
        return type_;
    }

    double get_energy(const Vector3R& velocity) const {
        return get_energy_particle(velocity, mass_, mpw_);
    }

    bool is_bound_injection() const {
        return type_ == "injection_bound";
    }
};

// Deterministic low-discrepancy initial loader. Coordinate and velocity use
// disjoint prime bases so position-velocity correlations do not align.
class HaltonDistribution : public IDistribution {
   public:
    HaltonDistribution(const Vector3R& center, const Vector3R& half_extent, const Vector3R& mean,
                       const Vector3R& sigma, int count, double mass, double mpw)
        : IDistribution(nullptr, nullptr, "initial", mass, mpw), center_(center), half_extent_(half_extent),
          mean_(mean), sigma_(sigma), count_(count) {
    }

    int get_count_to_inject() override {
        return count_;
    }

    void sample_particle(std::size_t index, LehmerEngine&, LehmerEngine&, Vector3R& position,
                         Vector3R& velocity) override {
        static constexpr std::size_t coordinate_bases[3] = {2, 3, 5};
        static constexpr std::size_t radius_bases[3] = {7, 11, 13};
        static constexpr std::size_t phase_bases[3] = {17, 19, 23};
        const std::size_t n = index + 1;
        for (int d = 0; d < 3; ++d) {
            const double ux = van_der_corput(n, coordinate_bases[d]);
            position[d] = center_[d] + half_extent_[d] * (2.0 * ux - 1.0);
            const double ur = std::max(van_der_corput(n, radius_bases[d]), 1e-15);
            const double up = van_der_corput(n, phase_bases[d]);
            velocity[d] = mean_[d] + sigma_[d] * std::sqrt(-2.0 * std::log(ur)) * std::sin(2.0 * M_PI * up);
        }
    }

   private:
    Vector3R center_, half_extent_, mean_, sigma_;
    int count_ = 0;
};


// Для начального распределения (целое число частиц)

// Quiet-start loader: an even number of markers is split into antithetic
// pairs. Coordinates are Halton in the transverse plane and stratified along
// one axis; the two members of a pair share a coordinate and have opposite
// low-discrepancy Maxwellian velocities. This is the unperturbed xpic IAW
// quiet-start recipe.
class QuietStartDistribution : public IDistribution {
   public:
    QuietStartDistribution(const Vector3R& center, const Vector3R& half_extent, const Vector3R& mean,
                           const Vector3R& sigma, int count, int stratified_axis, double mass, double mpw)
        : IDistribution(nullptr, nullptr, "initial", mass, mpw), center_(center), half_extent_(half_extent),
          mean_(mean), sigma_(sigma), count_(count), stratified_axis_(stratified_axis) {
        if (count_ <= 0 || count_ % 2 != 0)
            throw std::invalid_argument("halton_quiet requires a positive even particle count");
        if (stratified_axis_ < 0 || stratified_axis_ > 2)
            throw std::invalid_argument("halton_quiet stratified axis must be x, y or z");
    }

    int get_count_to_inject() override {
        return count_;
    }

    void sample_particle(std::size_t index, LehmerEngine&, LehmerEngine&, Vector3R& position,
                         Vector3R& velocity) override {
        static constexpr std::size_t coordinate_bases[3] = {2, 3, 5};
        static constexpr std::size_t radius_bases[3] = {7, 11, 13};
        static constexpr std::size_t phase_bases[3] = {17, 19, 23};
        const std::size_t pair = index / 2;
        const double sign = (index % 2 == 0) ? 1.0 : -1.0;
        const std::size_t pairs = static_cast<std::size_t>(count_ / 2);
        for (int d = 0; d < 3; ++d) {
            const double u = (d == stratified_axis_)
                                 ? (static_cast<double>(pair) + 0.5) / static_cast<double>(pairs)
                                 : van_der_corput(pair + 1, coordinate_bases[d]);
            position[d] = center_[d] + half_extent_[d] * (2.0 * u - 1.0);
            const double radius = std::max(van_der_corput(pair + 1, radius_bases[d]), 1e-15);
            const double phase = van_der_corput(pair + 1, phase_bases[d]);
            velocity[d] = mean_[d] + sign * sigma_[d] * std::sqrt(-2.0 * std::log(radius)) *
                                         std::sin(2.0 * M_PI * phase);
        }
    }

   private:
    Vector3R center_, half_extent_, mean_, sigma_;
    int count_ = 0;
    int stratified_axis_ = 2;
};

class InitialDistribution : public IDistribution {
   private:
    int count_;

   public:
    InitialDistribution(std::shared_ptr<const IPositionDistribution> position,
                        std::shared_ptr<const IVelocityDistribution> velocity, const std::string& type, int count,
                        double mass, double mpw)
        : IDistribution(position, velocity, type, mass, mpw),   // Базовый класс
          count_(count) {
    }

    int get_count_to_inject() override {
        return count_;   // Всегда возвращаем полное количество
    }
};

// Для инжекционного распределения (дробная скорость)
class InjectionDistribution : public IDistribution {
   private:
    double rate_;          // частиц/шаг
    double accumulator_;   // накопленная дробная часть

   public:
    InjectionDistribution(std::shared_ptr<const IPositionDistribution> position,
                          std::shared_ptr<const IVelocityDistribution> velocity, const std::string& type, double rate,
                          double mass, double mpw)
        : IDistribution(position, velocity, type, mass, mpw), rate_(rate), accumulator_(0.0) {
    }

    int get_count_to_inject() override {
        accumulator_ += rate_;
        int to_inject = static_cast<int>(accumulator_);
        accumulator_ -= to_inject;
        return to_inject;
    }

    void reset_accumulator() {
        accumulator_ = 0.0;
    }
};

// Simulation owns the source processes and their fractional accumulators.
// Sources borrow a species only for one generation call; they never own it.
// Step-local Lehmer streams are reconstructed from species name and timestep;
// fractional accumulators retain the existing restart behavior (reset on input).
struct ParticleSources {
    std::vector<std::unique_ptr<IDistribution>> initial;
    std::vector<std::unique_ptr<IDistribution>> injection;
};

class ParticlesArray;

// Only particles accepted by the domain contribute to count and kinetic energy.
struct InjectedParticles {
    double energy = 0;
    int count = 0;
};

InjectedParticles add_particles_from_distribution(ParticlesArray& species, IDistribution& distribution,
                                                  LehmerEngine& space, LehmerEngine& momentum, const Domain& domain,
                                                  double dt, bool check_boundaries);
double distribute_initial_particles(ParticlesArray& species,
                                    const std::vector<std::unique_ptr<IDistribution>>& distributions,
                                    const Domain& domain, bool co_locate_species = true);
double inject_particles_step(ParticlesArray& species, std::vector<std::unique_ptr<IDistribution>>& distributions,
                             int timestep, const Domain& domain, double dt, bool co_locate_species = true);
