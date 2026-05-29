#ifndef PARTICLES_DISTRIBUTION_COLLECTION_H_
#define PARTICLES_DISTRIBUTION_COLLECTION_H_
#include <nlohmann/json.hpp>
#include <vector>

#include "Particle.h"
#include "random_generator.h"
#include "sgs.h"
#include "util.h"
#include "vector3.h"

// ===== утилитарный парсер (убираем дублирование) =====
namespace util {
static Vector3R parse_double3(const nlohmann::json& arr) {
    if (!arr.is_array() || arr.size() != 3) {
        throw std::runtime_error("Expected array of 3 numbers");
    }
    return Vector3R(arr[0].get<double>(), arr[1].get<double>(), arr[2].get<double>());
}
}   // namespace util

// ---------------- position types ----------------
struct IPositionDistribution {
    virtual ~IPositionDistribution() = default;
    virtual Vector3R sample(ThreadRandomGenerator& rng) const = 0;
    virtual double get_volume() const = 0;
};

struct IVelocityDistribution {
    virtual ~IVelocityDistribution() = default;
    virtual Vector3R sample(ThreadRandomGenerator& /*rng*/) const = 0;
    virtual Vector3R sample(const Vector3R& /*position*/, ThreadRandomGenerator& /*rng*/) const {
        throw std::runtime_error(
            "IVelocityDistribution::sample(position, rng) is pure virtual.\n"
            "Every concrete velocity distribution MUST override this method.\n"
            "Example: GaussianVelocity, RigidRotationVelocity, TangentialVelocityDistribution");
    }
};

using PosDistPtr = std::shared_ptr<IPositionDistribution>;
using VelDistPtr = std::shared_ptr<IVelocityDistribution>;

// RectangleDistribution - минимальные переименования для ясности
struct RectangleDistribution : public IPositionDistribution {
    Vector3R center;
    Vector3R half_extent;   // раньше "half"
    RectangleDistribution(const Vector3R& c, const Vector3R& half_) : center(c), half_extent(half_) {
    }
    double get_volume() const override {
        return 8.0 * half_extent.x() * half_extent.y() * half_extent.z();
    }
    Vector3R sample(ThreadRandomGenerator& rng) const override {
        double ux = rng.Uniform01(), uy = rng.Uniform01(), uz = rng.Uniform01();
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

    Vector3R sample(ThreadRandomGenerator& rng) const override {
        double c1, c2, l;
        do {
            c1 = (1.0 - 2.0 * rng.Uniform01()) * radius;
            c2 = (1.0 - 2.0 * rng.Uniform01()) * radius;
        } while (c1 * c1 + c2 * c2 > radius * radius);

        l = (1.0 - 2.0 * rng.Uniform01()) * half_length;

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
    Vector3R sample(ThreadRandomGenerator& rng) const override {
        double r = std::sqrt(rng.Uniform01() * (r2 * r2 - r1 * r1) + r1 * r1);
        double phi = 2.0 * M_PI * rng.Uniform01();
        double rx = r * std::cos(phi);
        double ry = r * std::sin(phi);
        double rz = (1.0 - 2.0 * rng.Uniform01()) * half_length_z;
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
    Vector3R sample(ThreadRandomGenerator& rng) const override {
        return Vector3R(mean.x() + rng.Gauss(sigma.x()), mean.y() + rng.Gauss(sigma.y()),
                        mean.z() + rng.Gauss(sigma.z()));
    }
    Vector3R sample(const Vector3R& /*position*/, ThreadRandomGenerator& rng) const override {
        return sample(rng);
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
    Vector3R sample(ThreadRandomGenerator&) const override {
        throw std::runtime_error(
            "TangentialVelocityDistribution requires position. Use sample(pos, "
            "rng).");
    }

    Vector3R sample(const Vector3R& pos, ThreadRandomGenerator& rng) const override {
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

        // Определяем линейную скорость вдоль касательной
        double v_tang = mean_speed;
        if (sigma_speed > 0.0) {
            v_tang = mean_speed + rng.Gauss(sigma_speed);
        }

        // Базовый вектор направленной скорости
        Vector3R v_base = e_phi * v_tang;

        // Добавляем тепловой шум
        v_base.x() += rng.Gauss(thermal_sigma.x());
        v_base.y() += rng.Gauss(thermal_sigma.y());
        v_base.z() += rng.Gauss(thermal_sigma.z());

        return v_base;
    }
};

struct RigidRotationVelocity : public IVelocityDistribution {
    Vector3R rotation_center;   // точка на оси вращения (используются только x,y)
    double omega;               // угловая скорость [рад/с]

    RigidRotationVelocity(const Vector3R& rot_center, double omega_val)
        : rotation_center(rot_center), omega(omega_val) {
    }

    Vector3R sample(ThreadRandomGenerator& /*rng*/) const override {
        throw std::runtime_error("RigidRotationVelocity requires position. Use sample(pos, rng).");
    }

    // Основная версия, зависящая от позиции
    Vector3R sample(const Vector3R& position, ThreadRandomGenerator& /*rng*/) const override {
        double rx = position.x() - rotation_center.x();
        double ry = position.y() - rotation_center.y();
        // v = ω × r  => v_x = -ω * ry, v_y = ω * rx, v_z = 0
        return Vector3R(-omega * ry, omega * rx, 0.0);
    }
};

// ====== фабрики (используем util::parse_double3 и value(...)) ======
class PositionDistributionFactory {
   public:
    static PosDistPtr create(const nlohmann::json& config) {
        std::string type = config.value("type", "");

        if (type == "rectangle") {
            Vector3R center = util::parse_double3(config.at("center"));
            Vector3R half_length = util::parse_double3(config.at("half_length"));
            return std::make_shared<RectangleDistribution>(center, half_length);
        } else if (type == "cylinder_z") {
            Vector3R center = util::parse_double3(config.at("center"));
            double radius = config.value("radius", 0.0);
            double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderAxisAlignedDistribution>(center, radius, half_length, Axis::Z);
        } else if (type == "cylinder_x") {
            Vector3R center = util::parse_double3(config.at("center"));
            double radius = config.value("radius", 0.0);
            double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderAxisAlignedDistribution>(center, radius, half_length, Axis::X);
        } else if (type == "cylinder_ring_z") {
            Vector3R center = util::parse_double3(config.at("center"));
            // TODO: check r1 < r2
            double r1 = config.value("r1", 0.0);
            double r2 = config.value("r2", 0.0);
            double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderRingZDistribution>(center, r1, r2, half_length);
        } else {
            throw std::runtime_error("Unknown position distribution type: " + type);
        }
    }
};

class VelocityDistributionFactory {
   public:
    static VelDistPtr create(const nlohmann::json& config, double mass = 1.0) {
        std::string type = config.value("type", "");

        if (type == "gaussian") {
            Vector3R mean = util::parse_double3(config.at("mean"));
            // in keV
            Vector3R sigma = util::parse_double3(config.at("sigma"));
            double sigmax = std::sqrt(sigma.x() / SGS::MC2 / mass);
            double sigmay = std::sqrt(sigma.y() / SGS::MC2 / mass);
            double sigmaz = std::sqrt(sigma.z() / SGS::MC2 / mass);
            sigma = Vector3R(sigmax, sigmay, sigmaz);
            return std::make_shared<GaussianVelocity>(mean, sigma);

        } else if (type == "rigid_rotation") {
            Vector3R rotation_center = util::parse_double3(config.at("rotation_center"));
            double omega = config.at("omega").get<double>();
            return std::make_shared<RigidRotationVelocity>(rotation_center, omega);

        } else if (type == "tangential") {
            Vector3R rotation_center = util::parse_double3(config.at("rotation_center"));
            double mean_speed = config.at("mean_speed").get<double>();
            double sigma_speed = config.value("sigma_speed", 0.0);

            // Парсинг thermal_sigma: число или тройка чисел
            Vector3R thermal_sigma(0.0, 0.0, 0.0);
            if (config.contains("thermal_sigma")) {
                if (config.at("thermal_sigma").is_number()) {
                    double val = config.at("thermal_sigma").get<double>();
                    thermal_sigma = Vector3R(val, val, val);
                } else {
                    thermal_sigma = util::parse_double3(config.at("thermal_sigma"));
                }
            }   // если отсутствует, остаются нули

            return std::make_shared<TangentialVelocityDistribution>(rotation_center, mean_speed, sigma_speed,
                                                                    thermal_sigma);

        } else {
            throw std::runtime_error("Unknown velocity distribution type: " + type);
        }
    }
};

struct Distribution {
    std::shared_ptr<const IPositionDistribution> position;
    std::shared_ptr<const IVelocityDistribution> velocity;
    std::string type;
    int count = 0;
    int get_count() const {
        return count;
    }
    const std::string& get_type() const {
        return type;
    }
};

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
    Vector3R sample_position(ThreadRandomGenerator& rng) {
        return position_->sample(rng);
    }

    Vector3R sample_velocity(const Vector3R& position, ThreadRandomGenerator& rng) {
        return velocity_->sample(position, rng);
    }

    Vector3R sample_velocity(ThreadRandomGenerator& rng) {
        return velocity_->sample(rng);
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

// Для начального распределения (целое число частиц)
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
class DistributionFactory {
   public:
    static std::unique_ptr<IDistribution> create(const nlohmann::json& config, const std::string& type,
                                                 double cell_volume, int num_part_per_cell, double mass, double mpw) {
        if (!config.contains("dist_space") || !config.contains("dist_pulse")) {
            throw std::runtime_error("Distribution config missing 'dist_space' or 'dist_pulse'");
        }

        // Создаем распределения позиций и скоростей
        auto position = PositionDistributionFactory::create(config["dist_space"]);
        auto velocity = VelocityDistributionFactory::create(config["dist_pulse"], mass);

        if (!config.contains("density")) {
            std::cout << "Warning: density not specified, density will set to 1.0 " << std::endl;
        }
        double dens = config.value("density", 1.0);
        double volume = position->get_volume();

        if (type == "initial") {
            // Для начального распределения - целое число частиц
            double calc = static_cast<double>(num_part_per_cell) * volume * dens / cell_volume;
            int count = std::max(0, static_cast<int>(std::lround(calc)));

            return std::make_unique<InitialDistribution>(position, velocity, type, count, mass, mpw);
        } else if (type == "injection" || type == "injection_bound") {
            // Для инжекции - вычисляем скорость инжекции
            double rate = static_cast<double>(num_part_per_cell) * volume * dens / cell_volume;

            return std::make_unique<InjectionDistribution>(position, velocity, type, rate, mass, mpw);
        } else {
            throw std::runtime_error("Unknown distribution type: " + type);
        }
    }

    // Альтернативный метод для создания из общей конфигурации
    static std::vector<std::unique_ptr<IDistribution>> createFromConfig(const nlohmann::json& config,
                                                                        double cell_volume, int num_part_per_cell,
                                                                        double mass, double mpw) {
        std::vector<std::unique_ptr<IDistribution>> distributions;

        if (!config.contains("distribution") || !config["distribution"].is_array()) {
            std::cout << "Expected array of distributions" << std::endl;
            return distributions;
        }

        for (const auto& dist_config : config["distribution"]) {
            if (!dist_config.contains("type")) {
                throw std::runtime_error("Distribution type not specified");
            }

            std::string type = dist_config.value("type", "");

            if (type == "initial" || type == "injection" || type == "injection_bound") {
                auto dist = create(dist_config, type, cell_volume, num_part_per_cell, mass, mpw);
                distributions.push_back(std::move(dist));
            } else if (type != "none") {
                throw std::runtime_error("Unknown distribution type: " + type);
            }
        }

        return distributions;
    }
};

#endif
