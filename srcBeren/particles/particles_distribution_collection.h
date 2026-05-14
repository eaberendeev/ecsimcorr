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
    virtual Vector3R sample(ThreadRandomGenerator& rng) const = 0;   // изменяет RNG -> не const
    virtual double get_volume() const = 0;
};

struct IVelocityDistribution {
    virtual ~IVelocityDistribution() = default;
    virtual Vector3R sample(ThreadRandomGenerator& rng) const = 0;
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

struct CylinderZDistribution : public IPositionDistribution {
    Vector3R center;
    double radius;
    double half_length_z;
    CylinderZDistribution(const Vector3R& c, double radius_, double halfZ_)
        : center(c), radius(radius_), half_length_z(halfZ_) {
    }
    double get_volume() const override {
        // height = 2 * half_length_z
        return M_PI * radius * radius * (2.0 * half_length_z);
    }
    Vector3R sample(ThreadRandomGenerator& rng) const override {
        double rx, ry;
        do {
            rx = (1.0 - 2.0 * rng.Uniform01()) * radius;
            ry = (1.0 - 2.0 * rng.Uniform01()) * radius;
        } while (rx * rx + ry * ry > radius * radius);
        double rz = (1.0 - 2.0 * rng.Uniform01()) * half_length_z;
        return Vector3R(center.x() + rx, center.y() + ry, center.z() + rz);
    }
};

struct CylinderXDistribution : public IPositionDistribution {
    Vector3R center;
    double radius;
    double half_length_x;
    CylinderXDistribution(const Vector3R& c, double radius_, double halfX_)
        : center(c), radius(radius_), half_length_x(halfX_) {
    }
    double get_volume() const override {
        if (radius > 0.0) {
            return M_PI * radius * radius * (2.0 * half_length_x);
        } else {
            return 0.0;
        }
    }
    Vector3R sample(ThreadRandomGenerator& rng) const override {
        double rz, ry;
        do {
            rz = (1.0 - 2.0 * rng.Uniform01()) * radius;
            ry = (1.0 - 2.0 * rng.Uniform01()) * radius;
        } while (rz * rz + ry * ry > radius * radius);
        double rx = (1.0 - 2.0 * rng.Uniform01()) * half_length_x;
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
            return std::make_shared<CylinderZDistribution>(center, radius, half_length);
        } else if (type == "cylinder_x") {
            Vector3R center = util::parse_double3(config.at("center"));
            double radius = config.value("radius", 0.0);
            double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderXDistribution>(center, radius, half_length);
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
