#ifndef PARTICLES_DISTRIBUTION_COLLECTION_H_
#define PARTICLES_DISTRIBUTION_COLLECTION_H_
#include <nlohmann/json.hpp>
#include <vector>

#include "Particle.h"
#include "Vec.h"
#include "random_generator.h"
#include "service.h"
#include "util.h"
#include "sgs.h"

// ===== утилитарный парсер (убираем дублирование) =====
namespace util {
static double3 parse_double3(const nlohmann::json& arr) {
    if (!arr.is_array() || arr.size() != 3) {
        throw std::runtime_error("Expected array of 3 numbers");
    }
    return double3(arr[0].get<double>(), arr[1].get<double>(),
                   arr[2].get<double>());
}
}   // namespace util

// ---------------- position types ----------------
struct IPositionDistribution {
    virtual ~IPositionDistribution() = default;
    virtual double3 sample(
        ThreadRandomGenerator& rng) const = 0;   // изменяет RNG -> не const
    virtual double get_volume() const = 0;
};

struct IVelocityDistribution {
    virtual ~IVelocityDistribution() = default;
    virtual double3 sample(ThreadRandomGenerator& rng) const = 0;
};

using PosDistPtr = std::shared_ptr<IPositionDistribution>;
using VelDistPtr = std::shared_ptr<IVelocityDistribution>;

// RectangleDistribution - минимальные переименования для ясности
struct RectangleDistribution : public IPositionDistribution {
    double3 center;
    double3 half_extent;   // раньше "half"
    RectangleDistribution(const double3& c, const double3& half_)
        : center(c), half_extent(half_) {}
    double get_volume() const override {
        return 8.0 * half_extent.x() * half_extent.y() * half_extent.z();
    }
    double3 sample(ThreadRandomGenerator& rng) const override {
        double ux = rng.Uniform01(), uy = rng.Uniform01(), uz = rng.Uniform01();
        return double3(center.x() + half_extent.x() * (1.0 - 2.0 * ux),
                       center.y() + half_extent.y() * (1.0 - 2.0 * uy),
                       center.z() + half_extent.z() * (1.0 - 2.0 * uz));
    }
};

struct CylinderZDistribution : public IPositionDistribution {
    double3 center;
    double radius;
    double half_length_z;
    CylinderZDistribution(const double3& c, double radius_, double halfZ_)
        : center(c), radius(radius_), half_length_z(halfZ_) {}
    double get_volume() const override {
        // height = 2 * half_length_z
        return M_PI * radius * radius * (2.0 * half_length_z);
    }
    double3 sample(ThreadRandomGenerator& rng) const override {
        double rx, ry;
        do {
            rx = (1.0 - 2.0 * rng.Uniform01()) * radius;
            ry = (1.0 - 2.0 * rng.Uniform01()) * radius;
        } while (rx * rx + ry * ry > radius * radius);
        double rz = (1.0 - 2.0 * rng.Uniform01()) * half_length_z;
        return double3(center.x() + rx, center.y() + ry, center.z() + rz);
    }
};

struct CylinderXDistribution : public IPositionDistribution {
    double3 center;
    double radius;
    double half_length_x;
    CylinderXDistribution(const double3& c, double radius_, double halfX_)
        : center(c), radius(radius_), half_length_x(halfX_) {}
    double get_volume() const override {
        if (radius > 0.0) {
            return M_PI * radius * radius * (2.0 * half_length_x);
        } else {
            return 0.0;
        }
    }
    double3 sample(ThreadRandomGenerator& rng) const override {
        double rz, ry;
        do {
            rz = (1.0 - 2.0 * rng.Uniform01()) * radius;
            ry = (1.0 - 2.0 * rng.Uniform01()) * radius;
        } while (rz * rz + ry * ry > radius * radius);
        double rx = (1.0 - 2.0 * rng.Uniform01()) * half_length_x;
        return double3(center.x() + rx, center.y() + ry, center.z() + rz);
    }
};

// ---------------- velocity ----------------
struct GaussianVelocity : public IVelocityDistribution {
    double3 mean;
    double3 sigma;
    // исправлено: принимаем sigma как const&
    GaussianVelocity(const double3& m, const double3& s) : mean(m), sigma(s) {}
    double3 sample(ThreadRandomGenerator& rng) const override {
        return double3(mean.x() + rng.Gauss(sigma.x()),
                       mean.y() + rng.Gauss(sigma.y()),
                       mean.z() + rng.Gauss(sigma.z()));
    }
};

// ====== фабрики (используем util::parse_double3 и value(...)) ======
class PositionDistributionFactory {
   public:
    static PosDistPtr create(const nlohmann::json& config) {
        std::string type = config.value("type", "");

        if (type == "rectangle") {
            double3 center = util::parse_double3(config.at("center"));
            double3 half_length = util::parse_double3(config.at("half_length"));
            return std::make_shared<RectangleDistribution>(center, half_length);
        } else if (type == "cylinder_z") {
            double3 center = util::parse_double3(config.at("center"));
            double radius = config.value("radius", 0.0);
            double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderZDistribution>(center, radius,
                                                           half_length);
        } else if (type == "cylinder_x") {
            double3 center = util::parse_double3(config.at("center"));
            double radius = config.value("radius", 0.0);
            double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderXDistribution>(center, radius,
                                                           half_length);
        } else {
            throw std::runtime_error("Unknown position distribution type: " +
                                     type);
        }
    }
};

class VelocityDistributionFactory {
   public:
    static VelDistPtr create(const nlohmann::json& config, double mass = 1.0) {
        std::string type = config.value("type", "");

        if (type == "gaussian") {
            double3 mean = util::parse_double3(config.at("mean"));
            // in keV
            double3 sigma = util::parse_double3(config.at("sigma"));
            double sigmax = std::sqrt(sigma.x() / SGS::MC2 / std::sqrt(mass));
            double sigmay = std::sqrt(sigma.y() / SGS::MC2 / std::sqrt(mass));
            double sigmaz = std::sqrt(sigma.z() / SGS::MC2 / std::sqrt(mass));
            sigma =
                double3(sigmax, sigmay, sigmaz);
            return std::make_shared<GaussianVelocity>(mean, sigma);
        } else {
            throw std::runtime_error("Unknown velocity distribution type: " +
                                     type);
        }
    }
};

// ===== структура Distribution: мелкие оптимизации =====
struct Distribution {
    std::shared_ptr<const IPositionDistribution> position;
    std::shared_ptr<const IVelocityDistribution> velocity;
    std::string type;
    int count = 0;
    int get_count() const { return count; }
    // возвращаем const& чтобы не копировать строку
    const std::string& get_type() const { return type; }
};

#endif
