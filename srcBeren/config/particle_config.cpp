#include "particle_config.h"

#include <cmath>
#include <stdexcept>

#include "config/json_settings.h"
#include "timer.h"

SpeciesParameters species_parameters_from_config(const nlohmann::json& config) {
    RECORD_TIMER;
    SpeciesParameters parameters{config.value("Name", ""),
                                 config.value("Charge", 1.0),
                                 config.value("Mass", 1.0),
                                 config.value("Density", 1.0),
                                 config.at("NumPartPerCell").get<int>(),
                                 0};
    if (parameters.particles_per_cell <= 0 || !std::isfinite(parameters.mass) || parameters.mass <= 0 ||
        !std::isfinite(parameters.density) || parameters.density < 0 || !std::isfinite(parameters.charge))
        throw std::invalid_argument(
            "species requires positive mass/NumPartPerCell and finite charge/nonnegative density");
    parameters.macro_weight = parameters.density / parameters.particles_per_cell;
    return parameters;
}

namespace {
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

class DistributionFactory {
   public:
    static std::unique_ptr<IDistribution> create(const nlohmann::json& config, const std::string& type,
                                                 double cell_volume, int num_part_per_cell, double mass, double mpw) {
        if (!config.contains("dist_space") || !config.contains("dist_pulse")) {
            throw std::runtime_error("Distribution config missing 'dist_space' or 'dist_pulse'");
        }

        // Создаем распределения позиций и скоростей
        PosDistPtr position = PositionDistributionFactory::create(config["dist_space"]);
        VelDistPtr velocity = VelocityDistributionFactory::create(config["dist_pulse"], mass);

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

        for (const nlohmann::json& dist_config : config["distribution"]) {
            if (!dist_config.contains("type")) {
                throw std::runtime_error("Distribution type not specified");
            }

            std::string type = dist_config.value("type", "");

            if (type == "initial" || type == "injection" || type == "injection_bound") {
                std::unique_ptr<IDistribution> dist =
                    create(dist_config, type, cell_volume, num_part_per_cell, mass, mpw);
                distributions.push_back(std::move(dist));
            } else if (type != "none") {
                throw std::runtime_error("Unknown distribution type: " + type);
            }
        }

        return distributions;
    }
};

}   // namespace

ParticleSources particle_sources_from_config(const nlohmann::json& config, const SpeciesParameters& parameters,
                                             const Domain& domain) {
    RECORD_TIMER;
    const double cell_volume = domain.cell_size().elements_product();
    std::vector<std::unique_ptr<IDistribution>> distributions = DistributionFactory::createFromConfig(
        config, cell_volume, parameters.particles_per_cell, parameters.mass, parameters.macro_weight);
    ParticleSources sources;
    for (std::unique_ptr<IDistribution>& distribution : distributions) {
        const std::string type = distribution->get_type();
        if (type == "initial")
            sources.initial.push_back(std::move(distribution));
        else if (type == "injection" || type == "injection_bound")
            sources.injection.push_back(std::move(distribution));
    }
    return sources;
}
