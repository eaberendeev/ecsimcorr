#include "particle_config.h"

#include <cmath>
#include <stdexcept>

#include "config/json_settings.h"
#include "runtime/profiler.h"

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
        const std::string type = config.value("type", "");

        if (type == "rectangle") {
            const Vector3R center = util::parse_double3(config.at("center"));
            const Vector3R half_length = util::parse_double3(config.at("half_length"));
            return std::make_shared<RectangleDistribution>(center, half_length);
        } else if (type == "cylinder_z") {
            const Vector3R center = util::parse_double3(config.at("center"));
            const double radius = config.value("radius", 0.0);
            const double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderAxisAlignedDistribution>(center, radius, half_length, Axis::Z);
        } else if (type == "cylinder_x") {
            const Vector3R center = util::parse_double3(config.at("center"));
            const double radius = config.value("radius", 0.0);
            const double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderAxisAlignedDistribution>(center, radius, half_length, Axis::X);
        } else if (type == "cylinder_ring_z") {
            const Vector3R center = util::parse_double3(config.at("center"));
            // TODO: check r1 < r2
            const double r1 = config.value("r1", 0.0);
            const double r2 = config.value("r2", 0.0);
            const double half_length = config.value("half_length", 0.0);
            return std::make_shared<CylinderRingZDistribution>(center, r1, r2, half_length);
        } else {
            throw std::runtime_error("Unknown position distribution type: " + type);
        }
    }
};

class VelocityDistributionFactory {
   public:
    static VelDistPtr create(const nlohmann::json& config, double mass = 1.0) {
        const std::string type = config.value("type", "");

        if (type == "gaussian") {
            const Vector3R mean = util::parse_double3(config.at("mean"));
            // in keV
            const Vector3R temperature_kev = util::parse_double3(config.at("sigma"));
            const double sigmax = std::sqrt(temperature_kev.x() / SGS::MC2 / mass);
            const double sigmay = std::sqrt(temperature_kev.y() / SGS::MC2 / mass);
            const double sigmaz = std::sqrt(temperature_kev.z() / SGS::MC2 / mass);
            const Vector3R sigma(sigmax, sigmay, sigmaz);
            return std::make_shared<GaussianVelocity>(mean, sigma);

        } else if (type == "rigid_rotation") {
            const Vector3R rotation_center = util::parse_double3(config.at("rotation_center"));
            const double omega = config.at("omega").get<double>();
            return std::make_shared<RigidRotationVelocity>(rotation_center, omega);

        } else if (type == "tangential") {
            const Vector3R rotation_center = util::parse_double3(config.at("rotation_center"));
            const double mean_speed = config.at("mean_speed").get<double>();
            const double sigma_speed = config.value("sigma_speed", 0.0);

            // A scalar applies to all components; absent thermal_sigma means zero noise.
            const Vector3R thermal_sigma = !config.contains("thermal_sigma") ? Vector3R(0.0)
                                           : config.at("thermal_sigma").is_number()
                                               ? Vector3R(config.at("thermal_sigma").get<double>())
                                               : util::parse_double3(config.at("thermal_sigma"));

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

        const std::string sampling = config.value("sampling", std::string());
        if (sampling == "halton_quiet") {
            if (type != "initial")
                throw std::runtime_error("halton_quiet currently supports only type=initial");
            const auto& space = config.at("dist_space");
            const auto& pulse = config.at("dist_pulse");
            if (space.value("type", std::string()) != "rectangle" ||
                pulse.value("type", std::string()) != "gaussian")
                throw std::runtime_error("halton_quiet requires rectangle dist_space and gaussian dist_pulse");
            const Vector3R center = util::parse_double3(space.at("center"));
            const Vector3R half_length = util::parse_double3(space.at("half_length"));
            const Vector3R mean = util::parse_double3(pulse.at("mean"));
            const Vector3R temperature_kev = util::parse_double3(pulse.at("sigma"));
            const Vector3R sigma(std::sqrt(temperature_kev.x() / SGS::MC2 / mass),
                                 std::sqrt(temperature_kev.y() / SGS::MC2 / mass),
                                 std::sqrt(temperature_kev.z() / SGS::MC2 / mass));
            const int axis = config.value("stratified_axis", 2);
            const double dens = config.value("density", 1.0);
            const double volume = 8.0 * half_length.x() * half_length.y() * half_length.z();
            const double calc = static_cast<double>(num_part_per_cell) * volume * dens / cell_volume;
            const int count = std::max(0, static_cast<int>(std::lround(calc)));
            if (count <= 0 || count % 2 != 0)
                throw std::runtime_error("halton_quiet requires an even positive particle count");
            SinePerturbation perturbation;
            if (config.contains("perturbation")) {
                const auto& p = config.at("perturbation");
                perturbation.amplitude = p.value("amplitude", 0.0);
                perturbation.mode = p.value("mode", 1);
                perturbation.phase = p.value("phase", 0.0);
                perturbation.flow = p.value("flow", 0.0);   // velocity, units of c
            }
            return std::make_unique<QuietStartDistribution>(center, half_length, mean, sigma, count, axis, mass, mpw,
                                                            perturbation);
        }

        if (sampling == "halton") {
            if (type != "initial")
                throw std::runtime_error("Halton sampling currently supports only type=initial");
            const auto& space = config.at("dist_space");
            const auto& pulse = config.at("dist_pulse");
            if (space.value("type", std::string()) != "rectangle" ||
                pulse.value("type", std::string()) != "gaussian")
                throw std::runtime_error("Halton sampling requires rectangle dist_space and gaussian dist_pulse");
            const Vector3R center = util::parse_double3(space.at("center"));
            const Vector3R half_length = util::parse_double3(space.at("half_length"));
            const Vector3R mean = util::parse_double3(pulse.at("mean"));
            const Vector3R temperature_kev = util::parse_double3(pulse.at("sigma"));
            const Vector3R sigma(std::sqrt(temperature_kev.x() / SGS::MC2 / mass),
                                 std::sqrt(temperature_kev.y() / SGS::MC2 / mass),
                                 std::sqrt(temperature_kev.z() / SGS::MC2 / mass));
            const double dens = config.value("density", 1.0);
            const double volume = 8.0 * half_length.x() * half_length.y() * half_length.z();
            const double calc = static_cast<double>(num_part_per_cell) * volume * dens / cell_volume;
            const int count = std::max(0, static_cast<int>(std::lround(calc)));
            return std::make_unique<HaltonDistribution>(center, half_length, mean, sigma, count, mass, mpw);
        }

        // Создаем распределения позиций и скоростей
        const PosDistPtr position = PositionDistributionFactory::create(config["dist_space"]);
        const VelDistPtr velocity = VelocityDistributionFactory::create(config["dist_pulse"], mass);

        if (!config.contains("density")) {
            std::cout << "Warning: density not specified, density will set to 1.0 " << std::endl;
        }
        const double dens = config.value("density", 1.0);
        const double volume = position->get_volume();

        if (type == "initial") {
            // Для начального распределения - целое число частиц
            const double calc = static_cast<double>(num_part_per_cell) * volume * dens / cell_volume;
            const int count = std::max(0, static_cast<int>(std::lround(calc)));

            return std::make_unique<InitialDistribution>(position, velocity, type, count, mass, mpw);
        } else if (type == "injection" || type == "injection_bound") {
            // Для инжекции - вычисляем скорость инжекции
            const double rate = static_cast<double>(num_part_per_cell) * volume * dens / cell_volume;

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

            const std::string type = dist_config.value("type", "");

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
