#ifndef COIL_H_
#define COIL_H_
#include "World.h"
#include "nlohmann/json.hpp"
using json = nlohmann::json;

struct Coil {
    double z0, R, I;
    Coil(double z0, double R, double I) : z0(z0), R(R), I(I) {}
};

struct CoilsArray {
    std::vector<Coil> coils;
    static const int N = 2000;
    const double hp = 2 * M_PI / N;
    alignas(64) double cs[N];
    CoilsArray(const json &config) {
        for (const auto &coil_json : config["Coils"]) {
            if (!coil_json.contains("z") || !coil_json.contains("R") ||
                !coil_json.contains("I")) {
                std::cerr << "Error: Invalid coil configuration. Missing z, R, "
                             "or I parameters.\n";
                return;
            }

            double z0 = coil_json["z"].get<double>();
            double R = coil_json["R"].get<double>();
            double I = coil_json["I"].get<double>();
            coils.push_back(Coil(z0, R, I));
        }

#pragma omp simd
        for (auto i = 0; i < N; i++) {
            cs[i] = cos(i * hp);
        }
    }
    ~CoilsArray() {}

    double get_Bz(double z, double r);
    double get_Br(double z, double r);
    double get_integ_z(double z, double r, double R);
    double get_integ_r(double z, double r, double R);
};
void set_coils(Field3d &fieldB, const Domain &domain, const json &config);

// // Базовый класс для конфигурации внешнего поля
// struct ExternalFieldConfig {
//     virtual ~ExternalFieldConfig() = default;
//     virtual void apply(SomeType &field, SomeType &domain) const = 0;
// };

// // Конфигурация для равномерно заряженного цилиндра
// struct UniformlyChargedCylinderConfig : ExternalFieldConfig {
//     double radius;
//     double value;

//     UniformlyChargedCylinderConfig(double r, double v) : radius(r), value(v)
//     {}

//     void apply(SomeType &field, SomeType &domain) const override {
//         set_uniformly_charged_cylinder(field, domain, radius, value);
//     }
// };

// // Конфигурация для однородного поля
// struct UniformFieldConfig : ExternalFieldConfig {
//     double value;
//     std::vector<double> direction;

//     UniformFieldConfig(double v, const std::vector<double> &dir = {0, 0, 1})
//         : value(v), direction(dir) {}

//     void apply(SomeType &field, SomeType &domain) const override {
//         set_uniform_field(field, domain, value, direction);
//     }
// };

// // Функция для создания конфигурации из JSON
// std::unique_ptr<ExternalFieldConfig> create_external_field_config(
//     const json &j) {
//     if (!j.contains("ExternalFieldE")) {
//         return nullptr;
//     }

//     auto &external_field = j["ExternalFieldE"];

//     if (external_field.contains("uniformly_charged_cylinder")) {
//         auto &config = external_field["uniformly_charged_cylinder"];

//         try {
//             double radius = config.at("radius").get<double>();
//             double value = config.at("value").get<double>();

//             return std::make_unique<UniformlyChargedCylinderConfig>(radius,
//                                                                     value);
//         } catch (const json::exception &e) {
//             throw std::runtime_error(
//                 "Invalid uniformly_charged_cylinder config: " +
//                 std::string(e.what()));
//         }
//     } else if (external_field.contains("uniform_field")) {
//         auto &config = external_field["uniform_field"];

//         try {
//             double value = config.at("value").get<double>();
//             std::vector<double> direction = {0, 0, 1};

//             if (config.contains("direction")) {
//                 direction =
//                 config.at("direction").get<std::vector<double>>();
//             }

//             return std::make_unique<UniformFieldConfig>(value, direction);
//         } catch (const json::exception &e) {
//             throw std::runtime_error("Invalid uniform_field config: " +
//                                      std::string(e.what()));
//         }
//     }

//     return nullptr;
// }

// // Использование
// try {
//     auto field_config = create_external_field_config(parameters);

//     if (field_config) {
//         field_config->apply(fieldE_external, domain);
//     } else {
//         std::cout << "No external field configuration found" << std::endl;
//     }
// } catch (const std::exception &e) {
//     std::cerr << "Error creating external field: " << e.what() << std::endl;
//     // Обработка ошибки
// }

#endif
