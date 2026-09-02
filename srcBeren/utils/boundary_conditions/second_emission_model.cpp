#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "timer.h"

constexpr double kPi = 3.14159265358979323846;

SecondaryEmissionModel::SecondaryEmissionModel(Face face, std::string product_species,
                                               std::vector<EmissionSourceRule> rules)
    : face_(face), product_(std::move(product_species)), rules_(std::move(rules)), eng_(17) {
    // CYLINDER поддерживается: параметры цилиндра (центр и радиус) доступны
    // через domain.geom в момент вызова emit().
}

ParticlesArray* SecondaryEmissionModel::find_product(
    std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) const {
    auto it = all_species.find(product_);
    if (it == all_species.end()) {
        return nullptr;
    }
    return it->second.get();
}

Vector3R SecondaryEmissionModel::reflect_inward(const Vector3R& coord, const Domain& domain) const {
    const auto& box_min = domain.geom.box_min;
    const auto& box_max = domain.geom.box_max;
    Vector3R res = coord;
    switch (face_) {
        case Face::XMIN:
            res.x() = 2.0 * box_min.x() - coord.x();
            break;
        case Face::XMAX:
            res.x() = 2.0 * box_max.x() - coord.x();
            break;
        case Face::YMIN:
            res.y() = 2.0 * box_min.y() - coord.y();
            break;
        case Face::YMAX:
            res.y() = 2.0 * box_max.y() - coord.y();
            break;
        case Face::ZMIN:
            res.z() = 2.0 * box_min.z() - coord.z();
            break;
        case Face::ZMAX:
            res.z() = 2.0 * box_max.z() - coord.z();
            break;
        case Face::CYLINDER: {
            // После проверки на этапе загрузки (load_from_json) эти ветки
            // недостижимы, но внутри omp parallel бросать нельзя — возвращаем
            // входную координату без изменений.
            if (!domain.geom.use_cylinder) {
                return coord;
            }
            const Vector3R& c = domain.geom.cyl_center;
            const Vector3R r(coord.x() - c.x(), coord.y() - c.y(), 0.0);
            const double R = r.norm();
            if (R <= 0.0) {
                return coord;
            }
            const double scale = (2.0 * domain.geom.cyl_radius - R) / R;
            res.x() = c.x() + r.x() * scale;
            res.y() = c.y() + r.y() * scale;
            break;
        }
        default:
            break;
    }
    return res;
}

Vector3R SecondaryEmissionModel::inward_normal(const Vector3R& coord, const Domain& domain) const {
    switch (face_) {
        case Face::XMIN:
            return Vector3R(1, 0, 0);
        case Face::XMAX:
            return Vector3R(-1, 0, 0);
        case Face::YMIN:
            return Vector3R(0, 1, 0);
        case Face::YMAX:
            return Vector3R(0, -1, 0);
        case Face::ZMIN:
            return Vector3R(0, 0, 1);
        case Face::ZMAX:
            return Vector3R(0, 0, -1);
        case Face::CYLINDER: {
            // После проверки на этапе загрузки (load_from_json) эти ветки
            // недостижимы, но внутри omp parallel бросать нельзя — возвращаем
            // нормаль внутрь домена по умолчанию (+X).
            if (!domain.geom.use_cylinder) {
                return Vector3R(1, 0, 0);
            }
            const Vector3R& c = domain.geom.cyl_center;
            const Vector3R r(coord.x() - c.x(), coord.y() - c.y(), 0.0);
            const double R = r.norm();
            if (R <= 0.0) {
                return Vector3R(1, 0, 0);
            }
            return Vector3R(-r.x() / R, -r.y() / R, 0.0);   // единичная нормаль внутрь домена
        }
        default:
            return Vector3R(1, 0, 0);
    }
}

void SecondaryEmissionModel::emit(const Particle& src, ParticlesArray& src_species,
                                  std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species,
                                  BoundaryEmitter& emitter, const Domain& domain) {
    // После validate_emissions() сорт-продукт гарантированно существует, но
    // внутри omp parallel бросать нельзя — при отсутствии сорта только
    // предупреждаем и ничего не эмитируем.
    ParticlesArray* product_species = find_product(all_species);
    if (product_species == nullptr) {
        std::cerr << "second_emission: product species \"" << product_ << "\" not found, emission skipped\n";
        return;
    }
    const double mpw_dst = product_species->mpw();
    const double mass_dst = product_species->mass();

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    std::uniform_real_distribution<double> uniform_phi(0.0, 2.0 * kPi);

    // Ламбертовское (косинусное) распределение направлений вокруг нормали n
    // внутрь домена: cosθ = sqrt(u1), φ = 2π·u2.
    const auto sample_lambertian = [&](const Vector3R& n, double speed) -> Vector3R {
        Vector3R t1;
        if (std::abs(n.x()) < 0.9) {
            t1 = Vector3R(1, 0, 0).cross(n);
        } else {
            t1 = Vector3R(0, 1, 0).cross(n);
        }
        t1 = t1.normalized();
        const Vector3R t2 = n.cross(t1);

        const double cos_theta = std::sqrt(uniform01(eng_));
        const double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        const double phi = uniform_phi(eng_);
        const Vector3R dir = cos_theta * n + sin_theta * (std::cos(phi) * t1 + std::sin(phi) * t2);
        return speed * dir;
    };

    for (const auto& rule : rules_) {
        if (rule.species != src_species.name() || rule.yield <= 0.0)
            continue;

        // Веса макрочастиц у сортов могут отличаться — пересчитываем
        // ожидаемое число вторичных в макрочастицах продукта.
        const double expected_macro = rule.yield * (src_species.mpw() / mpw_dst);
        int n = static_cast<int>(std::floor(expected_macro));
        if (uniform01(eng_) < (expected_macro - static_cast<double>(n)))
            n += 1;

        for (int i_p = 0; i_p < n; ++i_p) {
            const Vector3R coord = reflect_inward(src.coord, domain);
            const Vector3R nrm = inward_normal(coord, domain);
            Vector3R v;
            switch (rule.energy_type) {
                case EmissionSourceRule::EnergyType::Fixed: {
                    const double speed = convert_kev_to_velocity(rule.fixed_kev, mass_dst);
                    v = sample_lambertian(nrm, speed);
                    break;
                }
                case EmissionSourceRule::EnergyType::Temperature: {
                    GaussianVelocity gauss(rule.temperature_mean,
                                           convert_kev_to_sigma(rule.temperature_kev, mass_dst));
                    v = gauss.sample(eng_);
                    // Отражаем нормальную компоненту внутрь домена
                    const double vn = v.dot(nrm);
                    if (vn < 0.0)
                        v -= (2.0 * vn) * nrm;
                    break;
                }
                case EmissionSourceRule::EnergyType::Fraction: {
                    const double E_inc = 0.5 * src_species.mass() * src.velocity.dot(src.velocity);
                    const double E_em = rule.fraction * E_inc;
                    double speed = 0.0;
                    if (E_em > 0.0)
                        speed = std::sqrt(2.0 * E_em / mass_dst);
                    v = sample_lambertian(nrm, speed);
                    break;
                }
            }
            emitter.emit_to_species(product_, Particle{coord, v});
            product_species->diag.add_emitted(face_, get_energy_particle(v, mass_dst, mpw_dst));
        }
    }
}
