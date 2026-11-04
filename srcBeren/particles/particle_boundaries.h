#pragma once

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Particle.h"
#include "grid/domain.h"
#include "particles/particle_path.h"
#include "runtime/random.h"

class ParticlesArray;

class BoundaryEmitter {
   public:
    BoundaryEmitter() = default;

    // добавить частицу в текущий сорт
    void emit_current_species(const Particle& p) {
        current_species_buffer_.push_back(p);
    }

    // добавить частицу в другой сорт по имени
    void emit_to_species(const std::string& species_name, const Particle& p) {
        species_buffers_[species_name].push_back(p);
    }

    const std::vector<Particle>& current_species_particles() const {
        return current_species_buffer_;
    }

    const std::unordered_map<std::string, std::vector<Particle>>& other_species_particles() const {
        return species_buffers_;
    }

    void clear_current_species_buffer() {
        current_species_buffer_.clear();
    }
    void clear_other_species_buffers() {
        for (auto& kv : species_buffers_) {
            species_buffers_[kv.first].clear();
        }
    }

   private:
    std::vector<Particle> current_species_buffer_;
    std::unordered_map<std::string, std::vector<Particle>> species_buffers_;
};
class ParticleBoundaryCondition {
   public:
    virtual ~ParticleBoundaryCondition() = default;
    // hit_face specifies an already selected crossing; p.coord is then the hit
    // position. Without it, classify the endpoint for the ordinary ECSIM mover.
    virtual ParticleFateResult classify_particle(const Particle&, const ParticlesArray&, const Domain&,
                                                 Face = Face::UNAVAILABLE) const {
        return {};
    }
    ParticleFateResult apply_to_particle(const Particle& p, ParticlesArray& particles, BoundaryEmitter& emitter,
                                         const Domain& domain) {
        const auto result = classify_particle(p, particles, domain);
        if (result.fate == ParticleFate::Reflected) {
            Particle reflected = p;
            reflected.coord = domain.geom.reflect_from_face(result.face, p.coord);
            reflected.velocity.z() = -p.velocity.z();
            emitter.emit_current_species(reflected);
        }
        return result;
    }

    virtual void collect_removing_faces(std::vector<Face>& out) const = 0;
};
// Правило эмиссии для одного сорта-источника
struct EmissionSourceRule {
    std::string species;   // имя сорта-источника
    enum class YieldModel { Constant, Vaughan, Threshold };
    YieldModel yield_model = YieldModel::Constant;
    double yield = 0.0;   // Constant / Threshold: среднее число вторичных на одну физическую частицу источника
    double delta_max = 0.0;   // Vaughan: пиковый коэффициент вторичной эмиссии (обычно 1..3)
    double energy_max_kev = 0.0;   // Vaughan: энергия налетающей частицы в максимуме кривой [кэВ]
    double threshold_kev = 0.0;   // Vaughan/Threshold: ниже этой энергии эмиссия нулевая [кэВ]
    enum class EnergyType { Fixed, Temperature, Fraction };
    EnergyType energy_type = EnergyType::Fixed;
    double fixed_kev = 0.0;   // Fixed: моноэнергия [кэВ]
    // Temperature: температура kT по компонентам [кэВ] (не разброс энергий!).
    // Преобразуется через convert_kev_to_sigma: σ_v = sqrt(kT / (MC2·m)).
    Vector3R temperature_kev{0, 0, 0};
    // Temperature: средняя скорость (дрейфовая скорость) в единицах кода
    // (намеренная асимметрия относительно sigma).
    Vector3R temperature_mean{0, 0, 0};
    double fraction = 0.0;   // Fraction: доля кинетической энергии налетевшей частицы, [0..1]
};

class SecondaryEmissionModel {
   public:
    SecondaryEmissionModel(Face face, const std::string& product_species, const std::vector<EmissionSourceRule>& rules);
    Face face() const {
        return face_;
    }
    const std::string& product() const {
        return product_;
    }
    // Вызывается handler'ом только для частиц, ПОГЛОЩЁННЫХ на грани face_
    // (fate==Removed).
    // Потокобезопасность: вызывается только в ПОСЛЕДОВАТЕЛЬНОЙ фазе
    // apply_to_particles (после завершения параллельной классификации), поэтому
    // RNG eng_ и события шага (SpeciesStepState) здесь безопасно мутировать.
    void emit(const Particle& src, ParticlesArray& src_species,
              std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species, BoundaryEmitter& emitter,
              const Domain& domain);

   private:
    // Ищет сорт-продукт в all_species; возвращает nullptr, если сорта нет
    // (никуда не бросает). Вызывается только из emit (последовательная фаза).
    ParticlesArray* find_product(std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) const;
    Vector3R reflect_inward(const Vector3R& coord,
                            const Domain& domain) const;   // отражает координату внутрь домена относительно face_
    Vector3R inward_normal(const Vector3R& coord, const Domain& domain) const;   // единичная нормаль внутрь домена
    Face face_;
    std::string product_;
    std::vector<EmissionSourceRule> rules_;
    LehmerEngine eng_;
};

using BoundaryEmissionModels = std::array<std::vector<SecondaryEmissionModel>, static_cast<size_t>(Face::Count)>;
