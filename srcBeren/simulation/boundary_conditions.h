#pragma once

#include "fields/field_boundaries.h"
#include "particles/particle_boundaries.h"

class OpenBoundaryConditionArray : public FieldBoundaryCondition, public ParticleBoundaryCondition {
   public:
    static constexpr int maxBc = static_cast<int>(Face::Count);

    OpenBoundaryConditionArray(std::vector<Face> faces, int gap = 0, double eps = 1e-10)
        : FieldBoundaryCondition(Face::UNAVAILABLE), createdBc_(std::ssize(faces)) {
        if (std::ssize(faces) > maxBc || faces.empty()) {
            throw std::runtime_error("OpenBoundaryConditionArray: expected from 1 to " + std::to_string(maxBc) +
                                     " faces, got " + std::to_string(faces.size()));
        }
        for (int i = 0; i < createdBc_; ++i) {
            faces_[i] = faces[i];
            gaps_[i] = gap;
            epss_[i] = eps;
        }
    }

    OpenBoundaryConditionArray(Face face, double gap = 0.0) : FieldBoundaryCondition(face), createdBc_(1) {
        faces_[0] = face;
        gaps_[0] = gap;
        epss_[0] = 1e-10;
    }

    void apply_to_operator(Operator& mat, const Domain& domain) override;
    ParticleFateResult classify_particle(const Particle& p, const ParticlesArray& particles, const Domain& domain,
                                         Face hit_face = Face::UNAVAILABLE) const override;
    void collect_removing_faces(std::vector<Face>& out) const override {
        for (int i = 0; i < createdBc_; ++i) {
            out.push_back(faces_[i]);
        }
    }
    void apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) override;

   private:
    struct Unavailable {};

   public:
    std::array<double, maxBc> gaps_;
    std::array<double, maxBc> epss_;
    std::array<Face, maxBc> faces_;
    // avoid to usage face_ from base class
    Unavailable face_;
    int createdBc_;
};

class PeriodicBoundaryCondition : public FieldBoundaryCondition {
   public:
    PeriodicBoundaryCondition(Face face) : FieldBoundaryCondition(face) {
    }

    void apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) override;
    void apply_to_operator(Operator& mat, const Domain& domain) override;
};

class OpenBoundaryCondition : public FieldBoundaryCondition, public ParticleBoundaryCondition {
   public:
    OpenBoundaryCondition(Face face, double gap = 0.0) : FieldBoundaryCondition(face), gap_(gap), eps_(1.e-10) {
    }

    void apply_to_operator(Operator& mat, const Domain& domain) override;
    ParticleFateResult classify_particle(const Particle& p, const ParticlesArray& particles, const Domain& domain,
                                         Face hit_face = Face::UNAVAILABLE) const override;
    void collect_removing_faces(std::vector<Face>& out) const override {
        out.push_back(face_);
    }
    void apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) override;

    double gap_;
    double eps_;
};

// Отражение электронов на Z-границе: если электрон внутри центрального круга
// заданного радиуса и его z-кинетическая энергия ниже порога — отражается
// (vz = -vz, координата отражается). Наследует OpenBoundaryCondition, поэтому
// используется ВМЕСТО open для той же грани — порядок не имеет значения.
class ElectronReflectionCondition : public OpenBoundaryCondition {
   public:
    ElectronReflectionCondition(Face face, double radius, double energy_threshold, double gap = 0.0)
        : OpenBoundaryCondition(face, gap), radius_(radius), energy_threshold_(energy_threshold) {
    }
    ParticleFateResult classify_particle(const Particle& p, const ParticlesArray& particles, const Domain& domain,
                                         Face hit_face = Face::UNAVAILABLE) const override;

   private:
    bool is_inside_central_circle(const Vector3R& pos, const Domain& domain) const;
    double radius_;
    double energy_threshold_;
};

// Профиль радиального поля Er0 на кольце:
//   Step   — постоянное поле dPhi/L в кольце [inner_radius, inner_radius + L];
//   Cosine — Er = dPhi/L * (1 + cos(2*pi*(r - r0)/L)) на [r0 - L/2, r0 + L/2].
// Оба профиля дают одинаковый перепад потенциала dPhi.
enum class Er0Profile { Step, Cosine };
enum class Er0Mode { Overwrite, Solve };

// Условие Er0: задаёт радиальное электрическое поле в кольцевой области
// на граничном слое (k=ghost для ZMIN / k=size_z-ghost-1 для ZMAX) и работает
// как open на этой грани (используется ВМЕСТО open; дублировать "open" на той
// же грани нельзя).
// mode="overwrite" (по умолчанию): master-совместимое поведение — после решения
// базовый OpenBoundaryCondition зануляет слой, затем профиль кольца
// перезаписывается в динамическое поле.
// mode="solv": единая электрическая Дирихле-маска (нулевые DOF открытого слоя
// плюс предписанные Ex/Ey кольца). Lmat2 фильтруется ТОЛЬКО по строкам (столбцы
// сохраняются, чтобы свободные уравнения удерживали отклик тока на
// предписанное поле); Дирихле-строки один раз ставятся в IMmat, правая часть
// поправляется на каждом решении. В обоих режимах apply_to_fields затем
// перезаписывает кольцо.
class Er0Condition : public OpenBoundaryCondition {
   public:
    Er0Condition(Face face, Er0Profile profile, double inner_radius, double r0, double width, double potential_drop,
                 std::string mode = "overwrite", double gap = 0.0)
        : OpenBoundaryCondition(face, gap),
          profile_(profile),
          inner_radius_(inner_radius),
          r0_(r0),
          width_(width),
          potential_drop_(potential_drop),
          mode_(mode == "solv" ? Er0Mode::Solve : Er0Mode::Overwrite) {
    }
    void init_electric_field(Field3d& fieldE, const Domain& domain) const override;

    void apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) override {
        OpenBoundaryCondition::apply_to_fields(field, field_type, domain);
        if (field_type == FieldType::ELECTRIC)
            init_electric_field(field, domain);
    }

    void apply_to_operator(Operator& mat, const Domain& domain) override;

    void collect_electric_dirichlet(std::vector<ElectricDirichletEntry>& entries, const Domain& domain) const override {
        if (mode_ != Er0Mode::Solve)
            return;
        collect_entries(entries, domain);
    }

   private:
    // Единый источник геометрии: нули на всём открытом слое грани и профиль
    // кольца в предписанных Ex/Ey. Результат кешируется (домен неизменен).
    void collect_entries(std::vector<ElectricDirichletEntry>& out, const Domain& domain) const;

    double radial_field(double radius) const;

    Er0Profile profile_;
    double inner_radius_;     // Step: внутренний радиус кольца
    double r0_;               // Cosine: центр кольца
    double width_;            // L: ширина кольца
    double potential_drop_;   // dPhi (уже поделён на SGS::MC2 при разборе JSON)
    Er0Mode mode_;
    mutable std::vector<ElectricDirichletEntry> entries_;
    mutable bool entriesReady_ = false;
};

// Coordinates transport, sparse wall events and deterministic secondary emission.
class BoundaryConditionHandler : public FieldBoundaryConditions {
   public:
    BoundaryConditionHandler() = default;
    BoundaryConditionHandler(std::vector<std::unique_ptr<FieldBoundaryCondition>> conditions,
                             BoundaryEmissionModels emissions, std::array<bool, 3> periodic)
        : FieldBoundaryConditions(std::move(conditions), periodic), emissions_(std::move(emissions)) {
        std::vector<Face> faces;
        for (const std::unique_ptr<FieldBoundaryCondition>& condition : conditions_) {
            if (ParticleBoundaryCondition* particle = dynamic_cast<ParticleBoundaryCondition*>(condition.get())) {
                particle_conditions_.push_back(particle);
                faces.clear();
                particle->collect_removing_faces(faces);
                for (Face face : faces) pathPolicy_.enable(face);
            }
        }
    }

    // Применяет все условия к частицам
    ParticleFateResult classify_particle(const Particle& p, const ParticlesArray& particles, const Domain& domain,
                                         Face hit_face = Face::UNAVAILABLE) const;
    void apply_to_particles(ParticlesArray& particles,
                            std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species,
                            const Domain& domain, ParticleBoundaryEvents* captured = nullptr);

    void flush_species(std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species);
    // Проверяет, что сорта-продукты всех моделей вторичной эмиссии существуют
    // в all_species (вызывается ПОСЛЕ инициализации всех сортов, до основного
    // цикла). Бросает std::runtime_error с перечислением отсутствующих сортов.
    void validate_emissions(const std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) const;
    // Политика обрезки траекторий на удаляющих гранях (для депозита тока).
    PathPolicy path_policy(const ParticlesArray& particles, const Domain& domain,
                           ParticleBoundaryEvents* events = nullptr) const {
        PathPolicy result = pathPolicy_;
        result.events = events;
        result.classify = [this, &particles, &domain](const Particle& p, Face face) {
            return classify_particle(p, particles, domain, face);
        };
        return result;
    }
    const PathPolicy& path_policy() const {
        return pathPolicy_;
    }

   private:
    // Pointees belong to conditions_; moving the handler preserves their addresses.
    std::vector<ParticleBoundaryCondition*> particle_conditions_;
    BoundaryEmissionModels emissions_;
    BoundaryEmitter emitter;
    PathPolicy pathPolicy_;
};
