// boundary_condition.h
#pragma once

#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "Particle.h"
#include "World.h"
#include "config.h"
#include "containers.h"
#include "indexing.h"
#include "particles_distribution_collection.h"
#include "random.h"

// Результат обработки частицы на границе. Unhandled — условие не сработало,
// Removed — частица поглощена/удалена, Reflected — отражена обратно в домен.
enum class ParticleFate { Unhandled, Removed, Reflected };

struct ParticleFateResult {
    ParticleFate fate = ParticleFate::Unhandled;
    Face face = Face::UNAVAILABLE;
};

static inline Face string_to_face(const std::string& s) {
    std::string upper = s;
    std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) { return std::toupper(c); });

    if (upper == "XMIN")
        return Face::XMIN;
    if (upper == "XMAX")
        return Face::XMAX;
    if (upper == "YMIN")
        return Face::YMIN;
    if (upper == "YMAX")
        return Face::YMAX;
    if (upper == "ZMIN")
        return Face::ZMIN;
    if (upper == "ZMAX")
        return Face::ZMAX;
    if (upper == "CYLINDER")
        return Face::CYLINDER;

    // Если ничего не подошло — бросаем исключение
    throw std::invalid_argument("Unknown face string: " + s);
}

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
// -----------------------------------------------
// Базовый класс для всех дополнительных граничных условий
// -----------------------------------------------
class BoundaryCondition {
   public:
    virtual ~BoundaryCondition() = default;
    BoundaryCondition(Face face) : face_(face) {
    }

    // Модификация матричного оператора
    virtual void modify_curlE_stencil(int /*i*/, int /*j*/, int /*k*/, std::vector<Trip>& /*trips*/,
                                      const Domain& /*domain*/) {
        // по умолчанию ничего не делает
    }
    // Модификация матричного оператора
    virtual void modify_curlB_stencil(int /*i*/, int /*j*/, int /*k*/, std::vector<Trip>& /*trips*/,
                                      const Domain& /*domain*/) {
        // по умолчанию ничего не делает
    }

    // TODO: instead name use sort properties
    virtual ParticleFateResult apply_to_particle(const Particle& /*p*/, ParticlesArray& /*particles*/,
                                                 BoundaryEmitter& /*emitter*/, const Domain& /*domain*/) {
        return {};   // по умолчанию ничего не делаем
    }
    // Применение к полям (например, задать значение на границе)
    virtual void apply_to_fields(Field3d& /*fields*/, FieldType /*field*/, const Domain& /*domain*/) {
        // по умолчанию ничего
    }
    // Разовая инициализация электрического поля на границе (вызывается при старте)
    virtual void init_electric_field(Field3d& /*fieldE*/, const Domain& /*domain*/) const {
        // по умолчанию ничего
    }
    virtual void apply_to_operator(Operator& /*mat*/, const Domain& /*domain*/) {
        // по умолчанию ничего
    }
    // Грань, для которой задано условие (у OpenBoundaryConditionArray — UNAVAILABLE,
    // т.к. оно может покрывать сразу несколько граней).
    virtual Face face() const {
        return face_;
    }
    Face face_;
};

// -----------------------------------------------
// Конкретные реализации
// -----------------------------------------------

class OpenBoundaryConditionArray : public BoundaryCondition {
   public:
    static constexpr int maxBc = static_cast<int>(Face::Count);

    OpenBoundaryConditionArray(std::vector<Face> faces, int gap = 0, double eps = 1e-10)
        : BoundaryCondition(Face::UNAVAILABLE), createdBc_(std::ssize(faces)) {
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

    OpenBoundaryConditionArray(Face face, double gap = 0.0) : BoundaryCondition(face), createdBc_(1) {
        faces_[0] = face;
        gaps_[0] = gap;
        epss_[0] = 1e-10;
    }

    void apply_to_operator(Operator& mat, const Domain& domain) override;
    ParticleFateResult apply_to_particle(const Particle& particle, ParticlesArray& particles, BoundaryEmitter& emitter,
                                         const Domain& domain) override;
    void apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) override {
        RECORD_TIMER;

        auto set_zero = [this](double& value, int i, int j, int k, int d, const Domain& dom, FieldType ft) {
            auto pos = dom.get_node_position(i, j, k, ft, d);
            for (int numBc = 0; numBc < createdBc_; ++numBc) {
                if (dom.geom.is_outside_face(faces_[numBc], pos, -gaps_[numBc] + epss_[numBc])) {
                    value = 0.0;
                    break;
                }
            }
        };

        if (field_type == FieldType::CURRENT || field_type == FieldType::ELECTRIC) {
            auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); ++i)
                for (int j = 0; j < size.y(); ++j)
                    for (int k = 0; k < size.z(); ++k)
                        for (int d = 0; d < 3; ++d)
                            set_zero(field(i, j, k, d), i, j, k, d, domain, FieldType::ELECTRIC);
        } else if (field_type == FieldType::MAGNETIC) {
            auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); ++i)
                for (int j = 0; j < size.y(); ++j)
                    for (int k = 0; k < size.z(); ++k)
                        for (int d = 0; d < 3; ++d)
                            set_zero(field(i, j, k, d), i, j, k, d, domain, FieldType::MAGNETIC);
        } else if (field_type == FieldType::DENSITY) {
            auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); ++i)
                for (int j = 0; j < size.y(); ++j)
                    for (int k = 0; k < size.z(); ++k)
                        set_zero(field(i, j, k, 0), i, j, k, 0, domain, FieldType::DENSITY);
        } else {
            std::cerr << "OpenBoundaryCondition: unsupported field type\n";
            assert(false);
        }
    }

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

class PeriodicBoundaryCondition : public BoundaryCondition {
   public:
    PeriodicBoundaryCondition(Face face) : BoundaryCondition(face){};

    void apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) override;
    void apply_to_operator(Operator& mat, const Domain& domain) override;
};

class OpenBoundaryCondition : public BoundaryCondition {
   public:
    OpenBoundaryCondition(Face face, double gap = 0.0) : BoundaryCondition(face), gap_(gap), eps_(1.e-10) {
    }

    void apply_to_operator(Operator& mat, const Domain& domain) override;
    ParticleFateResult apply_to_particle(const Particle& particle, ParticlesArray& particles, BoundaryEmitter& emitter,
                                         const Domain& domain) override;
    void apply_to_fields(Field3d& field, FieldType field_type, const Domain& domain) override {
        RECORD_TIMER;

        auto set_zero = [gap = gap_, eps = eps_](double& value, int i, int j, int k, int d, Face face,
                                                 const Domain& dom, FieldType ft) {
            auto pos = dom.get_node_position(i, j, k, ft, d);
            if (dom.geom.is_outside_face(face, pos, -gap + eps))
                value = 0.0;
        };

        if (field_type == FieldType::CURRENT || field_type == FieldType::ELECTRIC) {
            auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); ++i)
                for (int j = 0; j < size.y(); ++j)
                    for (int k = 0; k < size.z(); ++k)
                        for (int d = 0; d < 3; ++d)
                            set_zero(field(i, j, k, d), i, j, k, d, face_, domain, FieldType::ELECTRIC);
        } else if (field_type == FieldType::MAGNETIC) {
            auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); ++i)
                for (int j = 0; j < size.y(); ++j)
                    for (int k = 0; k < size.z(); ++k)
                        for (int d = 0; d < 3; ++d)
                            set_zero(field(i, j, k, d), i, j, k, d, face_, domain, FieldType::MAGNETIC);
        } else if (field_type == FieldType::DENSITY) {
            auto size = field.sizes();
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); ++i)
                for (int j = 0; j < size.y(); ++j)
                    for (int k = 0; k < size.z(); ++k)
                        set_zero(field(i, j, k, 0), i, j, k, 0, face_, domain, FieldType::DENSITY);
        } else {
            std::cerr << "OpenBoundaryCondition: unsupported field type\n";
            assert(false);
        }
    }

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
    ParticleFateResult apply_to_particle(const Particle& p, ParticlesArray& particles, BoundaryEmitter& emitter,
                                         const Domain& domain) override;

   private:
    bool is_inside_central_circle(const Vector3R& pos, const Domain& domain) const;
    double radius_;
    double energy_threshold_;
};

// Условие Er0: задаёт радиальное электрическое поле в кольцевой области
// на граничном слое (k=ghost для ZMIN/ZMAX). Вызывается один раз при инициализации.
class Er0Condition : public BoundaryCondition {
   public:
    Er0Condition(Face face, double inner_radius, double width, double potential_drop)
        : BoundaryCondition(face), inner_radius_(inner_radius), width_(width), potential_drop_(potential_drop) {
    }
    void init_electric_field(Field3d& fieldE, const Domain& domain) const override;

   private:
    double inner_radius_;
    double width_;
    double potential_drop_;
};

// Условие Bphi: добавляет диагональ в оператор для магнитного поля на указанной
// границе
class BphiCondition : public OpenBoundaryCondition {
   public:
    BphiCondition(Face face, double gap, double radius, double electron_threshold_energy)
        : OpenBoundaryCondition(face, gap), radius_(radius), electron_threshold_energy_(electron_threshold_energy) {
    }
    ParticleFateResult apply_to_particle(const Particle& particle, ParticlesArray& particles, BoundaryEmitter& emitter,
                                         const Domain& domain) override;
    void modify_curlB_stencil(const int i, const int j, const int k, std::vector<Trip>& trips,
                              const Domain& domain) override {
        const auto cell_size = domain.cell_size();
        const int vindx = domain.vind(i, j, k, 0);
        const int vindy = domain.vind(i, j, k, 1);

        const int k_wall = (face_ == Face::ZMIN) ? 1 : domain.size().z() - 2;

        auto pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, 0);
        bool in_region = domain.geom.contains_ignoring_face(face_, pos);
        // (x)[i+1/2,j,k]
        if (k == k_wall && in_region) {
            // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
            double val = 1.0 / cell_size.y();
            trips.push_back(Trip(vindx, domain.vind(i, j, k, 2), val));
            trips.push_back(Trip(vindx, domain.vind(i, j - 1, k, 2), -val));
            // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
            val = -1.0 / cell_size.z();
            trips.push_back(Trip(vindx, domain.vind(i, j, k, 1), val));
            trips.push_back(Trip(vindx, domain.vind(i, j, k - 1, 1), -val));
        }

        pos = domain.get_node_position(i, j, k, FieldType::ELECTRIC, 1);
        in_region = domain.geom.contains_ignoring_face(face_, pos);
        // (y)[i,j+1/2,k]
        if (k == k_wall && in_region) {
            // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
            double val = -1.0 / cell_size.x();
            trips.push_back(Trip(vindy, domain.vind(i, j, k, 2), val));
            trips.push_back(Trip(vindy, domain.vind(i - 1, j, k, 2), -val));
            // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
            val = 1.0 / cell_size.z();
            trips.push_back(Trip(vindy, domain.vind(i, j, k, 0), val));
            trips.push_back(Trip(vindy, domain.vind(i, j, k - 1, 0), -val));
        }
    }

    void modify_curlE_stencil(const int i, const int j, const int k, std::vector<Trip>& trips,
                              const Domain& domain) override {
        const auto cell_size = domain.cell_size();
        const int vindz = domain.vind(i, j, k, 2);
        auto pos = domain.get_node_position(i, j, k, FieldType::MAGNETIC, Z);
        bool in_region = domain.geom.contains_ignoring_face(face_, pos);
        const int k_wall = (face_ == Face::ZMIN) ? 1 : domain.size().z() - 2;

        if (k == k_wall && in_region) {
            // (z)[i+1/2,j+1/2,k]
            // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
            double val = 1.0 / cell_size.x();
            trips.push_back(Trip(vindz, domain.vind(i + 1, j, k, 1), val));
            trips.push_back(Trip(vindz, domain.vind(i, j, k, 1), -val));
            // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
            val = -1.0 / cell_size.y();
            trips.push_back(Trip(vindz, domain.vind(i, j + 1, k, 0), val));
            trips.push_back(Trip(vindz, domain.vind(i, j, k, 0), -val));
        }
    }

   private:
    // Проверяет, находится ли точка (x,y) внутри центрального круга радиуса
    // radius_
    // Центр круга определяется как геометрический центр box’а по X и Y.
    bool is_inside_central_circle(const Vector3R& pos, const Domain& domain) const {
        double center_x = 0.5 * domain.num_cells().x() * domain.cell_size().x();
        double center_y = 0.5 * domain.num_cells().y() * domain.cell_size().y();
        double dx = pos.x() - center_x;
        double dy = pos.y() - center_y;
        return (dx * dx + dy * dy) <= radius_ * radius_;
    }
    // Для электронов: возвращает true, если частицу надо отразить (а не
    // потерять) Отражаем только если частица внутри центрального круга И её
    // кинетическая энергия по Z меньше пороговой величины.
    bool should_electron_reflect(bool inside_circle, double vz, double m) const {
        double kinetic_z = 0.5 * m * vz * vz;   // здесь vz уже в единицах скорости?
        // Если внутри круга и энергия мала – отражаем (true), иначе проникает
        // (false)
        return inside_circle && (kinetic_z <= electron_threshold_energy_);
    }

    double radius_;   // ограничение по радиусу
    double electron_threshold_energy_;
};

// Можно добавить другие: FieldMirrorCondition, AbsorbingCondition и т.д.

// -----------------------------------------------
// Вторичная эмиссия на границе
// -----------------------------------------------
// Правило эмиссии для одного сорта-источника
struct EmissionSourceRule {
    std::string species;   // имя сорта-источника
    double yield = 0.0;   // среднее число вторичных на ОДНУ физическую частицу источника (может быть дробным)
    enum class EnergyType { Fixed, Temperature, Fraction };
    EnergyType energy_type = EnergyType::Fixed;
    double fixed_kev = 0.0;            // Fixed: моноэнергия [кэВ]
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
    SecondaryEmissionModel(Face face, std::string product_species, std::vector<EmissionSourceRule> rules);
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
    // RNG eng_ и аккумуляторы diag (SpeciesDiagStats) здесь безопасно мутировать.
    void emit(const Particle& src, ParticlesArray& src_species,
              std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species, BoundaryEmitter& emitter,
              const Domain& domain);

   private:
    // Ищет сорт-продукт в all_species; возвращает nullptr, если сорта нет
    // (никуда не бросает). Вызывается только из emit (последовательная фаза).
    ParticlesArray* find_product(
        std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) const;
    Vector3R reflect_inward(const Vector3R& coord,
                            const Domain& domain) const;   // отражает координату внутрь домена относительно face_
    Vector3R inward_normal(const Vector3R& coord, const Domain& domain) const;   // единичная нормаль внутрь домена
    Face face_;
    std::string product_;
    std::vector<EmissionSourceRule> rules_;
    LehmerEngine eng_;
};

// -----------------------------------------------
// Основной обработчик – контейнер условий
// -----------------------------------------------
class BoundaryConditionHandler {
   public:
    // Загружает условия из JSON (ожидается массив объектов)
    void load_from_json(const nlohmann::json& sys_config, const Domain& domain);

    // Применяет все подходящие условия к стенсилу (для матрицы)
    void modify_curlB_stencil(int i, int j, int k, std::vector<Trip>& trips, const Domain& domain) const {
        RECORD_TIMER;
        for (const auto& cond : conditions_) {
            cond->modify_curlB_stencil(i, j, k, trips, domain);
        }
    }
    void modify_curlE_stencil(int i, int j, int k, std::vector<Trip>& trips, const Domain& domain) const {
        RECORD_TIMER;
        for (const auto& cond : conditions_) {
            cond->modify_curlE_stencil(i, j, k, trips, domain);
        }
    }
    // Применяет все условия к частицам
    void apply_to_particles(ParticlesArray& particles,
                            std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& species,
                            const Domain& domain);

    // Применяет все условия к полям
    void apply_to_fields(Field3d& fields, FieldType field_t, const Domain& domain) const {
        RECORD_TIMER;
        for (const auto& cond : conditions_) {
            cond->apply_to_fields(fields, field_t, domain);
        }
    }
    // Разовая инициализация электрического поля на границе (вызывается при старте)
    void init_electric_field(Field3d& fieldE, const Domain& domain) const {
        RECORD_TIMER;
        for (const auto& cond : conditions_) {
            cond->init_electric_field(fieldE, domain);
        }
    }
    void apply_to_operator(Operator& mat, const Domain& domain) const {
        RECORD_TIMER;
        for (const auto& cond : conditions_) {
            cond->apply_to_operator(mat, domain);
        }
    }
    void flush_species(std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species);
    // Проверяет, что сорта-продукты всех моделей вторичной эмиссии существуют
    // в all_species (вызывается ПОСЛЕ инициализации всех сортов, до основного
    // цикла). Бросает std::runtime_error с перечислением отсутствующих сортов.
    void validate_emissions(
        const std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) const;
    // Если нужно знать, есть ли активные условия (например, чтобы не
    // вызывать apply... без нужды)
    bool empty() const {
        return conditions_.empty() && emissions_.empty();
    }
    // Проверить, является ли ось периодической (0=X,1=Y,2=Z)
    bool is_periodic(int axis) const {
        return periodic_[axis];
    }

    IndexRange active_range(const Grid& grid) const {
        IndexRange range;
        range.start = Vector3I(0, 0, 0);
        range.end = Vector3I(grid.size().x(), grid.size().y(), grid.size().z());

        if (is_periodic(0)) {
            range.start.x() = grid.ghost_cells();
            range.end.x() -= (grid.ghost_cells() + 1);
        }
        if (is_periodic(1)) {
            range.start.y() = grid.ghost_cells();
            range.end.y() -= (grid.ghost_cells() + 1);
        }
        if (is_periodic(2)) {
            range.start.z() = grid.ghost_cells();
            range.end.z() -= (grid.ghost_cells() + 1);
        }
        return range;
    }
    Vector3R wrap_periodic(const Vector3R& coord, const Domain& domain) const {
        Vector3R res = coord;
        const auto& box_min = domain.geom.box_min;
        const auto& box_max = domain.geom.box_max;
        const double Lx = box_max.x() - box_min.x();
        const double Ly = box_max.y() - box_min.y();
        const double Lz = box_max.z() - box_min.z();

        if (periodic_[0]) {
            res.x() = box_min.x() + std::fmod(res.x() - box_min.x(), Lx);
            if (res.x() < box_min.x())
                res.x() += Lx;
        }
        if (periodic_[1]) {
            res.y() = box_min.y() + std::fmod(res.y() - box_min.y(), Ly);
            if (res.y() < box_min.y())
                res.y() += Ly;
        }
        if (periodic_[2]) {
            res.z() = box_min.z() + std::fmod(res.z() - box_min.z(), Lz);
            if (res.z() < box_min.z())
                res.z() += Lz;
        }
        return res;
    }
    /// Переносит индекс idx по оси axis (0=X, 1=Y, 2=Z) внутрь допустимого
    /// диапазона, если ось периодическая. Если ось не периодическая,
    /// возвращает idx без изменений.
    int wrap_index(int idx, int axis, const Domain& domain) const {
        if (!periodic_[axis])
            return idx;

        const auto& size = domain.grid.size();
        int all_cells = size[axis] - 1;
        if (idx < 0)
            return all_cells - 2 * domain.grid.ghost_cells() + idx;
        if (idx > size[axis] - 1)
            return idx - all_cells + 2 * domain.grid.ghost_cells();
        return idx;
    }

   private:
    // Добавляет одно условие типа type с параметрами params (один объект грани)
    void add_condition(const std::string& type, const nlohmann::json& params, const Domain& domain);

    std::vector<std::unique_ptr<BoundaryCondition>> conditions_;
    // Модели вторичной эмиссии по граням: key = Face.
    // В conditions_ не добавляются — это НЕ consuming-условия.
    std::unordered_map<Face, std::vector<SecondaryEmissionModel>> emissions_;
    BoundaryEmitter emitter;
    bool periodic_[3] = {false, false, false};
};
