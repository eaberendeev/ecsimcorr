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
#include "containers.h"
#include "particles_distribution_collection.h"
#include "config.h"
#include "indexing.h"

static inline Face string_to_face(const std::string& s) {
    std::string upper = s;
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   [](unsigned char c) { return std::toupper(c); });

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

static inline std::vector<Face> faces_except(Face excluded) {
    std::vector<Face> result;
    for (Face f : ALL_FACES) {
        if (f != excluded) {
            result.push_back(f);
        }
    }
    return result;
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

    const std::unordered_map<std::string, std::vector<Particle>>&
    other_species_particles() const {
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
    BoundaryCondition(Face face) : face_(face) {}

    // Модификация матричного оператора
    virtual void modify_curlE_stencil(int /*i*/, int /*j*/, int /*k*/,
                                      std::vector<Trip>& /*trips*/,
                                      const Domain& /*domain*/) const {
        // по умолчанию ничего не делает
    }
    // Модификация матричного оператора
    virtual void modify_curlB_stencil(int /*i*/, int /*j*/, int /*k*/,
                                      std::vector<Trip>& /*trips*/,
                                      const Domain& /*domain*/) const {
        // по умолчанию ничего не делает
    }

    // TODO: instead name use sort properties 
    virtual void apply_to_particle(const Particle& /*p*/,
                                   const std::string& /*species_name*/,
                                   BoundaryEmitter& /*emitter*/,
                                   const Domain& /*domain*/) {
        // по умолчанию ничего не делаем
    }
    // Применение к полям (например, задать значение на границе)
    virtual void apply_to_fields(Field3d& /*fields*/, FieldType /*field*/,
                                 const Domain& /*domain*/) const {
        // по умолчанию ничего
    }
    virtual void apply_to_operator(Operator& /*mat*/,
                                 const Domain& /*domain*/) const {
        // по умолчанию ничего
    }
    Face face_;
};

// -----------------------------------------------
// Конкретные реализации
// -----------------------------------------------

class PeriodicBoundaryCondition : public BoundaryCondition {
   public:
    PeriodicBoundaryCondition(Face face) : BoundaryCondition(face){};

    void apply_to_particle(const Particle& particle,
                           const std::string& species_name,
                           BoundaryEmitter& emitter,
                           const Domain& domain) override;
    void apply_to_fields(Field3d& field, FieldType field_type,
                         const Domain& domain) const override;
    void apply_to_operator(Operator& mat, const Domain& domain) const override;
};

class OpenBoundaryCondition : public BoundaryCondition {
   public:
    OpenBoundaryCondition(Face face) : BoundaryCondition(face){};
    void apply_to_operator(Operator& mat, const Domain& domain) const override;
    void apply_to_fields(Field3d& field, FieldType field_type,
                         const Domain& domain) const override {
        if (field_type == FieldType::CURRENT ||
            field_type == FieldType::ELECTRIC) {
            auto size = field.sizes();
            auto set_zero = [](double& value, int i, int j, int k, int d, Face face,
                                    const Domain& domain ) {
                auto pos =
                    domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
                bool setZero = domain.geom.is_outside_face(face, pos, 1.e-12);
                if (setZero) {
                    value = 0.;
                }
            };
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); i++) {
                for (int j = 0; j < size.y(); j++) {
                    for (int k = 0; k < size.z(); k++) {
                        for (int d = 0; d < 3; d++) {
                            set_zero(field(i,j,k,d), i, j, k, d, face_,domain);
                        }
                    }
                }
            }
        } else if (field_type == FieldType::DENSITY) {
            auto size = field.sizes();
            auto set_zero = [](double& value, int i, int j, int k, Face face,
                                    const Domain& domain ) {
                auto pos =
                    domain.get_node_position(i, j, k, FieldType::DENSITY, 0);
                bool setZero = domain.geom.is_outside_face(face, pos, 1.e-12);
                if (setZero) {
                    value = 0.;
                }
            };
#pragma omp parallel for schedule(dynamic, 32) collapse(3)
            for (int i = 0; i < size.x(); i++) {
                for (int j = 0; j < size.y(); j++) {
                    for (int k = 0; k < size.z(); k++) {
                        set_zero(field(i,j,k, 0), i, j, k, face_,domain);
                    }
                }
            }
        } else {
            std::cout << "OpenBoundaryCondition: unsupported field type" << std::endl;
        }
    }
};

class SecondEmissionCondition : public BoundaryCondition {
   public:
    SecondEmissionCondition(Face face, Vector3R mean, Vector3R sigma) : BoundaryCondition(face) {
        other_faces = faces_except(face);
        gauss_.set(mean, sigma);
    }
    void apply_to_particle(const Particle& particle,
                           const std::string& species_name,
                           BoundaryEmitter& emitter,
                           const Domain& domain) override;

   private:
    std::vector<Face> other_faces;
    bool is_outside_other_faces(const Vector3R& p, const Domain& domain) const {
        for (auto oface : other_faces)  {
            if (domain.geom.is_outside_face(oface, p))
                return true;
        }
        return false;
    }
    
    ThreadRandomGenerator pulse_gen_;
    GaussianVelocity gauss_;
};

// Условие Bphi: добавляет диагональ в оператор для магнитного поля на указанной
// границе
class BphiCondition : public BoundaryCondition {
   public:
    BphiCondition(Face face, double Jz = 0.0, double radius = 0.0)
        : BoundaryCondition(face), Jz_(Jz), radius_(radius) {}

   // void modify_curlB_stencil(int i, int j, int k, std::vector<Trip>& trips,
   //                           const Domain& domain) const override {
        // const auto cell_size = domain.cell_size();
        // const int vindx = domain.vind(i, j, k, 0);
        // const int vindy = domain.vind(i, j, k, 1);
        // const int vindz = domain.vind(i, j, k, 2);
        // bool in_region = domain.is_inside_cyl(i, j, k, FieldType::ELECTRIC, X);
        // // (x)[i+1/2,j,k]
        // // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
        // if (k == 1 && in_region) {
        //     // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
        //     if(j>1){
        //     double val = 1.0 / cell_size.y();
        //     trips.push_back(Trip(vindx, domain.vind(i, j, k, 2), val));
        //     trips.push_back(Trip(vindx, domain.vind(i, j - 1, k, 2), -val));
        //     }
        //     double val = -1.0 / cell_size.z();
        //         trips.push_back(Trip(vindx, domain.vind(i, j, k, 1), val));
        //     trips.push_back(Trip(vindx, domain.vind(i, j, k - 1, 1), -val));
        // }

        // in_region = domain.is_inside_cyl(i, j, k, FieldType::ELECTRIC, Y);
        // // (y)[i,j+1/2,k]
        // // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
        // if (k == 1 && in_region) {
        //     if(i>1){
        //         // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
        //         double val = -1.0 / cell_size.x();

        //         trips.push_back(Trip(vindy, domain.vind(i, j, k, 2), val));
        //         trips.push_back(Trip(vindy, domain.vind(i - 1, j, k, 2), -val));
        //     }
        //     double val = 1.0 / cell_size.z();
        //     trips.push_back(Trip(vindy, domain.vind(i, j, k, 0), val));
        //     if (k == 1)
        //         trips.push_back(Trip(vindy, domain.vind(i, j, k - 1, 0), -val));
        // }
   // }

    // void modify_curlE_stencil(int i, int j, int k, std::vector<Trip>& trips,
    //                           const Domain& domain) const override {
    //     const auto cell_size = domain.cell_size();
        //const int vindz = domain.vind(i, j, k, 2);

        // bool in_region = domain.is_inside_cyl(i, j, k, FieldType::MAGNETIC, Z);
        // if (k == 1 && in_region) {
        //     // (z)[i+1/2,j+1/2,k]
        //     // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
        //     double val = 1.0 / cell_size.x();
        //     trips.push_back(Trip(vindz, domain.vind(i + 1, j, k, 1), val));
        //     trips.push_back(Trip(vindz, domain.vind(i, j, k, 1), -val));
        //     // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
        //     val = -1.0 / cell_size.y();
        //     trips.push_back(Trip(vindz, domain.vind(i, j + 1, k, 0), val));
        //     trips.push_back(Trip(vindz, domain.vind(i, j, k, 0), -val));
        // }
   // }

    void apply_to_fields(Field3d& field, FieldType field_t,
                        const Domain& domain) const override{
                            std::cout << "APPLY BPHI " << "\n";

        if (field_t == FieldType::MAGNETIC) {
                int startz = 0;
                int endz = 1;
                set_Bphi(field, Jz_, radius_, startz, endz, domain);
        }
    }

    
   private:
    // Прототип локальной утилиты задания вихревого поля Bphi
    void set_Bphi(Field3d& fieldB, double Jz, double radius,
                  [[maybe_unused]] double startz, [[maybe_unused]] double endz,
                  const Domain& domain) const {
        const auto size_x = fieldB.sizes().x();
        const auto size_y = fieldB.sizes().y();
    //    const auto size_z = fieldB.sizes().z();
        const double dx = domain.cell_size().x();
        const double dy = domain.cell_size().y();
    //    const double dz = domain.cell_size().z();

        // const int sz = std::max(0, static_cast<int>(std::llround(startz / dz)));
        // const int ez = std::min(static_cast<int>(size_z),
        //                         static_cast<int>(std::llround(endz / dz)));

        const double center_x = 0.5 * (size_x - 3) * dx;
        const double center_y = 0.5 * (size_y - 3) * dy;
        const double eps_r = 1e-12;

        for (int k = 0; k < 1; ++k) {
            for (int i = 0; i < size_x; ++i) {
                for (int j = 0; j < size_y; ++j) {
                //     if (!(domain.is_inside(i, j, k, FieldType::MAGNETIC, 0) &&
                // domain.is_inside(i, j, k, FieldType::MAGNETIC, 1))) continue;
                    // X-face node
                    double xx = i * dx - center_x - dx * GHOST_CELLS;
                    double yy = (j + 0.5) * dy - center_y - dy * GHOST_CELLS;
                    double rr = std::hypot(xx, yy);
                    if (rr < radius) {
                        fieldB(i, j, k, Axis::X) = -0.5 * Jz * (yy);
                    } else {
                        rr = std::max(rr, eps_r);
                        fieldB(i, j, k, Axis::X) =
                            -0.5 * radius * radius * Jz * yy / (rr * rr);
                    }

                    // Y-face node
                    yy = j * dy - center_y - dy * GHOST_CELLS;
                    xx = (i + 0.5) * dx - center_x - dx * GHOST_CELLS;
                    rr = std::hypot(xx, yy);
                    if (rr < radius) {
                        fieldB(i, j, k, Axis::Y) = 0.5 * Jz * (xx);
                    } else {
                        rr = std::max(rr, eps_r);
                        fieldB(i, j, k, Axis::Y) =
                            0.5 * radius * radius * Jz * xx / (rr * rr);
                    }
                }
            }
        }
    }
    double Jz_;       // может пригодиться позже
    double radius_;   // ограничение по радиусу
};

// Можно добавить другие: FieldMirrorCondition, AbsorbingCondition и т.д.

// -----------------------------------------------
// Основной обработчик – контейнер условий
// -----------------------------------------------
class BoundaryConditionHandler {
   public:
    // Загружает условия из JSON (ожидается массив объектов)
    void load_from_json(const nlohmann::json& sys_config) {
        conditions_.clear();
        //Domain dom;
        //dom.set_domain(sys_config);
        if (!sys_config.contains("Boundary_conditions"))
            return;
        const auto config = sys_config["Boundary_conditions"];
        if (!config.is_array()) {
            std::cerr << "BoundaryConditionHandler: expected array in JSON\n";
            return;
        }

        for (const auto& item : config) {
            if (item.empty())
                continue;

            auto it = item.begin();
            std::string type = it.key();  
            const auto& params = it.value();
            std::string face_str = params.at("face");
            Face face = string_to_face(face_str);
            if (type == "bphi") {
                double Jz = params.at("Jz");
                double radius = params.at("radius");
                conditions_.push_back(
                    std::make_unique<BphiCondition>(face, Jz, radius));
            } else if (type == "second_emisson") {
                Vector3R mean = util::parse_double3(params.at("mean"));
                Vector3R temperature = util::parse_double3(params.at("sigma"));
                Vector3R sigma = convert_kev_to_sigma(temperature, 1.0);
                conditions_.push_back(
                    std::make_unique<SecondEmissionCondition>(face, mean, sigma));
            } else if (type == "periodic") {
                conditions_.push_back(
                    std::make_unique<PeriodicBoundaryCondition>(face));
                    if (face == Face::XMIN || face == Face::XMAX) {
                        periodic_[0] = true;
                    } else if (face == Face::YMIN || face == Face::YMAX) {
                        periodic_[1] = true;
                    } else if (face == Face::ZMIN || face == Face::ZMAX) {
                        periodic_[2] = true;
                    }
            } else if (type == "open") {
                conditions_.push_back(
                    std::make_unique<OpenBoundaryCondition>(face));
            } else {
                std::cerr << "Unknown boundary condition type: " << type
                          << std::endl;
            }
        }

    }

    // Применяет все подходящие условия к стенсилу (для матрицы)
    void modify_curlB_stencil(int i, int j, int k, std::vector<Trip>& trips,
                              const Domain& domain) const {
        for (const auto& cond : conditions_) {
            cond->modify_curlB_stencil(i, j, k, trips, domain);
        }
    }
    void modify_curlE_stencil(int i, int j, int k, std::vector<Trip>& trips,
                              const Domain& domain) const {
        for (const auto& cond : conditions_) {
            cond->modify_curlE_stencil(i, j, k, trips, domain);
        }
    }
    // Применяет все условия к частицам
    void apply_to_particles(
        ParticlesArray& particles,
        std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>&
            species,
        const Domain& domain);

    // Применяет все условия к полям
    void apply_to_fields(Field3d& fields, FieldType field_t,
                         const Domain& domain) const {
        for (const auto& cond : conditions_) {
            cond->apply_to_fields(fields, field_t, domain);
        }
    }
    void apply_to_operator(Operator& mat,
                         const Domain& domain) const {
        for (const auto& cond : conditions_) {
            cond->apply_to_operator(mat, domain);
        }
    }
    void flush_species(
        std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>&
            all_species);
    // Если нужно знать, есть ли активные условия (например, чтобы не
    // вызывать apply... без нужды)
    bool empty() const { return conditions_.empty(); }
    // Проверить, является ли ось периодической (0=X,1=Y,2=Z)
    bool is_periodic(int axis) const { return periodic_[axis]; }

    IndexRange active_range(
        const Grid& grid) const {
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

   private:
    std::vector<std::unique_ptr<BoundaryCondition>> conditions_;
    BoundaryEmitter emitter;
    bool periodic_[3] = {false, false, false};
};
