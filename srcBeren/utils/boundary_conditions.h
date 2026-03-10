// boundary_condition.h
#pragma once

#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "World.h"
#include "containers.h"
#include "Particle.h"

class ParticlesArray;

class BoundaryEmitter {
   public:
    BoundaryEmitter(size_t num_species) : species_buffers_(num_species) {}

    // добавить частицу в текущий сорт
    void emit_current_species(const Particle& p) {
        current_species_buffer_.push_back(p);
    }

    // добавить частицу в другой сорт
    void emit_to_species(size_t species, const Particle& p) {
        species_buffers_[species].push_back(p);
    }

    const std::vector<Particle>& current_species_particles() const {
        return current_species_buffer_;
    }

    const std::vector<std::vector<Particle>>& other_species_particles() const {
        return species_buffers_;
    }

   private:
    std::vector<Particle> current_species_buffer_;
    std::vector<std::vector<Particle>> species_buffers_;
};

// -----------------------------------------------
// Базовый класс для всех дополнительных граничных условий
// -----------------------------------------------
class BoundaryCondition {
   public:
    virtual ~BoundaryCondition() = default;

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

    virtual void apply_to_particle(
        const Particle& /*p*/, BoundaryEmitter& /*emitter*/,
        const Domain& /*domain*/) const {   
        // по умолчанию ничего не делаем, частица будет удалена если пересекла границу
    }
    // Применение к полям (например, задать значение на границе)
    virtual void apply_to_fields(Field3d& /*fields*/, FieldType /*field*/,
                                 const Domain& /*domain*/) const {
        // по умолчанию ничего
    }
};

// -----------------------------------------------
// Конкретные реализации
// -----------------------------------------------

class PeriodicParticlesCondition : public BoundaryCondition {
   public:
    void apply_to_particle(const Particle& particle, BoundaryEmitter& emitter,
                           const Domain& domain) const;
};

class SecondEmissionCondition : public BoundaryCondition {
   public:
    SecondEmissionCondition(const std::string& axis, const std::string& bound) {
        // преобразуем ось в число
        if (axis == "x")
            dim_ = 0;
        else if (axis == "y")
            dim_ = 1;
        else if (axis == "z")
            dim_ = 2;
        else
            throw std::runtime_error("SecondEmissionCondition: invalid axis");
        isLower_ = (bound == "lower");
    }
    void apply_to_particle(const Particle& particle, BoundaryEmitter& emitter,
                           const Domain& domain) const;

   private:
    int dim_;
    bool isLower_;
};

// Условие Bphi: добавляет диагональ в оператор для магнитного поля на указанной
// границе
class BphiCondition : public BoundaryCondition {
   public:
    BphiCondition(const std::string& axis, const std::string& bound,
                  double Jz = 0.0, double radius = 0.0)
        : Jz_(Jz), radius_(radius) {
        // преобразуем ось в число
        if (axis == "x")
            dim_ = 0;
        else if (axis == "y")
            dim_ = 1;
        else if (axis == "z")
            dim_ = 2;
        else
            throw std::runtime_error("BphiCondition: invalid axis");
        isLower_ = (bound == "lower");
    }

    // void modify_curlB_stencil(int i, int j, int k, std::vector<Trip>& trips,
    //                           const Domain& domain) const override {
    //     const auto cell_size = domain.cell_size();
    //     const int vindx = domain.vind(i, j, k, 0);
    //     const int vindy = domain.vind(i, j, k, 1);
    //     const int vindz = domain.vind(i, j, k, 2);
    //     bool in_region = domain.is_inside_side_open(i, X, FieldType::ELECTRIC,
    //                                                 X, Side::LOWER) &&
    //                      domain.is_inside_side_open(j, Y, FieldType::ELECTRIC,
    //                                                 X, Side::LOWER);
    //     // (x)[i+1/2,j,k]
    //     // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
    //     if (k <= 1 && in_region) {
    //         double val = -1.0 / cell_size.z();
    //         trips.push_back(Trip(vindx, domain.vind(i, j, k, 1), val));
    //         if (k == 1)
    //             trips.push_back(Trip(vindx, domain.vind(i, j, k - 1, 1), -val));
    //     }

    //     in_region = domain.is_inside_side_open(i, X, FieldType::ELECTRIC, Y,
    //                                            Side::LOWER) &&
    //                 domain.is_inside_side_open(j, Y, FieldType::ELECTRIC, Y,
    //                                            Side::LOWER);
    //     // (y)[i,j+1/2,k]
    //     // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
    //     if (k <= 1 && in_region) {
    //         double val = 1.0 / cell_size.z();
    //         trips.push_back(Trip(vindy, domain.vind(i, j, k, 0), val));
    //         if (k == 1)
    //             trips.push_back(Trip(vindy, domain.vind(i, j, k - 1, 0), -val));
    //     }
    //     in_region = domain.is_inside_side_open(i, X, FieldType::ELECTRIC, Z,
    //                                            Side::LOWER) &&
    //                 domain.is_inside_side_open(j, Y, FieldType::ELECTRIC, Z,
    //                                            Side::LOWER);
    //     // (z)[i,j,k+1/2]
    //     if (k <= 1 && in_region) {
    //         // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
    //         double val = 1.0 / cell_size.x();
    //         trips.push_back(Trip(vindz, domain.vind(i, j, k, 1), val));
    //         trips.push_back(Trip(vindz, domain.vind(i - 1, j, k, 1), -val));
    //         // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
    //         val = -1.0 / cell_size.y();
    //         trips.push_back(Trip(vindz, domain.vind(i, j, k, 0), val));
    //         trips.push_back(Trip(vindz, domain.vind(i, j - 1, k, 0), -val));
    //     }
    // }
    // void modify_curlE_stencil(int i, int j, int k, std::vector<Trip>& trips,
    //                           const Domain& domain) const override {
    //     const auto cell_size = domain.cell_size();
    //     const int vindz = domain.vind(i, j, k, 2);

    //     bool in_region = domain.is_inside_side_open(i, X, FieldType::MAGNETIC,
    //                                                 X, Side::LOWER) &&
    //                      domain.is_inside_side_open(j, Y, FieldType::MAGNETIC,
    //                                                 X, Side::LOWER);
    //     if (k == 1 && in_region) {
    //         // (z)[i+1/2,j+1/2,k]
    //         // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
    //         double val = 1.0 / cell_size.x();
    //         trips.push_back(Trip(vindz, domain.vind(i + 1, j, k, 1), val));
    //         trips.push_back(Trip(vindz, domain.vind(i, j, k, 1), -val));
    //         // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
    //         val = -1.0 / cell_size.y();
    //         trips.push_back(Trip(vindz, domain.vind(i, j, k, 0), val));
    //         trips.push_back(Trip(vindz, domain.vind(i, j, k, 0), -val));
    //     }
    // }

    void apply_to_field(Field3d& field, FieldType field_t,
                        const Domain& domain) {
        if (field_t == FieldType::MAGNETIC) {
            if (dim_ == 2 && isLower_) {
                int startz = 0;
                int endz = 1;
                set_Bphi(field, Jz_, radius_, startz, endz, domain);
            }
        }
    }

   private:
    // Прототип локальной утилиты задания вихревого поля Bphi
    void set_Bphi(Field3d& fieldB, double Jz, double radius, double startz,
                  double endz, const Domain& domain) {
        const auto size_x = fieldB.sizes().x();
        const auto size_y = fieldB.sizes().y();
        const auto size_z = fieldB.sizes().z();
        const double dx = domain.cell_size().x();
        const double dy = domain.cell_size().y();
        const double dz = domain.cell_size().z();

        const int sz = std::max(0, static_cast<int>(std::llround(startz / dz)));
        const int ez = std::min(static_cast<int>(size_z),
                                static_cast<int>(std::llround(endz / dz)));

        const double center_x = 0.5 * (size_x - 3) * dx;
        const double center_y = 0.5 * (size_y - 3) * dy;
        const double eps_r = 1e-12;

        for (int k = sz; k < ez; ++k) {
            for (int i = 0; i < size_x; ++i) {
                for (int j = 0; j < size_y; ++j) {
                    // X-face node
                    double xx = i * dx - center_x - dx * GHOST_CELLS;
                    double yy = (j + 0.5) * dy - center_y - dy * GHOST_CELLS;
                    double rr = std::hypot(xx, yy);
                    if (rr < radius) {
                        fieldB(i, j, k, Dim::X) = -0.5 * Jz * (yy);
                    } else {
                        rr = std::max(rr, eps_r);
                        fieldB(i, j, k, Dim::X) =
                            -0.5 * radius * radius * Jz * yy / (rr * rr);
                    }

                    // Y-face node
                    yy = j * dy - center_y - dy * GHOST_CELLS;
                    xx = (i + 0.5) * dx - center_x - dx * GHOST_CELLS;
                    rr = std::hypot(xx, yy);
                    if (rr < radius) {
                        fieldB(i, j, k, Dim::Y) = 0.5 * Jz * (xx);
                    } else {
                        rr = std::max(rr, eps_r);
                        fieldB(i, j, k, Dim::Y) =
                            0.5 * radius * radius * Jz * xx / (rr * rr);
                    }
                }
            }
        }
    }
    int dim_;
    bool isLower_;
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
        Domain dom;
        dom.set_domain(sys_config);
        if (!sys_config.contains("Boundary_conditions"))
            return;
        const auto config = sys_config["Boundary_conditions"];
        if (!config.is_array()) {
            std::cerr << "BoundaryConditionHandler: expected array in JSON\n";
            return;
        }

        for (const auto& item : config) {
            std::string type = item["type"].get<std::string>();
            if (type == "bphi") {
                std::string axis = item.at("axis");
                std::string bound = item.at("bound");
                double Jz = item.at("Jz");
                double radius = item.at("radius");
                conditions_.push_back(
                    std::make_unique<BphiCondition>(axis, bound, Jz, radius));
            } else if (type == "second_emisson") {
                std::string axis = item.at("axis");
                std::string bound = item.at("bound");
                conditions_.push_back(
                    std::make_unique<SecondEmissionCondition>(axis, bound));
            } else {
                std::cerr << "Unknown boundary condition type: " << type
                          << std::endl;
            }
        }
        if (dom.get_bounds().is_periodic(0) || dom.get_bounds().is_periodic(1) ||
        dom.get_bounds().is_periodic(2)) {
            conditions_.push_back(std::make_unique<PeriodicParticlesCondition>());
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
        std::vector<std::unique_ptr<ParticlesArray>>& species,
        const Domain& domain) const;

    // Применяет все условия к полям
    void apply_to_fields(Field3d& fields, FieldType field_t,
                         const Domain& domain) const {
        for (const auto& cond : conditions_) {
            cond->apply_to_fields(fields, field_t, domain);
        }
    }

    // Если нужно знать, есть ли активные условия (например, чтобы не
    // вызывать apply... без нужды)
    bool empty() const { return conditions_.empty(); }

   private:
    std::vector<std::unique_ptr<BoundaryCondition>> conditions_;
};
