#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "timer.h"

bool BphiCondition::apply_to_particle(const Particle& p, ParticlesArray& particles, BoundaryEmitter& emitter,
                                      const Domain& domain) {
    // Определяем, вышла ли частица через нужную Z-грань
    if (face_ != Face::ZMIN && face_ != Face::ZMAX)
        return false;   // на всякий случай

    const double eps = 1e-12;

    // Проверяем, что точка действительно вне области по Z (с учётом возможного
    // gap_)
    if (!domain.geom.is_outside_face(face_, p.coord, eps))
        return false;

    // Дополнительная проверка: точка должна быть внутри всех остальных граней
    // (т.е. не выходить также по X/Y/цилиндру)
    if (!domain.geom.contains_ignoring_face(face_, p.coord, eps))
        return false;

    double vz = p.velocity.z();

    bool inside_circle = is_inside_central_circle(p.coord, domain);
    double mass = particles.mass();
    bool reflect = should_electron_reflect(inside_circle, vz, mass);

    // Ветвление по сорту частиц (электроны/ионы) как в оригинале
    if (particles.name() == "Electrons" && reflect) {
        particles.diag.add_reflected(face_);
        Particle new_p = p;
        new_p.velocity.z() = -vz;
        new_p.coord = domain.geom.reflect_from_face(face_, p.coord);
        // Добавляем обратно в тот же сорт через эмиттер
        emitter.emit_current_species(new_p);
        return true;
    }

    // Частица потеряна на Z-грани — аккумулируем энергию
    double e_kin = get_energy_particle(p.velocity, mass, particles.mpw());
    particles.diag.add_loss(face_, e_kin);

    // Индексы ячейки для накопления тока
    auto idxs = domain.grid.get_field_node_index(p.coord, FieldType::CURRENT, Z);
    int i = idxs.x();
    int j = idxs.y();
    // Учитываем заряд и вес частицы
    double charge = particles.charge;
    double mpw = particles.mpw();
    double dJ = charge * mpw * vz;

    Jz_(i, j) += dJ;
    return true;
}

void BphiCondition::set_Bphi(Field3d& fieldB, const Array2D<double>& Jz, int k, const Domain& domain) {
    const auto size_x = fieldB.sizes().x();
    const auto size_y = fieldB.sizes().y();
    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();

    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            double Bx = 0;
            double By = 0;
            for (int i1 = 0; i1 < size_x; ++i1) {
                for (int j1 = 0; j1 < size_y; ++j1) {
                    double x = i * dx;
                    double y = (j + 0.5) * dy;
                    double x1 = i1 * dx;
                    double y1 = j1 * dy;

                    double rr2 = (x - x1) * (x - x1) + (y - y1) * (y - y1);
                    Bx += -Jz(i1, j1) * (y - y1) * dx * dy / (2 * M_PI * rr2);
                }
            }
            for (int i1 = 0; i1 < size_x; ++i1) {
                for (int j1 = 0; j1 < size_y; ++j1) {
                    double x = (i + 0.5) * dx;
                    double y = j * dy;
                    double x1 = i1 * dx;
                    double y1 = j1 * dy;
                    double rr2 = (x - x1) * (x - x1) + (y - y1) * (y - y1);
                    By += Jz(i1, j1) * (x - x1) * dx * dy / (2 * M_PI * rr2);
                }
            }
            fieldB(i, j, k, Axis::X) = Bx;
            fieldB(i, j, k, Axis::Y) = By;
        }
    }
}
