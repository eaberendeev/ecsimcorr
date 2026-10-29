#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "grid/Domain.h"
#include "timer.h"

ParticleFateResult BphiCondition::classify_particle(const Particle& p, const ParticlesArray& particles,
                                                    const Domain& domain, Face hit_face) const {
    // Определяем, вышла ли частица через нужную Z-грань
    if (face_ != Face::ZMIN && face_ != Face::ZMAX)
        return {};   // на всякий случай

    const double eps = 1e-12;

    // Проверяем, что точка действительно вне области по Z (с учётом возможного
    // gap_)
    if (hit_face != Face::UNAVAILABLE ? hit_face != face_ : !domain.geom.is_outside_face(face_, p.coord, eps))
        return {};

    // Дополнительная проверка: точка должна быть внутри всех остальных граней
    // (т.е. не выходить также по X/Y/цилиндру)
    if (hit_face == Face::UNAVAILABLE && !domain.geom.contains_ignoring_face(face_, p.coord, eps))
        return {};

    double vz = p.velocity.z();

    bool inside_circle = is_inside_central_circle(p.coord, domain);
    double mass = particles.mass();
    bool reflect = should_electron_reflect(inside_circle, vz, mass);

    // Ветвление по сорту частиц (электроны/ионы) как в оригинале
    if (particles.name() == "Electrons" && reflect) {
        return {ParticleFate::Reflected, face_};
    }

    // Частица потеряна на Z-грани — классифицируем как Removed; накопление
    // энергии/диагностики выполняет обработчик.
    return {ParticleFate::Removed, face_};
}
