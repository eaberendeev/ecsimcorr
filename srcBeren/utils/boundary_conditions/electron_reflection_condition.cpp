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

bool ElectronReflectionCondition::is_inside_central_circle(const Vector3R& pos, const Domain& domain) const {
    const double cx = 0.5 * domain.num_cells().x() * domain.cell_size().x();
    const double cy = 0.5 * domain.num_cells().y() * domain.cell_size().y();
    const double dx = pos.x() - cx;
    const double dy = pos.y() - cy;
    return (dx * dx + dy * dy) <= radius_ * radius_;
}

ParticleFateResult ElectronReflectionCondition::apply_to_particle(const Particle& p, ParticlesArray& particles,
                                                                  BoundaryEmitter& emitter, const Domain& domain) {
    if (face_ != Face::ZMIN && face_ != Face::ZMAX)
        return {};

    const double eps = 1e-12;
    if (!domain.geom.is_outside_face(face_, p.coord, eps))
        return {};
    if (!domain.geom.contains_ignoring_face(face_, p.coord, eps))
        return {};

    const double vz = p.velocity.z();
    const double mass = particles.mass();
    const double kinetic_z = 0.5 * mass * vz * vz;

    if (particles.name() == "Electrons" && is_inside_central_circle(p.coord, domain) &&
        kinetic_z <= energy_threshold_) {
        Particle new_p = p;
        new_p.velocity.z() = -vz;
        new_p.coord = domain.geom.reflect_from_face(face_, p.coord);
        emitter.emit_current_species(new_p);
        return {ParticleFate::Reflected, face_};
    }

    return {ParticleFate::Removed, face_};
}
