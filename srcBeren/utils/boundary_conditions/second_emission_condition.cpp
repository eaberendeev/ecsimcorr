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

bool SecondEmissionCondition::apply_to_particle(const Particle& p, ParticlesArray& particles, BoundaryEmitter& emitter,
                                                const Domain& domain) {
    if (particles.name() != "Ions") {
        return false;
    }

    if (!domain.geom.is_outside_face(face_, p.coord) && is_outside_other_faces(p.coord, domain))
        return false;

    // Track ion loss energy
    double e_kin_ion = get_energy_particle(p.velocity, particles.mass(), particles.mpw());
    particles.diag.add_loss(face_, e_kin_ion);

    Particle p_new = p;
    p_new.velocity = gauss_.sample(pulse_gen_);

    if (face_ == Face::ZMIN) {
        p_new.coord.z() = -p_new.coord.z();
        p_new.velocity.z() = std::abs(p_new.velocity.z());
    } else if (face_ == Face::ZMAX) {
        p_new.coord.z() = 2 * domain.cell_size().z() * domain.num_cells().z() - p_new.coord.z();
        p_new.velocity.z() = std::abs(p_new.velocity.z());
    } else {
        return false;
    }

    // Track emitted electron energy
    double e_kin_el = get_energy_particle(p_new.velocity, 1.0, particles.mpw());
    particles.diag.add_emitted(face_, e_kin_el);

    emitter.emit_to_species("Electrons", p_new);
    return true;
}
