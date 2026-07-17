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

void BoundaryConditionHandler::apply_to_particles(
    ParticlesArray& particles,
    [[maybe_unused]] std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species,
    const Domain& domain) {
    RECORD_TIMER;

    auto& data = particles.particlesData;
    int nx = data.size().x(), ny = data.size().y(), nz = data.size().z();

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                auto& cell = data(ix, iy, iz);

                auto it = std::remove_if(cell.begin(), cell.end(), [&](Particle& p) {
                    if (domain.contains(p.coord))
                        return false;

                    if (periodic_[0] || periodic_[1] || periodic_[2]) {
                        Particle p_new = p;
                        p_new.coord = wrap_periodic(p.coord, domain);

                        if (domain.contains(p_new.coord)) {
                            emitter.emit_current_species(p_new);
                            return true;
                        }
                    }

                    for (const auto& cond : conditions_) {
                        if (cond->apply_to_particle(p, particles, emitter, domain))
                            break;
                    }
                    return true;
                });

                cell.erase(it, cell.end());
            }
        }
    }
    int count = 0;
    // добавить новые частицы текущего сорта
    for (const auto& p : emitter.current_species_particles()) {
        particles.add_particle(p);
        count++;
    }
    emitter.clear_current_species_buffer();
}

void BoundaryConditionHandler::flush_species(
    std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) {
    // добавить частицы в другие сорта по имени
    const auto& other = emitter.other_species_particles();

    for (const auto& kv : other) {
        int count = 0;
        const std::string& name = kv.first;
        if (all_species.count(name) == 0) {
            throw std::runtime_error("Species " + name + " not found when adding from emitter");
        }
        auto& target_species = all_species[name];
        for (const auto& p : kv.second) {
            target_species->add_particle(p);
            count++;
        }
        // boundary particle counts written to boundary.txt
    }
    emitter.clear_other_species_buffers();
}

void BoundaryConditionHandler::load_from_json(const nlohmann::json& sys_config, const Domain& domain) {
    conditions_.clear();
    // Domain dom;
    // dom.set_domain(sys_config);
    if (!sys_config.contains("Boundary_conditions"))
        return;
    const auto& config = sys_config["Boundary_conditions"];
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
        std::string face_str;
        Face face;
        if (type != "open") {
            face_str = params.at("face");
            face = string_to_face(face_str);
        }

        if (face == Face::CYLINDER && !domain.geom.use_cylinder) {
            throw std::runtime_error(
                "Boundary condition uses face \"CYLINDER\", but \"CylinderDomain\" is not configured.\n"
                "Add to system_config:\n"
                "  \"CylinderDomain\": {\"radius\": <value>, \"center\": [<x>, <y>]}");
        }

        if (type == "bphi") {
            const double electron_threshold_energy = params.at("electron_threshold_energy_");
            const double radius = params.at("radius");
            const double gap = params.at("gap");
            conditions_.push_back(
                std::make_unique<BphiCondition>(face, domain, gap, radius, electron_threshold_energy));
        } else if (type == "second_emisson") {
            Vector3R mean = util::parse_double3(params.at("mean"));
            Vector3R temperature = util::parse_double3(params.at("sigma"));
            Vector3R sigma = convert_kev_to_sigma(temperature, 1.0);
            conditions_.push_back(std::make_unique<SecondEmissionCondition>(face, mean, sigma));
        } else if (type == "periodic") {
            conditions_.push_back(std::make_unique<PeriodicBoundaryCondition>(face));
            if (face == Face::XMIN || face == Face::XMAX) {
                periodic_[0] = true;
            } else if (face == Face::YMIN || face == Face::YMAX) {
                periodic_[1] = true;
            } else if (face == Face::ZMIN || face == Face::ZMAX) {
                periodic_[2] = true;
            }
        } else if (type == "open") {
            if (!params.is_array()) {
                throw std::runtime_error("Open boundary conditions acquire array, but obtained another type");
            }
            std::cout << "parsing params of OBC:\n" << params.dump(4) << std::endl;
            std::array<Face, OpenBoundaryConditionArray::maxBc> faces;
            const int bcCount = params.size();
            for (int i = 0; i < bcCount; ++i) {
                faces[i] = string_to_face(params[i].at("face"));
            }
            std::cout << "trying to create open BC" << std::endl;
            conditions_.push_back(std::make_unique<OpenBoundaryConditionArray>(bcCount, faces));
        } else if (type == "electron_reflection") {
            const double radius = params.value("radius", 5.0);
            const double energy_threshold = params.value("energy_threshold", 0.1) / SGS::MC2;
            conditions_.push_back(std::make_unique<ElectronReflectionCondition>(face, radius, energy_threshold));
        } else if (type == "er0") {
            const double inner_radius = params.value("inner_radius", 5.0);
            const double width = params.value("width", 2.0);
            const double potential_drop = params.value("potential_drop", 0.1) / SGS::MC2;
            conditions_.push_back(std::make_unique<Er0Condition>(face, inner_radius, width, potential_drop));
        } else {
            std::cerr << "Unknown boundary condition type: " << type << std::endl;
        }
    }
}
