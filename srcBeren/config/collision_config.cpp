#include "config/collision_config.h"

#include <iostream>
#include <nlohmann/json.hpp>

#include "collisions/coulomb.h"
#include "collisions/legacy_collisions.h"
#include "collisions/neutral_collision.h"
#include "config/json_settings.h"

namespace {
std::vector<std::unique_ptr<CollisionOperator>> configured_channels(const nlohmann::json& config, double n0) {
    std::vector<std::unique_ptr<CollisionOperator>> processes;
    if (!config.contains("Collisions")) {
        return processes;
    }

    const auto& collisions = config["Collisions"];
    if (!collisions.is_array()) {
        std::cerr << "CollisionManager: 'Collisions' must be an array\n";
        return processes;
    }

    for (const auto& entry : collisions) {
        if (!entry.contains("type")) {
            std::cerr << "CollisionManager: entry missing 'type'\n";
            continue;
        }

        const std::string type = entry["type"].get<std::string>();

        if (type == "coulomb") {
            std::string sp1 = entry.value("species1", "");
            std::string sp2 = entry.value("species2", "");
            double lnL = entry.value("coulomb_log", 15.0);
            if (sp1.empty() || sp2.empty()) {
                std::cerr << "CollisionManager: coulomb entry needs species1, species2\n";
                continue;
            }
            processes.push_back(std::make_unique<CoulombCollisionOperator>(sp1, sp2, n0, lnL));
        } else if (type == "neutral_ionization") {
            NeutralCollisionConfig nc;
            nc.charged_species = entry.value("charged_species", "");
            nc.neutral_species = entry.value("neutral_species", "Neutrals");
            nc.electron_product = entry.value("electron_product", "Electrons");
            nc.ion_product = entry.value("ion_product", "Ions");
            nc.scheme = entry.value("scheme", "physical_only");
            nc.electron_ionization = entry.value("electron_ionization", true);
            nc.proton_ionization = entry.value("proton_ionization", true);
            nc.proton_charge_exchange = entry.value("proton_charge_exchange", true);
            if (nc.charged_species.empty()) {
                std::cerr << "CollisionManager: neutral_ionization needs charged_species\n";
                continue;
            }
            processes.push_back(std::make_unique<NeutralCollisionOperator>(nc, n0));
        } else {
            std::cerr << "CollisionManager: unknown collision type '" << type << "'\n";
        }
    }

    return processes;
}
}   // namespace

CollisionManager collisions_from_config(const nlohmann::json& config, double n0) {
    auto processes = configured_channels(config, n0);
    if (!processes.empty()) {
        std::cout << "Collisions:\n";
        for (const auto& process : processes) std::cout << "  " << process->info() << "\n";
    } else {
        const auto name = get_checked<std::string>(config, "Collider");
        if (name != "None") {
            processes.push_back(std::make_unique<LegacyCollisionOperator>(name, n0));
            std::cout << "Collider: " << name << "\n";
        }
    }
    return CollisionManager(std::move(processes));
}
