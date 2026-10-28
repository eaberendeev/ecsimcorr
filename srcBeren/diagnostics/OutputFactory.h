#pragma once

#include <memory>
#include <nlohmann/json.hpp>
#include <vector>

#include "IDiagnosticOutput.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"

// === Factory: creates IDiagnosticOutput vector from JSON config ===
class OutputFactory {
   public:
    static std::vector<std::unique_ptr<IDiagnosticOutput>> create(const nlohmann::json& diag_config,
                                                                  const Domain& domain, const Species& species,
                                                                  const nlohmann::json& system_config,
                                                                  const Field3d& fieldEn, const Field3d& fieldBn,
                                                                  const Field3d& fieldE_external,
                                                                  const Field3d& fieldBInit);

   private:
    static std::vector<std::unique_ptr<IDiagnosticOutput>> create_from_list(
        const nlohmann::json& outputs, const Domain& domain, const Species& species,
        const nlohmann::json& system_config, const Field3d& fieldEn, const Field3d& fieldBn,
        const Field3d& fieldE_external, const Field3d& fieldBInit);
};
