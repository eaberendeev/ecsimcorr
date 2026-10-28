#include "StateOutputs.h"

#include <algorithm>

#include "DiagnosticMeasurements.h"
#include "output_util.h"
#include "recovery.h"
#include "util.h"

// ============================================================
//  RecoveryOutput
// ============================================================
RecoveryOutput::RecoveryOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Species& species, int interval)
    : fieldEn_(fieldEn), fieldBn_(fieldBn), species_(species), interval_(interval) {
    ensure_directory(".//Recovery");
    ensure_directory(".//Recovery//Fields");
    ensure_directory(".//Recovery//Particles");
    for (const auto& kv : species_) {
        ensure_directory(".//Recovery//Particles//" + kv.second->name());
    }
}

void RecoveryOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements&) {
    if (timestep % interval_ != 0)
        return;
    for (const auto& kv : species_) {
        write_particles_to_recovery(kv.second.get(), timestep, interval_);
    }
    write_fields_to_recovery(fieldEn_, fieldBn_, timestep, interval_);
}

// ============================================================
//  Fields3DOutput
// ============================================================
Fields3DOutput::Fields3DOutput(const std::vector<OutputField>& fields, const std::vector<int>& timesteps)
    : fields_(fields), timesteps_(timesteps) {
    ensure_directory(".//Fields");
    ensure_directory(".//Fields//Diag3D");
}

void Fields3DOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements&) {
    if (std::find(timesteps_.begin(), timesteps_.end(), timestep) == timesteps_.end())
        return;
    for (const auto& selected : fields_) {
        const std::string filename =
            ".//Fields//Diag3D//" + std::string(selected.basename) + "_" + std::to_string(timestep) + ".3d";
        if (selected.extra) {
            const Field3d total = *selected.field + *selected.extra;
            write_field_to_file(filename, total);
        } else {
            write_field_to_file(filename, *selected.field);
        }
    }
}

// ============================================================
//  Full3DAllOutput  — full dump at given timesteps
// ============================================================
Full3DAllOutput::Full3DAllOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                                 const Field3d& fieldBInit, const Species& species, const std::vector<int>& timesteps)
    : fieldEn_(fieldEn),
      fieldBn_(fieldBn),
      fieldE_external_(fieldE_external),
      fieldBInit_(fieldBInit),
      species_(species),
      timesteps_(timesteps) {
    ensure_directory(".//Full3D");
    ensure_directory(".//Full3D//Fields");
    ensure_directory(".//Full3D//Particles");
    for (const auto& kv : species_) ensure_directory(".//Full3D//Particles//" + kv.second->name());
}

void Full3DAllOutput::request_measurements(int timestep, DiagnosticMeasurements& measured) const {
    if (std::find(timesteps_.begin(), timesteps_.end(), timestep) != timesteps_.end())
        for (const auto& [name, species] : species_) measured.request_density(name);
}

void Full3DAllOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) {
    if (std::find(timesteps_.begin(), timesteps_.end(), timestep) == timesteps_.end())
        return;

    std::string tstr = std::to_string(timestep);

    auto dump = [&](const Field3d& f, const std::string& path) { write_field_to_file(path, f); };

    dump(fieldEn_, ".//Full3D//Fields//En_" + tstr + ".3d");
    dump(fieldBn_, ".//Full3D//Fields//Bn_" + tstr + ".3d");
    dump(fieldE_external_, ".//Full3D//Fields//E_ex_" + tstr + ".3d");
    dump(fieldBInit_, ".//Full3D//Fields//B_init_" + tstr + ".3d");

    Field3d eFull(fieldEn_.sizes(), fieldEn_.nd());
    eFull.data() = fieldEn_.data() + fieldE_external_.data();
    dump(eFull, ".//Full3D//Fields//E_" + tstr + ".3d");

    Field3d bFull(fieldBn_.sizes(), fieldBn_.nd());
    bFull.data() = fieldBn_.data() + fieldBInit_.data();
    dump(bFull, ".//Full3D//Fields//B_" + tstr + ".3d");

    for (const auto& kv : species_) {
        const auto& sp = *kv.second;
        std::string dir = ".//Full3D//Particles//" + sp.name() + "//";
        dump(measured.density(kv.first), dir + "Density_" + tstr + ".3d");
        dump(sp.currentOnGrid, dir + "Current_" + tstr + ".3d");
    }
}
