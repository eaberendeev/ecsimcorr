#include "StateOutputs.h"

#include <algorithm>

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

void RecoveryOutput::output(int timestep, const Diagnostics&) {
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
Fields3DOutput::Fields3DOutput(const Field3d& fieldEn, const Field3d& fieldBn, const Field3d& fieldE_external,
                               const Field3d& fieldBInit, const std::vector<std::string>& field_names,
                               const std::vector<int>& timesteps)
    : fieldEn_(fieldEn),
      fieldBn_(fieldBn),
      fieldE_external_(fieldE_external),
      fieldBInit_(fieldBInit),
      field_names_(field_names),
      timesteps_(timesteps) {
    ensure_directory(".//Fields");
    ensure_directory(".//Fields//Diag3D");
}

void Fields3DOutput::output(int timestep, const Diagnostics&) {
    if (std::find(timesteps_.begin(), timesteps_.end(), timestep) == timesteps_.end())
        return;

    std::string dir = ".//Fields//Diag3D//";
    auto write_field = [&](const Field3d& field, const std::string& name) {
        std::string fname = dir + name + "_" + std::to_string(timestep) + ".3d";
        write_field_to_file(fname, field);
    };

    for (const auto& fname : field_names_) {
        if (fname == "E") {
            Field3d fieldEFull = fieldEn_ + fieldE_external_;
            write_field(fieldEFull, "FieldE");
        } else if (fname == "B") {
            Field3d fieldBFull(fieldBn_.sizes(), fieldBn_.nd());
            fieldBFull.data() = fieldBn_.data() + fieldBInit_.data();
            write_field(fieldBFull, "FieldB");
        } else if (fname == "En") {
            write_field(fieldEn_, "FieldEn");
        } else if (fname == "Bn") {
            write_field(fieldBn_, "FieldBn");
        } else if (fname == "E_ex") {
            write_field(fieldE_external_, "FieldE_external");
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

void Full3DAllOutput::output(int timestep, const Diagnostics&) {
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
        dump(sp.densityOnGrid, dir + "Density_" + tstr + ".3d");
        dump(sp.currentOnGrid, dir + "Current_" + tstr + ".3d");
    }
}
