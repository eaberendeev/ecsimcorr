#include "EnergyMeasurements.h"

#include "ParticlesDiagnostic.h"
#include "timer.h"

namespace {
void collect_per_species_diagnostics(Diagnostics& diagnostic, const Species& species, double& kineticEnergy,
                                     double& kineticEnergyNew, double& totalLostEnergy) {
    RECORD_TIMER;
    static const Face all_faces[] = {Face::XMIN, Face::XMAX, Face::YMIN,    Face::YMAX,
                                     Face::ZMIN, Face::ZMAX, Face::CYLINDER};
    static const char* face_names[] = {"XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX", "CYLINDER"};
    static const PerFaceStats zero_stats;
    double totalInjectEnergy = 0;
    double totalEmitEnergy = 0;

    for (const Species::value_type& kv : species) {
        const ParticlesArray& sp = *kv.second;
        const ParticleMoments moments = sp.diagnostic_moments();
        diagnostic.addEnergy(sp.name(), moments.energy());
        diagnostic.addEnergy(sp.name() + "Particles", static_cast<double>(moments.count));
        diagnostic.addEnergy(sp.name() + "Inject", sp.diag.injection_energy);
        diagnostic.addEnergy(sp.name() + "Z", moments.kinetic.z());
        diagnostic.addEnergy(sp.name() + "XY", moments.kinetic.x() + moments.kinetic.y());
        kineticEnergy += sp.kineticEnergy;
        kineticEnergyNew += diagnostic.energy[sp.name()];
        totalInjectEnergy += sp.diag.injection_energy;

        diagnostic.addBoundary(sp.name() + "_InjectE", sp.diag.injection_energy);
        diagnostic.addBoundary(sp.name() + "_InjectN", static_cast<double>(sp.diag.injection_count));

        double sp_lost_energy = 0;
        for (int fi = 0; fi < 7; ++fi) {
            Face f = all_faces[fi];
            const std::unordered_map<Face, PerFaceStats>::const_iterator it = sp.diag.boundary.find(f);
            const PerFaceStats& s = (it != sp.diag.boundary.end()) ? it->second : zero_stats;
            std::string prefix = sp.name() + "_" + face_names[fi];
            diagnostic.addBoundary(prefix + "LostE", s.lost_energy);
            diagnostic.addBoundary(prefix + "LostN", static_cast<double>(s.lost_count));
            diagnostic.addBoundary(prefix + "ReflN", static_cast<double>(s.reflected_count));
            diagnostic.addBoundary(prefix + "EmitE", s.emitted_energy);
            diagnostic.addBoundary(prefix + "EmitN", static_cast<double>(s.emitted_count));
            sp_lost_energy += s.lost_energy;
            totalEmitEnergy += s.emitted_energy;
        }
        totalLostEnergy += sp_lost_energy;
    }
    diagnostic.addEnergy("totalInjectEnergy", totalInjectEnergy);
    diagnostic.addEnergy("totalEmitEnergy", totalEmitEnergy);
}

}   // namespace

void measure_field_energy_and_conservation(Diagnostics& diagnostic, const EnergyFields& fields,
                                           const IndexRange& irange, double dt, double kineticEnergy,
                                           double kineticEnergyNew, double totalLostEnergy, double totalInjectEnergy,
                                           double external_field_power, double dampingEnergy,
                                           double filtered_current_power) {
    RECORD_TIMER;
    diagnostic.addEnergy("totalLostEnergy", totalLostEnergy);
    diagnostic.addEnergy("energyDamping", dampingEnergy);
    // In code units, field energy is 1/2 sum_active |F|^2; irange excludes ghost copies.
    diagnostic.addEnergy("energyFieldE", 0.5 * dot_product_sum(fields.En, fields.En, irange));
    diagnostic.addEnergy("energyFieldB", 0.5 * dot_product_sum(fields.Bn, fields.Bn, irange));
    const Field3d total_magnetic = fields.Bn + fields.B_prescribed;
    diagnostic.addEnergy("energyFieldBFull", 0.5 * dot_product_sum(total_magnetic, total_magnetic, irange));
    double energyFieldEold = 0.5 * dot_product_sum(fields.E, fields.E, irange);
    double energyFieldBold = 0.5 * dot_product_sum(fields.B, fields.B, irange);

    double energyFieldDifference =
        diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] - energyFieldBold - energyFieldEold;
    // emitted (secondary) particles are an external source, like injection
    double totalEmitEnergy = diagnostic.energy["totalEmitEnergy"];

    // filtered_current_power_ = <E0, J_raw>, где J_raw — fields.Jp ДО apply_to_fields(CURRENT).
    // Ebar равно E0 только на кольце и нулю на остальном замаскированном слое,
    // а J_solver там равен нулю, поэтому это дискретная работа отфильтрованного
    // граничного тока dt*<Ebar, J_raw - J_solver>. Это НЕ физическая мощность
    // порта/Poynting и не параметр для подгонки баланса.
    const double rawResidual = kineticEnergyNew - kineticEnergy - totalInjectEnergy - totalEmitEnergy +
                               energyFieldDifference - dt * external_field_power + totalLostEnergy;
    const double discreteResidual = rawResidual - dt * filtered_current_power;
    diagnostic.addEnergy("energyConserve", std::abs(discreteResidual));
    diagnostic.addEnergy("energyResidualRaw", rawResidual);
    diagnostic.addEnergy("energyResidualDiscrete", discreteResidual);
}

void measure_ecsim_energy(Diagnostics& diagnostic, const Species& species, const EnergyFields& fields,
                          const Operator& response, const IndexRange& irange, double dt, double damping_energy,
                          double filtered_current_power) {
    RECORD_TIMER;
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    double external_field_power = 0;
    double particle_current_power = 0;
    double totalLostEnergy = 0;

    collect_per_species_diagnostics(diagnostic, species, kineticEnergy, kineticEnergyNew, totalLostEnergy);

    for (const Species::value_type& kv : species) {
        const ParticlesArray& sp = *kv.second;
        external_field_power += dot_product_sum(fields.E_external, sp.currentOnGrid, irange);
        particle_current_power += dot_product_sum(fields.Ep, sp.currentOnGrid, irange);
    }

    measure_field_energy_and_conservation(diagnostic, fields, irange, dt, kineticEnergy, kineticEnergyNew,
                                          totalLostEnergy, diagnostic.energy["totalInjectEnergy"], external_field_power,
                                          damping_energy, filtered_current_power);

    const Field3d response_current = fields.Jp + (1.0 / dt) * (response * (fields.E + fields.En));
    double response_current_power = dot_product_sum(fields.Ep, response_current, irange);
    double energyFieldDifference = diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] -
                                   0.5 * dot_product_sum(fields.E, fields.E, irange) -
                                   0.5 * dot_product_sum(fields.B, fields.B, irange);

    diagnostic.addEnergy("deltaKinetic", kineticEnergyNew - kineticEnergy);
    diagnostic.addEnergy("deltaField", energyFieldDifference);
    diagnostic.addEnergy("JePredWork", dt * response_current_power);
    diagnostic.addEnergy("JeRealWork", dt * particle_current_power);
    diagnostic.addEnergy("JeExWork", dt * external_field_power);
    diagnostic.addEnergy("energyRemovedJWork", dt * filtered_current_power);
}

void measure_corr_energy(Diagnostics& diagnostic, const Species& species, const EnergyFields& fields,
                         const IndexRange& irange, double dt, double damping_energy, double filtered_current_power) {
    RECORD_TIMER;
    const Field3d fieldEn12 = 0.5 * (fields.E + fields.En);
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    double external_field_power = 0;
    double particle_current_power = 0;
    double totalLostEnergy = 0;

    collect_per_species_diagnostics(diagnostic, species, kineticEnergy, kineticEnergyNew, totalLostEnergy);

    for (const Species::value_type& kv : species) {
        const ParticlesArray& sp = *kv.second;
        external_field_power += dot_product_sum(fields.E_external, sp.currentOnGrid, irange);
        particle_current_power += dot_product_sum(fieldEn12, sp.currentOnGrid, irange);
    }

    measure_field_energy_and_conservation(diagnostic, fields, irange, dt, kineticEnergy, kineticEnergyNew,
                                          totalLostEnergy, diagnostic.energy["totalInjectEnergy"], external_field_power,
                                          damping_energy, filtered_current_power);

    double energyFieldDifference = diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] -
                                   0.5 * dot_product_sum(fields.E, fields.E, irange) -
                                   0.5 * dot_product_sum(fields.B, fields.B, irange);

    diagnostic.addEnergy("deltaKinetic", kineticEnergyNew - kineticEnergy);
    diagnostic.addEnergy("deltaField", energyFieldDifference);
    diagnostic.addEnergy("JeCorrWork", dt * particle_current_power);
    diagnostic.addEnergy("JeExWork", dt * external_field_power);
    diagnostic.addEnergy("energyRemovedJWork", dt * filtered_current_power);
}
