#include "EnergySpectrumOutput.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "checked_output.h"
#include "logger.h"
#include "util.h"

// ============================================================
//  EnergySpectrumOutput
// ============================================================
EnergySpectrumOutput::EnergySpectrumOutput(const Species& species, const std::vector<std::string>& species_names,
                                           int interval)
    : species_(species), species_names_(species_names), interval_(interval) {
    ensure_directory(".//Particles");
    for (const auto& name : species_names_) ensure_directory(".//Particles//" + name);
}

void EnergySpectrumOutput::output(int timestep, const Diagnostics&) {
    if (timestep % interval_ != 0)
        return;

    for (const auto& name : species_names_) {
        const auto it = species_.find(name);
        if (it == species_.end()) {
            const std::string msg = "EnergySpectrumOutput: species '" + name + "' not found";
            logger::error(msg);
            throw std::runtime_error(msg);
        }

        const auto& sp = *it->second;
        EnergySpectrum spectrum = sp.calculate_energy_spectrum();

        std::ostringstream filename;
        filename << ".//Particles//" << name << "//energy_spectrum_" << std::setw(6) << std::setfill('0') << timestep
                 << ".txt";

        std::ofstream out(filename.str());
        if (!out.is_open()) {
            const std::string msg = "Cannot open " + filename.str() + " for writing";
            logger::error(msg);
            throw std::runtime_error(msg);
        }

        const int nbins = spectrum.spectrum.size();
        const double bin_width = (nbins > 0) ? (spectrum.maxEnergy - spectrum.minEnergy) / nbins : 0.0;
        for (int i = 0; i < nbins; ++i) {
            const double ec = spectrum.minEnergy + (i + 0.5) * bin_width;
            out << ec << " " << spectrum.spectrum[i] << "\n";
        }
        close_output(out, filename.str());
    }
}
