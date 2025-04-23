#include "Diagnostic.h"

#include "output_util.h"

void Diagnostics::make_folders() const {
    create_directory(".//Anime");
    create_directory(".//Fields");
    create_directory(".//Fields//Diag3D");
    create_directory(".//Fields//Diag2D");
    create_directory(".//Fields//Diag1D");
    create_directory(".//Recovery");
    create_directory(".//Recovery//Fields");
    create_directory(".//Recovery//Particles");
    create_directory(".//Particles");
    for (const auto& name : particleNames) {
        create_directory(".//Particles//" + name);
        create_directory(".//Particles//" + name + "//Diag3D");
        create_directory(".//Particles//" + name + "//Diag2D");
        create_directory(".//Particles//" + name + "//Diag1D");
        create_directory(".//Recovery//Particles//" + name);
    }
    std::cerr << "Folders for output has been created: SUCCESS\n";
}

void Diagnostics::output_energy_spectrum(const EnergySpectrum& spectrum,
                                         int timestep,
                                         const std::string& output_dir) const {
    const int num_bins = spectrum.spectrum.size();
    const double bin_width = (spectrum.maxEnergy - spectrum.minEnergy) / num_bins;

    // Формирование имени файла
    std::ostringstream filename;
    filename << output_dir << "/energy_spectrum_" << std::setw(6)
             << std::setfill('0') << timestep << ".txt";

    // Запись данных в файл
    std::ofstream outfile(filename.str());
    if (!outfile.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename.str());
    }

    for (int i = 0; i < num_bins; ++i) {
        double energy_center = spectrum.minEnergy + (i + 0.5) * bin_width;
        outfile << energy_center << " "
                << spectrum.spectrum[i] << "\n";
    }
}
