#include "ProbeOutput.h"

#include <iomanip>
#include <stdexcept>

#include "DiagnosticMeasurements.h"
#include "logger.h"
#include "output_util.h"
#include "util.h"

// ============================================================
//  ProbeOutput  — write prepared point samples to probes.txt
// ============================================================
ProbeOutput::ProbeOutput(const std::vector<std::string>& species_names, const std::vector<Vector3R>& points,
                         int interval)
    : species_names_(species_names), points_(points), interval_(interval), header_written_(false) {
    ensure_directory(".//Probes");
    fProbes_.open(".//Probes//probes.txt");
    if (!fProbes_.is_open()) {
        const std::string message = "Cannot open .//Probes//probes.txt for writing";
        logger::error(message);
        throw std::runtime_error(message);
    }
}

void ProbeOutput::request_measurements(int timestep, DiagnosticMeasurements& measured) const {
    if (timestep % interval_ == 0)
        for (const auto& point : points_) measured.request_probe(point, species_names_);
}

void ProbeOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) {
    if (timestep % interval_ != 0)
        return;
    if (points_.empty())
        return;

    if (!header_written_) {
        fProbes_ << "# Time x y z";
        fProbes_ << " Ex Ey Ez Bx By Bz";
        fProbes_ << " Enx Eny Enz Bnx Bny Bnz";
        fProbes_ << " E_ex_x E_ex_y E_ex_z B_init_x B_init_y B_init_z";
        for (const auto& name : species_names_) fProbes_ << " dens_" << name;
        fProbes_ << "\n";
        header_written_ = true;
    }

    for (const auto& pt : points_) {
        fProbes_ << timestep << " " << pt.x() << " " << pt.y() << " " << pt.z();

        const auto& sample = measured.probe(pt);
        fProbes_ << std::scientific << std::setprecision(6);
        for (const auto& value : sample.fields) fProbes_ << " " << value.x() << " " << value.y() << " " << value.z();
        for (const auto& name : species_names_) fProbes_ << " " << sample.densities.at(name);
        fProbes_ << "\n";
    }
    flush_output(fProbes_, "Probes/probes.txt");
}

void ProbeOutput::finalize() {
    if (fProbes_.is_open()) {
        close_output(fProbes_, "Probes/probes.txt");
    }
}
