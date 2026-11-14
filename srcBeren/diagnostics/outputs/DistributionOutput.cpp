#include "DistributionOutput.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "io/field_io.h"
#include "runtime/profiler.h"
#include "utils/util.h"

DistributionOutput::DistributionOutput(const std::vector<std::string>& species_names,
                                       const DistributionRequest& request, int interval, std::string basename)
    : species_names_(species_names), request_(request), interval_(interval), basename_(std::move(basename)) {
    ensure_directory(".//Particles");
    for (const auto& name : species_names_) {
        ensure_directory(".//Particles//" + name);
        ensure_directory(".//Particles//" + name + "//Distribution");
    }
}

void DistributionOutput::request_measurements(int timestep, DiagnosticMeasurements& measured) const {
    if (timestep % interval_ == 0)
        for (const auto& name : species_names_) measured.request_distribution(name, basename_, request_);
}

void DistributionOutput::output(int timestep, const Diagnostics&, const DiagnosticMeasurements& measured) {
    RECORD_TIMER;
    if (timestep % interval_ != 0)
        return;
    const auto& [x, y] = request_.axes;
    for (const auto& name : species_names_) {
        const Distribution& f = measured.distribution(name, basename_);
        std::ostringstream filename;
        filename << ".//Particles//" << name << "//Distribution//" << basename_ << "_" << std::setw(6)
                 << std::setfill('0') << timestep << ".txt";
        std::ofstream out(filename.str());
        if (!out.is_open())
            throw std::runtime_error("Cannot open " + filename.str() + " for writing");
        out << std::setprecision(10);
        out << "# species " << name << " step " << timestep << " projection " << request_.projection << "\n";
        out << "# axis_a " << x.min << " " << x.max << " " << x.bins << "\n";
        out << "# axis_b " << y.min << " " << y.max << " " << y.bins << "\n";
        out << "# region_min " << request_.region_min.x() << " " << request_.region_min.y() << " "
            << request_.region_min.z() << " region_max " << request_.region_max.x() << " " << request_.region_max.y()
            << " " << request_.region_max.z() << "\n";
        out << "# a_center b_center f (sum of weights per bin / (da db))\n";
        const double da = (x.max - x.min) / x.bins, db = (y.max - y.min) / y.bins;
        for (int i = 0; i < x.bins; ++i)
            for (int j = 0; j < y.bins; ++j)
                out << x.min + (i + 0.5) * da << " " << y.min + (j + 0.5) * db << " "
                    << f.values[static_cast<size_t>(i) * y.bins + j] << "\n";
        close_output(out, filename.str());
    }
}
