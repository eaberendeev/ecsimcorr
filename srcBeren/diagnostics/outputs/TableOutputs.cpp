#include "TableOutputs.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "logger.h"
#include "output_util.h"
#include "timer.h"

ChargeContinuityOutput::ChargeContinuityOutput(double dt, int interval)
    : file_("continuity.txt"), dt_(dt), interval_(interval) {
    RECORD_TIMER;
    if (!file_.is_open())
        throw std::runtime_error("Cannot open continuity.txt for writing");
    file_ << "Time L2 Linf\n";
}

void ChargeContinuityOutput::output(int timestep, const Diagnostics& values, const DiagnosticMeasurements&) {
    RECORD_TIMER;
    if (!checks_charge_continuity(timestep))
        return;
    if (!values.charge_continuity)
        throw std::logic_error("Charge continuity endpoints were not measured for this step");
    const auto& residual = *values.charge_continuity;
    file_ << std::setprecision(12) << timestep * dt_ << ' ' << residual.l2 << ' ' << residual.linf << '\n';
    file_.flush();
}

// ============================================================
//  TableOutput  (shared by energy.txt and boundary.txt)
// ============================================================
TableOutput::TableOutput(const std::string& filename, double dt, int interval)
    : filename_(filename), dt_(dt), interval_(interval), header_written_(false) {
    RECORD_TIMER;
    file_.open(filename_);
    if (!file_.is_open()) {
        const std::string msg = "Cannot open " + filename_ + " for writing";
        logger::error(msg);
        throw std::runtime_error(msg);
    }
}

void TableOutput::output(int timestep, const Diagnostics& diag, const DiagnosticMeasurements&) {
    RECORD_TIMER;
    if (timestep % interval_ != 0)
        return;

    std::stringstream ss;
    if (!header_written_) {
        ss << "Time ";
        for (const auto& key : column_order(diag)) {
            ss << key << " ";
        }
        header_written_ = true;
        ss << "\n";
    }
    // 12 digits: energy conservation values are ~1e-13 and need more than the
    // default 6 digits to be useful.
    ss << std::setprecision(12) << timestep * dt_ << " ";
    for (const auto& key : column_order(diag)) {
        ss << diag.value_or_zero(key) << " ";
    }
    ss << "\n";
    file_ << ss.str();
    flush_output(file_, filename_);
}

void TableOutput::finalize() {
    RECORD_TIMER;
    if (file_.is_open()) {
        close_output(file_, filename_);
    }
}

// ============================================================
//  ConsoleOutput
// ============================================================
ConsoleOutput::ConsoleOutput(double dt) : dt_(dt) {
}

void ConsoleOutput::output(int timestep, const Diagnostics& diag, const DiagnosticMeasurements&) {
    RECORD_TIMER;
    const double t = timestep * dt_;
    std::ostringstream line;
    line << std::setw(6) << timestep << "  " << std::fixed << std::setprecision(2) << std::setw(10) << t
         << std::scientific << std::setprecision(3) << "  Ef=" << std::setw(10) << diag.value_or_zero("energyFieldE")
         << " Bf=" << std::setw(10) << diag.value_or_zero("energyFieldB");

    for (const auto& name : diag.particleNames) {
        const int count = static_cast<int>(diag.value_or_zero(name + "Particles"));
        line << "  " << name[0] << ":" << std::setw(5) << count << " E=" << std::setw(10) << diag.value_or_zero(name);
    }

    line << "  dE=" << std::setw(9) << diag.value_or_zero("energyConserve");
    std::cout << line.str() << "\n";
}
