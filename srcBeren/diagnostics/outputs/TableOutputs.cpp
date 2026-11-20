#include "TableOutputs.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include "io/field_io.h"
#include "runtime/logger.h"
#include "runtime/profiler.h"

namespace {
// Decimal places of dt, so step * dt prints exactly without rounding noise.
int decimals_of(double dt) {
    int decimals = 0;
    for (double scaled = dt; decimals < 12 && std::abs(scaled - std::round(scaled)) > 1e-9 * std::max(1.0, scaled);
         scaled *= 10)
        ++decimals;
    return decimals;
}
}   // namespace

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
        ss << "Step Time ";
        for (const auto& key : column_order(diag)) {
            ss << key << " ";
        }
        header_written_ = true;
        ss << "\n";
    }
    // 12 digits: energy conservation values are ~1e-13 and need more than the
    // default 6 digits to be useful; counts below 10^12 stay exact.
    ss << std::setprecision(12) << timestep << " " << timestep * dt_ << " ";
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
ConsoleOutput::ConsoleOutput(double dt) : dt_(dt), time_decimals_(decimals_of(dt)) {
}

void ConsoleOutput::output(int timestep, const Diagnostics& diag, const DiagnosticMeasurements&) {
    RECORD_TIMER;
    const double t = timestep * dt_;
    std::ostringstream line;
    line << std::setw(6) << timestep << "  " << std::fixed << std::setprecision(time_decimals_) << std::setw(10) << t
         << std::scientific << std::setprecision(3) << "  Ef=" << std::setw(10) << diag.value_or_zero("energyFieldE")
         << " Bf=" << std::setw(10) << diag.value_or_zero("energyFieldB");

    for (const auto& name : diag.particleNames) {
        const int64_t count = static_cast<int64_t>(diag.value_or_zero(name + "Particles"));
        line << "  " << name[0] << ":" << std::setw(5) << count << " E=" << std::setw(10) << diag.value_or_zero(name);
    }

    // Symplectic PIC reports the signed Hamiltonian change since initialization;
    // ECSIM reports its existing discrete per-step work residual.
    const std::map<std::string, double>::const_iterator error = diag.values.find("energyError");
    line << "  dE=" << std::setw(9)
         << (error != diag.values.end() ? error->second : diag.value_or_zero("energyConserve"));
    std::cout << line.str() << "\n";
}
