#include "Diagnostic.h"

#include <filesystem>
#include <stdexcept>
#include <string>

#include "logger.h"
#include "output_util.h"

namespace {

// Ensure an output directory exists. Existing directories are fine; any other
// failure is fatal because subsequent writes would silently go nowhere.
void ensure_directory(const std::string& path) {
    std::error_code ec;
    const std::filesystem::path fsPath(path);
    if (std::filesystem::is_directory(fsPath, ec))
        return;
    if (std::filesystem::exists(fsPath, ec)) {
        const std::string msg = "Output path exists but is not a directory: " + path;
        logger::error(msg);
        throw std::runtime_error(msg);
    }
    if (!std::filesystem::create_directory(fsPath, ec) || ec) {
        const std::string msg =
            "Failed to create output directory " + path + (ec ? ": " + ec.message() : std::string());
        logger::error(msg);
        throw std::runtime_error(msg);
    }
}

}   // namespace

void Diagnostics::make_folders() const {
    ensure_directory(".//Fields");
    ensure_directory(".//Fields//Diag2D");
    ensure_directory(".//Fields//Diag3D");
    ensure_directory(".//Recovery");
    ensure_directory(".//Recovery//Fields");
    ensure_directory(".//Recovery//Particles");
    ensure_directory(".//Particles");
    for (const auto& name : particleNames) {
        ensure_directory(".//Particles//" + name);
        ensure_directory(".//Particles//" + name + "//Diag2D");
        ensure_directory(".//Recovery//Particles//" + name);
    }
}
