#pragma once

#include <filesystem>

// Configuration and filesystem failures are independent of numerical output contents.
void verify_diagnostic_output_failures(const std::filesystem::path& root);
