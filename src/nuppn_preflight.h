#ifndef IMNET_NUPPN_PREFLIGHT_H
#define IMNET_NUPPN_PREFLIGHT_H

#include <chrono>
#include <filesystem>

namespace imyann {
std::filesystem::path nuppn_probe_path();
void check_nuppn_startup(const std::filesystem::path &run_dir,
                         const std::filesystem::path &probe,
                         std::chrono::milliseconds timeout = std::chrono::seconds(120));
} // namespace imyann
#endif
