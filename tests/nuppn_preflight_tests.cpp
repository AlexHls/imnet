#include "nuppn_preflight.h"
#include <csignal>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>

int main(int argc, char **argv) {
  if (argc == 3) {
    const auto mode = std::filesystem::path(argv[1]).filename().string();
    if (mode == "success" || mode == "hang" || mode == "failed_with_marker") {
      std::ofstream ready(argv[2]);
      ready << "imnet-nuppn-ready\n";
    }
    if (mode == "hang") std::this_thread::sleep_for(std::chrono::seconds(10));
    if (mode == "signal") std::raise(SIGTERM);
    std::cerr << "native diagnostic: " << mode << '\n';
    return mode == "failure" || mode == "failed_with_marker" ? 1 : 0;
  }
  try {
    const auto executable = std::filesystem::absolute(argv[0]);
    imyann::check_nuppn_startup("success", executable, std::chrono::seconds(5));
    for (const char *mode : {"stop_zero", "failure", "failed_with_marker", "signal", "hang"}) {
      bool rejected = false;
      try {
        imyann::check_nuppn_startup(mode, executable,
            mode == std::string("hang") ? std::chrono::milliseconds(100) : std::chrono::seconds(5));
      } catch (const std::exception &e) {
        const std::string error(e.what());
        rejected = error.find("NuPPN startup check failed") != std::string::npos;
        if (mode == std::string("hang")) rejected &= error.find("timed out") != std::string::npos;
        if (mode == std::string("stop_zero")) rejected &= error.find("native diagnostic: stop_zero") != std::string::npos;
      }
      if (!rejected) throw std::runtime_error(std::string("Probe failure was missed: ") + mode);
    }
    bool missing = false;
    try { imyann::check_nuppn_startup("success", executable / "missing"); }
    catch (const std::exception &) { missing = true; }
    if (!missing) throw std::runtime_error("Missing probe was accepted");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
