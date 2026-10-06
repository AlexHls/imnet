// Separate executable: Fortran STOP, runtime errors and signals stay in this process.
#include "nuppn_c_api.h"
#include <filesystem>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
  try {
    if (argc != 3) return 2;
    const auto marker = std::filesystem::absolute(argv[2]);
    std::filesystem::current_path(argv[1]);
    const int status = nuppn_init();
    if (status != 0 || nuppn_num_species() <= 0) {
      std::cerr << "NuPPN initialization failed with status " << status << '\n';
      return 1;
    }
    // Exit status alone is insufficient: some native STOP paths return zero.
    std::ofstream ready(marker);
    ready << "imnet-nuppn-ready\n";
    ready.close();
    return ready ? 0 : 1;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
