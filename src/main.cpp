#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <stdexcept>
#include <vector>

#include "app_state.h"
#include "file_io.h"
#include "ui_main_window.h"

using namespace imyann;

namespace {

struct CliOptions {
#ifdef IMNET_USE_NUPPN
  std::string data_dir = IMNET_NUPPN_RUN_DIR;
#else
  std::string data_dir = "./data";
#endif
  bool headless_mode = false;
  std::string abundance_file;
  std::string trajectory_file;
  std::string output_file;
  std::string state_file;
  bool normalize = false;
  bool has_rho = false;
  bool has_temp = false;
  bool has_dt = false;
  bool has_final_time = false;
  bool has_dt_max = false;
  bool has_dt_factor = false;
  bool has_max_steps = false;
  double rho = 0.0;
  double temp = 0.0;
  double dt = 0.0;
  double final_time = 0.0;
  double dt_max = 0.0;
  double dt_factor = 0.0;
  int max_steps = 0;
};

bool parse_double_arg(const std::string &text, double &out) {
  try {
    size_t parsed = 0;
    out = std::stod(text, &parsed);
    return parsed == text.size() && std::isfinite(out);
  } catch (...) {
    return false;
  }
}

bool parse_int_arg(const std::string &text, int &out) {
  try {
    size_t parsed = 0;
    const long value = std::stol(text, &parsed);
    if (parsed != text.size() ||
        value < static_cast<long>(std::numeric_limits<int>::min()) ||
        value > static_cast<long>(std::numeric_limits<int>::max())) {
      return false;
    }
    out = static_cast<int>(value);
    return true;
  } catch (...) {
    return false;
  }
}

void write_headless_header(std::ostream &out,
                           const std::vector<std::string> &species_names) {
  out << "step,time,rho,temp,dt,dedt,substeps,status";
  for (const auto &name : species_names) {
    out << "," << name;
  }
  out << "\n";
}

void write_headless_row(std::ostream &out, size_t step, double time, double rho,
                        double temp, double dt, double dedt,
                        int substeps, const std::string &status,
                        const std::vector<double> &xnuc) {
  out << step << "," << time << "," << rho << "," << temp << "," << dt << ","
      << dedt << "," << substeps << "," << status;
  for (double value : xnuc) {
    out << "," << value;
  }
  out << "\n";
}

bool configure_initial_state(AppState &app_state, const CliOptions &options) {
  if (!options.abundance_file.empty() &&
      !app_state.load_abundances_from_file(options.abundance_file)) {
    return false;
  }

  auto &settings = app_state.integration_settings();
  if (options.has_rho) {
    settings.rho = options.rho;
  }
  if (options.has_temp) {
    settings.temp = options.temp;
  }
  if (options.has_dt) {
    settings.dt = options.dt;
  }
  if (options.has_final_time) {
    settings.final_time = options.final_time;
  }
  if (options.has_dt_max) {
    settings.dt_max = options.dt_max;
  } else if (options.has_final_time && options.has_dt) {
    if (settings.dt_max < options.dt) {
      settings.dt_max = options.dt;
    }
  }
  if (options.has_dt_factor) {
    settings.dt_factor = options.dt_factor;
  }
  if (options.has_max_steps) {
    settings.max_steps = options.max_steps;
  }
  if (options.normalize && !app_state.normalize_abundances()) {
    return false;
  }

  if (!options.trajectory_file.empty() &&
      !app_state.load_trajectory_file(options.trajectory_file)) return false;
  return true;
}

int run_headless(AppState &app_state, const CliOptions &options) {
  const bool has_history = !options.trajectory_file.empty() || options.has_final_time;
  double dedt = 0.0;
  std::string error;
  const bool success = !options.trajectory_file.empty()
      ? app_state.run_trajectory(error)
      : options.has_final_time ? app_state.run_to_time(false, error)
                               : app_state.integrate_single_step(false, dedt, error);
  if (!success) {
    std::cerr << error << std::endl;
    // Keep existing exports when no integration interval completed. Partial
    // histories remain exportable, with the failed row and a nonzero exit code.
    if (!has_history || app_state.trajectory_cache().size() <= 1) return 1;
  }

  std::vector<FileWrite> outputs;
  if (!options.output_file.empty()) {
    outputs.push_back({options.output_file, [&](std::ostream &out) {
      out << std::setprecision(17) << std::scientific;
      write_headless_header(out, app_state.get_species_names());
      if (has_history) {
        for (const auto &step : app_state.trajectory_cache()) {
          write_headless_row(out, step.index, step.time, step.rho, step.temp,
                             step.dt, step.dedt, step.substeps, step.status, step.xnuc);
        }
      } else {
        const auto &settings = app_state.integration_settings();
        write_headless_row(out, 0, app_state.current_time(), settings.rho, settings.temp,
                           settings.dt, dedt, app_state.get_network()->last_substeps(),
                           "ok", settings.xnuc);
      }
    }});
  }
  if (!options.state_file.empty()) {
    outputs.push_back({options.state_file, [&](std::ostream &out) {
      if (!app_state.write_state(out)) throw std::runtime_error("Cannot serialize state");
    }});
  }
  if (!write_files_transactional(outputs, error)) {
    std::cerr << "Cannot save results: " << error << std::endl;
    return 1;
  }
  return success ? 0 : 1;
}

} // namespace

void print_help() {
  std::cout
      << "imnet - Interactive " << Network::backend_name()
      << " Nuclear Network\n"
      << "Usage: imnet [options]\n"
      << "Options:\n"
      << "  --help                    Show this help message\n"
      << "  --version                 Show version information\n"
#ifdef IMNET_USE_NUPPN
      << "  --data-dir <dir>          NuPPN run directory (default: "
      << IMNET_NUPPN_RUN_DIR << ")\n"
#else
      << "  --data-dir <dir>          Specify data directory (default: "
         "./data)\n"
#endif
      << "  --headless                Run without GUI\n"
      << "  --abundances <file>       Load initial abundance file\n"
      << "  --trajectory <file>       Load trajectory file\n"
      << "  --output <file>           Write CSV result file\n"
      << "  --save-state <file>       Save portable JSON state\n"
      << "  --rho <value>             Set density\n"
      << "  --temp <value>            Set temperature\n"
      << "  --dt <value>              Set fixed or initial timestep\n"
      << "  --final-time <value>      Integrate to final time\n"
      << "  --dt-max <value>          Maximum final-time timestep\n"
      << "  --dt-factor <value>       Final-time timestep growth factor\n"
      << "  --max-steps <value>       Maximum final-time integration steps\n"
      << "  --normalize               Normalize initial abundances before "
         "running\n"
      << "  --run-trajectory          Legacy alias for --headless with "
         "--trajectory\n"
      << std::endl;
}

void print_version() {
  std::cout << "imnet version " << IMNET_VERSION << std::endl;
}

int main(int argc, char *argv[]) {
  CliOptions options;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    const bool requires_value =
        arg == "--data-dir" || arg == "--abundances" ||
        arg == "--trajectory" || arg == "--output" ||
        arg == "--save-state" || arg == "--rho" || arg == "--temp" ||
        arg == "--dt" || arg == "--final-time" || arg == "--dt-max" ||
        arg == "--dt-factor" || arg == "--max-steps";
    if (requires_value && i + 1 >= argc) {
      std::cerr << "Missing value for " << arg << std::endl;
      return 1;
    }

    if (arg == "--help") {
      print_help();
      return 0;
    } else if (arg == "--version") {
      print_version();
      return 0;
    } else if (arg == "--data-dir") {
      options.data_dir = argv[++i];
    } else if (arg == "--headless") {
      options.headless_mode = true;
    } else if (arg == "--abundances") {
      options.abundance_file = argv[++i];
    } else if (arg == "--trajectory") {
      options.trajectory_file = argv[++i];
    } else if (arg == "--output") {
      options.output_file = argv[++i];
    } else if (arg == "--save-state") {
      options.state_file = argv[++i];
    } else if (arg == "--rho") {
      options.has_rho = parse_double_arg(argv[++i], options.rho);
      if (!options.has_rho || options.rho <= 0.0) {
        std::cerr << "Invalid --rho value" << std::endl;
        return 1;
      }
    } else if (arg == "--temp") {
      options.has_temp = parse_double_arg(argv[++i], options.temp);
      if (!options.has_temp || options.temp <= 0.0) {
        std::cerr << "Invalid --temp value" << std::endl;
        return 1;
      }
    } else if (arg == "--dt") {
      options.has_dt = parse_double_arg(argv[++i], options.dt);
      if (!options.has_dt || options.dt < 0.0) {
        std::cerr << "Invalid --dt value" << std::endl;
        return 1;
      }
    } else if (arg == "--final-time") {
      options.has_final_time =
          parse_double_arg(argv[++i], options.final_time);
      if (!options.has_final_time || options.final_time <= 0.0) {
        std::cerr << "Invalid --final-time value" << std::endl;
        return 1;
      }
    } else if (arg == "--dt-max") {
      options.has_dt_max = parse_double_arg(argv[++i], options.dt_max);
      if (!options.has_dt_max || options.dt_max <= 0.0) {
        std::cerr << "Invalid --dt-max value" << std::endl;
        return 1;
      }
    } else if (arg == "--dt-factor") {
      options.has_dt_factor =
          parse_double_arg(argv[++i], options.dt_factor);
      if (!options.has_dt_factor || options.dt_factor <= 0.0) {
        std::cerr << "Invalid --dt-factor value" << std::endl;
        return 1;
      }
    } else if (arg == "--max-steps") {
      options.has_max_steps = parse_int_arg(argv[++i], options.max_steps);
      if (!options.has_max_steps || options.max_steps <= 0) {
        std::cerr << "Invalid --max-steps value" << std::endl;
        return 1;
      }
    } else if (arg == "--normalize") {
      options.normalize = true;
    } else if (arg == "--run-trajectory") {
      options.headless_mode = true;
    } else {
      std::cerr << "Unknown option: " << arg << std::endl;
      print_help();
      return 1;
    }
  }

  try {
    if (!options.headless_mode && (!options.output_file.empty() || !options.state_file.empty())) {
      throw std::runtime_error("--output and --save-state require --headless");
    }
    if (options.output_file == "-") {
      throw std::runtime_error("--output requires a file path; stdout output is not supported");
    }
    if (!options.trajectory_file.empty() && options.has_final_time) {
      throw std::runtime_error("--trajectory and --final-time are separate run modes");
    }
    if (options.headless_mode && (options.abundance_file.empty() ||
        (options.output_file.empty() && options.state_file.empty()))) {
      throw std::runtime_error("Headless mode requires --abundances and --output or --save-state");
    }
    const auto same_file = [](const std::string &a, const std::string &b) {
      if (a.empty() || b.empty()) return false;
      std::error_code error;
      return std::filesystem::equivalent(a, b, error) ||
             std::filesystem::weakly_canonical(a) == std::filesystem::weakly_canonical(b);
    };
    for (const auto &output : {options.output_file, options.state_file}) {
      if (same_file(output, options.abundance_file) || same_file(output, options.trajectory_file))
        throw std::runtime_error("Output file must not overwrite an input file");
    }
    if (same_file(options.output_file, options.state_file))
      throw std::runtime_error("CSV and JSON output paths must be different");
    std::cout << "=== imnet: Interactive " << Network::backend_name()
              << " Nuclear Network ===" << std::endl;
    std::cout << std::endl;

    // Create application state
    AppState app_state;

    std::cout << "Data directory: " << options.data_dir << std::endl;
    std::cout << std::endl;

#ifdef IMNET_USE_NUPPN
    const std::string species_file = options.data_dir;
    const std::string reaclib_file;
    const std::string partition_file;
    const std::string mass_file;
    const std::string weak_file;
#else
    // Construct file paths
    const std::filesystem::path data_dir = options.data_dir;
    const std::string species_file = (data_dir / "species.txt").string();
    const std::string reaclib_file = (data_dir / "jinareaclib.dat").string();
    const std::string partition_file = (data_dir / "part.txt").string();
    const std::string mass_file = (data_dir / "mass.txt").string();
    const std::string weak_file = (data_dir / "lmp_weak_rates.txt").string();
#endif

    // Initialize network
    std::cout << "Loading nuclear network data..." << std::endl;
    if (!app_state.initialize_network(species_file, reaclib_file,
                                      partition_file, mass_file, weak_file)) {
      std::cerr << "Failed to initialize network" << std::endl;
      return 1;
    }

    std::cout << "Network loaded successfully" << std::endl;
    std::cout << std::endl;

    if (!configure_initial_state(app_state, options)) return 1;
    if (options.headless_mode) {
      return run_headless(app_state, options);
    }

    // GUI mode: create window and event loop
    try {
      std::cout << "Launching imgui window..." << std::endl;
      std::cout << std::endl;

      MainWindow window(1600, 900, std::string("imnet - Interactive ") +
                                       Network::backend_name() +
                                       " Nuclear Network");
      window.set_app_state(&app_state, options.has_final_time ? 1 :
                            app_state.has_trajectory() ? 2 : 0);

      std::cout << "Main window ready. Use File > Exit or close the window to exit."
                << std::endl;

      // Main event loop
      while (!window.should_close()) {
        window.process_frame();
      }
    } catch (const std::exception &e) {
      std::cerr << "GUI initialization failed (headless environment?): "
                << e.what() << std::endl;
      return 1;
    }

    std::cout << "Application closing..." << std::endl;
    return 0;

  } catch (const std::exception &e) {
    std::cerr << "Fatal error: " << e.what() << std::endl;
    return 1;
  }
}
