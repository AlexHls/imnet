#include "file_io.h"
#include "app_state.h"
#include "network_wrapper.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition) {
  if (!condition) {
    throw std::runtime_error("app-layer self-check failed");
  }
}

bool near(double actual, double expected) {
  return std::abs(actual - expected) < 1e-12;
}

} // namespace

int main() {
  using namespace imyann;

  const std::vector<Species> species = {
      {"p", 1, 1, 0},
      {"he4", 2, 4, 2},
  };
  const std::vector<double> xnuc = {0.5, 0.5};
  require(near(calculate_abar(xnuc, species), 1.6));
  require(near(calculate_zbar(xnuc, species), 1.2));
  require(near(calculate_ye(xnuc, species), 0.75));

  const auto directory =
      std::filesystem::temp_directory_path() / "imnet_app_layer_tests";
  std::filesystem::create_directories(directory);
  const auto trajectory = directory / "trajectory.txt";
  {
    std::ofstream file(trajectory);
    file << "# time rho temp\n"
            "TUNIT T9\n"
            "0 1e9 5\n";
  }

  std::vector<double> times;
  std::vector<double> rhos;
  std::vector<double> temps;
  require(load_trajectory(trajectory.string(), times, rhos, temps));
  require(times.size() == 1 && near(temps[0], 5e9));

  const auto invalid = directory / "invalid.txt";
  {
    std::ofstream file(invalid);
    file << "not a trajectory\n";
  }
  require(!load_trajectory(invalid.string(), times, rhos, temps));
  require(times.size() == 1 && near(temps[0], 5e9));

  const auto abundances = directory / "abundances.txt";
  {
    std::ofstream file(abundances);
    file << "  # leading-space comment\n"
            "p 0.5 # inline comment\n"
            "he4 0.5\n";
  }
  std::vector<double> loaded;
  require(load_abundances(abundances.string(), {"p", "he4"}, loaded));
  require(loaded == xnuc);
  const auto saved_abundances = directory / "saved_abundances.txt";
  require(save_abundances(saved_abundances.string(), {"p", "he4"}, xnuc));
  loaded.clear();
  require(load_abundances(saved_abundances.string(), {"p", "he4"}, loaded));
  require(loaded == xnuc);
  require(!save_abundances((directory / "bad.txt").string(), {"p"}, xnuc));

  const auto data = std::filesystem::path(IMNET_TEST_DATA_DIR);
  const auto empty_weak_rates = directory / "empty_weak_rates.txt";
  std::ofstream{empty_weak_rates};
  Network reaclib_only(
      (data / "species.txt").string(), (data / "jinareaclib.dat").string(),
      (data / "part.txt").string(), (data / "mass.txt").string(),
      empty_weak_rates.string());
  std::vector<double> neutron_abundance(reaclib_only.num_species(), 0.0);
  neutron_abundance[0] = 1.0;
  const auto reaclib_fluxes = reaclib_only.get_reaction_fluxes(
      1e9, 5e9, neutron_abundance, 0.0, 0,
      std::numeric_limits<size_t>::max());
  require(std::any_of(reaclib_fluxes.begin(), reaclib_fluxes.end(),
                      [](const ReactionFlux &flux) {
                        return flux.source_index == 0 &&
                               flux.target_index == 1 && flux.weak;
                      }));

  AppState state;
  require(state.initialize_network(
      (data / "species.txt").string(), (data / "jinareaclib.dat").string(),
      (data / "part.txt").string(), (data / "mass.txt").string(),
      (data / "lmp_weak_rates.txt").string()));
  auto &settings = state.integration_settings();
  settings.xnuc.assign(state.num_species(), 0.0);
  settings.xnuc[1] = 0.5;
  settings.xnuc[2] = 0.5;
  settings.rho = 1e9;
  settings.temp = 5e9;
  settings.dt = 0.0;
  double flux_state_dedt = 0.0;
  std::string flux_state_error;
  require(state.integrate_single_step(false, flux_state_dedt,
                                      flux_state_error));
  require(state.has_reaction_flux_state());
  settings.temp = 4e9;
  require(near(state.reaction_flux_state_temp(), 5e9));
  const auto fluxes =
      state.get_reaction_fluxes(0.0, 0,
                                std::numeric_limits<size_t>::max());
  require(std::any_of(fluxes.begin(), fluxes.end(),
                      [](const ReactionFlux &flux) {
                        return flux.source_index == 1 &&
                               flux.target_index == 0 && flux.weak;
                      }));
  const auto regular_fluxes =
      state.get_reaction_fluxes(0.0, 0, 1, true, false);
  const auto weak_fluxes =
      state.get_reaction_fluxes(0.0, 0, 1, false, true);
  const auto mixed_fluxes =
      state.get_reaction_fluxes(0.0, 0, 1, true, true);
  const auto rate_fluxes =
      state.get_reaction_fluxes(0.0, 2, 1, true, false);
  require(regular_fluxes.size() == 1 && !regular_fluxes[0].weak);
  require(weak_fluxes.size() == 1 && weak_fluxes[0].weak);
  require(mixed_fluxes.size() == 2 &&
          std::any_of(mixed_fluxes.begin(), mixed_fluxes.end(),
                      [](const ReactionFlux &flux) { return flux.weak; }) &&
          std::any_of(mixed_fluxes.begin(), mixed_fluxes.end(),
                      [](const ReactionFlux &flux) { return !flux.weak; }));
  require(rate_fluxes.size() == 1 && rate_fluxes[0].strength > 0.0 &&
          rate_fluxes[0].strength == rate_fluxes[0].strength_rate);
  settings.rho = 1e9;
  settings.temp = 5e9;
  settings.dt = 1e-12;
  settings.final_time = 2e-12;
  settings.dt_max = 1e-12;
  settings.dt_factor = 1.0;
  settings.max_steps = 4;
  std::string integration_error;
  require(state.run_to_time(false, integration_error));
  require(state.trajectory_size() == 3);
  require(state.trajectory_cache().size() == 3);
  require(state.get_network()->last_substeps() == 2);
  require(std::abs(state.trajectory_cache().back().time - 2e-12) < 1e-24);
  require(state.trajectory_cache()[1].substeps == 1);
  require(state.trajectory_cache().back().substeps == 1);
  const auto saved_state = directory / "state.json";
  require(state.save_state_to_file(saved_state.string()));
  std::ifstream saved(saved_state);
  std::stringstream json;
  json << saved.rdbuf();
  require(json.str().find("\"format\": \"imnet-state-v1\"") !=
          std::string::npos);
  require(json.str().find("\"strength_dxdt\"") != std::string::npos);
  require(json.str().find("\"strength_rate\"") != std::string::npos);

  std::filesystem::remove_all(directory);
  return 0;
}
