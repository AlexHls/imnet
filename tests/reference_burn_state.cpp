// Exercise the state/cache and flux APIs used by the GUI, without a window.
#include "app_state.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

static void require(bool ok, const std::string &message) {
  if (!ok) throw std::runtime_error(message);
}

int main(int argc, char **argv) {
  try {
    require(argc == 5, "Usage: imnet_reference_state data initial trajectory output");
    imyann::AppState app;
    const std::filesystem::path data(argv[1]);
#ifdef IMNET_USE_NUPPN
    require(app.initialize_network(data.string(), "", "", "", ""), "Network initialization failed");
#else
    require(app.initialize_network((data / "species.txt").string(),
        (data / "jinareaclib.dat").string(), (data / "part.txt").string(),
        (data / "mass.txt").string(), (data / "lmp_weak_rates.txt").string()),
        "Network initialization failed");
#endif
    require(app.load_abundances_from_file(argv[2]), "Cannot load composition");
    require(app.load_trajectory_file(argv[3]), "Cannot load trajectory");
    std::string error;
    const bool success = app.run_trajectory(error);
    require(success, error);
    const auto &cache = app.trajectory_cache();
    require(cache.size() == app.trajectory_size(), "Incomplete cache");
    // Browse backwards: fluxes must follow the selected cached composition,
    // rather than retaining the final integration state.
    for (size_t row : {cache.size() - 1, cache.size() / 2, size_t(1), size_t(0)}) {
      require(app.set_trajectory_step(row), "Cannot select cached step");
      require(app.integration_settings().xnuc == cache[row].xnuc, "Stale GUI composition");
      for (int metric = 0; metric < 3; ++metric) {
        const auto actual = app.get_reaction_fluxes(0.0, metric, 1000, true, true);
        const auto expected = app.get_network()->get_reaction_fluxes(
            cache[row].rho, cache[row].temp, cache[row].xnuc, 0.0, metric, 1000, true, true);
        require(!actual.empty() && actual.size() == expected.size(), "Missing GUI fluxes");
        for (size_t i = 0; i < actual.size(); ++i) {
          const auto &a = actual[i];
          const auto &b = expected[i];
          require(a.source_index == b.source_index && a.target_index == b.target_index &&
                  a.weak == b.weak && a.strength == b.strength, "Stale GUI fluxes");
          const double selected = metric == 0 ? a.strength_dydt :
                                  metric == 1 ? a.strength_dxdt : a.strength_rate;
          require(std::isfinite(selected) && selected > 0 && selected == a.strength,
                  "Incorrect arrow metric");
          if (i) require(actual[i-1].strength >= selected, "Unsorted arrows");
        }
        for (bool weak : {false, true}) {
          const auto filtered = app.get_reaction_fluxes(0.0, metric, 1, !weak, weak);
          require(filtered.size() <= 1, "Arrow limit ignored");
          for (const auto &f : filtered) require(f.weak == weak, "Arrow type filter ignored");
        }
        require(app.get_reaction_fluxes(actual.front().strength * 2, metric, 1000,
                                       true, true).empty(), "Arrow threshold ignored");
      }
    }
    require(app.save_state_to_file(argv[4]), "Cannot export GUI state");
    std::cout << "Cached-step browsing and all three flux metrics passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
