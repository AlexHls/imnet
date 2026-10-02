// Exercise the state/cache and flux APIs used by the GUI, without a window.
#include "app_state.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

static void require(bool ok, const std::string &message) {
  if (!ok) throw std::runtime_error(message);
}

int main(int argc, char **argv) {
  try {
    // Endothermic energy is valid; NaN/Inf and invalid history fields are not.
    imyann::IntegrationStepSnapshot valid{0, 0, 1e4, 2e8, 0, -1, 0, {1}};
    imyann::validate_integration_snapshot(valid, 1);
    for (int field = 0; field < 8; ++field) {
      auto invalid = valid;
      switch (field) {
        case 0: invalid.dedt = std::numeric_limits<double>::quiet_NaN(); break;
        case 1: invalid.dedt = std::numeric_limits<double>::infinity(); break;
        case 2: invalid.xnuc[0] = std::numeric_limits<double>::quiet_NaN(); break;
        case 3: invalid.xnuc[0] = -1; break;
        case 4: invalid.rho = 0; break;
        case 5: invalid.dt = -1; break;
        case 6: invalid.substeps = -1; break;
        case 7: invalid.xnuc.clear(); break;
      }
      bool rejected = false;
      try { imyann::validate_integration_snapshot(invalid, 1); }
      catch (const std::exception &) { rejected = true; }
      require(rejected, "Invalid solver history accepted");
    }
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
    const auto saved_size = cache.size();
    const auto saved_final = cache.back().xnuc;
    const auto saved_x = app.integration_settings().xnuc;
    const auto saved_time = app.current_time();
    app.integration_settings().max_steps = std::numeric_limits<int>::max() - 1;
    require(!app.run_to_time(true, error) && error.find("256 MiB") != std::string::npos,
            "Oversized history was accepted");
    require(cache.size() == saved_size && cache.back().xnuc == saved_final &&
            app.integration_settings().xnuc == saved_x && app.current_time() == saved_time,
            "Rejected run destroyed existing results");
    auto direct_x = saved_x;
    bool rejected = false;
    try {
      app.get_network()->integrate_to_time(1e4, 2e8, direct_x, 1, 1, 1, 1,
                                           std::numeric_limits<int>::max() - 1);
    } catch (const std::exception &e) {
      rejected = std::string(e.what()).find("256 MiB") != std::string::npos;
    }
    require(rejected && direct_x == saved_x, "Direct backend bypassed history budget");
    require(app.save_state_to_file(argv[4]), "Cannot export GUI state");
    // Both native histories must still capture successful and partial runs.
    auto &settings = app.integration_settings();
    settings.dt = settings.dt_max = 1e-6;
    settings.final_time = 2e-6;
    settings.dt_factor = 1;
    settings.max_steps = 2;
    require(app.run_to_time(true, error), "Small final-time run failed: " + error);
    require(app.trajectory_cache().size() == 3 && app.current_time() == 2e-6,
            "Final-time history incomplete");
    settings.max_steps = 1;
    require(!app.run_to_time(false, error), "Step limit was ignored");
    require(app.trajectory_cache().size() == 2 && app.current_time() == 1e-6 &&
            !app.trajectory_cache().back().success, "Partial results were lost");
    std::cout << "Cached-step browsing and all three flux metrics passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
