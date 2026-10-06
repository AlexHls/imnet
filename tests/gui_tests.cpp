// Run manually in a desktop session: imnet_gui_tests network initial
// trajectory.
#include "GLFW/glfw3.h"
#include "imgui_internal.h"
#include "ui_main_window.h"
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

static void require(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}

namespace imyann {
class MainWindowTest {
public:
  static void run(AppState &state) {
    const auto initial_composition = state.integration_settings().xnuc;
    MainWindow window(1600, 900, "imnet GUI regression");
    glfwHideWindow(window.native_window());
    window.set_app_state(&state, 2);
    ImGui::GetIO().IniFilename = nullptr;
    auto finish = [&] {
      for (int i = 0; i < 10000 && window.trajectory_job_.active; ++i)
        window.process_trajectory_integration_job();
      require(!window.trajectory_job_.failed, "Trajectory scheduler failed");
      require(window.trajectory_job_.complete,
              "Trajectory scheduler did not complete");
    };
    state.view_settings().max_cached_trajectory_steps = 10;
    window.start_trajectory_integration_job();
    finish();
    require(state.trajectory_cache().size() == 10,
            "Initial cache window incomplete");
    window.trajectory_step_ui_ = 30;
    window.retarget_trajectory_integration_job();
    finish();
    require(state.get_trajectory_cache_step(30),
            "Moving cache omitted selected step");
    window.trajectory_step_ui_ = 3;
    window.retarget_trajectory_integration_job();
    finish();
    require(state.get_trajectory_cache_step(3),
            "Backwards cache window incomplete");
    state.view_settings().max_cached_trajectory_steps = 40;
    window.retarget_trajectory_integration_job();
    finish();
    require(state.trajectory_cache().size() == 40,
            "Expanded cache window incomplete");
    window.cancel_trajectory_integration_job();
    window.trajectory_step_ui_ = 60;
    window.retarget_trajectory_integration_job();
    require(!window.trajectory_job_.active &&
                window.trajectory_job_.initial_xnuc.empty(),
            "Cancelled trajectory restarted while browsing");
    require(state.set_trajectory_step(0), "Cannot return to initial state");
    window.trajectory_step_ui_ = 0;
    window.start_trajectory_integration_job();
    state.view_settings().max_cached_trajectory_steps =
        static_cast<int>(state.trajectory_size());
    window.retarget_trajectory_integration_job();
    finish();
    require(state.trajectory_cache().size() == state.trajectory_size(),
            "Full trajectory incomplete");
    const auto expected_final = state.trajectory_cache().back().xnuc;
    require(state.set_trajectory_step(0),
            "Cannot rewind for headless comparison");
    std::string error;
    require(state.run_trajectory(error), "Headless comparison failed");
    require(state.trajectory_cache().back().xnuc == expected_final,
            "Incremental GUI integration differs from headless integration");
    state.view_settings().show_fluxes = true;
    require(state.set_trajectory_step(state.trajectory_size() / 2),
            "Cannot browse midpoint");
    // Exercise the actual render path, including layout, flux arrows and plots.
    for (int i = 0; i < 3; ++i)
      window.process_frame();
    state.view_settings().flux_threshold = std::numeric_limits<double>::quiet_NaN();
    window.process_frame();
    require(!window.nuclide_chart_.flux_error().empty(),
            "Chart hid a flux failure as an empty arrow list");
    state.view_settings().flux_threshold = 0;
    window.process_frame();
    require(window.nuclide_chart_.flux_error().empty(),
            "Chart flux error did not clear after recovery");
    window.show_trajectory_plot_ = true;
    window.trajectory_plot_log_x_ = true;
    window.show_trajectory_editor_ = true;
    for (int i = 0; i < 3; ++i)
      window.process_frame();
    window.refresh_isotope_diagnostics(0);
    state.integration_settings().rho = -1;
    window.refresh_isotope_diagnostics(0);
    require(!window.isotope_info_error_.empty(),
            "Diagnostic error was hidden as zero");

    // Final-time GUI runs yield between steps and match the headless integrator.
    auto &settings = state.integration_settings();
    settings.xnuc = initial_composition;
    settings.rho = 1e4;
    settings.temp = 2e8;
    settings.dt = 0.001;
    settings.dt_max = 0.008;
    settings.dt_factor = 2;
    settings.final_time = 0.021;
    settings.max_steps = 32;
    require(state.run_to_time(false, error), "Headless final-time comparison failed");
    const auto expected_history = state.trajectory_cache();
    settings.xnuc = initial_composition;
    window.normalize_before_integrate_ = false;
    window.run_mode_ = 1;
    state.view_settings().max_cached_trajectory_steps = 2;
    window.start_final_time_integration_job();
    require(window.trajectory_job_.active && state.trajectory_cache().size() == 1,
            "Starting a final-time run performed solver work synchronously");
    finish();
    require(state.trajectory_cache().size() == expected_history.size(),
            "Final-time history was limited by the trajectory cache window");
    for (size_t row = 0; row < expected_history.size(); ++row) {
      const auto &actual = state.trajectory_cache()[row];
      const auto &expected = expected_history[row];
      require(std::abs(actual.time - expected.time) < 1e-14, "Final-time grid differs");
      for (size_t i = 0; i < actual.xnuc.size(); ++i)
        require(std::abs(actual.xnuc[i] - expected.xnuc[i]) < 1e-12 + 1e-10 * expected.xnuc[i],
                "GUI final-time composition differs from headless");
    }
    require(settings.dt == 0.001 && state.current_time() == settings.final_time,
            "Completed final-time run changed initial dt or selected the wrong row");
    settings.dt = settings.dt_max = 1.0 / 1024;
    settings.dt_factor = 1;
    settings.final_time = 0.125;
    settings.max_steps = 128;
    window.start_final_time_integration_job();
    window.process_frame();
    require(window.trajectory_job_.active && window.trajectory_job_.current_step > 0 &&
            window.trajectory_job_.current_step <= 32,
            "Final-time run did not respect the frame batch limit");
    const auto partial = state.trajectory_cache();
    const auto partial_time = state.current_time();
    window.cancel_trajectory_integration_job();
    window.process_frame();
    require(!window.trajectory_job_.active && state.trajectory_cache().size() == partial.size() &&
            state.current_time() == partial_time && settings.xnuc == partial.back().xnuc,
            "Cancelled final-time run lost results or continued running");
    settings.max_steps = 1;
    window.start_final_time_integration_job();
    require(!window.trajectory_job_.active && state.current_time() == partial_time &&
            state.trajectory_cache().size() == partial.size(),
            "Invalid final-time plan destroyed previous results");

  }
};
} // namespace imyann

int main(int argc, char **argv) {
  try {
    require(argc == 4, "Usage: imnet_gui_tests network initial trajectory");
    const std::filesystem::path data(argv[1]);
    imyann::AppState state;
#ifdef IMNET_USE_NUPPN
    require(state.initialize_network(data.string(), "", "", "", ""),
            "Network init failed");
#else
    require(state.initialize_network((data / "species.txt").string(),
                                     (data / "jinareaclib.dat").string(),
                                     (data / "part.txt").string(),
                                     (data / "mass.txt").string(),
                                     (data / "lmp_weak_rates.txt").string()),
            "Network init failed");
#endif
    require(state.load_abundances_from_file(argv[2]),
            "Cannot load composition");
    require(state.load_trajectory_file(argv[3]), "Cannot load trajectory");
    require(state.trajectory_size() > 60,
            "GUI test needs at least 61 trajectory rows");
    imyann::MainWindowTest::run(state);
    std::cout << "GUI scheduler, cancellation, cache browsing, and diagnostics "
                 "passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
