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
