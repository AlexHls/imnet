#ifndef IMYANN_UI_MAIN_WINDOW_H
#define IMYANN_UI_MAIN_WINDOW_H

#include "app_state.h"
#include "ui_nuclide_chart.h"
#include <string>
#include <vector>

struct GLFWwindow;

namespace imyann {

/**
 * @class MainWindow
 * @brief Main application window and event loop
 *
 * Manages the GLFW window, imgui context, and main render loop.
 * Coordinates updates between the application state and UI components.
 */
class MainWindow {
public:
  /**
   * @brief Create and initialize the main window
   * @param width Initial window width (pixels)
   * @param height Initial window height (pixels)
   * @param title Window title
   */
  MainWindow(int width = 1200, int height = 800,
             const std::string &title = "imnet");

  ~MainWindow();

  // Deleted copy/move
  MainWindow(const MainWindow &) = delete;
  MainWindow &operator=(const MainWindow &) = delete;

  /**
   * @brief Set the application state
   */
  void set_app_state(AppState *state);

  /**
   * @brief Main event loop
   * @return true if window should close, false to continue
   */
  bool should_close() const;

  /**
   * @brief Process one frame
   * Handles input, updates UI, renders
   */
  void process_frame();

  /**
   * @brief Request window close
   */
  void request_close();

  /**
   * @brief Get window dimensions
   */
  int width() const { return width_; }
  int height() const { return height_; }

  /**
   * @brief Get raw GLFW window pointer (for advanced use)
   */
  GLFWwindow *native_window() { return window_; }

private:
  GLFWwindow *window_;
  int width_, height_;
  bool should_close_;
  AppState *app_state_;

  // UI layout state
  bool show_nuclide_chart_;
  bool show_integration_panel_;
  bool show_isotope_info_;
  bool show_trajectory_plot_;
  bool show_abundance_plot_;
  bool show_trajectory_editor_;
  bool show_settings_;
  bool trajectory_plot_log_x_;
  bool trajectory_plot_log_y_;
  bool abundance_plot_log_x_;
  bool abundance_plot_log_y_;
  bool abundance_plot_legend_;
  bool normalize_before_integrate_;
  bool show_all_abundances_;
  bool show_nonzero_only_;
  int run_mode_;
  int composition_preset_index_;
  double last_dedt_;
  std::string status_message_;
  char isotope_filter_[64];
  char abundance_plot_filter_[64];
  double abundance_plot_log_floor_;
  std::vector<int> abundance_plot_isotopes_;
  int isotope_info_index_;
  double isotope_info_rate_;
  std::vector<ReactionDiagnostic> isotope_reactions_;
  std::string isotope_info_error_;

  NuclideChart nuclide_chart_;
  bool chart_initialized_;

  bool show_load_abund_popup_;
  bool show_save_abund_popup_;
  bool show_save_state_popup_;
  bool show_load_traj_popup_;
  bool show_open_species_popup_;
#ifndef IMNET_USE_NUPPN
  bool show_open_rates_popup_;
#endif
  int trajectory_step_ui_;
  int trajectory_editor_selected_row_;
  char species_path_[512];
  char reaclib_path_[512];
  char partition_path_[512];
  char mass_path_[512];
  char weak_path_[512];
  char load_abund_path_[512];
  char save_abund_path_[512];
  char save_state_path_[512];
  char load_traj_path_[512];
  char save_traj_path_[512];
#ifdef IMNET_USE_NUPPN
  char nuppn_frame_input_[8192];
  char nuppn_physics_input_[65536];
  char nuppn_solver_input_[8192];
  bool nuppn_inputs_loaded_;
#endif

  struct TrajectoryIntegrationJob {
    bool active = false;
    bool failed = false;
    bool complete = false;
    size_t current_step = 0;
    size_t target_end_step = 0;
    size_t cache_start_step = 0;
    size_t cache_end_step = 0;
    std::vector<double> initial_xnuc;
    std::vector<double> xnuc;
    std::vector<TrajectoryStepCache> cache;
    std::string error;
  };
  TrajectoryIntegrationJob trajectory_job_;

  // File browser state
  std::string current_browse_dir_;
  std::vector<std::string> browse_files_;
  std::vector<std::string> browse_dirs_;
  char file_name_input_[256];

  /**
   * @brief Initialize imgui context
   */
  void init_imgui();

  /**
   * @brief Cleanup imgui context
   */
  void shutdown_imgui();

  /**
   * @brief Render main menu bar
   */
  void render_menu_bar();

  /**
   * @brief Render the main docking layout
   */
  void render_docking_layout();

  /**
   * @brief Render modal file operation dialogs
   */
  void render_file_dialogs();

  void sync_network_path_inputs();
  void reset_ui_after_network_reload();
#ifdef IMNET_USE_NUPPN
  bool load_nuppn_input_files();
  bool save_nuppn_input_files();
#endif

  /**
   * @brief Helper: refresh file browser for given directory
   */
  void refresh_file_browser(const std::string &dir);

  /**
   * @brief Render file browser component
   */
  bool render_file_browser(char *path_buffer, size_t buffer_size);

  /**
   * @brief Render selected isotope abundance editor
   */
  void render_selected_abundance_editor();

  /**
   * @brief Render single-step integration controls
   */
  void render_single_step_controls();

  /**
   * @brief Render editable composition table
   */
  void render_composition_editor();

  /**
   * @brief Render trajectory integration controls
   */
  void render_trajectory_controls();

  /**
   * @brief Start non-blocking trajectory integration
   */
  void start_trajectory_integration_job();

  /**
   * @brief Stop non-blocking trajectory integration
   */
  void cancel_trajectory_integration_job();

  /**
   * @brief Advance non-blocking trajectory integration
   */
  void process_trajectory_integration_job();

  /**
   * @brief Retarget non-blocking trajectory integration cache window
   */
  void retarget_trajectory_integration_job();

  /**
   * @brief Render loaded trajectory plots
   */
  void render_trajectory_plot_panel();

  /**
   * @brief Render cached abundance evolution plots
   */
  void render_abundance_plot_panel();

  /**
   * @brief Render editable trajectory table
   */
  void render_trajectory_editor_panel();

  /**
   * @brief Apply a simple initial composition preset
   */
  void apply_composition_preset(const std::vector<std::pair<std::string, double>>
                                    &composition);

  /**
   * @brief Select an isotope in app state and chart
   */
  void select_isotope(size_t index, bool append);

  /**
   * @brief Refresh cached diagnostics for selected isotope
   */
  void refresh_isotope_diagnostics(int index);

  void composition_changed();
};

} // namespace imyann

#endif // IMYANN_UI_MAIN_WINDOW_H
