#include "ui_main_window.h"
#include "ui_theme.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "GLFW/glfw3.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
#include "imgui.h"
#include "implot.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>

namespace imyann {

// ============================================================================
// GLFW Callbacks (static, forward to instance)
// ============================================================================

static MainWindow *g_main_window = nullptr;

static TrajectoryStepCache make_trajectory_step_cache(
    size_t index, const std::vector<double> &times,
    const std::vector<double> &rhos, const std::vector<double> &temps,
    const std::vector<double> &xnuc, double dedt, int substeps, bool success,
    const std::string &status, const std::string &error) {
  TrajectoryStepCache step;
  step.index = index;
  step.time = times[index];
  step.rho = rhos[index];
  step.temp = temps[index];
  step.dt =
      (index + 1 < times.size()) ? std::max(times[index + 1] - times[index],
                                            0.0)
                                 : 0.0;
  step.dedt = dedt;
  step.substeps = substeps;
  step.success = success;
  step.status = status;
  step.error = error;
  step.xnuc = xnuc;
  return step;
}

static void upsert_trajectory_cache_step(
    std::vector<TrajectoryStepCache> &cache, const TrajectoryStepCache &step) {
  auto it = std::find_if(cache.begin(), cache.end(),
                         [&step](const TrajectoryStepCache &cached) {
                           return cached.index == step.index;
                         });
  if (it == cache.end()) {
    cache.push_back(step);
  } else {
    *it = step;
  }
  std::sort(cache.begin(), cache.end(),
            [](const TrajectoryStepCache &a, const TrajectoryStepCache &b) {
              return a.index < b.index;
            });
}

static void prune_trajectory_cache(std::vector<TrajectoryStepCache> &cache,
                                   size_t start, size_t end) {
  cache.erase(std::remove_if(cache.begin(), cache.end(),
                             [start, end](const TrajectoryStepCache &step) {
                               return step.index < start || step.index > end;
                             }),
              cache.end());
}

static void glfw_error_callback(int error, const char *description) {
  std::cerr << "GLFW Error " << error << ": " << description << std::endl;
}

static void glfw_key_callback(GLFWwindow *window, int key, int scancode,
                              int action, int mods) {
  (void)window;
  (void)scancode;
  (void)mods;
  if (g_main_window && key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    g_main_window->request_close();
  }
}

// ============================================================================
// MainWindow Implementation
// ============================================================================

MainWindow::MainWindow(int width, int height, const std::string &title)
    : window_(nullptr), width_(width), height_(height), should_close_(false),
      app_state_(nullptr), show_nuclide_chart_(true),
      show_integration_panel_(true), show_isotope_info_(true),
      show_trajectory_plot_(false), show_abundance_plot_(false),
      show_trajectory_editor_(false), show_settings_(false),
      trajectory_plot_log_x_(false), trajectory_plot_log_y_(false),
      abundance_plot_log_x_(false), abundance_plot_log_y_(false),
      abundance_plot_legend_(true), normalize_before_integrate_(false),
      show_all_abundances_(true), show_nonzero_only_(false), run_mode_(0),
      composition_preset_index_(0), last_dedt_(0.0),
      status_message_("Ready"), isotope_info_index_(-1),
      isotope_info_rate_(0.0), chart_initialized_(false),
      show_load_abund_popup_(false), show_save_abund_popup_(false),
      show_save_state_popup_(false), show_load_traj_popup_(false),
      show_open_species_popup_(false),
#ifndef IMNET_USE_NUPPN
      show_open_rates_popup_(false), trajectory_step_ui_(0),
#else
      trajectory_step_ui_(0),
#endif
      trajectory_editor_selected_row_(0)
#ifdef IMNET_USE_NUPPN
      ,
      nuppn_inputs_loaded_(false)
#endif
{

  std::strncpy(species_path_, "./data/species.txt", sizeof(species_path_));
  std::strncpy(reaclib_path_, "./data/jinareaclib.dat", sizeof(reaclib_path_));
  std::strncpy(partition_path_, "./data/part.txt", sizeof(partition_path_));
  std::strncpy(mass_path_, "./data/mass.txt", sizeof(mass_path_));
  std::strncpy(weak_path_, "./data/lmp_weak_rates.txt", sizeof(weak_path_));
  std::strncpy(load_abund_path_, "./abundances.txt", sizeof(load_abund_path_));
  std::strncpy(save_abund_path_, "./abundances.txt", sizeof(save_abund_path_));
  std::strncpy(save_state_path_, "./imnet_state.json",
               sizeof(save_state_path_));
  std::strncpy(load_traj_path_, "./trajectory.txt", sizeof(load_traj_path_));
  std::strncpy(save_traj_path_, "./trajectory_modified.txt",
               sizeof(save_traj_path_));
#ifdef IMNET_USE_NUPPN
  std::strncpy(species_path_, IMNET_NUPPN_RUN_DIR, sizeof(species_path_));
  std::memset(nuppn_frame_input_, 0, sizeof(nuppn_frame_input_));
  std::memset(nuppn_physics_input_, 0, sizeof(nuppn_physics_input_));
  std::memset(nuppn_solver_input_, 0, sizeof(nuppn_solver_input_));
#endif
  species_path_[sizeof(species_path_) - 1] = '\0';
  reaclib_path_[sizeof(reaclib_path_) - 1] = '\0';
  partition_path_[sizeof(partition_path_) - 1] = '\0';
  mass_path_[sizeof(mass_path_) - 1] = '\0';
  weak_path_[sizeof(weak_path_) - 1] = '\0';
  load_abund_path_[sizeof(load_abund_path_) - 1] = '\0';
  save_abund_path_[sizeof(save_abund_path_) - 1] = '\0';
  save_state_path_[sizeof(save_state_path_) - 1] = '\0';
  load_traj_path_[sizeof(load_traj_path_) - 1] = '\0';
  save_traj_path_[sizeof(save_traj_path_) - 1] = '\0';
  std::memset(isotope_filter_, 0, sizeof(isotope_filter_));
  std::memset(abundance_plot_filter_, 0, sizeof(abundance_plot_filter_));
  abundance_plot_log_floor_ = 1e-99;
  std::memset(file_name_input_, 0, sizeof(file_name_input_));

  // Initialize file browser
  current_browse_dir_ = std::filesystem::current_path().string();

  // Initialize GLFW
  glfwSetErrorCallback(glfw_error_callback);

  if (!glfwInit()) {
    throw std::runtime_error("Failed to initialize GLFW");
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  // Create window
  window_ = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
  if (!window_) {
    glfwTerminate();
    throw std::runtime_error("Failed to create GLFW window");
  }

  glfwMakeContextCurrent(window_);
  glfwSwapInterval(1); // Enable vsync

  g_main_window = this;
  glfwSetKeyCallback(window_, glfw_key_callback);

  try {
    init_imgui();
  } catch (...) {
    g_main_window = nullptr;
    glfwDestroyWindow(window_);
    window_ = nullptr;
    glfwTerminate();
    throw;
  }

  std::cout << "MainWindow initialized: " << width << "x" << height
            << std::endl;
}

MainWindow::~MainWindow() {
  shutdown_imgui();

  if (window_) {
    glfwDestroyWindow(window_);
  }
  glfwTerminate();

  g_main_window = nullptr;
}

void MainWindow::init_imgui() {
  const char *glsl_version = "#version 330";

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImPlot::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

  // Setup Dear ImGui style
  ui_theme::apply();

  // Setup Platform/Renderer backends
  if (!ImGui_ImplGlfw_InitForOpenGL(window_, true)) {
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    throw std::runtime_error("Failed to initialize ImGui GLFW backend");
  }
  if (!ImGui_ImplOpenGL3_Init(glsl_version)) {
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    throw std::runtime_error("Failed to initialize ImGui OpenGL backend");
  }
}

void MainWindow::shutdown_imgui() {
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImPlot::DestroyContext();
  ImGui::DestroyContext();
}

bool MainWindow::should_close() const {
  return should_close_ || glfwWindowShouldClose(window_);
}

void MainWindow::set_app_state(AppState *state) {
  app_state_ = state;
  chart_initialized_ = false;
  if (app_state_) {
    sync_network_path_inputs();
#ifdef IMNET_USE_NUPPN
    load_nuppn_input_files();
#endif
    nuclide_chart_.initialize(app_state_);
    chart_initialized_ = true;
  }
}

void MainWindow::request_close() { should_close_ = true; }

void MainWindow::process_frame() {
  // Poll events
  glfwPollEvents();

  // Check if window is still valid before rendering
  if (!window_ || glfwWindowShouldClose(window_)) {
    should_close_ = true;
    return;
  }

  process_trajectory_integration_job();

  // Start ImGui frame
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  // Render UI
  render_menu_bar();
  render_docking_layout();

  // Rendering
  ImGui::Render();

  // Get framebuffer size
  int display_w, display_h;
  glfwGetFramebufferSize(window_, &display_w, &display_h);

  // Only render if we have valid dimensions
  if (display_w > 0 && display_h > 0) {
    glViewport(0, 0, display_w, display_h);
    const ImVec4 clear_color = ui_theme::background();
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);

    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  }

  glfwSwapBuffers(window_);
}

void MainWindow::render_menu_bar() {
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
#ifdef IMNET_USE_NUPPN
      if (ImGui::MenuItem("Open NuPPN Run Directory", "Ctrl+O")) {
        sync_network_path_inputs();
        show_open_species_popup_ = true;
      }
#else
      if (ImGui::MenuItem("Open Rate Files", "Ctrl+O")) {
        sync_network_path_inputs();
        show_open_rates_popup_ = true;
      }
      if (ImGui::MenuItem("Open Species", "Ctrl+Shift+O")) {
        sync_network_path_inputs();
        show_open_species_popup_ = true;
      }
#endif
      ImGui::Separator();
      if (ImGui::MenuItem("Load Abundances")) {
        show_load_abund_popup_ = true;
      }
      if (ImGui::MenuItem("Save Abundances")) {
        show_save_abund_popup_ = true;
      }
      if (ImGui::MenuItem("Save State", "Ctrl+S")) {
        show_save_state_popup_ = true;
      }
      if (ImGui::MenuItem("Load Trajectory")) {
        show_load_traj_popup_ = true;
      }
      ImGui::Separator();
      if (ImGui::MenuItem("Exit", "Ctrl+Q")) {
        request_close();
      }
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("View")) {
      ImGui::MenuItem("Nuclide Chart", nullptr, &show_nuclide_chart_);
      ImGui::MenuItem("Workflow", nullptr, &show_integration_panel_);
      ImGui::MenuItem("Isotope Info", nullptr, &show_isotope_info_);
      ImGui::MenuItem("Trajectory Plot", nullptr, &show_trajectory_plot_);
      ImGui::MenuItem("Abundance Plot", nullptr, &show_abundance_plot_);
      ImGui::MenuItem("Trajectory Editor", nullptr, &show_trajectory_editor_);
      ImGui::MenuItem("Settings", nullptr, &show_settings_);
      ImGui::EndMenu();
    }

    if (ImGui::BeginMenu("Help")) {
      if (ImGui::MenuItem("About imnet")) {
        show_settings_ = true;
      }
      ImGui::EndMenu();
    }

    ImGui::EndMainMenuBar();
  }

  render_file_dialogs();
}

void MainWindow::render_docking_layout() {
  // Simplified layout without docking (docking requires imgui config)
  // Create a side-by-side layout

  ImGuiViewport *viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x, viewport->WorkPos.y + 20),
                          ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(740, 620), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSizeConstraints(ImVec2(350, 300),
                                      ImVec2(FLT_MAX, FLT_MAX));

  // Render panels
  if (show_nuclide_chart_) {
    ImGui::Begin("Nuclide Chart", &show_nuclide_chart_,
                 ImGuiWindowFlags_HorizontalScrollbar);
    if (app_state_) {
      ImGui::Text("Network: %d species loaded", app_state_->num_species());
      ImGui::Separator();
      if (!chart_initialized_) {
        nuclide_chart_.initialize(app_state_);
        chart_initialized_ = true;
      }
      nuclide_chart_.render();

      const auto &selected = nuclide_chart_.get_selected();
      auto &state_selected = app_state_->integration_settings().selected;
      if (selected.size() == state_selected.size()) {
        state_selected = selected;
      }
    } else {
      ImGui::TextDisabled("(No network loaded)");
    }
    ImGui::End();
  }

  // Position the integration panel next to the chart, but don't fix its size
  ImGui::SetNextWindowPos(
      ImVec2(viewport->WorkPos.x + 760, viewport->WorkPos.y + 20),
      ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(440, 700), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSizeConstraints(ImVec2(320, 360),
                                      ImVec2(FLT_MAX, FLT_MAX));

  if (show_integration_panel_) {
    ImGui::Begin("Workflow", &show_integration_panel_,
                 ImGuiWindowFlags_HorizontalScrollbar);

    if (app_state_) {
      ImGui::Text("Backend: %s", Network::backend_name());
      ImGui::SameLine();
      ImGui::TextDisabled("%d species", app_state_->num_species());
      ImGui::Separator();

      render_single_step_controls();

      ImGui::Separator();
      render_trajectory_controls();

      ImGui::Separator();
      if (ImGui::CollapsingHeader("Chart View",
                                  ImGuiTreeNodeFlags_DefaultOpen)) {
        auto &view = app_state_->view_settings();
        ImGui::Checkbox("Color by abundance", &view.show_abundance);
        ImGui::Checkbox("Highlight unconnected", &view.show_unconnected);
      }

      ImGui::Separator();
      ImGui::TextWrapped("Status: %s", status_message_.c_str());
    } else {
      ImGui::TextDisabled("No network loaded.");
    }

    // Result popup
    if (ImGui::BeginPopupModal("Integration Result", nullptr,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
      ImGui::TextWrapped("Integration completed!");
      ImGui::Text("dE/dt = %.6e erg/g/s", last_dedt_);
      ImGui::Text("substeps = %d",
                  app_state_->get_network()
                      ? app_state_->get_network()->last_substeps()
                      : 0);
      if (ImGui::Button("OK", ImVec2(120, 0))) {
        ImGui::CloseCurrentPopup();
      }
      ImGui::EndPopup();
    }

    ImGui::End();
  }

  if (show_isotope_info_) {
    ImGui::SetNextWindowPos(
        ImVec2(viewport->WorkPos.x, viewport->WorkPos.y + 450),
        ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(300, 100), ImVec2(FLT_MAX, 300));
    ImGui::Begin("Isotope Information", &show_isotope_info_);
    if (!app_state_ || !chart_initialized_) {
      ImGui::TextDisabled("No isotope selected.");
    } else {
      int idx = nuclide_chart_.get_single_selected_isotope();

      if (idx >= 0 && idx < app_state_->num_species()) {
        if (idx != isotope_info_index_) {
          refresh_isotope_diagnostics(idx);
        }
        const auto &names = app_state_->get_species_names();
        const auto &settings = app_state_->integration_settings();
        const double x = (idx < static_cast<int>(settings.xnuc.size()))
                             ? settings.xnuc[idx]
                             : 0.0;
        const bool fixed = (idx < static_cast<int>(settings.fixed.size()))
                               ? settings.fixed[idx]
                               : false;
        const bool active = (idx < static_cast<int>(settings.active.size()))
                                ? settings.active[idx]
                                : false;
        const bool selected = (idx < static_cast<int>(settings.selected.size()))
                                  ? settings.selected[idx]
                                  : false;

        ImGui::Text("Species: %s", names[idx].c_str());
        ImGui::Text("Index: %d", idx);
        ImGui::Text("X: %.6e", x);
        ImGui::Text("Net dX/dt: %.6e 1/s", isotope_info_rate_);
        ImGui::Text("rho = %.6e g/cm^3", settings.rho);
        ImGui::Text("T   = %.6e K", settings.temp);
        ImGui::Text("dt  = %.6e s", settings.dt);
        ImGui::Text("Flags: active=%s, fixed=%s, selected=%s",
                    active ? "true" : "false", fixed ? "true" : "false",
                    selected ? "true" : "false");
        ImGui::Separator();
        if (ImGui::Button("Refresh Reaction Rates")) {
          refresh_isotope_diagnostics(idx);
        }
        if (!isotope_info_error_.empty()) {
          ImGui::TextColored(ui_theme::danger(), "%s",
                             isotope_info_error_.c_str());
        } else if (isotope_reactions_.empty()) {
          ImGui::TextDisabled("No reaction contributors found.");
        } else if (ImGui::BeginTable("SelectedReactionRates", 5,
                                     ImGuiTableFlags_Borders |
                                         ImGuiTableFlags_RowBg |
                                         ImGuiTableFlags_Resizable |
                                         ImGuiTableFlags_Sortable |
                                         ImGuiTableFlags_SortMulti)) {
          ImGui::TableSetupColumn("Reaction",
                                  ImGuiTableColumnFlags_DefaultSort);
          ImGui::TableSetupColumn("Q");
          ImGui::TableSetupColumn("Rate");
          ImGui::TableSetupColumn("dY/dt");
          ImGui::TableSetupColumn("dX/dt",
                                  ImGuiTableColumnFlags_PreferSortDescending);
          ImGui::TableHeadersRow();
          std::vector<ReactionDiagnostic> sorted_reactions = isotope_reactions_;
          if (ImGuiTableSortSpecs *sort_specs = ImGui::TableGetSortSpecs()) {
            if (sort_specs->SpecsCount > 0) {
              std::sort(
                  sorted_reactions.begin(), sorted_reactions.end(),
                  [sort_specs](const ReactionDiagnostic &a,
                               const ReactionDiagnostic &b) {
                    for (int n = 0; n < sort_specs->SpecsCount; ++n) {
                      const ImGuiTableColumnSortSpecs &spec =
                          sort_specs->Specs[n];
                      int cmp = 0;
                      switch (spec.ColumnIndex) {
                      case 0:
                        cmp = a.equation.compare(b.equation);
                        break;
                      case 1:
                        cmp = (a.q_value > b.q_value) - (a.q_value < b.q_value);
                        break;
                      case 2:
                        cmp = (a.rate > b.rate) - (a.rate < b.rate);
                        break;
                      case 3:
                        cmp = (a.contribution_dYdt > b.contribution_dYdt) -
                              (a.contribution_dYdt < b.contribution_dYdt);
                        break;
                      case 4:
                        cmp = (a.contribution_dXdt > b.contribution_dXdt) -
                              (a.contribution_dXdt < b.contribution_dXdt);
                        break;
                      default:
                        break;
                      }
                      if (cmp != 0) {
                        return spec.SortDirection ==
                                       ImGuiSortDirection_Ascending
                                   ? cmp < 0
                                   : cmp > 0;
                      }
                    }
                    return false;
                  });
              sort_specs->SpecsDirty = false;
            }
          }
          for (const auto &reaction : sorted_reactions) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextWrapped("%s%s", reaction.weak ? "[weak] " : "",
                               reaction.equation.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%.6e", reaction.q_value);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.6e", reaction.rate);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.6e", reaction.contribution_dYdt);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.6e", reaction.contribution_dXdt);
          }
          ImGui::EndTable();
        }
      } else {
        ImGui::TextDisabled("No isotope selected.");
      }
    }
    ImGui::End();
  }

  if (show_settings_) {
    ImGui::SetNextWindowPos(
        ImVec2(viewport->WorkPos.x + 100, viewport->WorkPos.y + 100),
        ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(400, 300),
                                        ImVec2(FLT_MAX, FLT_MAX));
    ImGui::Begin("About / Settings", &show_settings_);
    ImGui::TextWrapped("imnet v%s", IMNET_VERSION);
    ImGui::TextWrapped("Interactive %s Nuclear Network", Network::backend_name());
    ImGui::Separator();
    ImGui::TextWrapped("Built with:");
    ImGui::BulletText("%s network backend", Network::backend_name());
    ImGui::BulletText("C++17 application layer");
    ImGui::BulletText("imgui GUI framework");
    ImGui::BulletText("implot plotting");
    ImGui::BulletText("GLFW windowing");
    ImGui::BulletText("OpenGL 3.3+ rendering");
    ImGui::Separator();
    if (app_state_) {
#ifdef IMNET_USE_NUPPN
      ImGui::TextWrapped("NuPPN Run Directory");
      ImGui::TextWrapped("%s", app_state_->species_file().c_str());
      if (ImGui::Button("Reload Text##nuppn_inputs", ImVec2(130, 0))) {
        nuppn_inputs_loaded_ = load_nuppn_input_files();
      }
      ImGui::SameLine();
      if (ImGui::Button("Save Input Files##nuppn_inputs", ImVec2(170, 0))) {
        if (save_nuppn_input_files()) {
          status_message_ = "NuPPN input files saved";
        } else {
          status_message_ = "Could not save NuPPN input files";
        }
      }
      if (!nuppn_inputs_loaded_) {
        ImGui::TextColored(ui_theme::danger(),
                           "Could not load one or more NuPPN input files.");
      }
      ImGui::TextWrapped("ppn_frame.input");
      ImGui::InputTextMultiline("##nuppn_frame_input", nuppn_frame_input_,
                                sizeof(nuppn_frame_input_),
                                ImVec2(-1.0f, 120.0f));
      ImGui::TextWrapped("ppn_physics.input");
      ImGui::InputTextMultiline("##nuppn_physics_input", nuppn_physics_input_,
                                sizeof(nuppn_physics_input_),
                                ImVec2(-1.0f, 180.0f));
      ImGui::TextWrapped("ppn_solver.input");
      ImGui::InputTextMultiline("##nuppn_solver_input", nuppn_solver_input_,
                                sizeof(nuppn_solver_input_),
                                ImVec2(-1.0f, 120.0f));
      ImGui::Separator();
#endif
      auto &view = app_state_->view_settings();
      ImGui::TextWrapped("Trajectory Cache");
      ImGui::SetNextItemWidth(140.0f);
      if (ImGui::InputInt("Max cached timesteps",
                          &view.max_cached_trajectory_steps, 1, 10)) {
        view.max_cached_trajectory_steps =
            std::max(1, view.max_cached_trajectory_steps);
        retarget_trajectory_integration_job();
      }
    }
    ImGui::End();
  }

  render_trajectory_plot_panel();
  render_abundance_plot_panel();
  render_trajectory_editor_panel();
}

void MainWindow::render_single_step_controls() {
  auto &settings = app_state_->integration_settings();

  if (ImGui::CollapsingHeader("Run Setup",
                              ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::BeginTabBar("SingleStepTabs")) {
      if (ImGui::BeginTabItem("Conditions")) {
        constexpr float condition_input_width = 240.0f;
        float log_rho =
            static_cast<float>(std::log10(std::max(settings.rho, 1e-300)));
        float log_temp =
            static_cast<float>(std::log10(std::max(settings.temp, 1e-300)));
        float log_dt =
            static_cast<float>(std::log10(std::max(settings.dt, 1e-300)));

        ImGui::TextWrapped("Physical Conditions");
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::SliderFloat("log10(rho)##single", &log_rho, -2.0f, 12.0f,
                               "%.3f")) {
          settings.rho = std::pow(10.0, static_cast<double>(log_rho));
        }
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::InputDouble("rho (g/cm^3)", &settings.rho, 0.0, 0.0,
                               "%.9e")) {
          settings.rho = std::max(settings.rho, 0.0);
        }

        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::SliderFloat("log10(T)##single", &log_temp, 6.0f, 11.0f,
                               "%.3f")) {
          settings.temp = std::pow(10.0, static_cast<double>(log_temp));
        }
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::InputDouble("T (K)", &settings.temp, 0.0, 0.0, "%.9e")) {
          settings.temp = std::max(settings.temp, 0.0);
        }

        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::SliderFloat("log10(dt)##single", &log_dt, -12.0f, 4.0f,
                               "%.3f")) {
          settings.dt = std::pow(10.0, static_cast<double>(log_dt));
        }
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::InputDouble("dt (s)", &settings.dt, 0.0, 0.0, "%.9e")) {
          settings.dt = std::max(settings.dt, 0.0);
        }

        ImGui::Separator();
        ImGui::TextWrapped("Final-Time Mode");
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::InputDouble("final time (s)", &settings.final_time, 0.0,
                               0.0, "%.9e")) {
          settings.final_time = std::max(settings.final_time, 0.0);
        }
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::InputDouble("dt max (s)", &settings.dt_max, 0.0, 0.0,
                               "%.9e")) {
          settings.dt_max = std::max(settings.dt_max, 0.0);
        }
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::InputDouble("dt factor", &settings.dt_factor, 0.0, 0.0,
                               "%.6g")) {
          settings.dt_factor = std::max(settings.dt_factor, 0.0);
        }
        ImGui::SetNextItemWidth(condition_input_width);
        if (ImGui::InputInt("max steps", &settings.max_steps, 1, 10)) {
          settings.max_steps = std::max(settings.max_steps, 1);
        }

        ImGui::Separator();
        ImGui::Text("Abar %.9e", app_state_->calculate_abar());
        ImGui::Text("Zbar %.9e", app_state_->calculate_zbar());
        ImGui::Text("Ye   %.9e", app_state_->calculate_ye());

        std::string validation_error;
        if (app_state_->validate_integration_settings(validation_error)) {
          ImGui::TextColored(ui_theme::success(), "Current state is valid");
        } else {
          ImGui::TextColored(ui_theme::danger(), "%s",
                             validation_error.c_str());
        }

        ImGui::EndTabItem();
      }

      if (ImGui::BeginTabItem("Composition")) {
        render_composition_editor();
        ImGui::EndTabItem();
      }

      if (ImGui::BeginTabItem("Run")) {
        const char *run_modes[] = {"Fixed timestep", "Final time",
                                   "Loaded trajectory"};
        run_mode_ = std::clamp(run_mode_, 0, IM_ARRAYSIZE(run_modes) - 1);

        auto refresh_selected_isotope = [&]() {
          const int selected = nuclide_chart_.get_single_selected_isotope();
          if (selected >= 0) {
            refresh_isotope_diagnostics(selected);
          }
        };

        auto run_fixed_timestep = [&]() {
          std::string error;
          if (app_state_->integrate_single_step(normalize_before_integrate_,
                                                last_dedt_, error)) {
            cancel_trajectory_integration_job();
            status_message_ = "Single-step integration completed";
            refresh_selected_isotope();
            ImGui::OpenPopup("Integration Result");
          } else {
            status_message_ = "Integration failed: " + error;
          }
        };

        auto run_final_time = [&]() {
          std::string error;
          if (app_state_->run_to_time(normalize_before_integrate_, error)) {
            cancel_trajectory_integration_job();
            const auto &cache = app_state_->trajectory_cache();
            if (!cache.empty()) {
              last_dedt_ = cache.back().dedt;
            }
            trajectory_step_ui_ =
                static_cast<int>(app_state_->current_trajectory_step());
            const int substeps =
                app_state_->get_network()
                    ? app_state_->get_network()->last_substeps()
                    : 0;
            status_message_ =
                "Final-time integration completed, substeps: " +
                std::to_string(substeps);
            refresh_selected_isotope();
          } else {
            status_message_ = "Final-time integration failed: " + error;
          }
        };

        ImGui::TextWrapped("Run mode");
        ImGui::SetNextItemWidth(-1.0f);
        ImGui::Combo("##run_mode", &run_mode_, run_modes,
                     IM_ARRAYSIZE(run_modes));
        ImGui::Checkbox("Normalize before run", &normalize_before_integrate_);

        if (trajectory_job_.active) {
          if (ImGui::Button("Cancel Trajectory Run", ImVec2(-1.0f, 0))) {
            cancel_trajectory_integration_job();
          }
        } else {
          const bool can_run = run_mode_ != 2 || app_state_->has_trajectory();
          ImGui::BeginDisabled(!can_run);
          if (ImGui::Button("Run", ImVec2(-1.0f, 0))) {
            if (run_mode_ == 0) {
              run_fixed_timestep();
            } else if (run_mode_ == 1) {
              run_final_time();
            } else {
              start_trajectory_integration_job();
            }
          }
          ImGui::EndDisabled();
        }

        if (run_mode_ == 2 && !app_state_->has_trajectory()) {
          if (ImGui::Button("Load Trajectory...", ImVec2(-1.0f, 0))) {
            show_load_traj_popup_ = true;
          }
        }

        ImGui::Separator();
        ImGui::Text("Last dE/dt %.9e erg/g/s", last_dedt_);
        ImGui::Text("Last substeps %d",
                    app_state_->get_network()
                        ? app_state_->get_network()->last_substeps()
                        : 0);

        ImGui::Separator();
        ImGui::TextWrapped("Composition");
        const char *presets[] = {"he4", "c12/o16 50/50", "o16/ne20 50/50",
                                 "ni56"};
        composition_preset_index_ =
            std::clamp(composition_preset_index_, 0, IM_ARRAYSIZE(presets) - 1);
        ImGui::SetNextItemWidth(-88.0f);
        ImGui::Combo("##composition_preset", &composition_preset_index_,
                     presets, IM_ARRAYSIZE(presets));
        ImGui::SameLine();
        if (ImGui::Button("Apply", ImVec2(76, 0))) {
          switch (composition_preset_index_) {
          case 0:
            apply_composition_preset({{"he4", 1.0}});
            break;
          case 1:
            apply_composition_preset({{"c12", 0.5}, {"o16", 0.5}});
            break;
          case 2:
            apply_composition_preset({{"o16", 0.5}, {"ne20", 0.5}});
            break;
          default:
            apply_composition_preset({{"ni56", 1.0}});
            break;
          }
        }

        if (ImGui::Button("Compute NSE", ImVec2(-1.0f, 0))) {
          std::string error;
          if (app_state_->compute_nse(error)) {
            cancel_trajectory_integration_job();
            status_message_ = "NSE composition applied";
          } else {
            status_message_ = "NSE failed: " + error;
          }
        }

        const float half_button_width =
            (ImGui::GetContentRegionAvail().x -
             ImGui::GetStyle().ItemSpacing.x) *
            0.5f;
        if (ImGui::Button("Normalize", ImVec2(half_button_width, 0))) {
          if (app_state_->normalize_abundances()) {
            cancel_trajectory_integration_job();
            status_message_ = "Abundances normalized";
          } else {
            status_message_ = "Could not normalize abundances";
          }
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset", ImVec2(half_button_width, 0))) {
          app_state_->reset_abundances();
          cancel_trajectory_integration_job();
          status_message_ = "Abundances reset";
        }
        if (ImGui::Button("Clear X", ImVec2(-1.0f, 0))) {
          std::fill(settings.xnuc.begin(), settings.xnuc.end(), 0.0);
          composition_changed();
          status_message_ = "Abundances cleared";
        }

        ImGui::Separator();
        ImGui::Checkbox("Trajectory plot", &show_trajectory_plot_);
        ImGui::Checkbox("Abundance plot", &show_abundance_plot_);
        ImGui::Checkbox("Trajectory editor", &show_trajectory_editor_);
        ImGui::EndTabItem();
      }

      ImGui::EndTabBar();
    }
  }
}

void MainWindow::render_composition_editor() {
  auto &settings = app_state_->integration_settings();
  const auto &species_names = app_state_->get_species_names();

  ImGui::Checkbox("Show all isotopes", &show_all_abundances_);
  ImGui::SameLine();
  ImGui::Checkbox("Nonzero only", &show_nonzero_only_);
  ImGui::SetNextItemWidth(-1.0f);
  ImGui::InputTextWithHint("##isotope_filter", "Filter isotope name",
                           isotope_filter_, sizeof(isotope_filter_));

  if (!show_all_abundances_) {
    render_selected_abundance_editor();
    return;
  }

  std::string filter = isotope_filter_;
  std::transform(filter.begin(), filter.end(), filter.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  const ImGuiTableFlags flags =
      ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
      ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY;
  if (ImGui::BeginTable("CompositionTable", 5, flags, ImVec2(0, 300))) {
    ImGui::TableSetupColumn("Iso");
    ImGui::TableSetupColumn("X");
    ImGui::TableSetupColumn("Active");
    ImGui::TableSetupColumn("Fixed");
    ImGui::TableSetupColumn("Select");
    ImGui::TableHeadersRow();

    for (size_t i = 0; i < species_names.size(); ++i) {
      if (i >= settings.xnuc.size() || i >= settings.active.size() ||
          i >= settings.fixed.size() || i >= settings.selected.size()) {
        continue;
      }
      if (show_nonzero_only_ && settings.xnuc[i] <= 0.0) {
        continue;
      }
      std::string name = species_names[i];
      std::string name_lower = name;
      std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      if (!filter.empty() && name_lower.find(filter) == std::string::npos) {
        continue;
      }

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      if (ImGui::Selectable((name + "##row_" + std::to_string(i)).c_str(),
                            settings.selected[i])) {
        select_isotope(i, ImGui::GetIO().KeyCtrl);
      }

      ImGui::TableSetColumnIndex(1);
      ImGui::SetNextItemWidth(-1.0f);
      double x = settings.xnuc[i];
      if (ImGui::InputDouble(("##x_" + std::to_string(i)).c_str(), &x, 0.0, 0.0,
                             "%.6e")) {
        settings.xnuc[i] = std::isfinite(x) ? std::max(0.0, x) : 0.0;
        composition_changed();
      }

      ImGui::TableSetColumnIndex(2);
      bool active = settings.active[i];
      if (ImGui::Checkbox(("##active_" + std::to_string(i)).c_str(), &active)) {
        settings.active[i] = active;
      }

      ImGui::TableSetColumnIndex(3);
      bool fixed = settings.fixed[i];
      if (ImGui::Checkbox(("##fixed_" + std::to_string(i)).c_str(), &fixed)) {
        settings.fixed[i] = fixed;
      }

      ImGui::TableSetColumnIndex(4);
      if (ImGui::SmallButton(("Only##only_" + std::to_string(i)).c_str())) {
        select_isotope(i, false);
      }
    }

    ImGui::EndTable();
  }
}

void MainWindow::render_trajectory_controls() {
  if (!ImGui::CollapsingHeader("Cached Steps",
                               ImGuiTreeNodeFlags_DefaultOpen)) {
    return;
  }

  if (!app_state_->has_trajectory()) {
    ImGui::TextDisabled("No trajectory loaded.");
    if (ImGui::Button("Load Trajectory...", ImVec2(-1.0f, 0))) {
      show_load_traj_popup_ = true;
    }
    return;
  }

  const size_t loaded_steps = app_state_->trajectory_size();
  const auto &cache = app_state_->trajectory_cache();
  auto apply_step = [&]() {
    if (app_state_->set_trajectory_step(
            static_cast<size_t>(trajectory_step_ui_))) {
      status_message_ = "Applied trajectory step";
      retarget_trajectory_integration_job();
      const int selected = nuclide_chart_.get_single_selected_isotope();
      if (selected >= 0) {
        refresh_isotope_diagnostics(selected);
      }
    }
  };

  ImGui::Text("Rows %d", static_cast<int>(loaded_steps));
  ImGui::SameLine();
  ImGui::TextDisabled("cached %d / max %d", static_cast<int>(cache.size()),
                      app_state_->view_settings().max_cached_trajectory_steps);

  if (trajectory_job_.active || trajectory_job_.complete ||
      trajectory_job_.failed) {
    const float denominator =
        loaded_steps > 1 ? static_cast<float>(loaded_steps - 1) : 1.0f;
    const float progress =
        std::min(static_cast<float>(trajectory_job_.current_step) / denominator,
                 1.0f);
    ImGui::ProgressBar(progress, ImVec2(-1.0f, 0.0f));
    ImGui::Text("Integrated through step %d of %d",
                static_cast<int>(trajectory_job_.current_step),
                static_cast<int>(loaded_steps - 1));
    ImGui::Text("Cache window: %d - %d",
                static_cast<int>(trajectory_job_.cache_start_step),
                static_cast<int>(trajectory_job_.cache_end_step));
    if (trajectory_job_.failed) {
      ImGui::TextColored(ui_theme::danger(), "%s",
                         trajectory_job_.error.c_str());
    } else if (trajectory_job_.complete) {
      ImGui::TextColored(ui_theme::success(),
                         "Current cache window is ready.");
    }
  }

  if (loaded_steps == 0) {
    return;
  }

  const int max_available_step = static_cast<int>(loaded_steps - 1);
  if (trajectory_step_ui_ > max_available_step) {
    trajectory_step_ui_ = max_available_step;
  }

  ImGui::Separator();
  ImGui::TextWrapped("Step");
  ImGui::BeginDisabled(trajectory_step_ui_ <= 0);
  if (ImGui::Button("Prev") && trajectory_step_ui_ > 0) {
    --trajectory_step_ui_;
    apply_step();
  }
  ImGui::EndDisabled();
  ImGui::SameLine();
  ImGui::BeginDisabled(trajectory_step_ui_ >= max_available_step);
  if (ImGui::Button("Next") && trajectory_step_ui_ < max_available_step) {
    ++trajectory_step_ui_;
    apply_step();
  }
  ImGui::EndDisabled();

  ImGui::SetNextItemWidth(-1.0f);
  if (ImGui::SliderInt("##trajectory_step_slider", &trajectory_step_ui_, 0,
                       max_available_step)) {
    apply_step();
  }
  ImGui::SetNextItemWidth(120.0f);
  if (ImGui::InputInt("Step index", &trajectory_step_ui_, 1, 10)) {
    trajectory_step_ui_ =
        std::clamp(trajectory_step_ui_, 0, max_available_step);
    apply_step();
  }

  const auto *step = app_state_->get_trajectory_cache_step(
      static_cast<size_t>(trajectory_step_ui_));
  const auto &times = app_state_->trajectory_times();
  if (!step && trajectory_step_ui_ >= 0 &&
      trajectory_step_ui_ < static_cast<int>(times.size())) {
    ImGui::Text("t = %.6e s", times[static_cast<size_t>(trajectory_step_ui_)]);
    ImGui::TextDisabled("This row has not been cached yet.");
    return;
  }
  if (!step) {
    return;
  }

  ImGui::Separator();
  ImGui::Text("t      = %.9e s", step->time);
  ImGui::Text("rho    = %.9e g/cm^3", step->rho);
  ImGui::Text("T      = %.9e K", step->temp);
  ImGui::Text("dt     = %.9e s", step->dt);
  ImGui::Text("dE/dt  = %.9e erg/g/s", step->dedt);
  ImGui::Text("substeps = %d", step->substeps);
  ImGui::Text("status = %s", step->status.c_str());
  if (!step->success && !step->error.empty()) {
    ImGui::TextColored(ui_theme::danger(), "%s",
                       step->error.c_str());
  }

  if (!cache.empty() &&
      ImGui::BeginTable("TrajectoryCacheTable", 7,
                        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                            ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY,
                        ImVec2(0, 170))) {
    ImGui::TableSetupColumn("Step");
    ImGui::TableSetupColumn("t");
    ImGui::TableSetupColumn("rho");
    ImGui::TableSetupColumn("T");
    ImGui::TableSetupColumn("dE/dt");
    ImGui::TableSetupColumn("Substeps");
    ImGui::TableSetupColumn("Status");
    ImGui::TableHeadersRow();

    for (const auto &cached : cache) {
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      if (ImGui::Selectable((std::to_string(cached.index) + "##traj_row_" +
                             std::to_string(cached.index))
                                .c_str(),
                            cached.index ==
                                static_cast<size_t>(trajectory_step_ui_))) {
        trajectory_step_ui_ = static_cast<int>(cached.index);
        apply_step();
      }
      ImGui::TableSetColumnIndex(1);
      ImGui::Text("%.3e", cached.time);
      ImGui::TableSetColumnIndex(2);
      ImGui::Text("%.3e", cached.rho);
      ImGui::TableSetColumnIndex(3);
      ImGui::Text("%.3e", cached.temp);
      ImGui::TableSetColumnIndex(4);
      ImGui::Text("%.3e", cached.dedt);
      ImGui::TableSetColumnIndex(5);
      ImGui::Text("%d", cached.substeps);
      ImGui::TableSetColumnIndex(6);
      if (cached.success) {
        ImGui::Text("%s", cached.status.c_str());
      } else {
        ImGui::TextColored(ui_theme::danger(), "%s",
                           cached.status.c_str());
      }
    }

    ImGui::EndTable();
  }
}

void MainWindow::start_trajectory_integration_job() {
  if (!app_state_ || !app_state_->has_trajectory() ||
      !app_state_->has_network()) {
    status_message_ = "Cannot start trajectory integration";
    return;
  }

  std::string error;
  if (!app_state_->validate_integration_settings(error)) {
    status_message_ = "Trajectory integration cannot start: " + error;
    return;
  }

  trajectory_job_ = TrajectoryIntegrationJob();
  trajectory_job_.active = true;
  trajectory_job_.initial_xnuc = app_state_->integration_settings().xnuc;
  trajectory_job_.xnuc = trajectory_job_.initial_xnuc;
  retarget_trajectory_integration_job();

  status_message_ = "Background trajectory integration started";
}

void MainWindow::cancel_trajectory_integration_job() {
  if (trajectory_job_.active) {
    status_message_ = "Background trajectory integration cancelled";
  }
  trajectory_job_.active = false;
  trajectory_job_.complete = false;
}

void MainWindow::retarget_trajectory_integration_job() {
  if (!app_state_ || !app_state_->has_trajectory() ||
      trajectory_job_.initial_xnuc.empty() || trajectory_job_.failed) {
    return;
  }

  const size_t loaded_steps = app_state_->trajectory_size();
  if (loaded_steps == 0) {
    return;
  }

  const size_t max_cached_steps = static_cast<size_t>(
      std::max(1, app_state_->view_settings().max_cached_trajectory_steps));
  const size_t center = static_cast<size_t>(
      std::clamp(trajectory_step_ui_, 0, static_cast<int>(loaded_steps - 1)));
  const size_t half_window = max_cached_steps / 2;
  size_t start = center > half_window ? center - half_window : 0;
  size_t end = std::min(loaded_steps - 1, start + max_cached_steps - 1);
  if (end == loaded_steps - 1 && end + 1 >= max_cached_steps) {
    start = end + 1 - max_cached_steps;
  }

  const bool needs_restart =
      trajectory_job_.xnuc.empty() || end < trajectory_job_.current_step;
  if (needs_restart) {
    trajectory_job_.current_step = 0;
    trajectory_job_.xnuc = trajectory_job_.initial_xnuc;
    trajectory_job_.cache.clear();
    trajectory_job_.failed = false;
    trajectory_job_.error.clear();
  }

  trajectory_job_.cache_start_step = start;
  trajectory_job_.cache_end_step = end;
  trajectory_job_.target_end_step = end;
  prune_trajectory_cache(trajectory_job_.cache, start, end);

  const auto &times = app_state_->trajectory_times();
  const auto &rhos = app_state_->trajectory_rhos();
  const auto &temps = app_state_->trajectory_temps();
  if (trajectory_job_.current_step == 0 && start == 0) {
    upsert_trajectory_cache_step(
        trajectory_job_.cache,
        make_trajectory_step_cache(0, times, rhos, temps, trajectory_job_.xnuc,
                                   0.0, 0, true, "initial", ""));
  }

  trajectory_job_.complete =
      trajectory_job_.current_step >= trajectory_job_.target_end_step;
  trajectory_job_.active = !trajectory_job_.complete;
  app_state_->set_trajectory_cache(trajectory_job_.cache, center);
}

void MainWindow::process_trajectory_integration_job() {
  if (!trajectory_job_.active || !app_state_ || !app_state_->has_trajectory()) {
    return;
  }

  retarget_trajectory_integration_job();
  if (!trajectory_job_.active || trajectory_job_.failed) {
    return;
  }

  const auto &times = app_state_->trajectory_times();
  const auto &rhos = app_state_->trajectory_rhos();
  const auto &temps = app_state_->trajectory_temps();
  const size_t loaded_steps = times.size();
  if (trajectory_job_.current_step + 1 >= loaded_steps ||
      trajectory_job_.current_step >= trajectory_job_.target_end_step) {
    trajectory_job_.active = false;
    trajectory_job_.complete = true;
    app_state_->set_trajectory_cache(
        trajectory_job_.cache,
        static_cast<size_t>(std::clamp(
            trajectory_step_ui_, 0, static_cast<int>(loaded_steps - 1))));
    return;
  }

  const size_t i = trajectory_job_.current_step;
  const double dt = times[i + 1] - times[i];
  if (!std::isfinite(dt) || dt < 0.0) {
    trajectory_job_.failed = true;
    trajectory_job_.active = false;
    trajectory_job_.error =
        "Invalid timestep at trajectory interval " + std::to_string(i);
    status_message_ = "Trajectory integration failed: " + trajectory_job_.error;
    return;
  }

  double dedt = 0.0;
  int substeps = 0;
  try {
    std::vector<double> next_xnuc = trajectory_job_.xnuc;
    if (dt > 0.0) {
      dedt = app_state_->get_network()->integrate_interval(
          rhos[i], temps[i], rhos[i + 1], temps[i + 1], next_xnuc, dt);
      substeps = app_state_->get_network()->last_substeps();
    }
    for (size_t j = 0; j < next_xnuc.size(); ++j) {
      if (!std::isfinite(next_xnuc[j]) || next_xnuc[j] < 0.0) {
        throw std::runtime_error("Invalid abundance at species index " +
                                 std::to_string(j));
      }
    }
    trajectory_job_.xnuc = std::move(next_xnuc);
  } catch (const std::exception &e) {
    trajectory_job_.failed = true;
    trajectory_job_.active = false;
    trajectory_job_.error =
        "Step " + std::to_string(i + 1) + " failed: " + e.what();
    if (i + 1 >= trajectory_job_.cache_start_step &&
        i + 1 <= trajectory_job_.cache_end_step) {
      upsert_trajectory_cache_step(
          trajectory_job_.cache,
          make_trajectory_step_cache(i + 1, times, rhos, temps,
                                     trajectory_job_.xnuc, dedt, substeps,
                                     false,
                                     "failed", trajectory_job_.error));
    }
    status_message_ = "Trajectory integration failed: " + trajectory_job_.error;
    app_state_->set_trajectory_cache(
        trajectory_job_.cache,
        static_cast<size_t>(std::clamp(
            trajectory_step_ui_, 0, static_cast<int>(loaded_steps - 1))));
    return;
  }

  ++trajectory_job_.current_step;
  if (trajectory_job_.current_step >= trajectory_job_.cache_start_step &&
      trajectory_job_.current_step <= trajectory_job_.cache_end_step) {
    upsert_trajectory_cache_step(
        trajectory_job_.cache,
        make_trajectory_step_cache(trajectory_job_.current_step, times, rhos,
                                   temps, trajectory_job_.xnuc, dedt, substeps,
                                   true,
                                   "ok", ""));
  }
  prune_trajectory_cache(trajectory_job_.cache,
                         trajectory_job_.cache_start_step,
                         trajectory_job_.cache_end_step);

  const size_t selected_step = static_cast<size_t>(std::clamp(
      trajectory_step_ui_, 0, static_cast<int>(loaded_steps - 1)));
  app_state_->set_trajectory_cache(trajectory_job_.cache, selected_step);

  if (trajectory_job_.current_step >= trajectory_job_.target_end_step) {
    trajectory_job_.active = false;
    trajectory_job_.complete = true;
    status_message_ = "Trajectory cache window integrated";
  }
}

void MainWindow::render_trajectory_plot_panel() {
  if (!show_trajectory_plot_) {
    return;
  }

  ImGuiViewport *viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(
      ImVec2(viewport->WorkPos.x + 760, viewport->WorkPos.y + 20),
      ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(520, 520), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSizeConstraints(ImVec2(360, 320),
                                      ImVec2(FLT_MAX, FLT_MAX));

  ImGui::Begin("Trajectory Plot", &show_trajectory_plot_,
               ImGuiWindowFlags_HorizontalScrollbar);

  if (!app_state_ || !app_state_->has_trajectory()) {
    ImGui::TextDisabled("No trajectory loaded.");
    if (ImGui::Button("Load Trajectory...", ImVec2(-1.0f, 0))) {
      show_load_traj_popup_ = true;
    }
    ImGui::End();
    return;
  }

  const auto &times = app_state_->trajectory_times();
  const auto &rhos = app_state_->trajectory_rhos();
  const auto &temps = app_state_->trajectory_temps();
  const int count = static_cast<int>(
      std::min({times.size(), rhos.size(), temps.size()}));

  if (count <= 0) {
    ImGui::TextDisabled("No trajectory samples available.");
    ImGui::End();
    return;
  }

  const size_t current_index =
      std::min(app_state_->current_trajectory_step(), times.size() - 1);
  const double current_time = times[current_index];
  const bool can_use_log_x =
      std::all_of(times.begin(), times.end(), [](double t) { return t > 0.0; });

  ImGui::Text("Samples: %d", count);
  ImGui::SameLine();
  ImGui::Text("Current step: %d", static_cast<int>(current_index));
  ImGui::Checkbox("Log x", &trajectory_plot_log_x_);
  ImGui::SameLine();
  ImGui::Checkbox("Log y", &trajectory_plot_log_y_);
  if (trajectory_plot_log_x_ && !can_use_log_x) {
    ImGui::TextDisabled("Log x requires all trajectory times to be positive.");
  }

  if (ImPlot::BeginPlot("rho(t)##trajectory_rho", ImVec2(-1, 220))) {
    ImPlot::SetupAxes("t (s)", "rho (g/cm^3)", ImPlotAxisFlags_AutoFit,
                      ImPlotAxisFlags_AutoFit);
    if (trajectory_plot_log_x_ && can_use_log_x) {
      ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
    }
    if (trajectory_plot_log_y_) {
      ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
    }
    ImPlot::PlotLine("rho", times.data(), rhos.data(), count);
    ImPlot::PlotInfLines("current step", &current_time, 1);
    ImPlot::EndPlot();
  }

  if (ImPlot::BeginPlot("T(t)##trajectory_temp", ImVec2(-1, 220))) {
    ImPlot::SetupAxes("t (s)", "T (K)", ImPlotAxisFlags_AutoFit,
                      ImPlotAxisFlags_AutoFit);
    if (trajectory_plot_log_x_ && can_use_log_x) {
      ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
    }
    if (trajectory_plot_log_y_) {
      ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
    }
    ImPlot::PlotLine("T", times.data(), temps.data(), count);
    ImPlot::PlotInfLines("current step", &current_time, 1);
    ImPlot::EndPlot();
  }

  ImGui::End();
}

void MainWindow::render_abundance_plot_panel() {
  if (!show_abundance_plot_) {
    return;
  }

  ImGuiViewport *viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(
      ImVec2(viewport->WorkPos.x + 800, viewport->WorkPos.y + 80),
      ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(620, 520), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSizeConstraints(ImVec2(420, 320),
                                      ImVec2(FLT_MAX, FLT_MAX));

  ImGui::Begin("Abundance Plot", &show_abundance_plot_,
               ImGuiWindowFlags_HorizontalScrollbar);

  if (!app_state_) {
    ImGui::TextDisabled("No network loaded.");
    ImGui::End();
    return;
  }

  const auto &cache = app_state_->trajectory_cache();
  const auto &species_names = app_state_->get_species_names();
  abundance_plot_isotopes_.erase(
      std::remove_if(abundance_plot_isotopes_.begin(),
                     abundance_plot_isotopes_.end(),
                     [&species_names](int index) {
                       return index < 0 ||
                              index >= static_cast<int>(species_names.size());
                     }),
      abundance_plot_isotopes_.end());

  std::string filter = abundance_plot_filter_;
  std::transform(filter.begin(), filter.end(), filter.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  ImGui::SetNextItemWidth(-1.0f);
  ImGui::InputTextWithHint("##abundance_plot_filter", "Add isotope",
                           abundance_plot_filter_,
                           sizeof(abundance_plot_filter_));

  const int chart_index = nuclide_chart_.get_single_selected_isotope();
  if (chart_index >= 0 &&
      chart_index < static_cast<int>(species_names.size()) &&
      std::find(abundance_plot_isotopes_.begin(),
                abundance_plot_isotopes_.end(),
                chart_index) == abundance_plot_isotopes_.end()) {
    if (ImGui::Button(("Add " + species_names[chart_index]).c_str(),
                      ImVec2(140, 0))) {
      abundance_plot_isotopes_.push_back(chart_index);
    }
  }

  if (!filter.empty() &&
      ImGui::BeginTable("AbundancePlotSearch", 2,
                        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                            ImGuiTableFlags_ScrollY,
                        ImVec2(0, 120))) {
    ImGui::TableSetupColumn("Iso");
    ImGui::TableSetupColumn("Add", ImGuiTableColumnFlags_WidthFixed, 70.0f);
    ImGui::TableHeadersRow();

    int matches = 0;
    for (size_t i = 0; i < species_names.size() && matches < 24; ++i) {
      std::string name_lower = species_names[i];
      std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      if (name_lower.find(filter) == std::string::npos) {
        continue;
      }

      const int index = static_cast<int>(i);
      const bool already_selected =
          std::find(abundance_plot_isotopes_.begin(),
                    abundance_plot_isotopes_.end(),
                    index) != abundance_plot_isotopes_.end();
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      if (ImGui::Selectable((species_names[i] + "##abund_match_" +
                             std::to_string(i))
                                .c_str(),
                            false) &&
          !already_selected) {
        abundance_plot_isotopes_.push_back(index);
      }
      ImGui::TableSetColumnIndex(1);
      ImGui::BeginDisabled(already_selected);
      if (ImGui::SmallButton(("Add##abund_add_" + std::to_string(i)).c_str())) {
        abundance_plot_isotopes_.push_back(index);
      }
      ImGui::EndDisabled();
      ++matches;
    }
    ImGui::EndTable();
  }

  if (!abundance_plot_isotopes_.empty()) {
    ImGui::Text("Isotopes: %d",
                static_cast<int>(abundance_plot_isotopes_.size()));
    for (size_t i = 0; i < abundance_plot_isotopes_.size();) {
      const int index = abundance_plot_isotopes_[i];
      ImGui::PushID(static_cast<int>(i));
      ImGui::Text("%s", species_names[static_cast<size_t>(index)].c_str());
      ImGui::SameLine();
      if (ImGui::SmallButton("Remove")) {
        abundance_plot_isotopes_.erase(abundance_plot_isotopes_.begin() +
                                       static_cast<std::ptrdiff_t>(i));
        ImGui::PopID();
        continue;
      }
      ImGui::PopID();
      ++i;
    }
    if (ImGui::Button("Clear", ImVec2(100, 0))) {
      abundance_plot_isotopes_.clear();
    }
  }

  ImGui::Separator();
  ImGui::Checkbox("Log x", &abundance_plot_log_x_);
  ImGui::SameLine();
  ImGui::Checkbox("Log y", &abundance_plot_log_y_);
  ImGui::SameLine();
  ImGui::Checkbox("Legend", &abundance_plot_legend_);
  if (abundance_plot_log_y_) {
    ImGui::SetNextItemWidth(140.0f);
    if (ImGui::InputDouble("Log floor", &abundance_plot_log_floor_, 0.0, 0.0,
                           "%.3e")) {
      if (!std::isfinite(abundance_plot_log_floor_) ||
          abundance_plot_log_floor_ <= 0.0) {
        abundance_plot_log_floor_ = 1e-99;
      }
    }
  }

  if (cache.empty()) {
    ImGui::TextDisabled("No cached integration rows.");
    ImGui::End();
    return;
  }

  if (abundance_plot_isotopes_.empty()) {
    ImGui::TextDisabled("No isotopes selected.");
    ImGui::End();
    return;
  }

  double x_floor = 0.0;
  for (const auto &step : cache) {
    if (step.time > 0.0 && (x_floor <= 0.0 || step.time < x_floor)) {
      x_floor = step.time;
    }
  }
  if (x_floor > 0.0) {
    x_floor = std::max(x_floor * 0.1, 1e-99);
  } else {
    x_floor = 1e-99;
  }
  const double y_floor =
      (std::isfinite(abundance_plot_log_floor_) &&
       abundance_plot_log_floor_ > 0.0)
          ? abundance_plot_log_floor_
          : 1e-99;

  const ImPlotFlags plot_flags =
      abundance_plot_legend_ ? ImPlotFlags_None : ImPlotFlags_NoLegend;
  if (ImPlot::BeginPlot("X(t)##abundance_evolution", ImVec2(-1, 330),
                        plot_flags)) {
    ImPlot::SetupAxes("t (s)", "X", ImPlotAxisFlags_AutoFit,
                      ImPlotAxisFlags_AutoFit);
    if (abundance_plot_log_x_) {
      ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
    }
    if (abundance_plot_log_y_) {
      ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
    }
    if (abundance_plot_legend_) {
      ImPlot::SetupLegend(ImPlotLocation_NorthEast);
    }

    std::vector<double> xs;
    std::vector<double> ys;
    xs.reserve(cache.size());
    ys.reserve(cache.size());
    for (int isotope : abundance_plot_isotopes_) {
      xs.clear();
      ys.clear();
      for (const auto &step : cache) {
        if (isotope >= static_cast<int>(step.xnuc.size())) {
          continue;
        }
        const double raw_y = step.xnuc[static_cast<size_t>(isotope)];
        if (!std::isfinite(step.time) || !std::isfinite(raw_y)) {
          continue;
        }
        xs.push_back(abundance_plot_log_x_ ? std::max(step.time, x_floor)
                                           : step.time);
        ys.push_back(abundance_plot_log_y_ ? std::max(raw_y, y_floor)
                                           : std::max(raw_y, 0.0));
      }
      if (!xs.empty()) {
        ImPlot::PlotLine(species_names[static_cast<size_t>(isotope)].c_str(),
                         xs.data(), ys.data(), static_cast<int>(xs.size()));
      }
    }

    if (const auto *current =
            app_state_->get_trajectory_cache_step(
                app_state_->current_trajectory_step())) {
      double marker =
          abundance_plot_log_x_ ? std::max(current->time, x_floor)
                                : current->time;
      ImPlot::PlotInfLines("current", &marker, 1);
    }
    ImPlot::EndPlot();
  }

  ImGui::End();
}

void MainWindow::render_trajectory_editor_panel() {
  if (!show_trajectory_editor_) {
    return;
  }

  ImGuiViewport *viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(
      ImVec2(viewport->WorkPos.x + 720, viewport->WorkPos.y + 80),
      ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(620, 520), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSizeConstraints(ImVec2(420, 320),
                                      ImVec2(FLT_MAX, FLT_MAX));

  ImGui::Begin("Trajectory Editor", &show_trajectory_editor_,
               ImGuiWindowFlags_HorizontalScrollbar);

  if (!app_state_ || !app_state_->has_trajectory()) {
    ImGui::TextDisabled("No trajectory loaded.");
    if (ImGui::Button("Load Trajectory...", ImVec2(-1.0f, 0))) {
      show_load_traj_popup_ = true;
    }
    ImGui::End();
    return;
  }

  const auto &times = app_state_->trajectory_times();
  const auto &rhos = app_state_->trajectory_rhos();
  const auto &temps = app_state_->trajectory_temps();
  const int row_count =
      static_cast<int>(std::min({times.size(), rhos.size(), temps.size()}));
  if (row_count <= 0) {
    ImGui::TextDisabled("No trajectory rows available.");
    ImGui::End();
    return;
  }

  trajectory_editor_selected_row_ =
      std::clamp(trajectory_editor_selected_row_, 0, row_count - 1);

  ImGui::Text("Rows: %d", row_count);
  ImGui::SameLine();
  ImGui::Text("Cached rows: %d",
              static_cast<int>(app_state_->trajectory_cache().size()));
  if (!app_state_->trajectory_cache().empty()) {
    ImGui::SameLine();
    ImGui::TextDisabled("Editing clears cached integration results.");
  }

  const ImGuiTableFlags table_flags =
      ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
      ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY;
  if (ImGui::BeginTable("TrajectoryEditorTable", 4, table_flags,
                        ImVec2(0, 300))) {
    ImGui::TableSetupColumn("Row", ImGuiTableColumnFlags_WidthFixed, 70.0f);
    ImGui::TableSetupColumn("t (s)");
    ImGui::TableSetupColumn("rho (g/cm^3)");
    ImGui::TableSetupColumn("T (K)");
    ImGui::TableHeadersRow();

    for (int row = 0; row < row_count; ++row) {
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      if (ImGui::Selectable((std::to_string(row) + "##traj_edit_select_" +
                             std::to_string(row))
                                .c_str(),
                            row == trajectory_editor_selected_row_)) {
        trajectory_editor_selected_row_ = row;
      }

      double edited_time = times[static_cast<size_t>(row)];
      double edited_rho = rhos[static_cast<size_t>(row)];
      double edited_temp = temps[static_cast<size_t>(row)];

      auto apply_row_edit = [&]() {
        std::string error;
        if (app_state_->set_trajectory_row(static_cast<size_t>(row),
                                           edited_time, edited_rho,
                                           edited_temp, error)) {
          cancel_trajectory_integration_job();
          trajectory_step_ui_ =
              std::min(trajectory_step_ui_,
                       static_cast<int>(app_state_->trajectory_size() - 1));
          status_message_ = "Trajectory row updated";
        } else {
          status_message_ = "Trajectory edit rejected: " + error;
        }
      };

      ImGui::TableSetColumnIndex(1);
      ImGui::SetNextItemWidth(-1.0f);
      if (ImGui::InputDouble(("##traj_t_" + std::to_string(row)).c_str(),
                             &edited_time, 0.0, 0.0, "%.9e")) {
        apply_row_edit();
      }

      ImGui::TableSetColumnIndex(2);
      ImGui::SetNextItemWidth(-1.0f);
      if (ImGui::InputDouble(("##traj_rho_" + std::to_string(row)).c_str(),
                             &edited_rho, 0.0, 0.0, "%.9e")) {
        apply_row_edit();
      }

      ImGui::TableSetColumnIndex(3);
      ImGui::SetNextItemWidth(-1.0f);
      if (ImGui::InputDouble(("##traj_temp_" + std::to_string(row)).c_str(),
                             &edited_temp, 0.0, 0.0, "%.9e")) {
        apply_row_edit();
      }
    }

    ImGui::EndTable();
  }

  if (ImGui::Button("Add Row", ImVec2(110, 0))) {
    const double last_time = times.back();
    const double last_rho = rhos.back();
    const double last_temp = temps.back();
    const double dt =
        times.size() > 1 ? std::max(last_time - times[times.size() - 2], 0.0)
                         : 1.0;
    std::string error;
    if (app_state_->append_trajectory_row(last_time + dt, last_rho, last_temp,
                                          error)) {
      cancel_trajectory_integration_job();
      trajectory_editor_selected_row_ =
          static_cast<int>(app_state_->trajectory_size() - 1);
      status_message_ = "Trajectory row added";
    } else {
      status_message_ = "Could not add trajectory row: " + error;
    }
  }

  ImGui::SameLine();
  if (ImGui::Button("Remove Selected", ImVec2(150, 0))) {
    std::string error;
    if (app_state_->remove_trajectory_row(
            static_cast<size_t>(trajectory_editor_selected_row_), error)) {
      cancel_trajectory_integration_job();
      trajectory_editor_selected_row_ =
          std::min(trajectory_editor_selected_row_,
                   static_cast<int>(app_state_->trajectory_size() - 1));
      trajectory_step_ui_ =
          std::min(trajectory_step_ui_,
                   static_cast<int>(app_state_->trajectory_size() - 1));
      status_message_ = "Trajectory row removed";
    } else {
      status_message_ = "Could not remove trajectory row: " + error;
    }
  }

  ImGui::Separator();
  ImGui::TextWrapped("Save edited trajectory as:");
  ImGui::SetNextItemWidth(-120.0f);
  ImGui::InputText("##save_traj_path", save_traj_path_,
                   sizeof(save_traj_path_));
  ImGui::SameLine();
  if (ImGui::Button("Save As", ImVec2(100, 0))) {
    if (app_state_->save_trajectory_file(save_traj_path_)) {
      status_message_ =
          std::string("Saved trajectory: ") + save_traj_path_;
    } else {
      status_message_ =
          std::string("Failed to save trajectory: ") + save_traj_path_;
    }
  }

  ImGui::TextWrapped("Status: %s", status_message_.c_str());
  ImGui::End();
}

void MainWindow::apply_composition_preset(
    const std::vector<std::pair<std::string, double>> &composition) {
  auto &settings = app_state_->integration_settings();
  const auto &species_names = app_state_->get_species_names();

  std::fill(settings.xnuc.begin(), settings.xnuc.end(), 0.0);
  for (const auto &entry : composition) {
    auto it =
        std::find(species_names.begin(), species_names.end(), entry.first);
    if (it != species_names.end()) {
      settings
          .xnuc[static_cast<size_t>(std::distance(species_names.begin(), it))] =
          entry.second;
    }
  }
  composition_changed();
  status_message_ = "Composition preset applied";
}

void MainWindow::select_isotope(size_t index, bool append) {
  auto &selected = app_state_->integration_settings().selected;
  if (index >= selected.size()) {
    return;
  }
  if (!append) {
    std::fill(selected.begin(), selected.end(), false);
  }
  selected[index] = append ? !selected[index] : true;
  nuclide_chart_.set_selected(selected);
  if (selected[index]) {
    refresh_isotope_diagnostics(static_cast<int>(index));
  } else {
    isotope_info_index_ = -1;
    isotope_info_rate_ = 0.0;
    isotope_reactions_.clear();
    isotope_info_error_.clear();
  }
}

void MainWindow::refresh_isotope_diagnostics(int index) {
  isotope_info_index_ = index;
  isotope_info_rate_ = 0.0;
  isotope_reactions_.clear();
  isotope_info_error_.clear();

  if (!app_state_ || index < 0 || index >= app_state_->num_species()) {
    isotope_info_error_ = "No isotope selected";
    return;
  }

  try {
    isotope_info_rate_ = app_state_->get_species_rate(index);
    isotope_reactions_ = app_state_->get_species_reaction_diagnostics(index);
  } catch (const std::exception &e) {
    isotope_info_error_ = e.what();
  }
}

void MainWindow::refresh_file_browser(const std::string &dir) {
  browse_files_.clear();
  browse_dirs_.clear();

  try {
    namespace fs = std::filesystem;
    if (!fs::exists(dir) || !fs::is_directory(dir)) {
      return;
    }

    for (const auto &entry : fs::directory_iterator(dir)) {
      if (entry.is_directory()) {
        browse_dirs_.push_back(entry.path().filename().string());
      } else if (entry.is_regular_file()) {
        browse_files_.push_back(entry.path().filename().string());
      }
    }

    std::sort(browse_dirs_.begin(), browse_dirs_.end());
    std::sort(browse_files_.begin(), browse_files_.end());
  } catch (...) {
    // Handle errors silently
  }
}

bool MainWindow::render_file_browser(char *path_buffer, size_t buffer_size) {
  bool file_selected = false;

  ImGui::TextWrapped("Current directory: %s", current_browse_dir_.c_str());

  // Navigation buttons
  if (ImGui::Button("Up##dir")) {
    namespace fs = std::filesystem;
    auto parent = fs::path(current_browse_dir_).parent_path();
    if (parent != current_browse_dir_) {
      current_browse_dir_ = parent.string();
      refresh_file_browser(current_browse_dir_);
      std::memset(file_name_input_, 0, sizeof(file_name_input_));
    }
  }
  ImGui::SameLine();
  ImGui::TextWrapped("(navigate to parent)");

  ImGui::Separator();

  // Directory list
  if (ImGui::BeginChild("##dirs", ImVec2(0, 120), true)) {
    for (const auto &dir : browse_dirs_) {
      if (ImGui::Selectable(("[DIR] " + dir).c_str())) {
        current_browse_dir_ =
            (std::filesystem::path(current_browse_dir_) / dir).string();
        refresh_file_browser(current_browse_dir_);
        std::memset(file_name_input_, 0, sizeof(file_name_input_));
      }
    }
    ImGui::EndChild();
  }

  ImGui::TextWrapped("Files:");
  // File list
  if (ImGui::BeginChild("##files", ImVec2(0, 150), true)) {
    for (const auto &file : browse_files_) {
      if (ImGui::Selectable(file.c_str())) {
        std::strncpy(file_name_input_, file.c_str(),
                     sizeof(file_name_input_) - 1);
        file_name_input_[sizeof(file_name_input_) - 1] = '\0';
      }
    }
    ImGui::EndChild();
  }

  ImGui::Separator();

  // File name input
  ImGui::TextWrapped("File name:");
  ImGui::InputText("##filename", file_name_input_, sizeof(file_name_input_));

  ImGui::Separator();

  // Buttons
  if (ImGui::Button("Open##file", ImVec2(120, 0))) {
    namespace fs = std::filesystem;
    const std::string full_path =
        (fs::path(current_browse_dir_) / file_name_input_).string();
    std::error_code error;
    if (fs::is_regular_file(full_path, error)) {
      std::strncpy(path_buffer, full_path.c_str(), buffer_size - 1);
      path_buffer[buffer_size - 1] = '\0';
      file_selected = true;
    }
  }

  return file_selected;
}

void MainWindow::sync_network_path_inputs() {
  if (!app_state_) {
    return;
  }

  auto copy_path = [](char *buffer, size_t buffer_size,
                      const std::string &path) {
    std::strncpy(buffer, path.c_str(), buffer_size - 1);
    buffer[buffer_size - 1] = '\0';
  };

  copy_path(species_path_, sizeof(species_path_), app_state_->species_file());
  copy_path(reaclib_path_, sizeof(reaclib_path_), app_state_->reaclib_file());
  copy_path(partition_path_, sizeof(partition_path_),
            app_state_->partition_file());
  copy_path(mass_path_, sizeof(mass_path_), app_state_->mass_file());
  copy_path(weak_path_, sizeof(weak_path_), app_state_->weak_file());
}

void MainWindow::reset_ui_after_network_reload() {
  cancel_trajectory_integration_job();
  chart_initialized_ = false;
  nuclide_chart_.initialize(app_state_);
  chart_initialized_ = true;
  isotope_info_index_ = -1;
  isotope_info_rate_ = 0.0;
  isotope_reactions_.clear();
  isotope_info_error_.clear();
  trajectory_step_ui_ = 0;
  abundance_plot_isotopes_.clear();
}

void MainWindow::render_file_dialogs() {
#ifdef IMNET_USE_NUPPN
  if (show_open_species_popup_ &&
      !ImGui::IsPopupOpen("Open NuPPN Run Directory")) {
    ImGui::OpenPopup("Open NuPPN Run Directory");
  }
  if (show_open_species_popup_ &&
      ImGui::BeginPopupModal("Open NuPPN Run Directory", nullptr,
                             ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::TextWrapped("Run directory:");
    ImGui::InputText("##nuppn_run_dir", species_path_, sizeof(species_path_));
    ImGui::Separator();
    if (ImGui::Button("Load##nuppn_run_dir", ImVec2(120, 0))) {
      if (app_state_ && app_state_->reload_species_file(species_path_)) {
        sync_network_path_inputs();
        reset_ui_after_network_reload();
        nuppn_inputs_loaded_ = load_nuppn_input_files();
        status_message_ = std::string("Loaded NuPPN run: ") + species_path_;
      } else {
        status_message_ =
            std::string("Failed to load NuPPN run: ") + species_path_;
      }
      show_open_species_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel##nuppn_run_dir", ImVec2(120, 0))) {
      show_open_species_popup_ = false;
      sync_network_path_inputs();
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }
#else
  if (show_open_species_popup_ && !ImGui::IsPopupOpen("Open Species File")) {
    refresh_file_browser(current_browse_dir_);
    std::memset(file_name_input_, 0, sizeof(file_name_input_));
    ImGui::OpenPopup("Open Species File");
  }
  if (show_open_species_popup_ &&
      ImGui::BeginPopupModal("Open Species File", nullptr,
                             ImGuiWindowFlags_AlwaysAutoResize)) {
    if (render_file_browser(species_path_, sizeof(species_path_))) {
      if (app_state_ && app_state_->reload_species_file(species_path_)) {
        sync_network_path_inputs();
        reset_ui_after_network_reload();
        status_message_ = std::string("Reloaded species: ") + species_path_;
      } else {
        status_message_ =
            std::string("Failed to reload species: ") + species_path_;
      }
      show_open_species_popup_ = false;
      ImGui::CloseCurrentPopup();
    }

    ImGui::Separator();
    if (ImGui::Button("Cancel##open_species", ImVec2(120, 0))) {
      show_open_species_popup_ = false;
      sync_network_path_inputs();
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }

  if (show_open_rates_popup_ && !ImGui::IsPopupOpen("Open Rate Files")) {
    ImGui::OpenPopup("Open Rate Files");
  }
  if (show_open_rates_popup_ &&
      ImGui::BeginPopupModal("Open Rate Files", nullptr,
                             ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::TextWrapped("Reaclib rates:");
    ImGui::InputText("##reaclib_path", reaclib_path_, sizeof(reaclib_path_));
    ImGui::TextWrapped("Partition file:");
    ImGui::InputText("##partition_path", partition_path_,
                     sizeof(partition_path_));
    ImGui::TextWrapped("Mass file:");
    ImGui::InputText("##mass_path", mass_path_, sizeof(mass_path_));
    ImGui::TextWrapped("Weak rates:");
    ImGui::InputText("##weak_path", weak_path_, sizeof(weak_path_));

    ImGui::Separator();
    if (ImGui::Button("Reload##rates", ImVec2(120, 0))) {
      if (app_state_ && app_state_->reload_rate_files(
                            reaclib_path_, partition_path_, mass_path_,
                            weak_path_)) {
        sync_network_path_inputs();
        reset_ui_after_network_reload();
        status_message_ = "Reloaded rate files";
      } else {
        status_message_ = "Failed to reload rate files";
      }
      show_open_rates_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel##open_rates", ImVec2(120, 0))) {
      show_open_rates_popup_ = false;
      sync_network_path_inputs();
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }
#endif

  if (show_load_abund_popup_ && !ImGui::IsPopupOpen("Load Abundances File")) {
    refresh_file_browser(current_browse_dir_);
    std::memset(file_name_input_, 0, sizeof(file_name_input_));
    ImGui::OpenPopup("Load Abundances File");
  }
  if (show_load_abund_popup_ &&
      ImGui::BeginPopupModal("Load Abundances File", nullptr,
                             ImGuiWindowFlags_AlwaysAutoResize)) {
    if (render_file_browser(load_abund_path_, sizeof(load_abund_path_))) {
      if (app_state_ &&
          app_state_->load_abundances_from_file(load_abund_path_)) {
        cancel_trajectory_integration_job();
        status_message_ = std::string("Loaded abundances: ") + load_abund_path_;
      } else {
        status_message_ =
            std::string("Failed to load abundances: ") + load_abund_path_;
      }
      show_load_abund_popup_ = false;
      ImGui::CloseCurrentPopup();
    }

    ImGui::Separator();
    if (ImGui::Button("Cancel##load_abund", ImVec2(120, 0))) {
      show_load_abund_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }

  if (show_save_abund_popup_ && !ImGui::IsPopupOpen("Save Abundances File")) {
    ImGui::OpenPopup("Save Abundances File");
  }
  if (show_save_abund_popup_ &&
      ImGui::BeginPopupModal("Save Abundances File", nullptr,
                             ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::TextWrapped("Save abundances to:");
    ImGui::InputText("##save_abund_path", save_abund_path_,
                     sizeof(save_abund_path_));

    if (ImGui::Button("Save##abund", ImVec2(120, 0))) {
      if (app_state_ && app_state_->save_abundances_to_file(save_abund_path_)) {
        status_message_ = std::string("Saved abundances: ") + save_abund_path_;
      } else {
        status_message_ =
            std::string("Failed to save abundances: ") + save_abund_path_;
      }
      show_save_abund_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel##save_abund", ImVec2(120, 0))) {
      show_save_abund_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }

  if (show_save_state_popup_ && !ImGui::IsPopupOpen("Save State File")) {
    ImGui::OpenPopup("Save State File");
  }
  if (show_save_state_popup_ &&
      ImGui::BeginPopupModal("Save State File", nullptr,
                             ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::TextWrapped("Save portable analysis state to:");
    ImGui::InputText("##save_state_path", save_state_path_,
                     sizeof(save_state_path_));

    if (ImGui::Button("Save##state", ImVec2(120, 0))) {
      if (app_state_ && app_state_->save_state_to_file(save_state_path_)) {
        status_message_ = std::string("Saved state: ") + save_state_path_;
      } else {
        status_message_ = std::string("Failed to save state: ") +
                          save_state_path_;
      }
      show_save_state_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel##save_state", ImVec2(120, 0))) {
      show_save_state_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }

  if (show_load_traj_popup_ && !ImGui::IsPopupOpen("Load Trajectory File")) {
    refresh_file_browser(current_browse_dir_);
    std::memset(file_name_input_, 0, sizeof(file_name_input_));
    ImGui::OpenPopup("Load Trajectory File");
  }
  if (show_load_traj_popup_ &&
      ImGui::BeginPopupModal("Load Trajectory File", nullptr,
                             ImGuiWindowFlags_AlwaysAutoResize)) {
    if (render_file_browser(load_traj_path_, sizeof(load_traj_path_))) {
      if (app_state_ && app_state_->load_trajectory_file(load_traj_path_)) {
        cancel_trajectory_integration_job();
        trajectory_step_ui_ = 0;
        status_message_ = std::string("Loaded trajectory: ") + load_traj_path_;
      } else {
        status_message_ =
            std::string("Failed to load trajectory: ") + load_traj_path_;
      }
      show_load_traj_popup_ = false;
      ImGui::CloseCurrentPopup();
    }

    ImGui::Separator();
    if (ImGui::Button("Cancel##load_traj", ImVec2(120, 0))) {
      show_load_traj_popup_ = false;
      ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
  }
}

void MainWindow::render_selected_abundance_editor() {
  if (!app_state_) {
    return;
  }

  auto &settings = app_state_->integration_settings();
  const auto &species_names = app_state_->get_species_names();

  int selected_count = 0;
  for (bool s : settings.selected) {
    if (s) {
      ++selected_count;
    }
  }

  ImGui::Text("Selected isotopes: %d", selected_count);
  if (selected_count == 0) {
    ImGui::TextDisabled("No selected isotopes.");
    return;
  }

  ImGui::BeginChild("SelectedAbundanceEditor", ImVec2(0, 140), true);
  int rendered = 0;
  for (size_t i = 0; i < settings.selected.size(); ++i) {
    if (!settings.selected[i]) {
      continue;
    }

    if (i >= settings.xnuc.size() || i >= settings.fixed.size() ||
        i >= species_names.size()) {
      continue;
    }

    double x = settings.xnuc[i];
    std::string label = species_names[i] + "##xnuc_" + std::to_string(i);
    ImGui::PushItemWidth(180.0f);
    if (ImGui::InputDouble(label.c_str(), &x, 0.0, 0.0, "%.6e")) {
      settings.xnuc[i] = std::isfinite(x) ? std::max(0.0, x) : 0.0;
      composition_changed();
    }
    ImGui::PopItemWidth();

    ImGui::SameLine();
    std::string fixed_label = std::string("fixed##fixed_") + std::to_string(i);
    bool fixed_value = settings.fixed[i];
    if (ImGui::Checkbox(fixed_label.c_str(), &fixed_value)) {
      settings.fixed[i] = fixed_value;
    }

    ++rendered;
    if (rendered >= 12) {
      ImGui::TextDisabled("Showing first 12 selected isotopes...");
      break;
    }
  }
  ImGui::EndChild();
}

#ifdef IMNET_USE_NUPPN
bool MainWindow::load_nuppn_input_files() {
  if (!app_state_) {
    return false;
  }

  auto load_file = [](const std::filesystem::path &path, char *buffer,
                      size_t buffer_size) {
    std::ifstream file(path);
    if (!file.is_open()) {
      buffer[0] = '\0';
      return false;
    }
    const std::string text((std::istreambuf_iterator<char>(file)),
                           std::istreambuf_iterator<char>());
    const size_t copied = std::min(text.size(), buffer_size - 1);
    std::memcpy(buffer, text.data(), copied);
    buffer[copied] = '\0';
    return text.size() < buffer_size;
  };

  const std::filesystem::path run_dir = app_state_->species_file();
  return load_file(run_dir / "ppn_frame.input", nuppn_frame_input_,
                   sizeof(nuppn_frame_input_)) &&
         load_file(run_dir / "ppn_physics.input", nuppn_physics_input_,
                   sizeof(nuppn_physics_input_)) &&
         load_file(run_dir / "ppn_solver.input", nuppn_solver_input_,
                   sizeof(nuppn_solver_input_));
}

bool MainWindow::save_nuppn_input_files() {
  if (!app_state_) {
    return false;
  }

  auto save_file = [](const std::filesystem::path &path, const char *buffer) {
    std::ofstream file(path);
    if (!file.is_open()) {
      return false;
    }
    file << buffer;
    return static_cast<bool>(file);
  };

  const std::filesystem::path run_dir = app_state_->species_file();
  return save_file(run_dir / "ppn_frame.input", nuppn_frame_input_) &&
         save_file(run_dir / "ppn_physics.input", nuppn_physics_input_) &&
         save_file(run_dir / "ppn_solver.input", nuppn_solver_input_);
}
#endif

void MainWindow::composition_changed() {
  cancel_trajectory_integration_job();
  if (app_state_) {
    app_state_->clear_trajectory_cache();
  }
}

} // namespace imyann
