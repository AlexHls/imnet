#ifndef IMYANN_APP_STATE_H
#define IMYANN_APP_STATE_H

#include "network_wrapper.h"
#include <memory>
#include <string>
#include <vector>

namespace imyann {

/**
 * @struct IntegrationSettings
 * @brief Settings for network integration
 */
struct IntegrationSettings {
  double rho = 1e4;           ///< Density (g/cm³)
  double temp = 1e7;          ///< Temperature (K)
  double dt = 1.0;            ///< Fixed step or initial timestep (seconds)
  double final_time = 1.0;    ///< Final-time integration target (seconds)
  double dt_max = 1.0;        ///< Maximum final-time timestep (seconds)
  double dt_factor = 1.0;     ///< Final-time timestep growth factor
  int max_steps = 10000;      ///< Maximum final-time integration steps
  std::vector<double> xnuc;   ///< Species abundances
  std::vector<bool> fixed;    ///< Which species are fixed
  std::vector<bool> active;   ///< Which species are active in the network
  std::vector<bool> selected; ///< Which species are selected for display
};

/**
 * @struct ViewSettings
 * @brief UI view and display settings
 */
struct ViewSettings {
  float chart_zoom = 1.0f;        ///< Zoom level for nuclide chart
  bool show_unconnected = false;  ///< Highlight unconnected species
  bool show_abundance = false;    ///< Color-code by abundance
  bool show_fluxes = false;       ///< Draw reaction flux arrows
  bool show_regular_fluxes = true; ///< Draw regular reaction arrows
  bool show_weak_fluxes = true;    ///< Draw weak reaction arrows
  int flux_arrow_metric = 0;      ///< 0 = |dY/dt|, 1 = |dX/dt|, 2 = rate
  double flux_threshold = 1e-30;  ///< Minimum flux strength to display
  int max_flux_arrows = 80;       ///< Maximum displayed arrows per type
  int flux_color_mode = 0;        ///< 0 = flat, 1 = linear, 2 = log flux
  int flux_colormap = 0;          ///< Colormap for scaled flux arrows
  int flux_colorbar_position = 3; ///< Preset colorbar position
  float flux_color[4] = {1.00f, 0.22f, 0.64f, 0.95f};
  float weak_flux_color[4] = {0.50f, 1.00f, 0.36f, 0.95f};
  int max_cached_trajectory_steps = 100; ///< Max trajectory rows to cache
};

struct TrajectoryStepCache {
  size_t index = 0;
  double time = 0.0;
  double rho = 0.0;
  double temp = 0.0;
  double dt = 0.0;
  double dedt = 0.0;
  int substeps = 0;
  bool success = true;
  std::string status;
  std::string error;
  std::vector<double> xnuc;
};

/**
 * @class AppState
 * @brief Global application state container
 *
 * Holds the current network, integration settings, and UI state.
 * Provides a single source of truth for the application.
 */
class AppState {
public:
  AppState() = default;
  ~AppState() = default;

  // Deleted copy/move
  AppState(const AppState &) = delete;
  AppState &operator=(const AppState &) = delete;

  // ====================================================================
  // Network management
  // ====================================================================

  /**
   * @brief Initialize network from rate files
   * @return true if successful, false otherwise
   */
  bool initialize_network(const std::string &species_file,
                          const std::string &reaclib_file,
                          const std::string &partition_file,
                          const std::string &mass_file,
                          const std::string &weak_file);

  bool reload_species_file(const std::string &species_file);
  bool reload_rate_files(const std::string &reaclib_file,
                         const std::string &partition_file,
                         const std::string &mass_file,
                         const std::string &weak_file);

  const std::string &species_file() const { return species_file_; }
  const std::string &reaclib_file() const { return reaclib_file_; }
  const std::string &partition_file() const { return partition_file_; }
  const std::string &mass_file() const { return mass_file_; }
  const std::string &weak_file() const { return weak_file_; }

  /**
   * @brief Get current network (may be nullptr)
   */
  Network *get_network() { return network_.get(); }
  const Network *get_network() const { return network_.get(); }

  /**
   * @brief Check if network is initialized
   */
  bool has_network() const { return network_ != nullptr; }

  // ====================================================================
  // Species and abundance management
  // ====================================================================

  /**
   * @brief Get species names
   */
  const std::vector<std::string> &get_species_names() const {
    return species_names_;
  }

  /**
   * @brief Get number of species
   */
  int num_species() const { return static_cast<int>(species_names_.size()); }

  const std::vector<Species> &get_species() const { return species_data_; }

  /**
   * @brief Reset abundances to initial state
   */
  void reset_abundances();

  /**
   * @brief Normalize abundances (sum to 1)
   */
  bool normalize_abundances();

  // ====================================================================
  // Integration and computation
  // ====================================================================

  /**
   * @brief Perform network integration
   * @return Energy release rate (dE/dt in erg/g/s)
   */
  double integrate();

  /**
   * @brief Validate the current single-step integration state
   */
  bool validate_integration_settings(std::string &error) const;

  /**
   * @brief Perform one explicit network integration step
   */
  bool integrate_single_step(bool normalize_before, double &dedt,
                             std::string &error);

  /**
   * @brief Integrate at fixed rho/T until final_time
   */
  bool run_to_time(bool normalize_before, std::string &error);

  /**
   * @brief Compute NSE abundances
   */
  void compute_nse();

  /**
   * @brief Compute and apply NSE abundances
   */
  bool compute_nse(std::string &error);

  // ====================================================================
  // File I/O
  // ====================================================================

  /**
   * @brief Load abundances from file
   */
  bool load_abundances_from_file(const std::string &filename);

  /**
   * @brief Save abundances to file
   */
  bool save_abundances_to_file(const std::string &filename);

  /**
   * @brief Load trajectory (for v2 multi-timestep)
   */
  bool load_trajectory_file(const std::string &filename);

  /**
   * @brief Save current loaded trajectory
   */
  bool save_trajectory_file(const std::string &filename) const;

  /**
   * @brief Save the current scientific state for external analysis
   */
  bool save_state_to_file(const std::string &filename) const;

  /**
   * @brief Run full loaded trajectory and cache compositions
   * @return true if successful
   */
  bool run_trajectory();

  bool run_trajectory(std::string &error);

  /**
   * @brief Check whether a trajectory is loaded
   */
  bool has_trajectory() const { return !trajectory_times_.empty(); }

  /**
   * @brief Number of timesteps in loaded trajectory
   */
  size_t trajectory_size() const { return trajectory_times_.size(); }

  /**
   * @brief Current trajectory index
   */
  size_t current_trajectory_step() const { return trajectory_index_; }

  /**
   * @brief Select a cached trajectory step and apply it to settings
   */
  bool set_trajectory_step(size_t step);

  /**
   * @brief Edit a loaded trajectory row and invalidate cached results
   */
  bool set_trajectory_row(size_t step, double time, double rho, double temp,
                          std::string &error);

  /**
   * @brief Append a new trajectory row and invalidate cached results
   */
  bool append_trajectory_row(double time, double rho, double temp,
                             std::string &error);

  /**
   * @brief Remove a loaded trajectory row and invalidate cached results
   */
  bool remove_trajectory_row(size_t step, std::string &error);

  /**
   * @brief Get cached trajectory times
   */
  const std::vector<double> &trajectory_times() const {
    return trajectory_times_;
  }

  /**
   * @brief Get loaded trajectory densities
   */
  const std::vector<double> &trajectory_rhos() const {
    return trajectory_rhos_;
  }

  /**
   * @brief Get loaded trajectory temperatures
   */
  const std::vector<double> &trajectory_temps() const {
    return trajectory_temps_;
  }

  /**
   * @brief Last integrated trajectory dE/dt values
   */
  const std::vector<double> &trajectory_dedt() const {
    return trajectory_dedt_;
  }

  const std::vector<TrajectoryStepCache> &trajectory_cache() const {
    return trajectory_cache_;
  }

  const TrajectoryStepCache *get_trajectory_cache_step(size_t step) const;

  /**
   * @brief Replace cached trajectory rows
   */
  void set_trajectory_cache(std::vector<TrajectoryStepCache> cache,
                            size_t current_step);

  /**
   * @brief Clear cached trajectory integration results
   */
  void clear_trajectory_cache();

  // ====================================================================
  // Settings access (const-correct)
  // ====================================================================

  IntegrationSettings &integration_settings() { return integration_settings_; }
  const IntegrationSettings &integration_settings() const {
    return integration_settings_;
  }

  ViewSettings &view_settings() { return view_settings_; }
  const ViewSettings &view_settings() const { return view_settings_; }

  // ====================================================================
  // Computed properties
  // ====================================================================

  /**
   * @brief Calculate electron fraction (Ye = <Z>/<A>)
   */
  double calculate_ye() const;

  /**
   * @brief Calculate average mass number (<A>)
   */
  double calculate_abar() const;

  /**
   * @brief Calculate average proton number (<Z>)
   */
  double calculate_zbar() const;

  /**
   * @brief Estimate dX/dt for a species index at current state
   */
  double get_species_rate(int species_index) const;

  std::vector<ReactionDiagnostic>
  get_species_reaction_diagnostics(int species_index,
                                   size_t max_count = 12) const;

  std::vector<ReactionFlux> get_reaction_fluxes(double min_strength,
                                                int metric,
                                                size_t max_count = 100,
                                                bool include_regular = true,
                                                bool include_weak = true) const;
  bool has_reaction_flux_state() const { return flux_state_valid_; }
  double reaction_flux_state_rho() const { return flux_state_rho_; }
  double reaction_flux_state_temp() const { return flux_state_temp_; }

  /**
   * @brief Check network connectivity for a species against active mask
   */
  bool is_species_connected(int species_index) const;

private:
  std::unique_ptr<Network> network_;
  std::vector<std::string> species_names_;
  std::vector<Species> species_data_;
  std::vector<double> initial_xnuc_; ///< Backup of initial abundances
  std::string species_file_;
  std::string reaclib_file_;
  std::string partition_file_;
  std::string mass_file_;
  std::string weak_file_;

  // v2 trajectory data and cache
  std::vector<double> trajectory_times_;
  std::vector<double> trajectory_rhos_;
  std::vector<double> trajectory_temps_;
  std::vector<double> trajectory_dedt_;
  std::vector<TrajectoryStepCache> trajectory_cache_;
  size_t trajectory_index_ = 0;
  double last_dedt_ = 0.0;
  bool last_success_ = true;
  std::string last_status_ = "current";
  std::string last_error_;

  void invalidate_trajectory_cache();
  bool validate_composition(std::string &error) const;
  bool validate_trajectory_arrays(const std::vector<double> &times,
                                  const std::vector<double> &rhos,
                                  const std::vector<double> &temps,
                                  std::string &error) const;
  void record_reaction_flux_state(double rho, double temp,
                                  const std::vector<double> &xnuc);
  void invalidate_reaction_flux_state();

  IntegrationSettings integration_settings_;
  ViewSettings view_settings_;
  bool flux_state_valid_ = false;
  double flux_state_rho_ = 0.0;
  double flux_state_temp_ = 0.0;
  std::vector<double> flux_state_xnuc_;
};

} // namespace imyann

#endif // IMYANN_APP_STATE_H
