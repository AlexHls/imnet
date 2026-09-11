#include "app_state.h"
#include "file_io.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace imyann {
namespace {

void write_json_string(std::ostream &out, const std::string &value) {
  out << '"';
  for (unsigned char c : value) {
    switch (c) {
    case '"':
      out << "\\\"";
      break;
    case '\\':
      out << "\\\\";
      break;
    case '\b':
      out << "\\b";
      break;
    case '\f':
      out << "\\f";
      break;
    case '\n':
      out << "\\n";
      break;
    case '\r':
      out << "\\r";
      break;
    case '\t':
      out << "\\t";
      break;
    default:
      if (c < 0x20) {
        out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
            << static_cast<int>(c) << std::dec << std::setfill(' ');
      } else {
        out << static_cast<char>(c);
      }
    }
  }
  out << '"';
}

void validate_step_composition(const std::vector<double> &xnuc,
                               const std::vector<std::string> &species_names) {
  if (xnuc.size() != species_names.size()) {
    throw std::runtime_error(
        "Abundance vector size does not match species count");
  }
  for (size_t i = 0; i < xnuc.size(); ++i) {
    if (!std::isfinite(xnuc[i]) || xnuc[i] < 0.0) {
      const std::string name =
          i < species_names.size() ? species_names[i] : std::to_string(i);
      throw std::runtime_error("Invalid abundance for " + name);
    }
  }
}

} // namespace

bool AppState::initialize_network(const std::string &species_file,
                                  const std::string &reaclib_file,
                                  const std::string &partition_file,
                                  const std::string &mass_file,
                                  const std::string &weak_file) {
  try {
    const auto old_names = species_names_;
    const auto old_settings = integration_settings_;
    const auto old_initial_xnuc = initial_xnuc_;

    // Create network and derive species metadata from initialized C core
    auto new_network = std::make_unique<Network>(
        species_file, reaclib_file, partition_file, mass_file, weak_file);

    auto new_species_data = new_network->get_species();
    std::vector<std::string> new_species_names;
    new_species_names.reserve(new_species_data.size());
    for (const auto &s : new_species_data) {
      new_species_names.push_back(s.name);
    }

    if (new_species_names.empty()) {
      std::cerr << "Error: No species loaded from " << species_file
                << std::endl;
      return false;
    }

    std::cout << "Loaded " << new_species_names.size() << " species from "
              << species_file << std::endl;

    IntegrationSettings new_settings;
    new_settings.rho = old_settings.rho;
    new_settings.temp = old_settings.temp;
    new_settings.dt = old_settings.dt;
    new_settings.final_time = old_settings.final_time;
    new_settings.dt_max = old_settings.dt_max;
    new_settings.dt_factor = old_settings.dt_factor;
    new_settings.max_steps = old_settings.max_steps;
    new_settings.xnuc.assign(new_species_names.size(), 0.0);
    new_settings.fixed.assign(new_species_names.size(), false);
    new_settings.active.assign(new_species_names.size(), true);
    new_settings.selected.assign(new_species_names.size(), false);

    std::vector<double> new_initial_xnuc(new_species_names.size(), 0.0);
    std::unordered_map<std::string, size_t> old_species_indices;
    for (size_t i = 0; i < old_names.size(); ++i) {
      old_species_indices.emplace(old_names[i], i);
    }

    for (size_t i = 0; i < new_species_names.size(); ++i) {
      const auto it = old_species_indices.find(new_species_names[i]);
      if (it == old_species_indices.end()) {
        continue;
      }
      const size_t old_i = it->second;
      if (old_i < old_settings.xnuc.size()) {
        new_settings.xnuc[i] = old_settings.xnuc[old_i];
      }
      if (old_i < old_settings.fixed.size()) {
        new_settings.fixed[i] = old_settings.fixed[old_i];
      }
      if (old_i < old_settings.active.size()) {
        new_settings.active[i] = old_settings.active[old_i];
      }
      if (old_i < old_settings.selected.size()) {
        new_settings.selected[i] = old_settings.selected[old_i];
      }
      if (old_i < old_initial_xnuc.size()) {
        new_initial_xnuc[i] = old_initial_xnuc[old_i];
      }
    }

    network_ = std::move(new_network);
    species_data_ = std::move(new_species_data);
    species_names_ = std::move(new_species_names);
    integration_settings_ = std::move(new_settings);
    initial_xnuc_ = std::move(new_initial_xnuc);
    species_file_ = species_file;
    reaclib_file_ = reaclib_file;
    partition_file_ = partition_file;
    mass_file_ = mass_file;
    weak_file_ = weak_file;
    invalidate_trajectory_cache();

    std::cout << "Network initialized successfully" << std::endl;
    return true;

  } catch (const std::exception &e) {
    std::cerr << "Error initializing network: " << e.what() << std::endl;
    return false;
  }
}

bool AppState::reload_species_file(const std::string &species_file) {
  return initialize_network(species_file, reaclib_file_, partition_file_,
                            mass_file_, weak_file_);
}

bool AppState::reload_rate_files(const std::string &reaclib_file,
                                 const std::string &partition_file,
                                 const std::string &mass_file,
                                 const std::string &weak_file) {
  return initialize_network(species_file_, reaclib_file, partition_file,
                            mass_file, weak_file);
}

void AppState::reset_abundances() {
  integration_settings_.xnuc = initial_xnuc_;
  integration_settings_.fixed.assign(integration_settings_.xnuc.size(), false);
  invalidate_trajectory_cache();
}

bool AppState::normalize_abundances() {
  auto &xnuc = integration_settings_.xnuc;
  auto &fixed = integration_settings_.fixed;
  if (xnuc.size() != fixed.size()) {
    std::cerr << "Warning: Abundance and fixed-mask sizes do not match"
              << std::endl;
    return false;
  }

  double fixed_sum = 0.0;
  double unfixed_sum = 0.0;
  for (size_t i = 0; i < xnuc.size(); ++i) {
    if (!std::isfinite(xnuc[i]) || xnuc[i] < 0.0) {
      std::cerr << "Warning: Cannot normalize invalid abundances" << std::endl;
      return false;
    }
    if (fixed[i]) {
      fixed_sum += xnuc[i];
    } else {
      unfixed_sum += xnuc[i];
    }
  }

  constexpr double tolerance = 1e-12;
  if (fixed_sum > 1.0 + tolerance) {
    std::cerr << "Warning: Fixed abundances sum to " << fixed_sum
              << ", cannot normalize" << std::endl;
    return false;
  }

  if (unfixed_sum == 0.0 && fixed_sum < 1.0 - tolerance) {
    std::cerr << "Warning: No unfixed abundance is available to normalize"
              << std::endl;
    return false;
  }

  const double scale =
      unfixed_sum > 0.0 ? std::max(0.0, 1.0 - fixed_sum) / unfixed_sum : 0.0;
  for (size_t i = 0; i < xnuc.size(); ++i) {
    if (!fixed[i]) {
      xnuc[i] *= scale;
    }
  }
  invalidate_trajectory_cache();
  return true;
}

double AppState::integrate() {
  double dedt = 0.0;
  std::string error;
  if (integrate_single_step(true, dedt, error)) {
    return dedt;
  }

  std::cerr << "Error during integration: " << error << std::endl;
  return 0.0;
}

bool AppState::validate_integration_settings(std::string &error) const {
  if (!network_) {
    error = "Network not initialized";
    return false;
  }
  if (!validate_composition(error)) {
    return false;
  }
  if (integration_settings_.fixed.size() != species_names_.size() ||
      integration_settings_.active.size() != species_names_.size() ||
      integration_settings_.selected.size() != species_names_.size()) {
    error = "Species state vector size does not match species count";
    return false;
  }
  if (!std::isfinite(integration_settings_.rho) ||
      integration_settings_.rho <= 0.0) {
    error = "Density must be finite and positive";
    return false;
  }
  if (!std::isfinite(integration_settings_.temp) ||
      integration_settings_.temp <= 0.0) {
    error = "Temperature must be finite and positive";
    return false;
  }
  if (!std::isfinite(integration_settings_.dt) ||
      integration_settings_.dt < 0.0) {
    error = "Timestep must be finite and non-negative";
    return false;
  }

  error.clear();
  return true;
}

bool AppState::integrate_single_step(bool normalize_before, double &dedt,
                                     std::string &error) {
  dedt = 0.0;
  last_dedt_ = 0.0;
  last_success_ = false;
  last_status_ = "failed";
  last_error_.clear();
  if (!validate_integration_settings(error)) {
    last_error_ = error;
    return false;
  }

  try {
    if (normalize_before && !normalize_abundances()) {
      error = "Could not normalize abundances";
      last_error_ = error;
      return false;
    }
    if (!validate_integration_settings(error)) {
      last_error_ = error;
      return false;
    }
    dedt = network_->integrate(integration_settings_.rho,
                               integration_settings_.temp,
                               integration_settings_.xnuc,
                               integration_settings_.dt);
    last_dedt_ = dedt;
    last_success_ = true;
    last_status_ = "ok";
    invalidate_trajectory_cache();
    record_reaction_flux_state(integration_settings_.rho,
                               integration_settings_.temp,
                               integration_settings_.xnuc);
    std::cout << "Integration successful, dE/dt = " << dedt << std::endl;
    return true;
  } catch (const std::exception &e) {
    error = e.what();
    last_error_ = error;
    return false;
  }
}

bool AppState::run_to_time(bool normalize_before, std::string &error) {
  error.clear();
  last_dedt_ = 0.0;
  last_success_ = false;
  last_status_ = "failed";
  last_error_.clear();

  if (normalize_before && !normalize_abundances()) {
    error = "Could not normalize abundances";
    last_error_ = error;
    return false;
  }
  if (!validate_integration_settings(error)) {
    last_error_ = error;
    return false;
  }

  const double rho = integration_settings_.rho;
  const double temp = integration_settings_.temp;
  const double final_time = integration_settings_.final_time;
  const double dt_max = integration_settings_.dt_max;
  const double dt_factor = integration_settings_.dt_factor;
  const int max_steps = integration_settings_.max_steps;
  if (!std::isfinite(final_time) || final_time <= 0.0) {
    error = "Final time must be finite and positive";
  } else if (integration_settings_.dt <= 0.0) {
    error = "Initial timestep must be positive";
  } else if (!std::isfinite(dt_max) || dt_max <= 0.0) {
    error = "Maximum timestep must be finite and positive";
  } else if (!std::isfinite(dt_factor) || dt_factor <= 0.0) {
    error = "Timestep factor must be finite and positive";
  } else if (max_steps <= 0) {
    error = "Max steps must be positive";
  }
  if (!error.empty()) {
    last_error_ = error;
    return false;
  }

  const double initial_dt = integration_settings_.dt;
  const std::vector<double> initial_xnuc = integration_settings_.xnuc;
  trajectory_times_.clear();
  trajectory_rhos_.clear();
  trajectory_temps_.clear();
  trajectory_dedt_.clear();
  trajectory_cache_.clear();
  invalidate_reaction_flux_state();
  trajectory_dedt_.reserve(static_cast<size_t>(max_steps) + 1);
  trajectory_cache_.reserve(static_cast<size_t>(max_steps) + 1);

  auto append_step = [&](size_t index, double time, double row_rho,
                         double row_temp, double row_dt, double dedt,
                         int substeps, bool success, const std::string &status,
                         const std::string &step_error,
                         const std::vector<double> &xnuc) {
    trajectory_times_.push_back(time);
    trajectory_rhos_.push_back(row_rho);
    trajectory_temps_.push_back(row_temp);
    trajectory_dedt_.push_back(dedt);

    TrajectoryStepCache step;
    step.index = index;
    step.time = time;
    step.rho = row_rho;
    step.temp = row_temp;
    step.dt = row_dt;
    step.dedt = dedt;
    step.substeps = substeps;
    step.success = success;
    step.status = status;
    step.error = step_error;
    step.xnuc = xnuc;
    trajectory_cache_.push_back(std::move(step));
  };

  auto append_history = [&]() {
    for (const auto &snapshot : network_->last_step_history()) {
      append_step(snapshot.index, snapshot.time, snapshot.rho, snapshot.temp,
                  snapshot.dt, snapshot.dedt, snapshot.substeps, true,
                  snapshot.index == 0 ? "initial" : "ok", "", snapshot.xnuc);
    }
  };

  try {
    double dedt = network_->integrate_to_time(
        rho, temp, integration_settings_.xnuc, final_time, initial_dt, dt_max,
        dt_factor, max_steps);
    validate_step_composition(integration_settings_.xnuc, species_names_);
    trajectory_times_.clear();
    trajectory_rhos_.clear();
    trajectory_temps_.clear();
    trajectory_dedt_.clear();
    trajectory_cache_.clear();
    append_history();
    last_dedt_ = dedt;
    last_success_ = true;
    last_status_ = "ok";
    if (trajectory_cache_.empty()) {
      append_step(0, 0.0, rho, temp, initial_dt, 0.0, 0, true, "initial", "",
                  initial_xnuc);
      append_step(1, final_time, rho, temp, initial_dt, dedt,
                  network_->last_substeps(), true, "ok", "",
                  integration_settings_.xnuc);
    }
  } catch (const std::exception &e) {
    error = e.what();
    if (network_ && !network_->last_step_history().empty()) {
      trajectory_times_.clear();
      trajectory_rhos_.clear();
      trajectory_temps_.clear();
      trajectory_dedt_.clear();
      trajectory_cache_.clear();
      append_history();
      trajectory_cache_.back().success = false;
      trajectory_cache_.back().status = "failed";
      trajectory_cache_.back().error = error;
    } else {
      integration_settings_.xnuc = initial_xnuc;
      append_step(0, 0.0, rho, temp, initial_dt, 0.0, 0, true, "initial", "",
                  initial_xnuc);
      append_step(1, final_time, rho, temp, initial_dt, last_dedt_,
                  network_->last_substeps(), false, "failed", error,
                  integration_settings_.xnuc);
    }
    trajectory_index_ =
        trajectory_cache_.empty() ? 0 : trajectory_cache_.back().index;
    set_trajectory_step(trajectory_index_);
    integration_settings_.dt = initial_dt;
    last_error_ = error;
    return false;
  }

  trajectory_index_ =
      trajectory_cache_.empty() ? 0 : trajectory_cache_.back().index;
  set_trajectory_step(trajectory_index_);
  integration_settings_.dt = initial_dt;
  last_error_.clear();
  std::cout << "Final-time integration complete: "
            << trajectory_cache_.size() << " steps cached" << std::endl;
  error.clear();
  return true;
}

void AppState::compute_nse() {
  std::string error;
  if (!compute_nse(error)) {
    std::cerr << "Error computing NSE: " << error << std::endl;
  }
}

bool AppState::compute_nse(std::string &error) {
  if (!network_) {
    error = "Network not initialized";
    return false;
  }

  try {
    double ye = calculate_ye();
    if (!std::isfinite(ye) || ye < 0.0 || ye > 1.0) {
      error = "Current composition has invalid Ye";
      return false;
    }
    std::vector<double> nse_xnuc(integration_settings_.xnuc.size());

    network_->compute_nse(integration_settings_.rho, integration_settings_.temp,
                          ye, nse_xnuc);

    integration_settings_.xnuc = std::move(nse_xnuc);
    invalidate_trajectory_cache();
    std::cout << "NSE computation complete" << std::endl;
    error.clear();
    return true;
  } catch (const std::exception &e) {
    error = e.what();
    return false;
  }
}

bool AppState::load_abundances_from_file(const std::string &filename) {
  std::vector<double> xnuc;
  if (!load_abundances(filename, species_names_, xnuc)) {
    return false;
  }
  integration_settings_.xnuc = xnuc;
  initial_xnuc_ = std::move(xnuc);
  invalidate_trajectory_cache();
  return true;
}

bool AppState::save_abundances_to_file(const std::string &filename) {
  return save_abundances(filename, species_names_, integration_settings_.xnuc);
}

bool AppState::load_trajectory_file(const std::string &filename) {
  std::vector<double> times;
  std::vector<double> rhos;
  std::vector<double> temps;

  if (!load_trajectory(filename, times, rhos, temps)) {
    return false;
  }

  trajectory_times_ = std::move(times);
  trajectory_rhos_ = std::move(rhos);
  trajectory_temps_ = std::move(temps);
  invalidate_trajectory_cache();
  trajectory_index_ = 0;

  if (!trajectory_rhos_.empty()) {
    integration_settings_.rho = trajectory_rhos_[0];
  }
  if (!trajectory_temps_.empty()) {
    integration_settings_.temp = trajectory_temps_[0];
  }
  if (trajectory_times_.size() > 1) {
    integration_settings_.dt =
        std::max(trajectory_times_[1] - trajectory_times_[0], 0.0);
  } else {
    integration_settings_.dt = 0.0;
  }

  return true;
}

bool AppState::save_trajectory_file(const std::string &filename) const {
  return save_trajectory(filename, trajectory_times_, trajectory_rhos_,
                         trajectory_temps_);
}

bool AppState::save_state_to_file(const std::string &filename) const {
  std::string error;
  if (!network_ || !validate_composition(error) ||
      species_data_.size() != species_names_.size()) {
    std::cerr << "Error: Cannot save state: "
              << (error.empty() ? "network metadata is incomplete" : error)
              << std::endl;
    return false;
  }

  std::vector<TrajectoryStepCache> steps = trajectory_cache_;
  const bool current_is_cached =
      std::any_of(steps.begin(), steps.end(), [this](const auto &step) {
        return step.index == trajectory_index_;
      });
  if (steps.empty() || !current_is_cached) {
    TrajectoryStepCache current;
    current.index = trajectory_index_;
    current.time =
        trajectory_index_ < trajectory_times_.size()
            ? trajectory_times_[trajectory_index_]
            : 0.0;
    current.rho = integration_settings_.rho;
    current.temp = integration_settings_.temp;
    current.dt = integration_settings_.dt;
    current.dedt = last_dedt_;
    current.substeps = network_ ? network_->last_substeps() : 0;
    current.success = last_success_;
    current.status = last_status_;
    current.error = last_error_;
    current.xnuc = integration_settings_.xnuc;
    steps.push_back(std::move(current));
  }
  std::sort(steps.begin(), steps.end(),
            [](const auto &a, const auto &b) { return a.index < b.index; });

  std::ofstream out(filename);
  if (!out) {
    std::cerr << "Error: Could not open state output file: " << filename
              << std::endl;
    return false;
  }
  out << std::setprecision(17) << std::scientific;
  out << "{\n  \"format\": \"imnet-state-v1\",\n"
      << "  \"units\": {\"time\": \"s\", \"rho\": \"g/cm^3\", "
         "\"temp\": \"K\", \"dedt\": \"erg/g/s\", "
         "\"abundance\": \"mass_fraction\", "
         "\"flux_dydt\": \"1/s\", \"flux_dxdt\": \"1/s\"},\n"
      << "  \"selected_step\": " << trajectory_index_ << ",\n"
      << "  \"network_backend\": ";
  write_json_string(out, Network::backend_name());
  out << ",\n";
#ifdef IMNET_USE_NUPPN
  out << "  \"network_files\": {\"nuppn_run_dir\": ";
  write_json_string(out, species_file_);
  out << "},\n";
#else
  out
      << "  \"network_files\": {";
  const std::pair<const char *, const std::string *> files[] = {
      {"species", &species_file_},     {"reaclib", &reaclib_file_},
      {"partition", &partition_file_}, {"mass", &mass_file_},
      {"weak", &weak_file_},
  };
  for (size_t i = 0; i < std::size(files); ++i) {
    if (i) {
      out << ",";
    }
    out << "\n    \"" << files[i].first << "\": ";
    write_json_string(out, *files[i].second);
  }
  out << "\n  },\n";
#endif
  out << "  \"view\": {\"flux_metric\": ";
  write_json_string(out, view_settings_.flux_arrow_metric == 2
                             ? "rate"
                             : view_settings_.flux_arrow_metric == 1 ? "dXdt"
                                                                     : "dYdt");
  out << ", \"flux_threshold\": " << view_settings_.flux_threshold
      << ", \"max_flux_arrows\": " << view_settings_.max_flux_arrows
      << "},\n  \"species\": [";
  for (size_t i = 0; i < species_data_.size(); ++i) {
    const auto &species = species_data_[i];
    out << (i ? ",\n" : "\n") << "    {\"index\": " << i << ", \"name\": ";
    write_json_string(out, species.name);
    out << ", \"Z\": " << species.Z << ", \"A\": " << species.A
        << ", \"N\": " << species.N
        << ", \"fixed\": "
        << (i < integration_settings_.fixed.size() &&
                    integration_settings_.fixed[i]
                ? "true"
                : "false")
        << ", \"active\": "
        << (i < integration_settings_.active.size() &&
                    integration_settings_.active[i]
                ? "true"
                : "false")
        << ", \"selected\": "
        << (i < integration_settings_.selected.size() &&
                    integration_settings_.selected[i]
                ? "true"
                : "false")
        << "}";
  }
  out << "\n  ],\n  \"trajectory\": [";
  for (size_t i = 0; i < trajectory_times_.size(); ++i) {
    out << (i ? ",\n" : "\n") << "    {\"step\": " << i
        << ", \"time\": " << trajectory_times_[i]
        << ", \"rho\": " << trajectory_rhos_[i]
        << ", \"temp\": " << trajectory_temps_[i] << "}";
  }
  out << (trajectory_times_.empty() ? "" : "\n  ") << "],\n  \"steps\": [";

  try {
    for (size_t i = 0; i < steps.size(); ++i) {
      const auto &step = steps[i];
      if (step.xnuc.size() != species_data_.size()) {
        throw std::runtime_error("saved abundance vector size mismatch");
      }
      const auto fluxes = network_->get_reaction_fluxes(
          step.rho, step.temp, step.xnuc, 0.0, 0,
          std::numeric_limits<size_t>::max());
      out << (i ? ",\n" : "\n") << "    {\"step\": " << step.index
          << ", \"time\": " << step.time << ", \"rho\": " << step.rho
          << ", \"temp\": " << step.temp << ", \"dt\": " << step.dt
          << ", \"dedt\": " << step.dedt << ", \"substeps\": "
          << step.substeps
          << ", \"success\": " << (step.success ? "true" : "false")
          << ", \"status\": ";
      write_json_string(out, step.status);
      out << ", \"error\": ";
      write_json_string(out, step.error);
      out << ", \"abundances\": [";
      for (size_t j = 0; j < step.xnuc.size(); ++j) {
        out << (j ? ", " : "") << step.xnuc[j];
      }
      out << "], \"fluxes\": [";
      for (size_t j = 0; j < fluxes.size(); ++j) {
        const auto &flux = fluxes[j];
        out << (j ? ", " : "") << "{\"source\": " << flux.source_index
            << ", \"target\": " << flux.target_index
            << ", \"strength_rate\": " << flux.strength_rate
            << ", \"strength_dydt\": " << flux.strength_dydt
            << ", \"strength_dxdt\": " << flux.strength_dxdt
            << ", \"weak\": " << (flux.weak ? "true" : "false")
            << ", \"equation\": ";
        write_json_string(out, flux.equation);
        out << "}";
      }
      out << "]}";
    }
  } catch (const std::exception &e) {
    std::cerr << "Error: Cannot save state: " << e.what() << std::endl;
    return false;
  }
  out << "\n  ]\n}\n";
  out.close();
  if (!out) {
    std::cerr << "Error: Failed while writing state file: " << filename
              << std::endl;
    return false;
  }
  std::cout << "Saved state to " << filename << std::endl;
  return true;
}

bool AppState::run_trajectory() {
  std::string error;
  if (run_trajectory(error)) {
    return true;
  }
  std::cerr << "Error running trajectory: " << error << std::endl;
  return false;
}

bool AppState::run_trajectory(std::string &error) {
  if (!network_) {
    error = "Network not initialized";
    return false;
  }
  if (!validate_trajectory_arrays(trajectory_times_, trajectory_rhos_,
                                  trajectory_temps_, error)) {
    return false;
  }

  if (!validate_composition(error)) {
    return false;
  }

  std::vector<double> xnuc = integration_settings_.xnuc;
  trajectory_dedt_.clear();
  trajectory_cache_.clear();
  invalidate_reaction_flux_state();
  trajectory_dedt_.reserve(trajectory_times_.size());
  trajectory_cache_.reserve(trajectory_times_.size());

  auto make_step = [&](size_t i, double dedt, int substeps, bool success,
                       const std::string &status,
                       const std::string &step_error) {
    TrajectoryStepCache step;
    step.index = i;
    step.time = trajectory_times_[i];
    step.rho = trajectory_rhos_[i];
    step.temp = trajectory_temps_[i];
    step.dt = (i + 1 < trajectory_times_.size())
                  ? std::max(trajectory_times_[i + 1] - trajectory_times_[i],
                             0.0)
                  : 0.0;
    step.dedt = dedt;
    step.substeps = substeps;
    step.success = success;
    step.status = status;
    step.error = step_error;
    step.xnuc = xnuc;
    trajectory_cache_.push_back(step);
    trajectory_dedt_.push_back(step.dedt);
  };

  make_step(0, 0.0, 0, true, "initial", "");

  for (size_t i = 0; i + 1 < trajectory_times_.size(); ++i) {
    const double dt = trajectory_times_[i + 1] - trajectory_times_[i];
    if (!std::isfinite(dt) || dt < 0.0) {
      error = "Invalid timestep at trajectory interval " + std::to_string(i);
      make_step(i + 1, 0.0, 0, false, "failed", error);
      trajectory_index_ = i + 1;
      set_trajectory_step(trajectory_index_);
      return false;
    }

    double dedt = 0.0;
    try {
      std::vector<double> next_xnuc = xnuc;
      if (dt > 0.0) {
        dedt = network_->integrate_interval(
            trajectory_rhos_[i], trajectory_temps_[i], trajectory_rhos_[i + 1],
            trajectory_temps_[i + 1], next_xnuc, dt);
      }
      validate_step_composition(next_xnuc, species_names_);
      xnuc = std::move(next_xnuc);
      make_step(i + 1, dedt, network_->last_substeps(), true, "ok", "");
    } catch (const std::exception &e) {
      error = "Step " + std::to_string(i + 1) + " failed: " + e.what();
      make_step(i + 1, dedt, network_->last_substeps(), false, "failed",
                error);
      trajectory_index_ = i + 1;
      set_trajectory_step(trajectory_index_);
      return false;
    }
  }

  trajectory_index_ = std::min(trajectory_index_, trajectory_cache_.size() - 1);
  set_trajectory_step(trajectory_index_);

  std::cout << "Trajectory integration complete: " << trajectory_cache_.size()
            << " steps cached" << std::endl;
  error.clear();
  return true;
}

bool AppState::set_trajectory_step(size_t step) {
  if (step >= trajectory_times_.size()) {
    return false;
  }
  trajectory_index_ = step;

  if (const auto *cached = get_trajectory_cache_step(step)) {
    integration_settings_.rho = cached->rho;
    integration_settings_.temp = cached->temp;
    integration_settings_.dt = cached->dt;
    integration_settings_.xnuc = cached->xnuc;
    if (cached->success) {
      record_reaction_flux_state(cached->rho, cached->temp, cached->xnuc);
    } else {
      invalidate_reaction_flux_state();
    }
  } else {
    integration_settings_.rho = trajectory_rhos_[step];
    integration_settings_.temp = trajectory_temps_[step];
    integration_settings_.dt =
        (step + 1 < trajectory_times_.size())
            ? std::max(trajectory_times_[step + 1] - trajectory_times_[step],
                       0.0)
            : 0.0;
    invalidate_reaction_flux_state();
  }

  return true;
}

bool AppState::set_trajectory_row(size_t step, double time, double rho,
                                  double temp, std::string &error) {
  if (step >= trajectory_times_.size()) {
    error = "Trajectory row index out of range";
    return false;
  }

  std::vector<double> candidate_times = trajectory_times_;
  std::vector<double> candidate_rhos = trajectory_rhos_;
  std::vector<double> candidate_temps = trajectory_temps_;
  candidate_times[step] = time;
  candidate_rhos[step] = rho;
  candidate_temps[step] = temp;

  if (!validate_trajectory_arrays(candidate_times, candidate_rhos,
                                  candidate_temps, error)) {
    return false;
  }

  trajectory_times_ = std::move(candidate_times);
  trajectory_rhos_ = std::move(candidate_rhos);
  trajectory_temps_ = std::move(candidate_temps);
  trajectory_index_ = std::min(trajectory_index_, trajectory_times_.size() - 1);
  invalidate_trajectory_cache();
  set_trajectory_step(trajectory_index_);
  error.clear();
  return true;
}

bool AppState::append_trajectory_row(double time, double rho, double temp,
                                     std::string &error) {
  std::vector<double> candidate_times = trajectory_times_;
  std::vector<double> candidate_rhos = trajectory_rhos_;
  std::vector<double> candidate_temps = trajectory_temps_;
  candidate_times.push_back(time);
  candidate_rhos.push_back(rho);
  candidate_temps.push_back(temp);

  if (!validate_trajectory_arrays(candidate_times, candidate_rhos,
                                  candidate_temps, error)) {
    return false;
  }

  trajectory_times_ = std::move(candidate_times);
  trajectory_rhos_ = std::move(candidate_rhos);
  trajectory_temps_ = std::move(candidate_temps);
  trajectory_index_ = std::min(trajectory_index_, trajectory_times_.size() - 1);
  invalidate_trajectory_cache();
  set_trajectory_step(trajectory_index_);
  error.clear();
  return true;
}

bool AppState::remove_trajectory_row(size_t step, std::string &error) {
  if (trajectory_times_.size() <= 1) {
    error = "Cannot remove the last trajectory row";
    return false;
  }
  if (step >= trajectory_times_.size()) {
    error = "Trajectory row index out of range";
    return false;
  }

  trajectory_times_.erase(trajectory_times_.begin() + step);
  trajectory_rhos_.erase(trajectory_rhos_.begin() + step);
  trajectory_temps_.erase(trajectory_temps_.begin() + step);
  trajectory_index_ = std::min(trajectory_index_, trajectory_times_.size() - 1);
  invalidate_trajectory_cache();
  set_trajectory_step(trajectory_index_);
  error.clear();
  return true;
}

const TrajectoryStepCache *
AppState::get_trajectory_cache_step(size_t step) const {
  const auto it =
      std::find_if(trajectory_cache_.begin(), trajectory_cache_.end(),
                   [step](const TrajectoryStepCache &cached) {
                     return cached.index == step;
                   });
  if (it == trajectory_cache_.end()) {
    return nullptr;
  }
  return &(*it);
}

void AppState::set_trajectory_cache(std::vector<TrajectoryStepCache> cache,
                                    size_t current_step) {
  trajectory_cache_ = std::move(cache);
  std::sort(trajectory_cache_.begin(), trajectory_cache_.end(),
            [](const TrajectoryStepCache &a, const TrajectoryStepCache &b) {
              return a.index < b.index;
            });

  trajectory_dedt_.clear();
  trajectory_dedt_.reserve(trajectory_cache_.size());
  for (const auto &step : trajectory_cache_) {
    trajectory_dedt_.push_back(step.dedt);
  }

  trajectory_index_ =
      trajectory_times_.empty()
          ? 0
          : std::min(current_step, trajectory_times_.size() - 1);
  set_trajectory_step(trajectory_index_);
}

void AppState::clear_trajectory_cache() { invalidate_trajectory_cache(); }

void AppState::invalidate_trajectory_cache() {
  trajectory_dedt_.clear();
  trajectory_cache_.clear();
  invalidate_reaction_flux_state();
}

bool AppState::validate_composition(std::string &error) const {
  if (integration_settings_.xnuc.size() != species_names_.size()) {
    error = "Abundance vector size does not match species count";
    return false;
  }
  for (size_t i = 0; i < integration_settings_.xnuc.size(); ++i) {
    const double x = integration_settings_.xnuc[i];
    if (!std::isfinite(x) || x < 0.0) {
      error = "Invalid abundance for " + species_names_[i];
      return false;
    }
  }
  error.clear();
  return true;
}

bool AppState::validate_trajectory_arrays(const std::vector<double> &times,
                                          const std::vector<double> &rhos,
                                          const std::vector<double> &temps,
                                          std::string &error) const {
  if (times.empty()) {
    error = "Trajectory must contain at least one row";
    return false;
  }
  if (times.size() != rhos.size() || times.size() != temps.size()) {
    error = "Trajectory column sizes do not match";
    return false;
  }

  for (size_t i = 0; i < times.size(); ++i) {
    if (!std::isfinite(times[i]) || !std::isfinite(rhos[i]) ||
        !std::isfinite(temps[i])) {
      error = "Trajectory row " + std::to_string(i) +
              " contains a non-finite value";
      return false;
    }
    if (rhos[i] <= 0.0 || temps[i] <= 0.0) {
      error = "Trajectory row " + std::to_string(i) +
              " must have positive rho and T";
      return false;
    }
    if (i > 0 && times[i] < times[i - 1]) {
      error = "Trajectory time decreases at row " + std::to_string(i);
      return false;
    }
  }

  error.clear();
  return true;
}

double AppState::calculate_ye() const {
  return imyann::calculate_ye(integration_settings_.xnuc, species_data_);
}

double AppState::calculate_abar() const {
  return imyann::calculate_abar(integration_settings_.xnuc, species_data_);
}

double AppState::calculate_zbar() const {
  return imyann::calculate_zbar(integration_settings_.xnuc, species_data_);
}

double AppState::get_species_rate(int species_index) const {
  if (!network_) {
    return 0.0;
  }
  if (species_index < 0 ||
      species_index >= static_cast<int>(integration_settings_.xnuc.size())) {
    return 0.0;
  }

  try {
    return network_->get_reaction_rate(species_index, integration_settings_.rho,
                                       integration_settings_.temp,
                                       integration_settings_.xnuc);
  } catch (...) {
    return 0.0;
  }
}

std::vector<ReactionDiagnostic>
AppState::get_species_reaction_diagnostics(int species_index,
                                           size_t max_count) const {
  if (!network_) {
    return {};
  }
  if (species_index < 0 ||
      species_index >= static_cast<int>(integration_settings_.xnuc.size())) {
    return {};
  }

  try {
    return network_->get_reaction_diagnostics(
        species_index, integration_settings_.rho, integration_settings_.temp,
        integration_settings_.xnuc, max_count);
  } catch (...) {
    return {};
  }
}

std::vector<ReactionFlux>
AppState::get_reaction_fluxes(double min_strength, int metric,
                              size_t max_count, bool include_regular,
                              bool include_weak) const {
  if (!network_ || !flux_state_valid_) {
    return {};
  }
  if (flux_state_xnuc_.size() != species_names_.size()) {
    return {};
  }

  try {
    return network_->get_reaction_fluxes(
        flux_state_rho_, flux_state_temp_, flux_state_xnuc_, min_strength,
        metric, max_count,
        include_regular, include_weak);
  } catch (...) {
    return {};
  }
}

void AppState::record_reaction_flux_state(double rho, double temp,
                                          const std::vector<double> &xnuc) {
  flux_state_rho_ = rho;
  flux_state_temp_ = temp;
  flux_state_xnuc_ = xnuc;
  flux_state_valid_ = true;
}

void AppState::invalidate_reaction_flux_state() {
  flux_state_valid_ = false;
  flux_state_xnuc_.clear();
}

bool AppState::is_species_connected(int species_index) const {
  if (!network_) {
    return false;
  }
  if (species_index < 0 ||
      species_index >= static_cast<int>(species_names_.size())) {
    return false;
  }

  try {
    return network_->is_connected(species_names_[species_index],
                                  integration_settings_.active);
  } catch (...) {
    return false;
  }
}

} // namespace imyann
