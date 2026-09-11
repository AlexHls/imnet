#include "network_wrapper.h"

#ifdef IMNET_USE_NUPPN

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <queue>
#include <stdexcept>
#include <unordered_map>

namespace imyann {

namespace {

class ScopedCurrentPath {
public:
  explicit ScopedCurrentPath(const std::string &path)
      : previous_(std::filesystem::current_path()) {
    std::filesystem::current_path(path);
  }

  ~ScopedCurrentPath() {
    std::error_code error;
    std::filesystem::current_path(previous_, error);
  }

private:
  std::filesystem::path previous_;
};

std::string normalize_nuppn_species_name(const char *name) {
  std::string value(name ? name : "");
  value.erase(std::remove_if(value.begin(), value.end(),
                             [](unsigned char c) { return std::isspace(c); }),
              value.end());
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (value == "prot") {
    return "p";
  }
  if (value == "neut") {
    return "n";
  }
  return value;
}

void validate_conditions(double rho, double temp) {
  if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(temp) ||
      temp <= 0.0) {
    throw std::runtime_error("Invalid density or temperature");
  }
}

void validate_composition(const std::vector<double> &xnuc,
                          size_t expected_size) {
  if (xnuc.size() != expected_size) {
    throw std::runtime_error("xnuc size mismatch with network");
  }
  if (std::any_of(xnuc.begin(), xnuc.end(),
                  [](double x) { return !std::isfinite(x) || x < 0.0; })) {
    throw std::runtime_error("Invalid abundance");
  }
}

void check_nuppn_status(int status, const char *operation) {
  if (status != 0) {
    throw std::runtime_error(std::string(operation) + " failed with status " +
                             std::to_string(status));
  }
}

std::string &initialized_nuppn_run_dir() {
  static std::string run_dir;
  return run_dir;
}

struct NuppnParticipant {
  int index = -1;
  int count = 0;
};

struct NuppnReaction {
  NuppnParticipant inputs[2];
  NuppnParticipant outputs[2];
  double rate = 0.0;
  double flow = 0.0;
  double q_value = 0.0;
  bool weak = false;
};

struct FluxAccumulator {
  double strength_rate = 0.0;
  double strength_dydt = 0.0;
  double strength_dxdt = 0.0;
  double strongest_reaction = 0.0;
  bool weak = false;
  std::string equation;
};

bool valid_participant(const NuppnParticipant &participant,
                       size_t species_count) {
  return participant.count > 0 && participant.index >= 0 &&
         participant.index < static_cast<int>(species_count);
}

std::string join_nuppn_participants(const NuppnParticipant participants[2],
                                    const std::vector<Species> &species) {
  std::string result;
  for (int i = 0; i < 2; ++i) {
    if (!valid_participant(participants[i], species.size())) {
      continue;
    }
    for (int n = 0; n < participants[i].count; ++n) {
      if (!result.empty()) {
        result += " + ";
      }
      result += species[static_cast<size_t>(participants[i].index)].name;
    }
  }
  return result.empty() ? "?" : result;
}

std::string reaction_equation(const NuppnReaction &reaction,
                              const std::vector<Species> &species) {
  return join_nuppn_participants(reaction.inputs, species) + " -> " +
         join_nuppn_participants(reaction.outputs, species);
}

NuppnReaction load_nuppn_reaction(int index) {
  NuppnReaction reaction;
  int weak = 0;
  check_nuppn_status(
      nuppn_get_reaction(index, &reaction.inputs[0].index,
                         &reaction.inputs[0].count,
                         &reaction.inputs[1].index,
                         &reaction.inputs[1].count,
                         &reaction.outputs[0].index,
                         &reaction.outputs[0].count,
                         &reaction.outputs[1].index,
                         &reaction.outputs[1].count, &reaction.rate,
                         &reaction.flow, &reaction.q_value, &weak),
      "nuppn_get_reaction");
  reaction.weak = weak != 0;
  return reaction;
}

std::vector<IntegrationStepSnapshot>
load_nuppn_step_history(const std::vector<Species> &species) {
  const int count = nuppn_history_size();
  if (count <= 0) {
    return {};
  }

  std::vector<IntegrationStepSnapshot> history;
  history.reserve(static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) {
    IntegrationStepSnapshot step;
    step.index = static_cast<size_t>(i);
    step.xnuc.assign(species.size(), 0.0);
    check_nuppn_status(
        nuppn_get_history_step(i, &step.time, &step.rho, &step.temp, &step.dt,
                               &step.dedt, &step.substeps, step.xnuc.data(),
                               static_cast<int>(step.xnuc.size())),
        "nuppn_get_history_step");
    history.push_back(std::move(step));
  }
  return history;
}

int participant_stoich(const NuppnReaction &reaction, int species_index) {
  int coefficient = 0;
  for (const auto &input : reaction.inputs) {
    if (input.index == species_index) {
      coefficient -= input.count;
    }
  }
  for (const auto &output : reaction.outputs) {
    if (output.index == species_index) {
      coefficient += output.count;
    }
  }
  return coefficient;
}

} // namespace

const char *Network::backend_name() { return "NuPPN"; }

Network::Network(const std::string &species_file,
                 const std::string &reaclib_file,
                 const std::string &partition_file,
                 const std::string &mass_file, const std::string &weak_file)
    : run_dir_(species_file), initialized_(false) {
  (void)reaclib_file;
  (void)partition_file;
  (void)mass_file;
  (void)weak_file;
  if (run_dir_.empty()) {
    run_dir_ = IMNET_NUPPN_RUN_DIR;
  }
  if (!std::filesystem::is_directory(run_dir_)) {
    throw std::runtime_error("NuPPN run directory not found: " + run_dir_);
  }
  run_dir_ = std::filesystem::weakly_canonical(run_dir_).string();
  auto &first_run_dir = initialized_nuppn_run_dir();
  if (!first_run_dir.empty() && first_run_dir != run_dir_) {
    throw std::runtime_error(
        "NuPPN backend is already initialized; restart to change run directory");
  }

  ScopedCurrentPath cwd(run_dir_);
  check_nuppn_status(nuppn_init(), "nuppn_init");
  if (first_run_dir.empty()) {
    first_run_dir = run_dir_;
  }

  const int count = nuppn_num_species();
  if (count <= 0) {
    throw std::runtime_error("NuPPN loaded no species");
  }

  species_.reserve(static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) {
    char name[16] = {};
    int z = 0;
    int a = 0;
    check_nuppn_status(nuppn_get_species(i, name, sizeof(name), &z, &a),
                       "nuppn_get_species");
    species_.push_back({normalize_nuppn_species_name(name), z, a, a - z});
  }
  initialized_ = true;
}

Network::~Network() = default;

double Network::integrate(double rho, double temp, std::vector<double> &xnuc,
                          double dt) {
  last_step_history_.clear();
  return integrate_interval(rho, temp, rho, temp, xnuc, dt);
}

double Network::integrate_interval(double rho0, double temp0, double rho1,
                                   double temp1, std::vector<double> &xnuc,
                                   double dt) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  validate_conditions(rho0, temp0);
  validate_conditions(rho1, temp1);
  validate_composition(xnuc, species_.size());
  if (!std::isfinite(dt) || dt < 0.0) {
    throw std::runtime_error("Invalid timestep");
  }

  last_step_history_.clear();
  double dedt = 0.0;
  ScopedCurrentPath cwd(run_dir_);
  check_nuppn_status(
      nuppn_integrate(rho0, temp0, rho1, temp1, xnuc.data(),
                      static_cast<int>(xnuc.size()), dt, &dedt),
      "nuppn_integrate");
  last_substeps_ = nuppn_last_substeps();
  validate_composition(xnuc, species_.size());
  return dedt;
}

double Network::integrate_to_time(double rho, double temp,
                                  std::vector<double> &xnuc,
                                  double final_time, double initial_dt,
                                  double max_dt, double dt_factor,
                                  int max_steps) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  last_step_history_.clear();
  validate_conditions(rho, temp);
  validate_composition(xnuc, species_.size());
  if (!std::isfinite(final_time) || final_time <= 0.0 ||
      !std::isfinite(initial_dt) || initial_dt <= 0.0 ||
      !std::isfinite(max_dt) || max_dt <= 0.0 ||
      !std::isfinite(dt_factor) || dt_factor <= 0.0 || max_steps <= 0) {
    throw std::runtime_error("Invalid final-time integration settings");
  }

  double dedt = 0.0;
  ScopedCurrentPath cwd(run_dir_);
  const int status =
      nuppn_integrate_to_time(rho, temp, xnuc.data(),
                              static_cast<int>(xnuc.size()), final_time,
                              initial_dt, max_dt, dt_factor, max_steps, &dedt);
  last_substeps_ = nuppn_last_substeps();
  last_step_history_ = load_nuppn_step_history(species_);
  check_nuppn_status(status, "nuppn_integrate_to_time");
  validate_composition(xnuc, species_.size());
  return dedt;
}

void Network::compute_nse(double rho, double temp, double ye,
                          std::vector<double> &xnuc_out) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc_out, species_.size());
  if (!std::isfinite(ye) || ye < 0.0 || ye > 1.0) {
    throw std::runtime_error("Invalid NSE parameters");
  }

  ScopedCurrentPath cwd(run_dir_);
  check_nuppn_status(
      nuppn_compute_nse(rho, temp, ye, xnuc_out.data(),
                        static_cast<int>(xnuc_out.size())),
      "nuppn_compute_nse");
  validate_composition(xnuc_out, species_.size());
}

double Network::get_reaction_rate(int species_index, double rho, double temp,
                                  const std::vector<double> &xnuc) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }
  if (species_index < 0 || species_index >= num_species()) {
    throw std::runtime_error("Species index out of range");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc, species_.size());

  std::vector<double> dxdt(xnuc.size(), 0.0);
  ScopedCurrentPath cwd(run_dir_);
  check_nuppn_status(
      nuppn_get_dxdt(rho, temp, xnuc.data(), static_cast<int>(xnuc.size()),
                     dxdt.data()),
      "nuppn_get_dxdt");
  return dxdt[static_cast<size_t>(species_index)];
}

std::vector<ReactionDiagnostic> Network::get_reaction_diagnostics(
    int species_index, double rho, double temp, const std::vector<double> &xnuc,
    size_t max_count) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }
  if (species_index < 0 || species_index >= num_species()) {
    throw std::runtime_error("Species index out of range");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc, species_.size());

  ScopedCurrentPath cwd(run_dir_);
  check_nuppn_status(
      nuppn_evaluate_reactions(rho, temp, xnuc.data(),
                               static_cast<int>(xnuc.size())),
      "nuppn_evaluate_reactions");

  const int reaction_count = nuppn_num_reactions();
  std::vector<ReactionDiagnostic> diagnostics;
  diagnostics.reserve(static_cast<size_t>(reaction_count));
  for (int i = 0; i < reaction_count; ++i) {
    const NuppnReaction reaction = load_nuppn_reaction(i);
    const int coefficient = participant_stoich(reaction, species_index);
    if (coefficient == 0) {
      continue;
    }

    ReactionDiagnostic item;
    item.equation = reaction_equation(reaction, species_);
    item.q_value = reaction.q_value;
    item.rate = reaction.rate;
    item.abundance_weighted_rate = reaction.flow;
    item.contribution_dYdt = static_cast<double>(coefficient) * reaction.flow;
    item.contribution_dXdt =
        item.contribution_dYdt *
        static_cast<double>(species_[static_cast<size_t>(species_index)].A);
    item.weak = reaction.weak;
    diagnostics.push_back(std::move(item));
  }

  std::sort(diagnostics.begin(), diagnostics.end(),
            [](const ReactionDiagnostic &a, const ReactionDiagnostic &b) {
              return std::abs(a.contribution_dXdt) >
                     std::abs(b.contribution_dXdt);
            });
  if (diagnostics.size() > max_count) {
    diagnostics.resize(max_count);
  }
  return diagnostics;
}

std::vector<ReactionFlux> Network::get_reaction_fluxes(
    double rho, double temp, const std::vector<double> &xnuc,
    double min_strength, int metric, size_t max_count, bool include_regular,
    bool include_weak) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc, species_.size());
  if (!std::isfinite(min_strength) || min_strength < 0.0) {
    throw std::runtime_error("Invalid minimum flux strength");
  }
  if (metric < 0 || metric > 2) {
    throw std::runtime_error("Invalid flux metric");
  }

  ScopedCurrentPath cwd(run_dir_);
  check_nuppn_status(
      nuppn_evaluate_reactions(rho, temp, xnuc.data(),
                               static_cast<int>(xnuc.size())),
      "nuppn_evaluate_reactions");

  std::unordered_map<long long, FluxAccumulator> flux_by_pair;
  const auto add_flux = [&](int source, int target, double strength_rate,
                            double strength_dydt, double strength_dxdt,
                            bool weak, const std::string &equation) {
    const double selected_strength = metric == 2   ? strength_rate
                                     : metric == 1 ? strength_dxdt
                                                   : strength_dydt;
    if (!std::isfinite(strength_rate) || !std::isfinite(strength_dydt) ||
        !std::isfinite(strength_dxdt) || selected_strength <= 0.0 ||
        source < 0 || target < 0 || source == target ||
        source >= num_species() || target >= num_species()) {
      return;
    }
    const long long pair =
        static_cast<long long>(source) * static_cast<long long>(num_species()) +
        static_cast<long long>(target);
    const long long key = pair * 2 + (weak ? 1 : 0);
    auto &acc = flux_by_pair[key];
    acc.strength_rate += strength_rate;
    acc.strength_dydt += strength_dydt;
    acc.strength_dxdt += strength_dxdt;
    acc.weak = weak;
    if (selected_strength >= acc.strongest_reaction || acc.equation.empty()) {
      acc.strongest_reaction = selected_strength;
      acc.equation = equation;
    }
  };

  const int reaction_count = nuppn_num_reactions();
  for (int r = 0; r < reaction_count; ++r) {
    const NuppnReaction reaction = load_nuppn_reaction(r);
    const double reaction_rate = std::abs(reaction.rate);
    const double event_flow = std::abs(reaction.flow);
    const std::string equation = reaction_equation(reaction, species_);
    for (const auto &input : reaction.inputs) {
      if (!valid_participant(input, species_.size())) {
        continue;
      }
      for (const auto &output : reaction.outputs) {
        if (!valid_participant(output, species_.size())) {
          continue;
        }
        const double dydt_flow =
            event_flow * static_cast<double>(std::max(input.count,
                                                      output.count));
        const double dxdt_flow =
            event_flow *
            std::max(static_cast<double>(
                         input.count * species_[input.index].A),
                     static_cast<double>(
                         output.count * species_[output.index].A));
        add_flux(input.index, output.index, reaction_rate, dydt_flow,
                 dxdt_flow, reaction.weak, equation);
      }
    }
  }

  std::vector<ReactionFlux> fluxes;
  fluxes.reserve(flux_by_pair.size());
  for (const auto &[key, acc] : flux_by_pair) {
    if ((acc.weak && !include_weak) || (!acc.weak && !include_regular)) {
      continue;
    }
    const long long pair = key / 2;
    const double strength = metric == 2   ? acc.strength_rate
                            : metric == 1 ? acc.strength_dxdt
                                          : acc.strength_dydt;
    if (strength < min_strength) {
      continue;
    }
    ReactionFlux flux;
    flux.source_index =
        static_cast<int>(pair / static_cast<long long>(num_species()));
    flux.target_index =
        static_cast<int>(pair % static_cast<long long>(num_species()));
    flux.strength = strength;
    flux.strength_rate = acc.strength_rate;
    flux.strength_dydt = acc.strength_dydt;
    flux.strength_dxdt = acc.strength_dxdt;
    flux.weak = acc.weak;
    flux.equation = acc.equation;
    fluxes.push_back(std::move(flux));
  }

  std::sort(fluxes.begin(), fluxes.end(),
            [](const ReactionFlux &a, const ReactionFlux &b) {
              return a.strength > b.strength;
            });
  size_t regular_count = 0;
  size_t weak_count = 0;
  fluxes.erase(
      std::remove_if(fluxes.begin(), fluxes.end(),
                     [&](const ReactionFlux &flux) {
                       size_t &count = flux.weak ? weak_count : regular_count;
                       return count++ >= max_count;
                     }),
      fluxes.end());
  return fluxes;
}

bool Network::is_connected(const std::string &species_name,
                           const std::vector<bool> &active_mask) {
  if (active_mask.size() != species_.size()) {
    throw std::runtime_error("active_mask size mismatch");
  }
  const int target = find_species_index(species_name);
  if (target < 0) {
    return false;
  }
  if (active_mask[static_cast<size_t>(target)]) {
    return true;
  }

  const auto graph = build_connectivity_graph();
  std::vector<bool> visited(species_.size(), false);
  std::queue<int> q;
  for (size_t i = 0; i < active_mask.size(); ++i) {
    if (active_mask[i]) {
      visited[i] = true;
      q.push(static_cast<int>(i));
    }
  }

  while (!q.empty()) {
    const int u = q.front();
    q.pop();
    if (u == target) {
      return true;
    }
    for (int v : graph[static_cast<size_t>(u)]) {
      if (!visited[static_cast<size_t>(v)]) {
        visited[static_cast<size_t>(v)] = true;
        q.push(v);
      }
    }
  }
  return false;
}

int Network::find_species_index(const std::string &species_name) const {
  for (size_t i = 0; i < species_.size(); ++i) {
    if (species_name == species_[i].name) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

std::vector<std::vector<int>> Network::build_connectivity_graph() const {
  std::vector<std::vector<int>> graph(species_.size());
  const int reaction_count = nuppn_num_reactions();
  for (int r = 0; r < reaction_count; ++r) {
    const NuppnReaction reaction = load_nuppn_reaction(r);
    std::vector<int> nodes;
    for (const auto &input : reaction.inputs) {
      if (valid_participant(input, species_.size())) {
        nodes.push_back(input.index);
      }
    }
    for (const auto &output : reaction.outputs) {
      if (valid_participant(output, species_.size())) {
        nodes.push_back(output.index);
      }
    }
    for (size_t i = 0; i < nodes.size(); ++i) {
      for (size_t j = i + 1; j < nodes.size(); ++j) {
        const int a = nodes[i];
        const int b = nodes[j];
        if (a != b) {
          graph[static_cast<size_t>(a)].push_back(b);
          graph[static_cast<size_t>(b)].push_back(a);
        }
      }
    }
  }
  return graph;
}

} // namespace imyann

#else

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <queue>
#include <stdexcept>
#include <unordered_map>

namespace imyann {

namespace {

std::string trim_species_name(const char *name) {
  std::string value(name ? name : "");
  const auto first =
      std::find_if_not(value.begin(), value.end(), [](unsigned char c) {
        return std::isspace(c);
      });
  const auto last =
      std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) {
        return std::isspace(c);
      }).base();
  if (first >= last) {
    return "";
  }
  return std::string(first, last);
}

std::string join_species(const int indices[], int count,
                         const struct network_data &nd) {
  std::string result;
  for (int i = 0; i < count; ++i) {
    const int index = indices[i];
    if (i > 0) {
      result += " + ";
    }
    if (index >= 0 && index < static_cast<int>(nd.nuc_count)) {
      result += trim_species_name(nd.nucdata[index].name);
    } else {
      result += "?";
    }
  }
  return result;
}

struct FluxAccumulator {
  double strength_rate = 0.0;
  double strength_dydt = 0.0;
  double strength_dxdt = 0.0;
  double strongest_reaction = 0.0;
  bool weak = false;
  std::string equation;
};

void validate_conditions(double rho, double temp) {
  if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(temp) ||
      temp <= 0.0) {
    throw std::runtime_error("Invalid density or temperature");
  }
}

void validate_composition(const std::vector<double> &xnuc,
                          size_t expected_size) {
  if (xnuc.size() != expected_size) {
    throw std::runtime_error("xnuc size mismatch with network");
  }
  if (std::any_of(xnuc.begin(), xnuc.end(),
                  [](double x) { return !std::isfinite(x) || x < 0.0; })) {
    throw std::runtime_error("Invalid abundance");
  }
}

} // namespace

// ============================================================================
// Network Implementation
// ============================================================================

const char *Network::backend_name() { return "YANN"; }

Network::Network(const std::string &species_file,
                 const std::string &reaclib_file,
                 const std::string &partition_file,
                 const std::string &mass_file, const std::string &weak_file)
    : initialized_(false) {
  std::memset(&nd_, 0, sizeof(nd_));
  std::memset(&nw_, 0, sizeof(nw_));

  const int rc =
      network_init(const_cast<char *>(species_file.c_str()),
                   const_cast<char *>(reaclib_file.c_str()),
                   const_cast<char *>(partition_file.c_str()),
                   const_cast<char *>(mass_file.c_str()),
                   const_cast<char *>(weak_file.c_str()), nullptr,
                   static_cast<eos_mode>(EOS_IONIZED | EOS_RADIATION |
                                         EOS_DEGENERATE | EOS_COULOMB),
                   &nd_);

  if (rc != 0 || !nd_.initialized) {
    throw std::runtime_error("network_init failed");
  }

  if (network_workspace_init(&nd_, &nw_) != 0 || !nw_.initialized) {
    network_deinit(&nd_);
    throw std::runtime_error("network_workspace_init failed");
  }

  try {
    species_.reserve(nd_.nuc_count);
    for (size_t i = 0; i < nd_.nuc_count; ++i) {
      species_.push_back({trim_species_name(nd_.nucdata[i].name),
                          nd_.nucdata[i].nz, nd_.nucdata[i].na,
                          nd_.nucdata[i].nn});
    }
  } catch (...) {
    network_workspace_deinit(&nw_);
    network_deinit(&nd_);
    throw;
  }
  initialized_ = true;
}

Network::~Network() {
  if (initialized_) {
    if (nw_.initialized) {
      network_workspace_deinit(&nw_);
    }
    if (nd_.initialized) {
      network_deinit(&nd_);
    }
    initialized_ = false;
  }
}

double Network::integrate(double rho, double temp, std::vector<double> &xnuc,
                          double dt) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  last_step_history_.clear();
  validate_conditions(rho, temp);
  validate_composition(xnuc, nd_.nuc_count);
  if (!std::isfinite(dt) || dt < 0.0) {
    throw std::runtime_error("Invalid timestep");
  }

  last_step_history_.clear();
  std::vector<double> dx(nd_.nuc_count, 0.0);
  double dedt = 0.0;
  double drhodt = 0.0;

  const int rc = network_integrate(temp, rho, xnuc.data(), dx.data(), dt, &dedt,
                                   &drhodt, &nd_, &nw_);

  if (rc != 0) {
    throw std::runtime_error("network_integrate failed");
  }
  if (!std::isfinite(dedt)) {
    throw std::runtime_error("network_integrate returned non-finite energy");
  }

#ifdef NETWORK_ABSOLUTE_RESULT
  std::vector<double> result = std::move(dx);
#else
  std::vector<double> result = xnuc;
  for (size_t i = 0; i < xnuc.size(); ++i) {
    result[i] += dx[i] * dt;
    if (result[i] < 0.0) {
      result[i] = 0.0;
    }
  }
#endif
  validate_composition(result, nd_.nuc_count);
  xnuc = std::move(result);
  last_substeps_ = dt > 0.0 ? 1 : 0;
  return dedt;
}

double Network::integrate_interval(double rho0, double temp0, double rho1,
                                   double temp1, std::vector<double> &xnuc,
                                   double dt) {
  validate_conditions(rho1, temp1);
  return integrate(rho0, temp0, xnuc, dt);
}

double Network::integrate_to_time(double rho, double temp,
                                  std::vector<double> &xnuc,
                                  double final_time, double initial_dt,
                                  double max_dt, double dt_factor,
                                  int max_steps) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc, nd_.nuc_count);
  if (!std::isfinite(final_time) || final_time <= 0.0 ||
      !std::isfinite(initial_dt) || initial_dt <= 0.0 ||
      !std::isfinite(max_dt) || max_dt <= 0.0 ||
      !std::isfinite(dt_factor) || dt_factor <= 0.0 || max_steps <= 0) {
    throw std::runtime_error("Invalid final-time integration settings");
  }

  double time = 0.0;
  double next_dt = std::min(initial_dt, max_dt);
  double dedt = 0.0;
  int steps = 0;
  std::vector<IntegrationStepSnapshot> history;
  history.push_back({0, 0.0, rho, temp, next_dt, 0.0, 0, xnuc});
  for (; time < final_time && steps < max_steps; ++steps) {
    const double step_dt = std::min(next_dt, final_time - time);
    if (!std::isfinite(step_dt) || step_dt <= 0.0) {
      throw std::runtime_error("Invalid final-time timestep");
    }
    dedt = integrate(rho, temp, xnuc, step_dt);
    time += step_dt;
    history.push_back({static_cast<size_t>(steps + 1), time, rho, temp,
                       step_dt, dedt, last_substeps_, xnuc});
    const double grown = step_dt * dt_factor;
    next_dt =
        std::min(std::isfinite(grown) && grown > 0.0 ? grown : max_dt,
                 max_dt);
  }
  if (time < final_time) {
    last_substeps_ = steps;
    last_step_history_ = std::move(history);
    throw std::runtime_error("Final-time integration exceeded max steps");
  }

  last_substeps_ = steps;
  last_step_history_ = std::move(history);
  return dedt;
}

void Network::compute_nse(double rho, double temp, double ye,
                          std::vector<double> &xnuc_out) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  validate_conditions(rho, temp);
  if (xnuc_out.size() != nd_.nuc_count) {
    throw std::runtime_error("xnuc_out size mismatch with network");
  }
  if (!std::isfinite(ye) || ye < 0.0 || ye > 1.0) {
    throw std::runtime_error("Invalid NSE parameters");
  }

  std::vector<double> result(nd_.nuc_count);
  const int rc = network_nse(temp, rho, ye, result.data(), &nd_, &nw_);
  if (rc != 0) {
    throw std::runtime_error("network_nse failed");
  }
  validate_composition(result, nd_.nuc_count);
  xnuc_out = std::move(result);
}

double Network::get_reaction_rate(int species_index, double rho, double temp,
                                  const std::vector<double> &xnuc) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  if (species_index < 0 || species_index >= static_cast<int>(nd_.nuc_count)) {
    throw std::runtime_error("Species index out of range");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc, nd_.nuc_count);

  std::vector<double> y(nd_.nuc_count, 0.0);
  for (size_t i = 0; i < nd_.nuc_count; ++i) {
    y[i] = xnuc[i] / static_cast<double>(nd_.nucdata[i].na);
  }

  std::vector<double> rhs(nd_.n_matrix, 0.0);
  std::vector<network_var> deriv(nd_.n_matrix);
  const int rc = network_getrhs(rho, temp, y.data(), 0, &nd_, &nw_, rhs.data(),
                                deriv.data());

  if (rc != 0) {
    throw std::runtime_error("network_getrhs failed");
  }

  const double rate =
      rhs[species_index] * static_cast<double>(nd_.nucdata[species_index].na);
  if (!std::isfinite(rate)) {
    throw std::runtime_error("network_getrhs returned a non-finite rate");
  }
  return rate;
}

std::vector<ReactionDiagnostic> Network::get_reaction_diagnostics(
    int species_index, double rho, double temp, const std::vector<double> &xnuc,
    size_t max_count) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  if (species_index < 0 || species_index >= static_cast<int>(nd_.nuc_count)) {
    throw std::runtime_error("Species index out of range");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc, nd_.nuc_count);

  std::vector<double> y(nd_.nuc_count, 0.0);
  for (size_t i = 0; i < nd_.nuc_count; ++i) {
    y[i] = xnuc[i] / static_cast<double>(nd_.nucdata[i].na);
  }

  std::vector<double> rhs(nd_.n_matrix, 0.0);
  std::vector<network_var> deriv(nd_.n_matrix);
  const int rc = network_getrhs(rho, temp, y.data(), 0, &nd_, &nw_, rhs.data(),
                                deriv.data());
  if (rc != 0) {
    throw std::runtime_error("network_getrhs failed");
  }

  const auto &nuc = nd_.nucdata[species_index];
  std::vector<ReactionDiagnostic> diagnostics;
  diagnostics.reserve(nuc.nrates + nuc.nweakrates);

  for (size_t i = 0; i < nuc.nrates; ++i) {
    const size_t rate_index = nuc.irates[i];
    if (rate_index >= nd_.rate_count) {
      continue;
    }
    const auto &rate = nd_.rates[rate_index];
    ReactionDiagnostic item;
    item.equation = join_species(rate.input, rate.ninput, nd_) + " -> " +
                    join_species(rate.output, rate.noutput, nd_);
    item.q_value = rate.q;
    item.rate = nw_.rates[rate_index].v;
    item.abundance_weighted_rate = nw_.yrates[rate_index].v;
    item.contribution_dYdt = nuc.w[i] * item.abundance_weighted_rate;
    item.contribution_dXdt =
        item.contribution_dYdt * static_cast<double>(nuc.na);
    item.weak = false;
    diagnostics.push_back(std::move(item));
  }

  for (size_t i = 0; i < nuc.nweakrates; ++i) {
    const size_t rate_index = nuc.iweakrates[i];
    if (rate_index >= nd_.weakrate_count) {
      continue;
    }
    const auto &rate = nd_.weakrates[rate_index];
    int input[1] = {rate.input};
    int output[1] = {rate.output};
    ReactionDiagnostic item;
    item.equation =
        join_species(input, 1, nd_) + " -> " + join_species(output, 1, nd_);
    item.q_value = rate.q1;
    item.rate = nw_.weakrates[rate_index].v;
    item.abundance_weighted_rate = nw_.yweakrates[rate_index].v;
    item.contribution_dYdt = nuc.wweak[i] * item.abundance_weighted_rate;
    item.contribution_dXdt =
        item.contribution_dYdt * static_cast<double>(nuc.na);
    item.weak = true;
    diagnostics.push_back(std::move(item));
  }

  std::sort(diagnostics.begin(), diagnostics.end(),
            [](const ReactionDiagnostic &a, const ReactionDiagnostic &b) {
              return std::abs(a.contribution_dXdt) >
                     std::abs(b.contribution_dXdt);
            });
  if (diagnostics.size() > max_count) {
    diagnostics.resize(max_count);
  }
  return diagnostics;
}

std::vector<ReactionFlux> Network::get_reaction_fluxes(
    double rho, double temp, const std::vector<double> &xnuc,
    double min_strength, int metric, size_t max_count,
    bool include_regular, bool include_weak) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  validate_conditions(rho, temp);
  validate_composition(xnuc, nd_.nuc_count);
  if (!std::isfinite(min_strength) || min_strength < 0.0) {
    throw std::runtime_error("Invalid minimum flux strength");
  }
  if (metric < 0 || metric > 2) {
    throw std::runtime_error("Invalid flux metric");
  }

  std::vector<double> y(nd_.nuc_count, 0.0);
  for (size_t i = 0; i < nd_.nuc_count; ++i) {
    y[i] = xnuc[i] / static_cast<double>(nd_.nucdata[i].na);
  }

  std::vector<double> rhs(nd_.n_matrix, 0.0);
  std::vector<network_var> deriv(nd_.n_matrix);
  const int rc = network_getrhs(rho, temp, y.data(), 0, &nd_, &nw_, rhs.data(),
                                deriv.data());
  if (rc != 0) {
    throw std::runtime_error("network_getrhs failed");
  }

  std::unordered_map<long long, FluxAccumulator> flux_by_pair;
  const auto add_flux = [&](int source, int target, double strength_rate,
                            double strength_dydt, double strength_dxdt,
                            bool weak, const std::string &equation) {
    const double selected_strength = metric == 2   ? strength_rate
                                     : metric == 1 ? strength_dxdt
                                                   : strength_dydt;
    if (!std::isfinite(strength_rate) || !std::isfinite(strength_dydt) ||
        !std::isfinite(strength_dxdt) || selected_strength <= 0.0 ||
        source < 0 || target < 0 || source == target ||
        source >= static_cast<int>(nd_.nuc_count) ||
        target >= static_cast<int>(nd_.nuc_count)) {
      return;
    }
    const long long pair =
        static_cast<long long>(source) * static_cast<long long>(nd_.nuc_count) +
        static_cast<long long>(target);
    const long long key = pair * 2 + (weak ? 1 : 0);
    auto &acc = flux_by_pair[key];
    acc.strength_rate += strength_rate;
    acc.strength_dydt += strength_dydt;
    acc.strength_dxdt += strength_dxdt;
    acc.weak = weak;
    if (selected_strength >= acc.strongest_reaction || acc.equation.empty()) {
      acc.strongest_reaction = selected_strength;
      acc.equation = equation;
    }
  };

  for (size_t r = 0; r < nd_.rate_count; ++r) {
    const auto &rate = nd_.rates[r];
    const double reaction_rate = std::abs(nw_.rates[r].v);
    const double event_flow = std::abs(nw_.yrates[r].v);
    const std::string equation = join_species(rate.input, rate.ninput, nd_) +
                                 " -> " +
                                 join_species(rate.output, rate.noutput, nd_);
    for (int i = 0; i < rate.ninput; ++i) {
      bool repeated_input = false;
      for (int k = 0; k < i; ++k) {
        repeated_input = repeated_input || rate.input[k] == rate.input[i];
      }
      if (repeated_input) {
        continue;
      }
      for (int j = 0; j < rate.noutput; ++j) {
        bool repeated_output = false;
        for (int k = 0; k < j; ++k) {
          repeated_output = repeated_output || rate.output[k] == rate.output[j];
        }
        if (repeated_output) {
          continue;
        }
        int input_multiplicity = 0;
        for (int k = 0; k < rate.ninput; ++k) {
          if (rate.input[k] == rate.input[i]) {
            ++input_multiplicity;
          }
        }
        int output_multiplicity = 0;
        for (int k = 0; k < rate.noutput; ++k) {
          if (rate.output[k] == rate.output[j]) {
            ++output_multiplicity;
          }
        }
        const double dydt_flow =
            event_flow *
            static_cast<double>(std::max(input_multiplicity,
                                         output_multiplicity));
        const double dxdt_flow =
            event_flow *
            std::max(static_cast<double>(input_multiplicity *
                                         nd_.nucdata[rate.input[i]].na),
                     static_cast<double>(output_multiplicity *
                                         nd_.nucdata[rate.output[j]].na));
        add_flux(rate.input[i], rate.output[j], reaction_rate, dydt_flow,
                 dxdt_flow, rate.isWeak != 0, equation);
      }
    }
  }

  for (size_t r = 0; r < nd_.weakrate_count; ++r) {
    const auto &rate = nd_.weakrates[r];
    const double reaction_rate = std::abs(nw_.weakrates[r].v);
    const double dydt_flow = std::abs(nw_.yweakrates[r].v);
    const double dxdt_flow =
        dydt_flow *
        static_cast<double>(std::max(nd_.nucdata[rate.input].na,
                                     nd_.nucdata[rate.output].na));
    int input[1] = {rate.input};
    int output[1] = {rate.output};
    const std::string equation =
        join_species(input, 1, nd_) + " -> " + join_species(output, 1, nd_);
    add_flux(rate.input, rate.output, reaction_rate, dydt_flow, dxdt_flow, true,
             equation);
  }

  std::vector<ReactionFlux> fluxes;
  fluxes.reserve(flux_by_pair.size());
  for (const auto &[key, acc] : flux_by_pair) {
    if ((acc.weak && !include_weak) || (!acc.weak && !include_regular)) {
      continue;
    }
    const long long pair = key / 2;
    const double strength = metric == 2   ? acc.strength_rate
                            : metric == 1 ? acc.strength_dxdt
                                          : acc.strength_dydt;
    if (strength < min_strength) {
      continue;
    }
    ReactionFlux flux;
    flux.source_index =
        static_cast<int>(pair / static_cast<long long>(nd_.nuc_count));
    flux.target_index =
        static_cast<int>(pair % static_cast<long long>(nd_.nuc_count));
    flux.strength = strength;
    flux.strength_rate = acc.strength_rate;
    flux.strength_dydt = acc.strength_dydt;
    flux.strength_dxdt = acc.strength_dxdt;
    flux.weak = acc.weak;
    flux.equation = acc.equation;
    fluxes.push_back(std::move(flux));
  }

  std::sort(fluxes.begin(), fluxes.end(),
            [](const ReactionFlux &a, const ReactionFlux &b) {
              return a.strength > b.strength;
            });
  size_t regular_count = 0;
  size_t weak_count = 0;
  fluxes.erase(
      std::remove_if(fluxes.begin(), fluxes.end(),
                     [&](const ReactionFlux &flux) {
                       size_t &count = flux.weak ? weak_count : regular_count;
                       return count++ >= max_count;
                     }),
      fluxes.end());
  return fluxes;
}

bool Network::is_connected(const std::string &species_name,
                           const std::vector<bool> &active_mask) {
  if (!initialized_) {
    throw std::runtime_error("Network not initialized");
  }

  if (active_mask.size() != nd_.nuc_count) {
    throw std::runtime_error("active_mask size mismatch");
  }

  const int target = find_species_index(species_name);
  if (target < 0) {
    return false;
  }

  if (active_mask[target]) {
    return true;
  }

  const auto graph = build_connectivity_graph();
  std::vector<bool> visited(nd_.nuc_count, false);
  std::queue<int> q;

  for (size_t i = 0; i < active_mask.size(); ++i) {
    if (active_mask[i]) {
      visited[i] = true;
      q.push(static_cast<int>(i));
    }
  }

  while (!q.empty()) {
    const int u = q.front();
    q.pop();
    if (u == target) {
      return true;
    }
    for (int v : graph[u]) {
      if (!visited[v]) {
        visited[v] = true;
        q.push(v);
      }
    }
  }

  return false;
}

int Network::find_species_index(const std::string &species_name) const {
  for (size_t i = 0; i < nd_.nuc_count; ++i) {
    if (species_name == trim_species_name(nd_.nucdata[i].name)) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

std::vector<std::vector<int>> Network::build_connectivity_graph() const {
  std::vector<std::vector<int>> g(nd_.nuc_count);

  for (size_t r = 0; r < nd_.rate_count; ++r) {
    const auto &rate = nd_.rates[r];
    std::vector<int> nodes;
    nodes.reserve(static_cast<size_t>(rate.ninput + rate.noutput));
    for (int i = 0; i < rate.ninput; ++i) {
      nodes.push_back(rate.input[i]);
    }
    for (int i = 0; i < rate.noutput; ++i) {
      nodes.push_back(rate.output[i]);
    }
    for (size_t i = 0; i < nodes.size(); ++i) {
      for (size_t j = i + 1; j < nodes.size(); ++j) {
        const int a = nodes[i];
        const int b = nodes[j];
        if (a >= 0 && b >= 0 && a < static_cast<int>(nd_.nuc_count) &&
            b < static_cast<int>(nd_.nuc_count)) {
          g[a].push_back(b);
          g[b].push_back(a);
        }
      }
    }
  }

  return g;
}

} // namespace imyann

#endif
