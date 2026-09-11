#ifndef IMYANN_NETWORK_WRAPPER_H
#define IMYANN_NETWORK_WRAPPER_H

#ifdef IMNET_USE_NUPPN
#include "nuppn_c_api.h"
#else
extern "C" {
#include "helm_eos.h"
#include "network.h"
#include "network_integrate.h"
#include "network_nse.h"
#include "utilities.h"
}
#endif

#include <algorithm>
#include <string>
#include <vector>

/**
 * @brief Minimal C++ wrapper around the selected nuclear network backend
 *
 * This header provides type-safe access to the C network functions.
 * It exposes the C API without adding unnecessary abstraction layers,
 * allowing direct calls to the underlying C code.
 */

namespace imyann {

/**
 * @struct Species
 * @brief Represents a single nuclear species (isotope)
 */
struct Species {
  std::string name;
  int Z; // proton number
  int A; // mass number
  int N; // neutron number (A - Z)
};

struct ReactionDiagnostic {
  std::string equation;
  double q_value = 0.0;
  double rate = 0.0;
  double abundance_weighted_rate = 0.0;
  double contribution_dYdt = 0.0;
  double contribution_dXdt = 0.0;
  bool weak = false;
};

struct ReactionFlux {
  int source_index = -1;
  int target_index = -1;
  double strength = 0.0;
  double strength_rate = 0.0;
  double strength_dydt = 0.0;
  double strength_dxdt = 0.0;
  bool weak = false;
  std::string equation;
};

struct IntegrationStepSnapshot {
  size_t index = 0;
  double time = 0.0;
  double rho = 0.0;
  double temp = 0.0;
  double dt = 0.0;
  double dedt = 0.0;
  int substeps = 0;
  std::vector<double> xnuc;
};

/**
 * @class Network
 * @brief Wrapper around the selected nuclear network backend
 *
 * Manages C network state and provides convenient C++ access to
 * integration, NSE, and rate computation functions.
 */
class Network {
public:
  /**
   * @brief Initialize network with rate files
   * @param num_species Number of species in the network
   * @param reaclib_file Path to REACLIB rate file
   * @param partition_file Path to partition function file
   * @param mass_file Path to atomic mass file
   * @param weak_file Path to weak rate file (optional, can be empty or
   * /dev/null)
   */
  Network(const std::string &species_file, const std::string &reaclib_file,
          const std::string &partition_file, const std::string &mass_file,
          const std::string &weak_file);

  ~Network();

  // Deleted copy/move to simplify lifetime management
  Network(const Network &) = delete;
  Network &operator=(const Network &) = delete;
  Network(Network &&) = delete;
  Network &operator=(Network &&) = delete;

  /**
   * @brief Integrate the network forward in time
   * @param rho Density in g/cm³
   * @param temp Temperature in K
   * @param xnuc Species abundances (input/output)
   * @param dt Integration time in seconds
   * @return Energy release rate (dE/dt in erg/g/s)
   */
  double integrate(double rho, double temp, std::vector<double> &xnuc,
                   double dt);

  /**
   * @brief Integrate from one thermodynamic state to another over dt
   */
  double integrate_interval(double rho0, double temp0, double rho1,
                            double temp1, std::vector<double> &xnuc,
                            double dt);

  double integrate_to_time(double rho, double temp, std::vector<double> &xnuc,
                           double final_time, double initial_dt,
                           double max_dt, double dt_factor, int max_steps);

  int last_substeps() const { return last_substeps_; }
  const std::vector<IntegrationStepSnapshot> &last_step_history() const {
    return last_step_history_;
  }

  /**
   * @brief Compute NSE abundances at given conditions
   * @param rho Density in g/cm³
   * @param temp Temperature in K
   * @param ye Electron fraction (Ye = <Z>/<A>)
   * @param xnuc_out Output NSE abundances
   */
  void compute_nse(double rho, double temp, double ye,
                   std::vector<double> &xnuc_out);

  /**
   * @brief Get number of species
   */
  int num_species() const { return static_cast<int>(species_.size()); }

  static const char *backend_name();

  /**
   * @brief Get species metadata from initialized network
   */
  const std::vector<Species> &get_species() const { return species_; }

  /**
   * @brief Get reaction rate for a specific isotope at given conditions
   * @param species_index Index of the species
   * @param rho Density in g/cm³
   * @param temp Temperature in K
   * @param xnuc Current abundances
   * @return Reaction rate value
   */
  double get_reaction_rate(int species_index, double rho, double temp,
                           const std::vector<double> &xnuc);

  std::vector<ReactionDiagnostic>
  get_reaction_diagnostics(int species_index, double rho, double temp,
                           const std::vector<double> &xnuc,
                           size_t max_count = 12);

  std::vector<ReactionFlux> get_reaction_fluxes(double rho, double temp,
                                                const std::vector<double> &xnuc,
                                                double min_strength,
                                                int metric,
                                                size_t max_count = 100,
                                                bool include_regular = true,
                                                bool include_weak = true);

  /**
   * @brief Check if a species is connected (reachable) from active species
   * @param species_name Name of the species
   * @param active_mask Boolean mask of active species
   * @return true if reachable, false otherwise
   */
  bool is_connected(const std::string &species_name,
                    const std::vector<bool> &active_mask);

private:
#ifdef IMNET_USE_NUPPN
  std::string run_dir_;
#else
  struct network_data nd_;
  struct network_workspace nw_;
#endif
  bool initialized_;
  int last_substeps_ = 0;
  std::vector<IntegrationStepSnapshot> last_step_history_;
  std::vector<Species> species_;

  int find_species_index(const std::string &species_name) const;
  std::vector<std::vector<int>> build_connectivity_graph() const;
};

/**
 * @brief Utility: Calculate average mass number from abundances
 * @param xnuc Species abundances
 * @param species_array Network species array
 * @return Average mass number <A>
 */
inline double calculate_abar(const std::vector<double> &xnuc,
                             const std::vector<Species> &species) {
  double mass_fraction_sum = 0.0;
  double molar_abundance_sum = 0.0;
  const size_t count = std::min(xnuc.size(), species.size());
  for (size_t i = 0; i < count; ++i) {
    if (species[i].A > 0) {
      mass_fraction_sum += xnuc[i];
      molar_abundance_sum += xnuc[i] / static_cast<double>(species[i].A);
    }
  }
  return molar_abundance_sum > 0.0
             ? mass_fraction_sum / molar_abundance_sum
             : 0.0;
}

/**
 * @brief Utility: Calculate average proton number from abundances
 * @param xnuc Species abundances
 * @param species_array Network species array
 * @return Average proton number <Z>
 */
inline double calculate_zbar(const std::vector<double> &xnuc,
                             const std::vector<Species> &species) {
  double molar_abundance_sum = 0.0;
  double proton_abundance_sum = 0.0;
  const size_t count = std::min(xnuc.size(), species.size());
  for (size_t i = 0; i < count; ++i) {
    if (species[i].A > 0) {
      const double abundance = xnuc[i] / static_cast<double>(species[i].A);
      molar_abundance_sum += abundance;
      proton_abundance_sum += abundance * species[i].Z;
    }
  }
  return molar_abundance_sum > 0.0
             ? proton_abundance_sum / molar_abundance_sum
             : 0.0;
}

/**
 * @brief Utility: Calculate electron fraction Ye = <Z>/<A>
 */
inline double calculate_ye(const std::vector<double> &xnuc,
                           const std::vector<Species> &species) {
  double mass_fraction_sum = 0.0;
  double proton_abundance_sum = 0.0;
  const size_t count = std::min(xnuc.size(), species.size());
  for (size_t i = 0; i < count; ++i) {
    if (species[i].A > 0) {
      mass_fraction_sum += xnuc[i];
      proton_abundance_sum +=
          xnuc[i] * species[i].Z / static_cast<double>(species[i].A);
    }
  }
  return mass_fraction_sum > 0.0
             ? proton_abundance_sum / mass_fraction_sum
             : 0.0;
}

} // namespace imyann

#endif // IMYANN_NETWORK_WRAPPER_H
