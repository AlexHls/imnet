#ifndef IMYANN_FILE_IO_H
#define IMYANN_FILE_IO_H

#include <string>
#include <vector>

namespace imyann {

/**
 * @brief Save abundances to a file
 * @param filename Output file path
 * @param species_names Vector of species names
 * @param xnuc Abundance vector
 * @return true if successful, false otherwise
 */
bool save_abundances(const std::string &filename,
                     const std::vector<std::string> &species_names,
                     const std::vector<double> &xnuc);

/**
 * @brief Load abundances from a file
 * @param filename Input file path
 * @param species_names Vector of species names for reference
 * @param xnuc_out Output: abundance vector
 * @return true if successful, false otherwise
 */
bool load_abundances(const std::string &filename,
                     const std::vector<std::string> &species_names,
                     std::vector<double> &xnuc_out);

/**
 * @brief Load a (rho, T) trajectory file for multi-timestep integration
 * @param filename Input file path
 * @param times Output: vector of times (seconds)
 * @param rhos Output: vector of densities (g/cm³)
 * @param temps Output: vector of temperatures (K)
 * @return true if successful, false otherwise
 */
bool load_trajectory(const std::string &filename, std::vector<double> &times,
                     std::vector<double> &rhos, std::vector<double> &temps);

/**
 * @brief Save a (time, rho, T) trajectory file
 * @param filename Output file path
 * @param times Vector of times (seconds)
 * @param rhos Vector of densities (g/cm³)
 * @param temps Vector of temperatures (K)
 * @return true if successful, false otherwise
 */
bool save_trajectory(const std::string &filename,
                     const std::vector<double> &times,
                     const std::vector<double> &rhos,
                     const std::vector<double> &temps);

} // namespace imyann

#endif // IMYANN_FILE_IO_H
