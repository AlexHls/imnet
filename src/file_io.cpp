#include "file_io.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <utility>

namespace imyann {
namespace {

enum class TrajectoryColumnOrder {
  TimeRhoTemp,
  TimeTempRho,
};

std::string trim_copy(const std::string &value) {
  auto begin = std::find_if_not(value.begin(), value.end(), [](char c) {
    return std::isspace(static_cast<unsigned char>(c));
  });
  auto end = std::find_if_not(value.rbegin(), value.rend(), [](char c) {
               return std::isspace(static_cast<unsigned char>(c));
             }).base();

  if (begin >= end) {
    return {};
  }
  return std::string(begin, end);
}

std::string uppercase_copy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](char c) {
    return static_cast<char>(
        std::toupper(static_cast<unsigned char>(c)));
  });
  return value;
}

std::string lowercase_alnum_copy(const std::string &value) {
  std::string normalized;
  for (char c : value) {
    if (std::isalnum(static_cast<unsigned char>(c))) {
      normalized.push_back(
          static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
  }
  return normalized;
}

bool starts_with(const std::string &value, const std::string &prefix) {
  return value.size() >= prefix.size() &&
         std::equal(prefix.begin(), prefix.end(), value.begin());
}

bool has_only_comment_remaining(std::istringstream &stream) {
  std::string extra;
  return !(stream >> extra) || (!extra.empty() && extra[0] == '#');
}

bool is_time_header_token(const std::string &token) {
  return starts_with(token, "time") || starts_with(token, "age");
}

bool is_temp_header_token(const std::string &token) {
  return token == "t" || starts_with(token, "temp") ||
         starts_with(token, "temperature") || starts_with(token, "t9");
}

bool is_rho_header_token(const std::string &token) {
  return starts_with(token, "rho") || starts_with(token, "density");
}

bool parse_trajectory_column_header(const std::string &line,
                                    TrajectoryColumnOrder &order) {
  std::istringstream stream(line);
  std::vector<std::string> columns;
  std::string token;

  while (stream >> token) {
    token = lowercase_alnum_copy(token);
    if (token.empty()) {
      continue;
    }

    if (is_time_header_token(token)) {
      columns.push_back("time");
    } else if (is_temp_header_token(token)) {
      columns.push_back("temp");
    } else if (is_rho_header_token(token)) {
      columns.push_back("rho");
    }
  }

  for (size_t i = 0; i + 2 < columns.size(); ++i) {
    if (columns[i] != "time") {
      continue;
    }
    if (columns[i + 1] == "rho" && columns[i + 2] == "temp") {
      order = TrajectoryColumnOrder::TimeRhoTemp;
      return true;
    }
    if (columns[i + 1] == "temp" && columns[i + 2] == "rho") {
      order = TrajectoryColumnOrder::TimeTempRho;
      return true;
    }
  }

  return false;
}

bool parse_trajectory_metadata(const std::string &line,
                               TrajectoryColumnOrder &order,
                               bool has_column_header,
                               double &temperature_scale,
                               std::string &error_message) {
  std::istringstream stream(line);
  std::string key;
  stream >> key;
  key = uppercase_copy(key);

  if (key == "ID") {
    return true;
  }

  if (key == "AGEUNIT") {
    std::string unit;
    if (!(stream >> unit)) {
      error_message = "AGEUNIT metadata is missing a unit";
      return false;
    }
    unit = uppercase_copy(unit);
    if (unit != "SEC" && unit != "S" && unit != "SECOND" &&
        unit != "SECONDS") {
      error_message =
          "Unsupported trajectory time unit '" + unit + "'; expected seconds";
      return false;
    }
    return true;
  }

  if (key == "RHOUNIT") {
    std::string unit;
    if (!(stream >> unit)) {
      error_message = "RHOUNIT metadata is missing a unit";
      return false;
    }
    unit = uppercase_copy(unit);
    if (unit != "CGS" && unit != "G/CM3" && unit != "G/CM^3") {
      error_message = "Unsupported trajectory density unit '" + unit +
                      "'; expected CGS g/cm^3";
      return false;
    }
    return true;
  }

  if (key == "TUNIT") {
    std::string unit;
    if (!(stream >> unit)) {
      error_message = "TUNIT metadata is missing a unit";
      return false;
    }
    unit = uppercase_copy(unit);
    if (unit == "T9K" || unit == "T9" || unit == "GK" ||
        unit == "GIGAKELVIN") {
      temperature_scale = 1.0e9;
    } else if (unit == "K" || unit == "KELVIN") {
      temperature_scale = 1.0;
    } else {
      error_message = "Unsupported trajectory temperature unit '" + unit +
                      "'; expected K or T9/GK";
      return false;
    }

    if (!has_column_header) {
      order = TrajectoryColumnOrder::TimeTempRho;
    }
    return true;
  }

  return false;
}

} // namespace

bool save_abundances(const std::string &filename,
                     const std::vector<std::string> &species_names,
                     const std::vector<double> &xnuc) {
  if (species_names.size() != xnuc.size()) {
    std::cerr << "Error: Cannot save mismatched species and abundance arrays"
              << std::endl;
    return false;
  }
  for (size_t i = 0; i < xnuc.size(); ++i) {
    if (!std::isfinite(xnuc[i]) || xnuc[i] < 0.0) {
      std::cerr << "Error: Cannot save invalid abundance for "
                << species_names[i] << std::endl;
      return false;
    }
  }

  std::ofstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open output file: " << filename << std::endl;
    return false;
  }

  file << "# species mass_fraction\n";
  file << std::setprecision(17) << std::scientific;
  for (size_t i = 0; i < xnuc.size(); ++i) {
    file << species_names[i] << " " << xnuc[i] << "\n";
  }

  file.close();
  if (!file) {
    std::cerr << "Error: Failed while writing abundance file: " << filename
              << std::endl;
    return false;
  }
  std::cout << "Saved abundances to " << filename << std::endl;
  return true;
}

bool load_abundances(const std::string &filename,
                     const std::vector<std::string> &species_names,
                     std::vector<double> &xnuc_out) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open abundance file: " << filename
              << std::endl;
    return false;
  }

  std::unordered_map<std::string, size_t> species_indices;
  species_indices.reserve(species_names.size());
  for (size_t i = 0; i < species_names.size(); ++i) {
    species_indices.emplace(species_names[i], i);
  }

  std::vector<double> parsed_xnuc(species_names.size(), 0.0);
  std::vector<bool> seen(species_names.size(), false);
  size_t matched_species_count = 0;

  std::string line;
  size_t line_number = 0;
  while (std::getline(file, line)) {
    ++line_number;
    line = trim_copy(line);
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::string key;
    double value;

    if (!(iss >> key >> value) || !has_only_comment_remaining(iss)) {
      std::cerr << "Error: Could not parse abundance line " << line_number
                << ": " << line << std::endl;
      return false;
    }

    if (!std::isfinite(value)) {
      std::cerr << "Error: Non-finite abundance file value on line: " << line
                << std::endl;
      return false;
    }

    if (key == "dt" || key == "rho" || key == "temp" || key == "T") {
      std::cerr << "Error: Abundance files contain only species mass "
                   "fractions; found thermodynamic key '"
                << key << "' on line " << line_number << " of " << filename
                << std::endl;
      return false;
    }

    const auto it = species_indices.find(key);
    if (it != species_indices.end()) {
      if (value < 0.0) {
        std::cerr << "Error: Negative abundance for species " << key << " in "
                  << filename << std::endl;
        return false;
      }
      const size_t index = it->second;
      if (seen[index]) {
        std::cerr << "Warning: Duplicate abundance for species " << key
                  << " in "
                  << filename << std::endl;
      } else {
        seen[index] = true;
        ++matched_species_count;
      }
      parsed_xnuc[index] = value;
    } else {
      std::cerr << "Warning: Unknown species in abundance file: " << key
                << std::endl;
    }
  }

  if (matched_species_count == 0) {
    std::cerr << "Error: No known species abundances found in " << filename
              << std::endl;
    return false;
  }

  xnuc_out = std::move(parsed_xnuc);
  std::cout << "Loaded abundances from " << filename << std::endl;
  return true;
}

bool load_trajectory(const std::string &filename, std::vector<double> &times,
                     std::vector<double> &rhos, std::vector<double> &temps) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open trajectory file: " << filename
              << std::endl;
    return false;
  }

  std::vector<double> parsed_times;
  std::vector<double> parsed_rhos;
  std::vector<double> parsed_temps;

  std::string line;
  size_t line_number = 0;
  TrajectoryColumnOrder column_order = TrajectoryColumnOrder::TimeRhoTemp;
  bool has_column_header = false;
  double temperature_scale = 1.0;
  while (std::getline(file, line)) {
    ++line_number;
    line = trim_copy(line);

    if (line.empty()) {
      continue;
    }

    if (line[0] == '#') {
      TrajectoryColumnOrder parsed_order = column_order;
      if (parse_trajectory_column_header(line.substr(1), parsed_order)) {
        column_order = parsed_order;
        has_column_header = true;
      }
      continue;
    }

    std::string metadata_error;
    if (parse_trajectory_metadata(line, column_order, has_column_header,
                                  temperature_scale, metadata_error)) {
      continue;
    }
    if (!metadata_error.empty()) {
      std::cerr << "Error: " << metadata_error << " on trajectory line "
                << line_number << std::endl;
      return false;
    }

    std::istringstream iss(line);
    double t, second, third;

    if (!(iss >> t >> second >> third) || !has_only_comment_remaining(iss)) {
      std::cerr << "Error: Could not parse trajectory line " << line_number
                << ": " << line << std::endl;
      return false;
    }

    double rho = second;
    double temp = third;
    if (column_order == TrajectoryColumnOrder::TimeTempRho) {
      temp = second;
      rho = third;
    }
    temp *= temperature_scale;

    if (!std::isfinite(t) || !std::isfinite(rho) || !std::isfinite(temp)) {
      std::cerr << "Error: Non-finite value on trajectory line " << line_number
                << std::endl;
      return false;
    }
    if (rho <= 0.0 || temp <= 0.0) {
      std::cerr << "Error: rho and temp must be positive on trajectory line "
                << line_number << std::endl;
      return false;
    }
    if (!parsed_times.empty() && t < parsed_times.back()) {
      std::cerr << "Error: Trajectory time decreases on line " << line_number
                << std::endl;
      return false;
    }

    parsed_times.push_back(t);
    parsed_rhos.push_back(rho);
    parsed_temps.push_back(temp);
  }

  if (parsed_times.empty()) {
    std::cerr << "Error: No trajectory data loaded from " << filename
              << std::endl;
    return false;
  }

  times = std::move(parsed_times);
  rhos = std::move(parsed_rhos);
  temps = std::move(parsed_temps);
  std::cout << "Loaded trajectory with " << times.size()
            << " timesteps from " << filename << std::endl;
  return true;
}

bool save_trajectory(const std::string &filename,
                     const std::vector<double> &times,
                     const std::vector<double> &rhos,
                     const std::vector<double> &temps) {
  if (times.empty() || times.size() != rhos.size() ||
      times.size() != temps.size()) {
    std::cerr << "Error: Cannot save invalid trajectory arrays" << std::endl;
    return false;
  }

  for (size_t i = 0; i < times.size(); ++i) {
    if (!std::isfinite(times[i]) || !std::isfinite(rhos[i]) ||
        !std::isfinite(temps[i])) {
      std::cerr << "Error: Cannot save trajectory with non-finite value at row "
                << i << std::endl;
      return false;
    }
    if (rhos[i] <= 0.0 || temps[i] <= 0.0) {
      std::cerr << "Error: Cannot save trajectory with non-positive rho/T at "
                   "row "
                << i << std::endl;
      return false;
    }
    if (i > 0 && times[i] < times[i - 1]) {
      std::cerr << "Error: Cannot save trajectory with decreasing time at row "
                << i << std::endl;
      return false;
    }
  }

  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open trajectory output file: " << filename
              << std::endl;
    return false;
  }

  file << "# time(s) rho(g/cm^3) temp(K)\n";
  file << std::setprecision(17) << std::scientific;
  for (size_t i = 0; i < times.size(); ++i) {
    file << times[i] << " " << rhos[i] << " " << temps[i] << "\n";
  }

  file.close();
  if (!file) {
    std::cerr << "Error: Failed while writing trajectory file: " << filename
              << std::endl;
    return false;
  }
  std::cout << "Saved trajectory to " << filename << std::endl;
  return true;
}

} // namespace imyann
