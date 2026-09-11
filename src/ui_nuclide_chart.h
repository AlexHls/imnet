#ifndef IMYANN_UI_NUCLIDE_CHART_H
#define IMYANN_UI_NUCLIDE_CHART_H

#include "app_state.h"
#include "imgui.h"
#include <map>
#include <string>
#include <vector>

namespace imyann {

/**
 * @struct IsotopeButton
 * @brief Represents a single isotope button in the chart
 */
struct IsotopeButton {
  int index;          ///< Index in species array
  int Z;              ///< Proton number
  int N;              ///< Neutron number (A - Z)
  std::string name;   ///< Species name
};

/**
 * @class NuclideChart
 * @brief Renders the chart of nuclides (Z vs N grid)
 *
 * Manages the interactive grid of isotope buttons, color-coding,
 * and interaction handlers.
 */
class NuclideChart {
public:
  NuclideChart();

  /**
   * @brief Initialize chart from app state
   */
  void initialize(AppState *app_state);

  /**
   * @brief Render the chart
   */
  void render();

  /**
   * @brief Get hovered isotope index (-1 if none)
   */
  int get_hovered_isotope() const { return hovered_isotope_; }

  /**
   * @brief Get selected isotopes
   */
  const std::vector<bool> &get_selected() const { return selected_; }

  /**
   * @brief Replace selected isotope mask
   */
  void set_selected(const std::vector<bool> &selected);

  /**
   * @brief Get one selected isotope index (-1 if none)
   */
  int get_single_selected_isotope() const;

private:
  AppState *app_state_;
  std::vector<IsotopeButton> buttons_;
  std::map<int, IsotopeButton> button_map_; // index -> button

  int hovered_isotope_;
  std::vector<bool> selected_;

  int z_min_, z_max_;
  int n_min_, n_max_;

  float zoom_;

  /**
   * @brief Build button grid from app state
   */
  void build_grid();

  /**
   * @brief Compute color for isotope based on abundance
   */
  ImVec4 get_abundance_color(int index) const;
};

} // namespace imyann

#endif // IMYANN_UI_NUCLIDE_CHART_H
