#ifndef IMNET_UI_PLOT_HELPERS_H
#define IMNET_UI_PLOT_HELPERS_H

#include "imgui.h"
#include "implot.h"
#include <cmath>
#include <limits>
#include <unordered_map>

namespace imyann {

// Only publish numeric edits when Enter is pressed or focus leaves the field.
// Keeping drafts separate also allows validation to reject an edit safely.
inline bool input_double_committed(const char *label, double *value,
                                   double step = 0.0, double step_fast = 0.0,
                                   const char *format = "%.6e") {
  static std::unordered_map<ImGuiID, double> drafts;
  const auto id = ImGui::GetID(label);
  auto &draft = drafts.try_emplace(id, *value).first->second;
  ImGui::InputDouble(label, &draft, step, step_fast, format);
  const bool commit = ImGui::IsItemDeactivatedAfterEdit();
  const bool changed = commit && std::isfinite(draft) && draft != *value;
  if (changed) *value = draft;
  if (!ImGui::IsItemActive()) drafts.erase(id);
  return changed;
}

// Keep partial keyboard input separate from the floor used by the plot.
struct LogPlotFloor {
  double value;
  double draft;

  bool draw(const char *label) {
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble(label, &draft, 0.0, 0.0, "%.3e");
    if (!ImGui::IsItemDeactivatedAfterEdit()) return false;
    if (!std::isfinite(draft) || draft <= 0.0) {
      draft = value;
      return false;
    }
    const bool changed = draft != value;
    value = draft;
    return changed;
  }
};

inline void setup_log_plot_axis(ImAxis axis, double floor) {
  ImPlot::SetupAxisScale(axis, ImPlotScale_Log10);
  // ImPlot fit padding is linear and can cross zero even for positive data.
  ImPlot::SetupAxisLimitsConstraints(axis, floor,
                                    std::numeric_limits<double>::max());
}

inline void plot_current_time(const char *label, double time) {
  // InfLines spans the visible Y range. A cursor must never enlarge a fit.
  ImPlot::PlotInfLines(label, &time, 1,
                      {ImPlotProp_Flags, ImPlotItemFlags_NoFit});
}

} // namespace imyann
#endif
