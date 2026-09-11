#include "ui_nuclide_chart.h"
#include "ui_theme.h"
#include "imgui.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace imyann {

namespace {

struct DrawnFlux {
  ReactionFlux flux;
  ImVec2 start;
  ImVec2 end;
  float relative_strength = 0.0f;
  float color_strength = 0.0f;
  float thickness = 1.0f;
};

ImVec4 lerp_color(const ImVec4 &a, const ImVec4 &b, float t) {
  return ImVec4(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                a.z + (b.z - a.z) * t, a.w + (b.w - a.w) * t);
}

ImVec4 flux_colormap_color(int colormap, float t, float alpha) {
  t = std::clamp(t, 0.0f, 1.0f);
  ImVec4 color;
  switch (colormap) {
  case 1:
    color = t < 0.5f
                ? lerp_color(ui_theme::rgb(64, 50, 118, alpha),
                             ui_theme::rgb(219, 75, 118, alpha), t / 0.5f)
                : lerp_color(ui_theme::rgb(219, 75, 118, alpha),
                             ui_theme::warning(alpha),
                             (t - 0.5f) / 0.5f);
    break;
  case 2:
    color = t < 0.5f
                ? lerp_color(ui_theme::rgb(18, 42, 50, alpha),
                             ui_theme::primary(alpha), t / 0.5f)
                : lerp_color(ui_theme::primary(alpha),
                             ui_theme::rgb(178, 218, 92, alpha),
                             (t - 0.5f) / 0.5f);
    break;
  default:
    color = lerp_color(ui_theme::primary(alpha), ui_theme::danger(alpha), t);
    break;
  }
  color.w = alpha;
  return color;
}

ImU32 isotope_text_color(const ImVec4 &fill) {
  const float luminance = 0.2126f * fill.x + 0.7152f * fill.y +
                          0.0722f * fill.z;
  return ImGui::GetColorU32(luminance > 0.55f ? ui_theme::rgb(7, 12, 14)
                                              : ui_theme::text());
}

} // namespace

NuclideChart::NuclideChart()
    : app_state_(nullptr), hovered_isotope_(-1), z_min_(0), z_max_(118),
      n_min_(0), n_max_(178), zoom_(1.0f) {}

void NuclideChart::initialize(AppState *app_state) {
  app_state_ = app_state;
  if (app_state_) {
    build_grid();
  }
}

void NuclideChart::set_selected(const std::vector<bool> &selected) {
  if (selected.size() == selected_.size()) {
    selected_ = selected;
  }
}

void NuclideChart::build_grid() {
  if (!app_state_)
    return;

  buttons_.clear();
  button_map_.clear();
  selected_.assign(app_state_->num_species(), false);

  // Build Z and N ranges and isotope buttons
  z_min_ = 999;
  z_max_ = -1;
  n_min_ = 999;
  n_max_ = -1;

  const auto &species = app_state_->get_species();
  for (size_t i = 0; i < species.size(); ++i) {
    const int z = species[i].Z;
    const int n = species[i].N;
    if (z < 0 || n < 0) {
      continue;
    }

    z_min_ = std::min(z_min_, z);
    z_max_ = std::max(z_max_, z);
    n_min_ = std::min(n_min_, n);
    n_max_ = std::max(n_max_, n);

    IsotopeButton btn;
    btn.index = static_cast<int>(i);
    btn.Z = z;
    btn.N = n;
    btn.name = species[i].name;
    buttons_.push_back(btn);
    button_map_[btn.index] = btn;
  }

  if (buttons_.empty()) {
    z_min_ = 0;
    z_max_ = 1;
    n_min_ = 0;
    n_max_ = 1;
  }
}

void NuclideChart::render() {
  if (!app_state_) {
    ImGui::TextDisabled("Chart not initialized");
    return;
  }

  const auto &state_selected = app_state_->integration_settings().selected;
  if (state_selected.size() == selected_.size()) {
    selected_ = state_selected;
  }

  bool reset_view = false;
  ImGui::SetNextItemWidth(150.0f);
  ImGui::SliderFloat("Zoom", &zoom_, 0.5f, 2.5f, "%.2fx");
  ImGui::SameLine();
  if (ImGui::Button("Reset View")) {
    zoom_ = 1.0f;
    reset_view = true;
  }

  auto &view_settings = app_state_->view_settings();
  ImGui::Checkbox("Flux arrows", &view_settings.show_fluxes);
  if (view_settings.show_fluxes) {
    ImGui::SameLine();
    ImGui::Checkbox("Regular", &view_settings.show_regular_fluxes);
    ImGui::SameLine();
    ImGui::Checkbox("Weak", &view_settings.show_weak_fluxes);

    const char *arrow_metrics[] = {"|dY/dt|", "|dX/dt|", "Rate"};
    const char *color_modes[] = {"Flat", "Linear", "Log"};
    if (ImGui::BeginTable("FluxControls", 3,
                          ImGuiTableFlags_SizingFixedFit)) {
      ImGui::TableSetupColumn("Metric", ImGuiTableColumnFlags_WidthFixed,
                              120.0f);
      ImGui::TableSetupColumn("Minimum", ImGuiTableColumnFlags_WidthFixed,
                              120.0f);
      ImGui::TableSetupColumn("Maximum", ImGuiTableColumnFlags_WidthFixed,
                              128.0f);
      ImGui::TableNextRow();
      ImGui::TableNextColumn();
      ImGui::TextUnformatted("Metric");
      ImGui::SetNextItemWidth(112.0f);
      ImGui::Combo("##flux_metric", &view_settings.flux_arrow_metric,
                   arrow_metrics, IM_ARRAYSIZE(arrow_metrics));
      ImGui::TableNextColumn();
      ImGui::TextUnformatted("Minimum");
      ImGui::SetNextItemWidth(112.0f);
      ImGui::InputDouble("##flux_min", &view_settings.flux_threshold, 0.0, 0.0,
                         "%.1e");
      ImGui::TableNextColumn();
      ImGui::TextUnformatted("Arrows/type");
      ImGui::SetNextItemWidth(120.0f);
      ImGui::SliderInt("##flux_max", &view_settings.max_flux_arrows, 1, 250);

      ImGui::TableNextRow();
      ImGui::TableNextColumn();
      ImGui::TextUnformatted("Color scale");
      ImGui::SetNextItemWidth(112.0f);
      ImGui::Combo("##flux_color_mode", &view_settings.flux_color_mode,
                   color_modes, IM_ARRAYSIZE(color_modes));
      ImGui::TableNextColumn();
      if (view_settings.flux_color_mode == 0) {
        ImGui::TextUnformatted("Regular");
        ImGui::ColorEdit4("##flux_regular", view_settings.flux_color,
                          ImGuiColorEditFlags_NoInputs |
                              ImGuiColorEditFlags_AlphaBar);
      } else {
        const char *colormaps[] = {"Blue-red", "Plasma", "Viridis"};
        ImGui::TextUnformatted("Map");
        ImGui::SetNextItemWidth(112.0f);
        ImGui::Combo("##flux_colormap", &view_settings.flux_colormap,
                     colormaps, IM_ARRAYSIZE(colormaps));
      }
      ImGui::TableNextColumn();
      if (view_settings.flux_color_mode == 0) {
        ImGui::TextUnformatted("Weak");
        ImGui::ColorEdit4("##flux_weak", view_settings.weak_flux_color,
                          ImGuiColorEditFlags_NoInputs |
                              ImGuiColorEditFlags_AlphaBar);
      } else {
        const char *positions[] = {"Top left", "Top right", "Bottom left",
                                   "Bottom right"};
        ImGui::TextUnformatted("Bar");
        ImGui::SetNextItemWidth(120.0f);
        ImGui::Combo("##flux_colorbar",
                     &view_settings.flux_colorbar_position, positions,
                     IM_ARRAYSIZE(positions));
      }
      ImGui::EndTable();
    }
    if (view_settings.flux_threshold < 0.0) {
      view_settings.flux_threshold = 0.0;
    }
    if (!app_state_->has_reaction_flux_state()) {
      ImGui::TextDisabled("Fluxes update after a run or cached-step selection.");
    } else {
      ImGui::TextDisabled("Flux state: rho %.3e, T %.3e",
                          app_state_->reaction_flux_state_rho(),
                          app_state_->reaction_flux_state_temp());
    }
  }
  ImGui::Separator();

  ImGui::BeginChild("ChartArea", ImVec2(0, 0), true,
                    ImGuiWindowFlags_HorizontalScrollbar);
  if (reset_view) {
    ImGui::SetScrollX(0.0f);
    ImGui::SetScrollY(0.0f);
  }

  hovered_isotope_ = -1;

  const float cell = 22.0f * zoom_;
  const float left_margin = 64.0f;
  const float top_margin = 56.0f;
  const float right_margin = 20.0f;
  const float bottom_margin = 20.0f;
  const ImVec2 canvas_pos = ImGui::GetCursorScreenPos();
  const ImVec2 virtual_size(
      left_margin + (n_max_ - n_min_ + 1) * cell + right_margin,
      top_margin + (z_max_ - z_min_ + 1) * cell + bottom_margin);
  const ImVec2 origin(canvas_pos.x + left_margin, canvas_pos.y + top_margin);
  const ImVec2 mouse = ImGui::GetIO().MousePos;
  const bool chart_hovered = ImGui::IsWindowHovered();
  ImDrawList *draw_list = ImGui::GetWindowDrawList();

  ImU32 grid_color = ImGui::GetColorU32(ui_theme::border_light(0.44f));
  const ImU32 axis_color = ImGui::GetColorU32(ui_theme::text_muted());

  for (int n = n_min_; n <= n_max_; ++n) {
    float x = origin.x + (n - n_min_) * cell;
    draw_list->AddLine(ImVec2(x, origin.y),
                       ImVec2(x, origin.y + (z_max_ - z_min_ + 1) * cell),
                       grid_color);
    if ((n - n_min_) % 5 == 0) {
      draw_list->AddText(ImVec2(x + 2.0f, canvas_pos.y + 30.0f),
                         axis_color,
                         std::to_string(n).c_str());
    }
  }

  for (int z = z_min_; z <= z_max_; ++z) {
    float y = origin.y + (z - z_min_) * cell;
    draw_list->AddLine(ImVec2(origin.x, y),
                       ImVec2(origin.x + (n_max_ - n_min_ + 1) * cell, y),
                       grid_color);
    if ((z - z_min_) % 5 == 0) {
      draw_list->AddText(ImVec2(canvas_pos.x + 6.0f, y + 2.0f),
                         axis_color,
                         std::to_string(z).c_str());
    }
  }

  draw_list->AddText(ImVec2(canvas_pos.x + left_margin, canvas_pos.y + 8.0f),
                     axis_color, "N (neutrons) ->");
  draw_list->AddText(ImVec2(canvas_pos.x + 6.0f,
                            canvas_pos.y + top_margin - 22.0f),
                     axis_color, "Z (protons)");

  std::vector<DrawnFlux> drawn_fluxes;
  const int selected_isotope = get_single_selected_isotope();
  double flux_min_log = 0.0;
  double flux_max_log = 0.0;
  double flux_min_value = 0.0;
  double flux_max_value = 0.0;

  if (view_settings.show_fluxes) {
    const auto fluxes = app_state_->get_reaction_fluxes(
        view_settings.flux_threshold, view_settings.flux_arrow_metric,
        static_cast<size_t>(std::max(1, view_settings.max_flux_arrows)),
        view_settings.show_regular_fluxes, view_settings.show_weak_fluxes);

    double strongest_flux = 0.0;
    for (const auto &flux : fluxes) {
      strongest_flux = std::max(strongest_flux, flux.strength);
    }

    const double min_log =
        std::log10(std::max(view_settings.flux_threshold, 1e-300));
    const double max_log = std::log10(std::max(strongest_flux, 1e-300));
    const double log_span = std::max(max_log - min_log, 1.0);
    flux_min_log = min_log;
    flux_max_log = max_log;
    flux_min_value = view_settings.flux_threshold;
    flux_max_value = strongest_flux;
    const double linear_span =
        std::max(strongest_flux - view_settings.flux_threshold, 1e-300);

    for (const auto &flux : fluxes) {
      const auto source_it = button_map_.find(flux.source_index);
      const auto target_it = button_map_.find(flux.target_index);
      if (source_it == button_map_.end() || target_it == button_map_.end()) {
        continue;
      }

      const auto &source = source_it->second;
      const auto &target = target_it->second;
      ImVec2 p0(origin.x + (source.N - n_min_ + 0.5f) * cell,
                origin.y + (source.Z - z_min_ + 0.5f) * cell);
      ImVec2 p1(origin.x + (target.N - n_min_ + 0.5f) * cell,
                origin.y + (target.Z - z_min_ + 0.5f) * cell);
      const float dx = p1.x - p0.x;
      const float dy = p1.y - p0.y;
      const float len = std::sqrt(dx * dx + dy * dy);
      if (len < 1.0f) {
        continue;
      }

      const float ux = dx / len;
      const float uy = dy / len;
      p0.x += ux * cell * 0.35f;
      p0.y += uy * cell * 0.35f;
      p1.x -= ux * cell * 0.35f;
      p1.y -= uy * cell * 0.35f;

      const double flux_log = std::log10(std::max(flux.strength, 1e-300));
      const float t =
          static_cast<float>(std::clamp((flux_log - min_log) / log_span, 0.0,
                                        1.0));
      const float color_t =
          view_settings.flux_color_mode == 1
              ? static_cast<float>(
                    std::clamp((flux.strength - view_settings.flux_threshold) /
                                   linear_span,
                               0.0, 1.0))
              : t;
      DrawnFlux drawn;
      drawn.flux = flux;
      drawn.start = p0;
      drawn.end = p1;
      drawn.relative_strength = t;
      drawn.color_strength = color_t;
      drawn.thickness = 1.0f + 5.0f * t;
      drawn_fluxes.push_back(std::move(drawn));
    }
  }

  for (const auto &btn : buttons_) {
    const float x0 = origin.x + (btn.N - n_min_) * cell;
    const float y0 = origin.y + (btn.Z - z_min_) * cell;
    const ImVec2 p0(x0 + 1.0f, y0 + 1.0f);
    const ImVec2 p1(x0 + cell - 1.0f, y0 + cell - 1.0f);

    const bool contains_mouse = (mouse.x >= p0.x && mouse.x <= p1.x &&
                                 mouse.y >= p0.y && mouse.y <= p1.y);
    if (contains_mouse && chart_hovered) {
      hovered_isotope_ = btn.index;
      if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        const bool multi_select = ImGui::GetIO().KeyCtrl;
        if (!multi_select) {
          std::fill(selected_.begin(), selected_.end(), false);
        }
        if (btn.index >= 0 && btn.index < static_cast<int>(selected_.size())) {
          selected_[btn.index] = !selected_[btn.index];
        }
      }
    }

    const ImVec4 c = get_abundance_color(btn.index);
    const bool is_selected =
        btn.index >= 0 && btn.index < static_cast<int>(selected_.size())
            ? selected_[btn.index]
            : false;
    const auto &settings = app_state_->integration_settings();
    const bool active =
        btn.index >= 0 && btn.index < static_cast<int>(settings.active.size())
            ? settings.active[btn.index]
            : true;
    const bool show_unconnected = app_state_->view_settings().show_unconnected;
    const bool connected =
        show_unconnected ? app_state_->is_species_connected(btn.index) : true;

    ImVec4 fill = c;
    if (show_unconnected && !connected) {
      fill = ui_theme::rgb(62, 26, 25);
    } else if (!active) {
      fill = lerp_color(ui_theme::panel(), fill, 0.35f);
    }

    const ImU32 fill_col = ImGui::GetColorU32(fill);
    const ImU32 border_col =
        is_selected
            ? ImGui::GetColorU32(ui_theme::primary())
            : (show_unconnected && !connected
                   ? ImGui::GetColorU32(ui_theme::danger())
                   : ImGui::GetColorU32(ui_theme::background()));
    draw_list->AddRectFilled(p0, p1, fill_col, 2.0f);
    draw_list->AddRect(p0, p1, border_col, 2.0f, 0, is_selected ? 2.0f : 1.0f);

    if (zoom_ >= 1.0f) {
      draw_list->AddText(ImVec2(p0.x + 2.0f, p0.y + 2.0f),
                         isotope_text_color(fill), btn.name.c_str());
    }
  }

  for (const auto &drawn : drawn_fluxes) {
    const float dx = drawn.end.x - drawn.start.x;
    const float dy = drawn.end.y - drawn.start.y;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1.0f) {
      continue;
    }

    const float ux = dx / len;
    const float uy = dy / len;
    const bool connected_to_selection =
        selected_isotope >= 0 &&
        (drawn.flux.source_index == selected_isotope ||
         drawn.flux.target_index == selected_isotope);
    const bool dim_for_selection =
        selected_isotope >= 0 && !connected_to_selection;
    const float thickness =
        dim_for_selection
            ? std::max(1.0f, drawn.thickness * 0.55f)
            : (connected_to_selection ? drawn.thickness + 2.5f
                                      : drawn.thickness);
    const float head_len =
        std::max(5.0f, 7.0f + 5.0f * drawn.relative_strength) * zoom_;
    const float head_half_width = std::max(3.0f, thickness * 1.4f);
    ImVec4 color;
    if (view_settings.flux_color_mode != 0) {
      color = flux_colormap_color(view_settings.flux_colormap,
                                  drawn.color_strength, 0.90f);
      if (drawn.flux.weak) {
        color = lerp_color(color, ui_theme::rgb(185, 134, 252, color.w),
                           0.35f);
      }
    } else {
      const float *source_color = drawn.flux.weak ? view_settings.weak_flux_color
                                                  : view_settings.flux_color;
      color = ImVec4(source_color[0], source_color[1], source_color[2],
                     source_color[3]);
    }
    if (dim_for_selection) {
      color.w *= 0.24f;
    } else if (connected_to_selection) {
      color.w = std::min(color.w + 0.12f, 1.0f);
    }
    const ImU32 arrow_color = ImGui::GetColorU32(color);

    draw_list->AddLine(drawn.start, drawn.end,
                       ImGui::GetColorU32(ui_theme::rgb(10, 14, 17, 0.82f)),
                       thickness + 2.2f);
    draw_list->AddLine(drawn.start, drawn.end, arrow_color, thickness);
    const ImVec2 back(drawn.end.x - ux * head_len, drawn.end.y - uy * head_len);
    const ImVec2 normal(-uy, ux);
    draw_list->AddTriangleFilled(
        drawn.end,
        ImVec2(back.x + normal.x * head_half_width,
               back.y + normal.y * head_half_width),
        ImVec2(back.x - normal.x * head_half_width,
               back.y - normal.y * head_half_width),
        arrow_color);
  }

  if (view_settings.show_fluxes && view_settings.flux_color_mode != 0 &&
      !drawn_fluxes.empty()) {
    const float bar_width = 180.0f;
    const float bar_height = 12.0f;
    const float pad = 14.0f;
    const ImVec2 child_pos = ImGui::GetWindowPos();
    const ImVec2 child_size = ImGui::GetWindowSize();
    const bool right = view_settings.flux_colorbar_position == 1 ||
                       view_settings.flux_colorbar_position == 3;
    const bool bottom = view_settings.flux_colorbar_position == 2 ||
                        view_settings.flux_colorbar_position == 3;
    const float bar_x = right ? child_pos.x + child_size.x - bar_width - pad
                              : child_pos.x + pad;
    const float bar_y = bottom ? child_pos.y + child_size.y - 44.0f
                               : child_pos.y + 26.0f;
    const ImVec2 bar_min(bar_x, bar_y);
    const ImVec2 bar_max(bar_min.x + bar_width, bar_min.y + bar_height);
    const int segments = 48;
    for (int i = 0; i < segments; ++i) {
      const float t0 = static_cast<float>(i) / static_cast<float>(segments);
      const float t1 = static_cast<float>(i + 1) / static_cast<float>(segments);
      const ImVec2 p0(bar_min.x + t0 * bar_width, bar_min.y);
      const ImVec2 p1(bar_min.x + t1 * bar_width + 1.0f, bar_max.y);
      draw_list->AddRectFilled(
          p0, p1,
          ImGui::GetColorU32(
              flux_colormap_color(view_settings.flux_colormap, t0, 0.95f)));
    }
    draw_list->AddRect(bar_min, bar_max, ImGui::GetColorU32(ImGuiCol_Text));
    const bool log_scale = view_settings.flux_color_mode == 2;
    const char *metric_label =
        view_settings.flux_arrow_metric == 2
            ? (log_scale ? "log10 Rate" : "Rate")
        : view_settings.flux_arrow_metric == 1
            ? (log_scale ? "log10 |dX/dt|" : "|dX/dt|")
            : (log_scale ? "log10 |dY/dt|" : "|dY/dt|");
    char min_label[32];
    char max_label[32];
    if (log_scale) {
      std::snprintf(min_label, sizeof(min_label), "%.2f", flux_min_log);
      std::snprintf(max_label, sizeof(max_label), "%.2f", flux_max_log);
    } else {
      std::snprintf(min_label, sizeof(min_label), "%.1e", flux_min_value);
      std::snprintf(max_label, sizeof(max_label), "%.1e", flux_max_value);
    }
    draw_list->AddText(ImVec2(bar_min.x, bar_min.y - 15.0f),
                       ImGui::GetColorU32(ImGuiCol_Text), metric_label);
    draw_list->AddText(ImVec2(bar_min.x, bar_max.y + 3.0f),
                       ImGui::GetColorU32(ImGuiCol_Text), min_label);
    const ImVec2 max_label_size = ImGui::CalcTextSize(max_label);
    draw_list->AddText(ImVec2(bar_max.x - max_label_size.x, bar_max.y + 3.0f),
                       ImGui::GetColorU32(ImGuiCol_Text), max_label);
  }

  if (hovered_isotope_ >= 0 && hovered_isotope_ < app_state_->num_species()) {
    const auto &names = app_state_->get_species_names();
    if (hovered_isotope_ < static_cast<int>(names.size())) {
      ImGui::SetTooltip("%s", names[hovered_isotope_].c_str());
    }
  }

  ImGui::Dummy(virtual_size);
  ImGui::EndChild();
}

ImVec4 NuclideChart::get_abundance_color(int index) const {
  if (!app_state_) {
    return ui_theme::border();
  }

  if (!app_state_->view_settings().show_abundance) {
    return ui_theme::rgb(58, 92, 102);
  }

  const auto &xnuc = app_state_->integration_settings().xnuc;
  if (index < 0 || index >= static_cast<int>(xnuc.size())) {
    return ui_theme::rgb(33, 41, 46);
  }

  const double x = xnuc[index];
  if (x <= 0.0) {
    return ui_theme::rgb(24, 30, 34);
  }

  const double lmin = -20.0;
  const double lmax = 0.0;
  const double lx = std::log10(x);
  const float t =
      static_cast<float>(std::clamp((lx - lmin) / (lmax - lmin), 0.0, 1.0));

  if (t < 0.5f) {
    const float u = t / 0.5f;
    return lerp_color(ui_theme::rgb(33, 43, 64),
                      ui_theme::rgb(47, 158, 180), u);
  }
  const float u = (t - 0.5f) / 0.5f;
  return lerp_color(ui_theme::rgb(47, 158, 180),
                    ui_theme::rgb(255, 218, 94), u);
}

int NuclideChart::get_single_selected_isotope() const {
  for (size_t i = 0; i < selected_.size(); ++i) {
    if (selected_[i]) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

} // namespace imyann
