#include "ui_theme.h"
#include "implot.h"

namespace imyann::ui_theme {

void apply() {
  ImGui::StyleColorsDark();
  ImGuiStyle &style = ImGui::GetStyle();
  style.WindowPadding = ImVec2(12.0f, 10.0f);
  style.FramePadding = ImVec2(8.0f, 4.0f);
  style.CellPadding = ImVec2(7.0f, 4.0f);
  style.ItemSpacing = ImVec2(8.0f, 6.0f);
  style.ItemInnerSpacing = ImVec2(6.0f, 4.0f);
  style.ScrollbarSize = 13.0f;
  style.GrabMinSize = 10.0f;
  style.WindowRounding = 5.0f;
  style.ChildRounding = 4.0f;
  style.FrameRounding = 3.0f;
  style.PopupRounding = 4.0f;
  style.ScrollbarRounding = 3.0f;
  style.GrabRounding = 3.0f;
  style.TabRounding = 3.0f;
  style.WindowBorderSize = 1.0f;
  style.ChildBorderSize = 1.0f;
  style.FrameBorderSize = 1.0f;
  style.PopupBorderSize = 1.0f;
  style.TabBorderSize = 1.0f;
  style.WindowTitleAlign = ImVec2(0.0f, 0.5f);

  ImVec4 *colors = style.Colors;
  colors[ImGuiCol_Text] = text();
  colors[ImGuiCol_TextDisabled] = text_muted();
  colors[ImGuiCol_WindowBg] = window();
  colors[ImGuiCol_ChildBg] = panel();
  colors[ImGuiCol_PopupBg] = rgb(13, 19, 23, 0.98f);
  colors[ImGuiCol_Border] = border();
  colors[ImGuiCol_BorderShadow] = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
  colors[ImGuiCol_FrameBg] = frame();
  colors[ImGuiCol_FrameBgHovered] = rgb(29, 43, 51);
  colors[ImGuiCol_FrameBgActive] = rgb(35, 56, 66);
  colors[ImGuiCol_TitleBg] = rgb(10, 16, 20);
  colors[ImGuiCol_TitleBgActive] = rgb(17, 33, 40);
  colors[ImGuiCol_TitleBgCollapsed] = rgb(9, 13, 16);
  colors[ImGuiCol_MenuBarBg] = rgb(8, 12, 15);
  colors[ImGuiCol_ScrollbarBg] = rgb(9, 14, 17);
  colors[ImGuiCol_ScrollbarGrab] = rgb(39, 54, 63);
  colors[ImGuiCol_ScrollbarGrabHovered] = rgb(50, 75, 86);
  colors[ImGuiCol_ScrollbarGrabActive] = rgb(63, 100, 113);
  colors[ImGuiCol_CheckMark] = primary();
  colors[ImGuiCol_SliderGrab] = primary(0.82f);
  colors[ImGuiCol_SliderGrabActive] = rgb(117, 226, 236);
  colors[ImGuiCol_Button] = rgb(24, 39, 47);
  colors[ImGuiCol_ButtonHovered] = rgb(31, 63, 74);
  colors[ImGuiCol_ButtonActive] = rgb(38, 92, 106);
  colors[ImGuiCol_Header] = rgb(24, 49, 58);
  colors[ImGuiCol_HeaderHovered] = rgb(31, 74, 86);
  colors[ImGuiCol_HeaderActive] = rgb(39, 97, 111);
  colors[ImGuiCol_Separator] = border();
  colors[ImGuiCol_SeparatorHovered] = primary(0.72f);
  colors[ImGuiCol_SeparatorActive] = primary();
  colors[ImGuiCol_ResizeGrip] = primary(0.18f);
  colors[ImGuiCol_ResizeGripHovered] = primary(0.52f);
  colors[ImGuiCol_ResizeGripActive] = primary(0.86f);
  colors[ImGuiCol_InputTextCursor] = primary();
  colors[ImGuiCol_Tab] = rgb(14, 23, 28);
  colors[ImGuiCol_TabHovered] = rgb(31, 72, 84);
  colors[ImGuiCol_TabSelected] = rgb(23, 46, 55);
  colors[ImGuiCol_TabSelectedOverline] = primary();
  colors[ImGuiCol_TabDimmed] = rgb(12, 18, 22);
  colors[ImGuiCol_TabDimmedSelected] = rgb(20, 32, 38);
  colors[ImGuiCol_TabDimmedSelectedOverline] = primary(0.45f);
  colors[ImGuiCol_DockingPreview] = primary(0.45f);
  colors[ImGuiCol_DockingEmptyBg] = background();
  colors[ImGuiCol_PlotLines] = primary();
  colors[ImGuiCol_PlotLinesHovered] = rgb(141, 234, 241);
  colors[ImGuiCol_PlotHistogram] = warning();
  colors[ImGuiCol_PlotHistogramHovered] = rgb(255, 207, 102);
  colors[ImGuiCol_TableHeaderBg] = rgb(18, 35, 42);
  colors[ImGuiCol_TableBorderStrong] = border_light();
  colors[ImGuiCol_TableBorderLight] = border();
  colors[ImGuiCol_TableRowBg] = rgb(12, 17, 21, 0.58f);
  colors[ImGuiCol_TableRowBgAlt] = rgb(20, 28, 34, 0.56f);
  colors[ImGuiCol_TextLink] = primary();
  colors[ImGuiCol_TextSelectedBg] = primary(0.30f);
  colors[ImGuiCol_TreeLines] = border_light(0.75f);
  colors[ImGuiCol_DragDropTarget] = warning(0.90f);
  colors[ImGuiCol_DragDropTargetBg] = warning(0.16f);
  colors[ImGuiCol_UnsavedMarker] = warning();
  colors[ImGuiCol_NavCursor] = primary();
  colors[ImGuiCol_NavWindowingHighlight] = primary(0.70f);
  colors[ImGuiCol_NavWindowingDimBg] = rgb(0, 0, 0, 0.45f);
  colors[ImGuiCol_ModalWindowDimBg] = rgb(0, 0, 0, 0.58f);

  ImPlot::StyleColorsDark();
  ImPlotStyle &plot = ImPlot::GetStyle();
  plot.PlotBorderSize = 1.0f;
  plot.MinorAlpha = 0.30f;
  plot.PlotPadding = ImVec2(10.0f, 8.0f);
  plot.LabelPadding = ImVec2(5.0f, 5.0f);
  plot.LegendPadding = ImVec2(10.0f, 10.0f);
  plot.LegendInnerPadding = ImVec2(7.0f, 5.0f);
  plot.FitPadding = ImVec2(0.03f, 0.08f);

  ImVec4 *plot_colors = plot.Colors;
  plot_colors[ImPlotCol_FrameBg] = frame();
  plot_colors[ImPlotCol_PlotBg] = panel();
  plot_colors[ImPlotCol_PlotBorder] = border();
  plot_colors[ImPlotCol_LegendBg] = rgb(13, 19, 23, 0.94f);
  plot_colors[ImPlotCol_LegendBorder] = border_light(0.85f);
  plot_colors[ImPlotCol_LegendText] = text();
  plot_colors[ImPlotCol_TitleText] = text();
  plot_colors[ImPlotCol_InlayText] = text_muted();
  plot_colors[ImPlotCol_AxisText] = text_muted();
  plot_colors[ImPlotCol_AxisGrid] = border_light(0.40f);
  plot_colors[ImPlotCol_AxisTick] = border_light(0.68f);
  plot_colors[ImPlotCol_AxisBg] = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
  plot_colors[ImPlotCol_AxisBgHovered] = primary(0.10f);
  plot_colors[ImPlotCol_AxisBgActive] = primary(0.18f);
  plot_colors[ImPlotCol_Selection] = warning(0.26f);
  plot_colors[ImPlotCol_Crosshairs] = primary(0.82f);

  const ImVec4 series_colors[] = {
      primary(),          success(),       warning(),          danger(),
      rgb(185, 134, 252), rgb(56, 189, 248), rgb(244, 114, 182),
      rgb(148, 163, 184)};
  plot.Colormap =
      ImPlot::AddColormap("imnet", series_colors, IM_ARRAYSIZE(series_colors));
}

} // namespace imyann::ui_theme
