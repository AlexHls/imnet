#include "ui_plot_helpers.h"
#include "implot_internal.h"
#include <iostream>
#include <stdexcept>

static void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}

int main() {
  try {
    ImGui::CreateContext();
    ImPlot::CreateContext();
    auto &io = ImGui::GetIO();
    io.DisplaySize = ImVec2(900, 700);
    io.DeltaTime = 1.0f / 60;
    io.IniFilename = nullptr;
    io.ConfigMacOSXBehaviors = false;
    unsigned char *pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);

    imyann::LogPlotFloor floor{1e-99, 1e-99};
    const auto edit_frame = [&](bool focus = false) {
      ImGui::NewFrame();
      ImGui::SetNextWindowSize(ImVec2(850, 650));
      ImGui::Begin("Test");
      if (focus) ImGui::SetKeyboardFocusHere();
      floor.draw("Log floor");
      ImGui::End();
      ImGui::Render();
    };
    edit_frame(true);
    edit_frame();
    io.AddKeyEvent(ImGuiMod_Ctrl, true);
    io.AddKeyEvent(ImGuiKey_A, true);
    edit_frame();
    io.AddKeyEvent(ImGuiKey_A, false);
    io.AddKeyEvent(ImGuiMod_Ctrl, false);
    io.AddKeyEvent(ImGuiKey_Backspace, true);
    edit_frame();
    require(floor.value == 1e-99, "Partial input changed the committed floor");
    io.AddKeyEvent(ImGuiKey_Backspace, false);
    io.AddInputCharactersUTF8("1e-8");
    edit_frame();
    io.AddKeyEvent(ImGuiKey_Enter, true);
    edit_frame();
    io.AddKeyEvent(ImGuiKey_Enter, false);
    edit_frame();
    require(floor.value == 1e-8, "Scientific-notation edit was not committed");

    double committed = 1.0;
    const auto number_frame = [&](bool focus = false) {
      ImGui::NewFrame();
      ImGui::Begin("Test");
      if (focus) ImGui::SetKeyboardFocusHere();
      imyann::input_double_committed("Density", &committed);
      ImGui::End();
      ImGui::Render();
    };
    number_frame(true);
    number_frame();
    io.AddKeyEvent(ImGuiMod_Ctrl, true);
    io.AddKeyEvent(ImGuiKey_A, true);
    number_frame();
    io.AddKeyEvent(ImGuiKey_A, false);
    io.AddKeyEvent(ImGuiMod_Ctrl, false);
    io.AddKeyEvent(ImGuiKey_Backspace, true);
    number_frame();
    require(committed == 1.0, "Partial numeric edit escaped into simulation state");
    io.AddKeyEvent(ImGuiKey_Backspace, false);
    io.AddInputCharactersUTF8("2e-8");
    number_frame();
    require(committed == 1.0, "Unfinished scientific notation changed the state");
    io.AddKeyEvent(ImGuiKey_Enter, true);
    number_frame();
    io.AddKeyEvent(ImGuiKey_Enter, false);
    number_frame();
    require(committed == 2e-8, "Numeric field did not commit scientific notation");

    // Run real ImPlot frames: fitting, a changed floor, manual ranges,
    // and off-screen cursors must not introduce near-zero axis limits.
    for (int frame = 0; frame < 6; ++frame) {
      const double x_floor = frame < 2 ? 1e-9 : 1e-4;
      const double y_floor = frame < 2 ? 1e-12 : 1e-6;
      ImGui::NewFrame();
      ImGui::Begin("Test");
      if (frame < 4) ImPlot::SetNextAxesToFit();
      require(ImPlot::BeginPlot("Abundances", ImVec2(700, 400)), "Plot did not open");
      ImPlot::SetupAxes("t", "X");
      imyann::setup_log_plot_axis(ImAxis_X1, x_floor);
      imyann::setup_log_plot_axis(ImAxis_Y1, y_floor);
      if (frame == 4) {
        ImPlot::SetupAxisLimits(ImAxis_X1, 1e-2, 1e2, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 1e-4, 0.5, ImPlotCond_Always);
      }
      const double xs[] = {x_floor, 1, 1e4};
      const double ys[] = {y_floor, 0.1, 1};
      ImPlot::PlotLine("p", xs, ys, 3);
      imyann::plot_current_time("current", frame % 2 ? 1e100 : 0);
      if (frame >= 4) {
        // The visible cursor must span only the plot rectangle, not Y=0.
        auto *draw = ImPlot::GetPlotDrawList();
        const int first_vertex = draw->VtxBuffer.Size;
        const auto pos = ImPlot::GetPlotPos();
        const auto size = ImPlot::GetPlotSize();
        imyann::plot_current_time("visible current", 1.0);
        require(draw->VtxBuffer.Size > first_vertex, "Visible cursor was not drawn");
        for (int i = first_vertex; i < draw->VtxBuffer.Size; ++i) {
          const float y = draw->VtxBuffer[i].pos.y;
          require(y >= pos.y - 2 && y <= pos.y + size.y + 2,
                  "Cursor extends beyond the visible Y range");
        }
      }
      auto *plot = ImPlot::GetCurrentPlot();
      ImPlot::EndPlot();
      const auto &x = plot->Axes[ImAxis_X1].Range;
      const auto &y = plot->Axes[ImAxis_Y1].Range;
      require(x.Min >= x_floor && x.Max < 1e5, "X fit escaped floor/data bounds");
      require(y.Min >= y_floor && y.Max < 2, "Y fit escaped floor/data bounds");
      if (frame >= 4) {
        require(x.Min == 1e-2 && x.Max == 1e2, "Manual X limits changed");
        require(y.Min == 1e-4 && y.Max == 0.5, "Manual Y limits changed");
      }
      ImGui::End();
      ImGui::Render();
    }
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    std::cout << "Log-floor editing, fit bounds, and cursor exclusion passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
