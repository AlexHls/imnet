#ifndef IMYANN_UI_THEME_H
#define IMYANN_UI_THEME_H

#include "imgui.h"

namespace imyann::ui_theme {

inline ImVec4 rgb(int r, int g, int b, float a = 1.0f) {
  return ImVec4(static_cast<float>(r) / 255.0f,
                static_cast<float>(g) / 255.0f,
                static_cast<float>(b) / 255.0f, a);
}

inline ImVec4 background() { return rgb(10, 14, 17); }
inline ImVec4 window() { return rgb(15, 21, 26); }
inline ImVec4 panel() { return rgb(11, 16, 20); }
inline ImVec4 frame() { return rgb(20, 29, 35); }
inline ImVec4 border() { return rgb(47, 64, 72); }
inline ImVec4 border_light(float alpha = 1.0f) {
  return rgb(70, 91, 101, alpha);
}
inline ImVec4 text() { return rgb(224, 234, 238); }
inline ImVec4 text_muted() { return rgb(130, 151, 160); }
inline ImVec4 primary(float alpha = 1.0f) { return rgb(75, 202, 216, alpha); }
inline ImVec4 success() { return rgb(118, 220, 154); }
inline ImVec4 warning(float alpha = 1.0f) { return rgb(239, 178, 71, alpha); }
inline ImVec4 danger(float alpha = 1.0f) { return rgb(238, 100, 86, alpha); }

void apply();

} // namespace imyann::ui_theme

#endif // IMYANN_UI_THEME_H
