#pragma once

#include "imgui.h"

/// Font set rebuilt whenever the DPI scale changes.
struct Fonts {
    ImFont* sans = nullptr;     // body text
    ImFont* heading = nullptr;  // section headers (medium weight, slightly larger)
    ImFont* mono = nullptr;     // log viewer / numeric readouts
};

namespace theme {

// Single accent for primary actions/highlights; everything else neutral.
inline constexpr ImVec4 kAccent{0.26f, 0.59f, 0.98f, 1.0f};
inline constexpr ImVec4 kAccentDim{0.20f, 0.41f, 0.68f, 1.0f};
inline constexpr ImVec4 kOk{0.35f, 0.75f, 0.45f, 1.0f};
inline constexpr ImVec4 kWarn{0.95f, 0.73f, 0.25f, 1.0f};
inline constexpr ImVec4 kError{0.93f, 0.36f, 0.36f, 1.0f};
inline constexpr ImVec4 kMuted{0.58f, 0.61f, 0.66f, 1.0f};
inline constexpr ImVec4 kInvalidBg{0.35f, 0.13f, 0.13f, 1.0f};

/// Applies the neutral-gray theme, scaled for the given DPI factor.
void apply_style(float dpi_scale);

/// Rebuilds the font atlas from embedded fonts at the given DPI scale.
/// Must be followed by a backend device-objects rebuild.
Fonts build_fonts(float dpi_scale);

}  // namespace theme
