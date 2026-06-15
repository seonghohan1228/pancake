#pragma once

#include "imgui.h"

/// Font set rebuilt whenever the DPI scale changes.
struct Fonts {
    ImFont* sans = nullptr;     // body text
    ImFont* heading = nullptr;  // section headers (medium weight, slightly larger)
    ImFont* mono = nullptr;     // log viewer / numeric readouts
};

namespace theme {

// Light theme: single blue accent for primary actions/highlights, neutral
// grays elsewhere, color reserved for meaning (state, warnings, data).
inline constexpr ImVec4 kAccent{0.15f, 0.42f, 0.80f, 1.0f};
inline constexpr ImVec4 kAccentDim{0.35f, 0.57f, 0.86f, 1.0f};
inline constexpr ImVec4 kOk{0.10f, 0.52f, 0.28f, 1.0f};
inline constexpr ImVec4 kWarn{0.66f, 0.45f, 0.03f, 1.0f};
inline constexpr ImVec4 kError{0.78f, 0.16f, 0.16f, 1.0f};
inline constexpr ImVec4 kMuted{0.44f, 0.47f, 0.52f, 1.0f};
inline constexpr ImVec4 kInvalidBg{0.99f, 0.89f, 0.89f, 1.0f};

/// Applies the neutral-gray theme, scaled for the given DPI factor.
void apply_style(float dpi_scale);

/// Rebuilds the font atlas from embedded fonts at the given DPI scale.
/// Must be followed by a backend device-objects rebuild.
Fonts build_fonts(float dpi_scale);

}  // namespace theme
