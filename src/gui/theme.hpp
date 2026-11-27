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
// Design-risk colours: verdicts and checks only (green = safe for the design).
inline constexpr ImVec4 kRiskLow = kOk;
inline constexpr ImVec4 kRiskMedium = kWarn;
inline constexpr ImVec4 kRiskHigh = kError;
// Model colours of the deck, only for curves that compare the two models.
inline constexpr ImVec4 kModelVc{0.89f, 0.29f, 0.29f, 1.0f};      // VC only (red)
inline constexpr ImVec4 kModelGcVc{0.18f, 0.72f, 0.42f, 1.0f};    // GC+VC (green)
inline constexpr ImVec4 kHighlight{0.58f, 0.36f, 0.86f, 1.0f};    // highlights (purple)
// Polar pressure fills relative to ambient (not risk colours).
inline constexpr ImVec4 kOverPressure{0.15f, 0.42f, 0.80f, 0.28f};
inline constexpr ImVec4 kUnderPressure{0.93f, 0.55f, 0.18f, 0.35f};

/// Applies the neutral-gray theme, scaled for the given DPI factor.
void apply_style(float dpi_scale);

/// Heading font of the current atlas (section titles, key values); the body
/// font before build_fonts ran.
ImFont* heading_font();

/// Rebuilds the font atlas from embedded fonts at the given DPI scale.
/// Must be followed by a backend device-objects rebuild.
Fonts build_fonts(float dpi_scale);

}  // namespace theme
