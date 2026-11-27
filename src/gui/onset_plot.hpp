#pragma once

// Gaseous-cavitation onset chart (deck p43) as an ImPlot widget: the pressure
// drop ratio Phi(S) of the gas-free runs on log axes, the gas-release region
// below the curve shaded in the high-risk colour, and case points coloured by
// their verdict. Used by the estimates box, the batch window and the report.

#include "estimates.hpp"
#include "imgui.h"

#include <string>
#include <vector>

namespace ui {

struct OnsetMarker {
    std::string label;
    double sommerfeld = 0.0;
    double margin_ratio = 0.0;  // (p_amb - p_sat) / P
    gui::OnsetVerdict verdict = gui::OnsetVerdict::Unknown;
};

/// Draws the chart in a region of `size` (use -1 for the available width).
/// `compact` hides the legend and axis titles (mini chart in the estimates box).
void onset_chart(const char* id, ImVec2 size, const std::vector<OnsetMarker>& markers,
                 bool compact = false);

}  // namespace ui
