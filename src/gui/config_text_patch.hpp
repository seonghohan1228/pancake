#pragma once

// Text-preserving save of a loaded config.
//
// A config written by hand or by a case driver carries comments, blank lines,
// CRLF line endings and paths (e.g. afml-r6 paths with Korean folder names)
// that a full rewrite would lose or re-encode. The GUI therefore saves a loaded
// file by patching its text: lines of keys whose value did not change stay
// byte-identical, changed keys are rewritten in place, keys that were absent
// are appended under one GUI comment, and `#@gui` metadata lines are updated.
// No toolkit dependency: unit-tested.

#include "case_model.hpp"
#include "config.hpp"

#include <string>

namespace gui {

/// Returns `original_text` changed so that it parses to `updated`.
/// `original` is the config parsed from `original_text`.
std::string patch_config_text(const std::string& original_text,
                              const SimulationConfig& original,
                              const SimulationConfig& updated,
                              const GuiMeta& meta);

/// Full text for a new case: the solver's writer plus the `#@gui` lines.
std::string new_config_text(const SimulationConfig& cfg, const GuiMeta& meta);

}  // namespace gui
