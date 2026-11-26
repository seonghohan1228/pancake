#pragma once

// Workflow-node status lights and their one-sentence reasons.
//
//   Green  - complete, no issue
//   Yellow - input needed (or an acknowledged range issue)
//   Red    - conflict or unsupported combination; the case will not run
//
// Run is enabled only when every setup node is green. Red reasons include the
// solver's own SimulationConfig::validate() messages, mapped to the node that
// owns the offending key. No toolkit dependency: unit-tested.

#include "case_model.hpp"
#include "estimates.hpp"
#include "fluid_library.hpp"

#include <array>
#include <string>
#include <vector>

namespace gui {

enum class Node { Geometry, Fluid, Operating, Boundaries, Mesh, Solution, Results };
constexpr int kNodeCount = 7;
const char* node_name(Node node);

enum class Light { Idle, Green, Yellow, Red };

struct NodeStatus {
    Light light = Light::Green;
    std::vector<std::string> reasons;  // first one is shown in the tree tooltip
};

using StatusArray = std::array<NodeStatus, kNodeCount>;

struct StatusInputs {
    const CaseModel* model = nullptr;
    const FluidPair* pair = nullptr;       // null when no pair is loaded
    std::string pair_error;                // load failure of the chosen pair
    const Estimates* estimates = nullptr;
    const SimulationConfig* config = nullptr;  // mapped config, tables attached
    bool allow_out_of_range = false;       // user accepted a state outside the fluid table
    bool load_table_found = true;          // transient: table file exists and parses
    std::string load_table_error;
    Light results_light = Light::Idle;     // set by the run outcome
    std::string results_reason;
};

/// Computes every node's light. Setup nodes only; Results is copied from inputs.
StatusArray compute_status(const StatusInputs& in);

/// True when all setup nodes (Geometry .. Solution) are green.
bool ready_to_run(const StatusArray& status);

/// Plain description of the solver form for the read-only Expanded line,
/// e.g. "FBNS (complementarity) with dissolved-gas pressure plateau; steady
/// load balance (Newton); isothermal".
std::string solver_form_text(const CaseModel& model, const FluidPair* pair);

/// Node that owns a solver validation message (by the config key it names).
Node node_for_message(const std::string& message);

}  // namespace gui
