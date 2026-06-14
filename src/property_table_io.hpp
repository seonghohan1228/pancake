#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "config.hpp"

namespace FluidProperties {

/// Load a rectangular `P_MPa, T_C, <value_header>` CSV into a PropertyTable2D,
/// converting to SI (Pa, K, value * value_scale). Columns are matched by header
/// name (order independent); a UTF-8 BOM, CRLF/LF line endings, blank lines, and
/// '#' comment lines are tolerated. On any structural problem `error_out` is set
/// to a message naming the file and offending line and an empty table is
/// returned. The grid resolution is taken from the file (no compiled-in size),
/// so refining/coarsening the CSV needs no recompile.
PropertyTable2D load_table_2d_csv(const std::filesystem::path& path,
                                  const std::string& value_header,
                                  double value_scale,
                                  std::string& error_out);

/// Resolve and load every configured `*_table_file` into cfg's 2-D tables.
/// Relative paths resolve against `cfg.config_dir` first, then the CWD. Any load
/// failure is appended to `errors` (the caller aborts startup, as for validate()).
void load_config_property_tables(SimulationConfig& cfg, std::vector<std::string>& errors);

}  // namespace FluidProperties
