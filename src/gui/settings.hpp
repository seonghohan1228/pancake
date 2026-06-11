#pragma once

#include <filesystem>
#include <string>

/// Locations for all runtime writes. Everything lives under
/// %LOCALAPPDATA%\PancakeGui so the exe can run from a read-only folder.
struct AppPaths {
    std::filesystem::path data_dir;       // %LOCALAPPDATA%\PancakeGui
    std::filesystem::path imgui_ini;      // data_dir\imgui.ini (dock layout)
    std::filesystem::path settings_file;  // data_dir\settings.txt

    /// Resolves the LocalAppData directory and creates data_dir if missing.
    static AppPaths resolve();
};

/// User settings persisted across sessions (not the simulation config).
struct UserSettings {
    std::string solver_path;  // optional override; empty = auto-discover
    std::string last_config;  // last opened/saved config file
    int mpi_ranks = 1;
    int window_width = 1480;
    int window_height = 900;
    bool window_maximized = false;

    void load(const std::filesystem::path& file);
    void save(const std::filesystem::path& file) const;
};
