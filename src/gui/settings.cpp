#include "settings.hpp"

#include "win_util.hpp"

#include <shlobj.h>

#include <fstream>
#include <sstream>

AppPaths AppPaths::resolve() {
    AppPaths paths;
    PWSTR local_app_data = nullptr;
    if (SUCCEEDED(SHGetKnownFolderPath(FOLDERID_LocalAppData, KF_FLAG_CREATE, nullptr,
                                       &local_app_data))) {
        paths.data_dir = std::filesystem::path(local_app_data) / L"PancakeGui";
        CoTaskMemFree(local_app_data);
    } else {
        // Extremely unlikely; fall back to the user temp dir rather than the
        // exe dir, which may be read-only.
        paths.data_dir = std::filesystem::temp_directory_path() / L"PancakeGui";
    }
    std::error_code ec;
    std::filesystem::create_directories(paths.data_dir, ec);
    paths.imgui_ini = paths.data_dir / L"imgui.ini";
    paths.settings_file = paths.data_dir / L"settings.txt";
    return paths;
}

void UserSettings::load(const std::filesystem::path& file) {
    std::ifstream input(file);
    if (!input.is_open()) return;
    std::string line;
    while (std::getline(input, line)) {
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = line.substr(0, eq);
        const std::string value = line.substr(eq + 1);
        if (key == "solver_path") solver_path = value;
        else if (key == "last_config") last_config = value;
        else if (key.rfind("unit.", 0) == 0) unit_choices[key.substr(5)] = std::atoi(value.c_str());
        else if (key == "mpi_ranks") mpi_ranks = std::max(1, std::atoi(value.c_str()));
        else if (key == "window_width") window_width = std::max(640, std::atoi(value.c_str()));
        else if (key == "window_height") window_height = std::max(480, std::atoi(value.c_str()));
        else if (key == "window_maximized") window_maximized = (value == "1");
    }
}

void UserSettings::save(const std::filesystem::path& file) const {
    std::ofstream output(file, std::ios::trunc);
    if (!output.is_open()) return;
    output << "solver_path=" << solver_path << "\n"
           << "last_config=" << last_config << "\n"
           << "mpi_ranks=" << mpi_ranks << "\n"
           << "window_width=" << window_width << "\n"
           << "window_height=" << window_height << "\n"
           << "window_maximized=" << (window_maximized ? 1 : 0) << "\n";
    for (const auto& [key, index] : unit_choices) {
        output << "unit." << key << "=" << index << "\n";
    }
}
