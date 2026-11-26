#include "config_text_patch.hpp"

#include <set>
#include <sstream>
#include <vector>

namespace gui {

namespace {

struct Line {
    std::string text;  // without the line ending
    std::string ending;
};

std::vector<Line> split_lines(const std::string& text) {
    std::vector<Line> lines;
    std::size_t start = 0;
    while (start < text.size()) {
        const std::size_t newline = text.find('\n', start);
        if (newline == std::string::npos) {
            lines.push_back({text.substr(start), ""});
            break;
        }
        std::size_t end = newline;
        std::string ending = "\n";
        if (end > start && text[end - 1] == '\r') {
            --end;
            ending = "\r\n";
        }
        lines.push_back({text.substr(start, end - start), ending});
        start = newline + 1;
    }
    return lines;
}

/// Key of a `key = value` line, or "" for comments and blank lines.
std::string line_key(const std::string& line) {
    const auto start = line.find_first_not_of(" \t");
    if (start == std::string::npos || line[start] == '#') return "";
    const auto eq = line.find('=');
    if (eq == std::string::npos) return "";
    std::string key = line.substr(start, eq - start);
    key.erase(key.find_last_not_of(" \t") + 1);
    return key;
}

bool is_meta_line(const std::string& line) {
    const auto start = line.find_first_not_of(" \t");
    return start != std::string::npos && line.compare(start, 5, "#@gui") == 0;
}

std::string meta_key(const std::string& line) {
    const auto start = line.find("#@gui") + 5;
    const auto eq = line.find('=', start);
    std::string key = line.substr(start, eq == std::string::npos ? std::string::npos : eq - start);
    const auto a = key.find_first_not_of(" \t");
    const auto b = key.find_last_not_of(" \t");
    return a == std::string::npos ? std::string() : key.substr(a, b - a + 1);
}

bool is_inlet_key(const std::string& key) {
    return key == "inlet_circular" || key == "inlet_groove" || key == "inlet_species";
}

/// Inlet lines in writer order (repeated keys), from the solver's writer.
std::vector<std::string> inlet_lines(const SimulationConfig& cfg) {
    std::vector<std::string> lines;
    std::istringstream stream(cfg.to_config_text());
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (is_inlet_key(line_key(line))) lines.push_back(line);
    }
    return lines;
}

}  // namespace

std::string patch_config_text(const std::string& original_text,
                              const SimulationConfig& original,
                              const SimulationConfig& updated,
                              const GuiMeta& meta) {
    const auto old_values = config_key_values(original.to_config_text());
    const auto new_values = config_key_values(updated.to_config_text());
    std::set<std::string> changed;
    bool inlets_changed = false;
    for (const auto& [key, value] : new_values) {
        const auto it = old_values.find(key);
        if (it != old_values.end() && it->second == value) continue;
        if (is_inlet_key(key)) inlets_changed = true;
        else changed.insert(key);
    }
    for (const auto& [key, value] : old_values) {
        if (is_inlet_key(key) && !new_values.count(key)) inlets_changed = true;
    }

    std::vector<Line> lines = split_lines(original_text);
    // Majority line ending, used for every new line.
    std::size_t crlf = 0;
    for (const Line& line : lines) crlf += line.ending == "\r\n";
    const std::string eol = crlf * 2 > lines.size() ? "\r\n" : "\n";
    if (!lines.empty() && lines.back().ending.empty()) lines.back().ending = eol;

    std::vector<Line> out;
    std::set<std::string> written;
    std::set<std::string> meta_written;
    bool inlets_placed = false;
    for (const Line& line : lines) {
        if (is_meta_line(line.text)) {
            const std::string key = meta_key(line.text);
            const auto it = meta.find(key);
            if (it == meta.end() || meta_written.count(key)) continue;  // drop stale metadata
            out.push_back({"#@gui " + key + " = " + it->second, line.ending});
            meta_written.insert(key);
            continue;
        }
        const std::string key = line_key(line.text);
        if (is_inlet_key(key) && inlets_changed) {
            if (!inlets_placed) {
                for (const std::string& inlet : inlet_lines(updated)) out.push_back({inlet, eol});
                inlets_placed = true;
            }
            continue;
        }
        if (!key.empty() && changed.count(key) && !written.count(key)) {
            out.push_back({key + " = " + new_values.at(key), line.ending});
            written.insert(key);
            continue;
        }
        out.push_back(line);
    }

    std::vector<Line> appended;
    for (const std::string& key : changed) {
        if (!written.count(key)) appended.push_back({key + " = " + new_values.at(key), eol});
    }
    if (inlets_changed && !inlets_placed) {
        for (const std::string& inlet : inlet_lines(updated)) appended.push_back({inlet, eol});
    }
    for (const auto& [key, value] : meta) {
        if (!meta_written.count(key)) appended.push_back({"#@gui " + key + " = " + value, eol});
    }
    if (!appended.empty()) {
        out.push_back({"", eol});
        out.push_back({"# Set by the pancake GUI", eol});
        out.insert(out.end(), appended.begin(), appended.end());
    }

    std::string text;
    for (const Line& line : out) text += line.text + line.ending;
    return text;
}

std::string new_config_text(const SimulationConfig& cfg, const GuiMeta& meta) {
    std::string text = "# pancake case written by the pancake GUI. Run: pancake -c <this file>\n";
    for (const auto& [key, value] : meta) text += "#@gui " + key + " = " + value + "\n";
    text += "\n" + cfg.to_config_text();
    return text;
}

}  // namespace gui
